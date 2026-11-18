"""Restart both continuation implementations from the same t=20 checkpoint."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import time

from compare_library import compare

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'legacy' / 'finite_difference_modified_impurity_Thomas_generalized_without_Matrix_hermitian_vel_contd_new_cluster_upwind_scheme.cpp'
PREVIOUS = ROOT / 'results_L500_v05_library_comparison'
CHECKPOINT = PREVIOUS / 'library/configurations.tar.gz'
CHECKPOINT_MEMBER = 'phi_2000.dat'
OUTPUT = ROOT / 'results_L500_v05_continuation_comparison'
PARAMETERS = dict(length=500.0, dx=0.001, dt=0.01, t_initial=20.0, t_final=25.0,
                  V_0=0.5, mass_ratio=0.333, G=1.0, gamma=0.1,
                  picard_tolerance=1e-8, max_picard_iterations=10,
                  config_save_every=50, velocity_save_every=10)
TOTAL_STEPS = int(round((PARAMETERS['t_final']-PARAMETERS['t_initial'])/PARAMETERS['dt']))
BUFFER_SIZE = 1024 * 1024
PRECISION = 17
READINESS_POLL_SECONDS = 0.25
PROGRESS_POLL_SECONDS = 30


def checkpoint_velocity():
    manifest = json.loads((PREVIOUS/'manifest.json').read_text())
    for key, value in PARAMETERS.items():
        if key not in ('t_initial', 't_final'):
            assert manifest['parameters'][key] == value, key
    rows = [list(map(float, line.split())) for line in (PREVIOUS/'library/velocity.dat').read_text().splitlines()
            if line and not line.startswith('#')]
    matches = [row for row in rows if row[0] == PARAMETERS['t_initial']]
    assert len(matches) == 1
    return matches[0][1:]


def prepare(eigen, archive_prefix):
    velocity = checkpoint_velocity()
    OUTPUT.mkdir(exist_ok=False)
    (OUTPUT/'reference').mkdir()
    # Only the baseline needs a plain file, as required by its original reader.
    digest = hashlib.sha256()
    with tarfile.open(CHECKPOINT, 'r|gz') as archive:
        for member in archive:
            if member.name == CHECKPOINT_MEMBER:
                with archive.extractfile(member) as stream, (OUTPUT/'reference/checkpoint.dat').open('wb') as target:
                    while chunk := stream.read(BUFFER_SIZE):
                        digest.update(chunk)
                        target.write(chunk)
                break
        else:
            raise RuntimeError('t=20 checkpoint is missing')
    original = SOURCE.read_text()
    reference = original
    replacements = {
        'double m_tilde=1.0/8.0, G_tilde=-1.5, gamma=0.1;':
            f"double m_tilde={PARAMETERS['mass_ratio']}, G_tilde={PARAMETERS['G']}, gamma={PARAMETERS['gamma']};",
        'double delta_t=0.001, delta_x = (1.0/10000.0),t_init=142.5;':
            f"double delta_t={PARAMETERS['dt']}, delta_x={PARAMETERS['dx']}, t_init={PARAMETERS['t_initial']};",
        'complex<double> V_0(2.0,0.0),V;': f"complex<double> V_0({PARAMETERS['V_0']},0.0),V;",
        'int Length=2000,i;': f"int Length={int(PARAMETERS['length'])},i;",
        't_fin =20000.0': f"t_fin ={PARAMETERS['t_final']}",
        'int k, K_fin =12;': f"int k, K_fin ={PARAMETERS['max_picard_iterations']};",
        'double err_iter = 10e-9, err_config = 10e-5;':
            f"double err_iter={PARAMETERS['picard_tolerance']}, err_config=10e-5;",
        'V=0.274999+I*0.0;': f'V=Complex({velocity[0]!r},{velocity[1]!r});',
        'for(t=t_init;t<=t_fin;t=t+delta_t)':
            'for(t=t_init;count<=TOTAL_STEPS;t=t_init+(count-1)*delta_t)',
        'myfile_velocity<<': 'myfile_velocity<<setprecision(OUTPUT_PRECISION)<<',
        'myfile_velocity.open(filename_velocity, ios::app);':
            'myfile_velocity.open("velocity.dat", ios::app);',
    }
    for before, after in replacements.items():
        assert before in reference, before
        reference = reference.replace(before, after)
    filename_start = reference.index('    sprintf(filename_input_Phi,')
    filename_end = reference.index('\n', filename_start)
    reference = reference[:filename_start] + '    sprintf(filename_input_Phi,"checkpoint.dat");' + reference[filename_end:]
    # Keep the original field reader and the complete numerical update.
    # Replace legacy stopping/output only, retaining full precision checkpoints.
    constants = f'''#include "snapshot_archive.hpp"
constexpr int TOTAL_STEPS={TOTAL_STEPS};
constexpr int VELOCITY_EVERY={PARAMETERS['velocity_save_every']};
constexpr int CONFIG_EVERY={PARAMETERS['config_save_every']};
constexpr int OUTPUT_PRECISION={PRECISION};
'''
    reference = constants + reference
    instrumentation = '''
    input_Phi.close();
    SnapshotArchive archive("configurations.tar.gz");
    auto snapshot = [&](int step, const VectorXcd& field) {
        ostringstream text;
        text<<setprecision(OUTPUT_PRECISION);
        for (int j=0;j<space_points;++j)
            text<<j*delta_x-Length/2.0<<" "<<real(field(j))<<" "<<imag(field(j))<<"\\n";
        archive.add("phi_"+to_string(step)+".dat",text.str());
    };
    snapshot(0,Phi_n);
    { ofstream ready("checkpoint_loaded"); ready<<"loaded\\n"; }
'''
    reference = reference.replace('    int count = 1;', instrumentation+'\n    int count = 1;')
    reference = reference[:reference.index('         //Density writing')] + '''
         if (k>K_fin || !Phi_n1.allFinite() || !isfinite(real(V)) || !isfinite(imag(V)))
             throw runtime_error("Original continuation did not converge or became nonfinite");
         if (count % VELOCITY_EVERY == 0) {
             myfile_velocity.open("velocity.dat",ios::app);
             myfile_velocity<<setprecision(OUTPUT_PRECISION)<<t_init+count*delta_t<<" "<<real(V)<<" "<<imag(V)<<endl;
             myfile_velocity.close();
         }
         if (count % CONFIG_EVERY == 0) snapshot(count,Phi_n1);
         Phi_n=Phi_n1;
         count++;
     }
    archive.close();
    ofstream final_field("final_field.bin",ios::binary);
    final_field.write(reinterpret_cast<const char*>(Phi_n.data()),Phi_n.size()*sizeof(Complex));
    return 0;
}
'''
    (OUTPUT/'reference.cpp').write_text(reference)
    assignments = '\n'.join(f'    p.{key} = {value!r};' for key, value in PARAMETERS.items())
    (OUTPUT/'library.cpp').write_text('''#include <gpe/GPE_solver_contd.hpp>
#include <fstream>
int main(int argc, char** argv) {
    if (argc!=2) return 1;
    gpe::SimulationParameters p;
''' + assignments + f'''
    p.initial_velocity = std::complex<double>({velocity[0]!r},{velocity[1]!r});
    p.output_directory = "library";
    const auto result = gpe::GPE_solver_contd(p, argv[1], "{CHECKPOINT_MEMBER}");
    std::ofstream field("library/final_field.bin",std::ios::binary);
    field.write(reinterpret_cast<const char*>(result.final_field.data()),
                result.final_field.size()*sizeof(std::complex<double>));
}}
''')
    manifest = dict(parameters=PARAMETERS, checkpoint_velocity=velocity,
                    checkpoint_file=str(CHECKPOINT), checkpoint_member=CHECKPOINT_MEMBER,
                    checkpoint_sha256=digest.hexdigest(),
                    source_sha256={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                   for name in (str(SOURCE.relative_to(ROOT)), 'src/solver.cpp', 'src/continuation.cpp', 'src/numerics.hpp')},
                    eigen=str(eigen), libarchive_prefix=str(archive_prefix),
                    compiler=subprocess.check_output(['g++','--version'],text=True),
                    reference_adjustments='Parameters, checkpoint filename/velocity, integer step count, output and stopping only. Original reader, numerical routines and update retained. Temporary plain input removed after loading. Library reads original compressed archive directly.')
    (OUTPUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    common = ['g++','-std=c++17','-O3','-DNDEBUG','-I'+str(eigen)]
    subprocess.run(common+['-I'+str(ROOT/'src'),str(OUTPUT/'reference.cpp'),'-lz','-o',str(OUTPUT/'reference_solver')],check=True)
    subprocess.run(common+['-I'+str(ROOT/'include'),'-I'+str(archive_prefix/'include'),
                          str(OUTPUT/'library.cpp'),str(ROOT/'src/solver.cpp'),str(ROOT/'src/continuation.cpp'),
                          '-L'+str(archive_prefix/'lib'),'-lz','-larchive','-o',str(OUTPUT/'library_solver')],check=True)
    print(f'Prepared t=20 checkpoint, velocity={velocity}; output: {OUTPUT}',flush=True)


def run(archive_prefix):
    env = dict(os.environ)
    env['LD_LIBRARY_PATH'] = str(archive_prefix/'lib') + os.pathsep + env.get('LD_LIBRARY_PATH','')
    jobs = []
    for name, cwd, arguments in [('reference', OUTPUT/'reference', []), ('library', OUTPUT, [str(CHECKPOINT)])]:
        log = (OUTPUT/f'{name}.log').open('w')
        process = subprocess.Popen([str(OUTPUT/f'{name}_solver'),*arguments],cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
        jobs.append((name,process,log))
    # The original has finished reading before publishing its initial snapshot.
    while not (OUTPUT/'reference/checkpoint_loaded').exists() and jobs[0][1].poll() is None:
        time.sleep(READINESS_POLL_SECONDS)
    if (OUTPUT/'reference/checkpoint_loaded').exists():
        (OUTPUT/'reference/checkpoint.dat').unlink()
    while any(process.poll() is None for _,process,_ in jobs):
        time.sleep(PROGRESS_POLL_SECONDS)
    failed = []
    for name, process, log in jobs:
        log.close()
        if process.returncode:
            failed.append((name,process.returncode))
    if failed:
        raise RuntimeError(f'Solver failures: {failed}; partial data retained')


def report():
    compare(output=OUTPUT, parameters=PARAMETERS,
            initial_description='Initial state: saved phi_2000.dat at t=20; periodic boundaries.',
            reference_description='Reference: original continuation reader, numerical routines and update.',
            library_label='Library GPE_solver_contd')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eigen',type=Path)
    parser.add_argument('--libarchive-prefix',type=Path)
    parser.add_argument('--compare-only',action='store_true')
    args = parser.parse_args()
    if not args.compare_only:
        if args.eigen is None or not (args.eigen/'Eigen/Core').is_file():
            parser.error('--eigen must contain Eigen/Core')
        if args.libarchive_prefix is None or not (args.libarchive_prefix/'include/archive.h').is_file():
            parser.error('--libarchive-prefix must contain include/archive.h')
        prepare(args.eigen,args.libarchive_prefix)
        started = time.monotonic()
        run(args.libarchive_prefix)
        print(f'Both continuations finished in {time.monotonic()-started:.1f} seconds',flush=True)
    report()
