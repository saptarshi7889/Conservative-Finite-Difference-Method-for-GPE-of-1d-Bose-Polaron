"""Compare the original fresh-run evolution with GPE_solver every 10 steps.

Uses existing compiler, Eigen, zlib, NumPy and Matplotlib; installs nothing.
Original research sources and library sources are never modified.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import time

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'legacy' / 'finite_difference_modified_impurity_Thomas_generalized_without_Matrix_hermitian_vel_new_cluster_upwind_scheme_corrected.cpp'
OUTPUT = ROOT / 'results_L500_v05_library_comparison'
PARAMETERS = dict(length=500.0, dx=0.001, dt=0.01, t_final=25.0,
                  V_0=0.5, mass_ratio=0.333, G=1.0, gamma=0.1,
                  picard_tolerance=1e-8, max_picard_iterations=10,
                  config_save_every=50, velocity_save_every=10)
FIGURE_SIZE = (9, 7)
DPI = 160


def prepare(eigen):
    OUTPUT.mkdir(exist_ok=False)
    original = SOURCE.read_text()
    reference = original
    replacements = {
        'double m_tilde=1.0/8.0, G_tilde=-1.5, gamma=0.1;': 'double m_tilde=0.333, G_tilde=1.0, gamma=0.1;',
        'double delta_t=0.001, delta_x = (1.0/10000.0);': 'double delta_t=0.01, delta_x=0.001;',
        'complex<double> V_0(2.0,0.0),V;': 'complex<double> V_0(0.5,0.0),V;',
        'int Length=2000,i;': 'int Length=500,i;',
        't_fin =20000.0': 't_fin =25.0',
        'int k, K_fin =12;': 'int k, K_fin =10;',
        'for(t=t_init;t<=t_fin;t=t+delta_t)': 'for(t=t_init;count<=2500;t=t_init+(count-1)*delta_t)',
        'myfile_velocity<<': 'myfile_velocity<<setprecision(17)<<',
    }
    for before, after in replacements.items():
        assert before in reference, before
        reference = reference.replace(before, after)
    # Retain the original complete numerical step, and replace only legacy
    # equilibrium/CPU stopping and output scheduling with the requested cadence.
    reference = '#include "snapshot_archive.hpp"\n' + reference
    initial_record = '''
    SnapshotArchive archive("configurations.tar.gz");
    auto snapshot = [&](int step, const VectorXcd& field) {
        ostringstream text;
        text<<setprecision(17);
        for (int j=0;j<space_points;++j)
            text<<j*delta_x-Length/2.0<<" "<<real(field(j))<<" "<<imag(field(j))<<"\\n";
        archive.add("phi_"+to_string(step)+".dat",text.str());
    };
    snapshot(0,Phi_n);
'''
    reference = reference.replace('    int count = 1;',initial_record+'\n    int count = 1;')
    cut = reference.index('         //Density writing')
    reference = reference[:cut] + '''
         if (k>K_fin || !Phi_n1.allFinite() || !isfinite(real(V)) || !isfinite(imag(V)))
             throw runtime_error("Reference did not converge or became nonfinite");
         if (count % 10 == 0) {
             myfile_velocity.open("velocity.dat",ios::app);
             myfile_velocity<<setprecision(17)<<count*delta_t<<" "<<real(V)<<" "<<imag(V)<<endl;
             myfile_velocity.close();
         }
         if (count % 50 == 0) snapshot(count,Phi_n1);
         Phi_n=Phi_n1;
         count++;
     }
    archive.close();
    ofstream final_field("final_field.bin",ios::binary);
    final_field.write(reinterpret_cast<const char*>(Phi_n.data()),Phi_n.size()*sizeof(Complex));
    return 0;
}
'''
    # Initial velocity must go into the same fixed filename as later rows.
    reference = reference.replace('myfile_velocity.open(filename_velocity, ios::app);',
                                  'myfile_velocity.open("velocity.dat", ios::app);')
    (OUTPUT/'reference.cpp').write_text(reference)
    assignments = '\n'.join(f'    p.{key} = {value!r};' for key, value in PARAMETERS.items())
    (OUTPUT/'library.cpp').write_text('''#include <gpe/solver.hpp>
#include <fstream>
int main() {
    gpe::SimulationParameters p;
''' + assignments + '''
    p.output_directory = "library";
    const auto result = gpe::GPE_solver(p);
    std::ofstream field("library/final_field.bin",std::ios::binary);
    field.write(reinterpret_cast<const char*>(result.final_field.data()),
                result.final_field.size()*sizeof(std::complex<double>));
}
''')
    manifest = dict(parameters=PARAMETERS, initial_field='phi(x)=1',
                    original_sha256=hashlib.sha256(original.encode()).hexdigest(),
                    library_sha256=hashlib.sha256((ROOT/'src/solver.cpp').read_bytes()).hexdigest(),
                    eigen=str(eigen), compiler=subprocess.check_output(['g++','--version'],text=True),
                    reference_adjustments='Parameters; integer step count; full precision velocity every 10 steps; compressed snapshots every 50 steps; fail on nonconvergence; legacy output and stopping replaced. Numerical update retained.')
    (OUTPUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    common = ['g++','-std=c++17','-O3','-DNDEBUG','-I'+str(eigen)]
    subprocess.run(common+['-I'+str(ROOT/'src'),str(OUTPUT/'reference.cpp'),'-lz','-o',str(OUTPUT/'reference_solver')],check=True)
    subprocess.run(common+['-I'+str(ROOT/'include'),str(OUTPUT/'library.cpp'),
                          str(ROOT/'src/solver.cpp'),'-lz','-o',str(OUTPUT/'library_solver')],check=True)
    (OUTPUT/'reference').mkdir()


def run():
    jobs = []
    for name, cwd in [('reference',OUTPUT/'reference'),('library',OUTPUT)]:
        log = (OUTPUT/f'{name}.log').open('w')
        process = subprocess.Popen([str(OUTPUT/f'{name}_solver')],cwd=cwd,stdout=log,stderr=subprocess.STDOUT)
        jobs.append((name,process,log))
    failed = []
    for name, process, log in jobs:
        code = process.wait()
        log.close()
        if code:
            failed.append((name,code))
    if failed:
        raise RuntimeError(f'Solver failures: {failed}; partial data retained')


def compare(output=OUTPUT, parameters=PARAMETERS,
            initial_description='Initial field: phi(x)=1; periodic boundaries.',
            reference_description='Reference: original fresh-run numerical update and routines.',
            library_label='Library GPE_solver'):
    t_initial = parameters.get('t_initial', 0.0)
    total_steps = int(round((parameters['t_final']-t_initial)/parameters['dt']))
    expected_rows = total_steps//parameters['velocity_save_every']+1
    expected_snapshots = total_steps//parameters['config_save_every']+1
    os.environ.setdefault('MPLCONFIGDIR','/tmp/gpe-comparison-mpl')
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    reference = np.loadtxt(output/'reference/velocity.dat')
    library = np.loadtxt(output/'library/velocity.dat')
    assert reference.shape == library.shape == (expected_rows,3)
    assert np.isfinite(reference).all() and np.isfinite(library).all()
    expected_times = t_initial + np.arange(expected_rows)*parameters['velocity_save_every']*parameters['dt']
    assert np.array_equal(reference[:,0],expected_times)
    assert np.array_equal(library[:,0],expected_times)
    difference = library[:,1:]-reference[:,1:]
    exact = np.all(difference==0,axis=1)
    final_reference = np.fromfile(output/'reference/final_field.bin',dtype=np.complex128)
    final_library = np.fromfile(output/'library/final_field.bin',dtype=np.complex128)
    assert final_reference.shape == final_library.shape == (int(parameters['length']/parameters['dx']),)
    stats = dict(time_points=expected_rows, identical_velocity_points=int(exact.sum()),
                 max_real_velocity_difference=float(np.max(abs(difference[:,0]))),
                 max_imaginary_velocity_difference=float(np.max(abs(difference[:,1]))),
                 final_reference_velocity=float(reference[-1,1]), final_library_velocity=float(library[-1,1]),
                 final_field_exact=bool(np.array_equal(final_reference,final_library)),
                 max_final_field_difference=float(np.max(abs(final_reference-final_library))))
    # Read every member to EOF to check complete tar contents and gzip integrity.
    expected_names = [f'phi_{step}.dat' for step in range(0,total_steps+1,parameters['config_save_every'])]
    for name in ('reference','library'):
        members = []
        with tarfile.open(output/name/'configurations.tar.gz','r|gz') as archive:
            for member in archive:
                members.append(member.name)
                stream = archive.extractfile(member)
                while stream.read(1024*1024):
                    pass
        assert members == expected_names and len(members) == expected_snapshots
        subprocess.run(['gzip','-t',str(output/name/'configurations.tar.gz')],check=True)
        stats[name+'_compressed_configurations'] = len(members)
    np.savetxt(output/'comparison.csv',np.column_stack((expected_times,reference[:,1:],library[:,1:],difference,exact)),
               delimiter=',',header='time,reference_real,reference_imag,library_real,library_imag,delta_real,delta_imag,exact',comments='',fmt='%.17g')
    (output/'summary.json').write_text(json.dumps(stats,indent=2)+'\n')
    fig, axes = plt.subplots(2,1,figsize=FIGURE_SIZE,sharex=True)
    axes[0].plot(expected_times,reference[:,1],label='Original research update')
    axes[0].plot(expected_times,library[:,1],'--',label=library_label)
    axes[0].set(ylabel='Impurity velocity',title='Identical parameters: V₀=0.5, L=500, Δx=0.001, Δt=0.01')
    axes[0].legend()
    axes[1].plot(expected_times,difference[:,0],label='Real velocity difference')
    axes[1].plot(expected_times,difference[:,1],'--',label='Imaginary velocity difference')
    axes[1].set(xlabel='Time',ylabel='Library − original',title=f'{exact.sum()} / {expected_rows} time points exactly equal')
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(output/'velocity_comparison.png',dpi=DPI)
    with PdfPages(output/'comparison_report.pdf') as pdf:
        pdf.savefig(fig)
        page = plt.figure(figsize=FIGURE_SIZE)
        details = ('Velocity regression comparison\n\n'
                   + '\n'.join(f'{key}: {value}' for key,value in parameters.items())
                   + '\n'+initial_description+'\n\n'
                   + '\n'.join(f'{key}: {value}' for key,value in stats.items())
                   + '\n\n'+reference_description+'\n'
                   'Only parameters, stopping and output handling were adapted.\n'
                   'Both compiled using the same g++ -O3 -DNDEBUG and Eigen.\n'
                   f'All {expected_rows} saved times checked at round-trip double precision.\n'
                   'This checks equivalence for this run, not physical convergence\n'
                   'or conservation. No literature or novelty assessment is involved.\n'
                   'Full pointwise evidence: comparison.csv; provenance: manifest.json.')
        page.text(0.08,0.95,details,va='top',fontsize=10)
        pdf.savefig(page)
        plt.close(page)
    plt.close(fig)
    print(json.dumps(stats,indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eigen',type=Path)
    parser.add_argument('--compare-only',action='store_true')
    args = parser.parse_args()
    if not args.compare_only:
        if args.eigen is None or not (args.eigen/'Eigen/Core').is_file():
            parser.error('--eigen must contain Eigen/Core')
        prepare(args.eigen)
        start = time.monotonic()
        run()
        print(f'Both runs completed in {time.monotonic()-start:.1f} seconds',flush=True)
    compare()
