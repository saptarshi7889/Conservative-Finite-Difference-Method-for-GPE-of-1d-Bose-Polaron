"""Use the named research solver's routines with the verified CLI update driver."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESEARCH = ROOT / 'legacy' / 'finite_difference_modified_impurity_Thomas_generalized_without_Matrix_hermitian_vel_contd_new_cluster_upwind_scheme.cpp'
OUTPUT = ROOT / 'build/L500_continuation.cpp'


def prepare():
    research = RESEARCH.read_text()
    prefix = research[:research.index('//start of the main code')]
    driver = (ROOT / 'src/main.cpp').read_text().replace('#include "numerics.hpp"', '')
    driver = driver.replace('    if (!filesystem::create_directory(output)) throw runtime_error("Output directory already exists.");', '''    const bool append = args.count("--restart") && filesystem::is_directory(output);
    if (!append && !filesystem::create_directory(output)) throw runtime_error("Output directory already exists.");
    const int step_offset = static_cast<int>(round(t_init/delta_t));
    if (append) {
        ifstream previous(output/"velocity.dat");
        string line;
        double last_time=-1, last_velocity=0, last_imaginary=0;
        while (getline(previous,line)) {
            if (line.empty() || line[0]=='#') continue;
            istringstream row(line);
            if (!(row>>last_time>>last_velocity>>last_imaginary)) throw runtime_error("Invalid existing velocity row.");
        }
        if (abs(last_time-t_init)>GRID_TOLERANCE || abs(last_velocity-real(V))>GRID_TOLERANCE || abs(last_imaginary)>GRID_TOLERANCE)
            throw runtime_error("Velocity history does not match restart checkpoint.");
        for (int step=1;step<=total_steps;++step) {
            if (step%save_every!=0 && step!=total_steps) continue;
            auto path=output/("phi_"+to_string(step_offset+step)+".dat");
            if (filesystem::exists(path) || filesystem::exists(path.string()+".xz") || filesystem::exists(path.string()+".tmp"))
                throw runtime_error("Continuation snapshot already exists.");
        }
    }''')
    original = '    ofstream velocity(output/"velocity.dat"), norm(output/"norm.dat"), convergence(output/"convergence.dat"), parameters(output/"parameters.txt");'
    assert original in driver
    driver = driver.replace(original, '''    const auto mode=append ? ios::app : ios::out;
    ofstream velocity(output/"velocity.dat",mode), norm(output/"norm.dat",mode), convergence(output/"convergence.dat",mode);
    const auto parameter_path=output/(append ? "parameters_continuation_"+to_string(step_offset)+".txt" : "parameters.txt");
    if (append && filesystem::exists(parameter_path)) throw runtime_error("Continuation metadata already exists.");
    ofstream parameters(parameter_path);''')
    driver = driver.replace('    velocity<<"# time real_velocity imaginary_velocity\\n";', '    if (!append) velocity<<"# time real_velocity imaginary_velocity\\n";')
    driver = driver.replace('    norm<<"# time euclidean_norm dx_sum_abs_phi_squared\\n";', '    if (!append) norm<<"# time euclidean_norm dx_sum_abs_phi_squared\\n";')
    driver = driver.replace('    convergence<<"# time iterations iterate_difference_norm\\n";', '    if (!append) convergence<<"# time iterations iterate_difference_norm\\n";')
    driver = driver.replace('    record_velocity(t_init);\n    record(t_init,0);', '    if (!append) { record_velocity(t_init); record(t_init,step_offset); }')
    driver = driver.replace('record(t_init+count*delta_t,count)', 'record(t_init+count*delta_t,step_offset+count)')
    # Use the requested continuation source's update block verbatim, removing
    # only its legacy convergence-file writes (the driver records every step).
    update = research[research.index('     mult3 ='):research.index('         Phi_n1 = Phi_hat;')]
    log_start = update.index('               if (count % 50 == 0)')
    log_end = update.index('               break;', log_start)
    update = update[:log_start] + update[log_end:]
    begin = driver.index('     mult3 =')
    end = driver.index('         if (k > K_fin', begin)
    driver = driver[:begin] + update + driver[end:]
    OUTPUT.write_text(prefix + '\n// Restart and output driver adapted from src/main.cpp.\n' + driver)
    return OUTPUT


if __name__ == '__main__':
    print(prepare())
