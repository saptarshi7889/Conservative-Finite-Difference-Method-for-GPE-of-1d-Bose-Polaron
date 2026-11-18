#include <gpe/GPE_solver_contd.hpp>
#include <iostream>

namespace {
constexpr double CONTINUATION_END_TIME = 0.04;
}

int main() {
    gpe::SimulationParameters parameters;
    // Small demonstration, not a converged research simulation.
    parameters.output_directory = "library-results";
    try {
        const auto result = gpe::GPE_solver(parameters);
        std::cout << "Configurations: " << result.configurations_archive << '\n'
                  << "Velocity: " << result.velocity_file << '\n';
        parameters.t_initial = result.final_time;
        parameters.t_final = CONTINUATION_END_TIME;
        parameters.initial_velocity = result.final_velocity;
        parameters.output_directory = "continued-library-results";
        const auto continued = gpe::GPE_solver_contd(
            parameters, result.configurations_archive, "phi_20.dat");
        std::cout << "Continued velocity: " << continued.velocity_file << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
