#include <gpe/GPE_solver_driven_contd.hpp>
#include <iostream>

namespace {
constexpr double FORCE = 0.5;
constexpr double COUPLING = 0.5;
constexpr double GAMMA = 1.0 / 400.0;
constexpr double MASS_RATIO = 1.0 / 3.0;
constexpr double CONTINUATION_END_TIME = 0.04;
}

int main() {
    gpe::DrivenSimulationParameters parameters;
    // Small demonstration, not a converged research simulation.
    parameters.F = FORCE;
    parameters.G = COUPLING;
    parameters.gamma = GAMMA;
    parameters.mass_ratio = MASS_RATIO;
    parameters.output_directory = "driven-results";
    try {
        const auto result = gpe::GPE_solver_driven(parameters);
        std::cout << "Configurations: " << result.configurations_archive << '\n'
                  << "Velocity: " << result.velocity_file << '\n';
        parameters.t_initial = result.final_time;
        parameters.t_final = CONTINUATION_END_TIME;
        parameters.initial_velocity = result.final_velocity;
        parameters.output_directory = "continued-driven-results";
        const auto continued = gpe::GPE_solver_driven_contd(
            parameters, result.configurations_archive, "phi_20.dat");
        std::cout << "Continued velocity: " << continued.velocity_file << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
