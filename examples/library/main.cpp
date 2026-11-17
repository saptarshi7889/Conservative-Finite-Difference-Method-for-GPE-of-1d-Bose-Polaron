#include <gpe/solver.hpp>
#include <iostream>

int main() {
    gpe::SimulationParameters parameters;
    // Small demonstration, not a converged research simulation.
    parameters.output_directory = "library-results";
    try {
        const auto result = gpe::GPE_solver(parameters);
        std::cout << "Configurations: " << result.configurations_archive << '\n'
                  << "Velocity: " << result.velocity_file << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
