#pragma once
#include <gpe/solver.hpp>

namespace gpe {
struct DrivenSimulationParameters {
    double t_final = 0.02;
    double length = 16.0;
    double F = 0.5; // Dimensionless constant driving force; initial drive velocity is zero.
    double mass_ratio = 0.125; // m/M
    double G = -1.5;
    double gamma = 0.1;
    double dx = 0.125;
    double dt = 0.001;
    double picard_tolerance = 1e-8;
    int max_picard_iterations = 12;
    int config_save_every = 10;
    int velocity_save_every = 10;
    std::filesystem::path output_directory = "results";
    Eigen::VectorXcd initial_field; // Empty means phi(x) = 1.
    double t_initial = 0.0;
    std::optional<std::complex<double>> initial_velocity;
    bool write_plain_snapshots = false; // CLI checkpoint compatibility.
};
// Constant-force evolution with drive velocity (m/M)*F*t at absolute time t.
// The initial field defaults to uniform and the fresh-run velocity to zero.
// Throws on Picard nonconvergence before accepting or recording the failed step.
// Output, convergence and in-memory restart conventions match GPE_solver.
SimulationResult GPE_solver_driven(const DrivenSimulationParameters& parameters);
}
