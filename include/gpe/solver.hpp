#pragma once
#include <Eigen/Core>
#include <complex>
#include <filesystem>
#include <optional>

namespace gpe {
struct SimulationParameters {
    double t_final = 0.02;
    double length = 16.0;
    double V_0 = 2.0;
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
struct SimulationResult {
    std::filesystem::path configurations_archive;
    std::filesystem::path velocity_file;
    Eigen::VectorXcd final_field;
    std::complex<double> final_velocity;
    double final_time;
};
// Writes and flushes each velocity record directly to disk; retains no velocity history.
// Throws std::exception on invalid input, nonconvergence, or output failure.
SimulationResult GPE_solver(const SimulationParameters& parameters);
}
