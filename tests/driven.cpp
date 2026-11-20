#include <gpe/GPE_solver_driven_contd.hpp>
#include <cmath>
#include <stdexcept>

namespace {
constexpr double FORCES[] = {0.5, -0.5};
constexpr double COUPLINGS[] = {0.5, -0.5};
constexpr double CHECKPOINT_TIME = 0.01;
constexpr double FINAL_TIME = 0.02;
constexpr double TOLERANCE = 2e-13;
constexpr double ZERO_COUPLING = 0.0;
constexpr int SAVE_INTERVAL = 10;
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    const std::filesystem::path root = argv[1];
    for (double force : FORCES) {
        const auto force_root = root / (force > 0 ? "positive" : "negative");
        for (double coupling : COUPLINGS) {
            const auto run_root = force_root / (coupling > 0 ? "repulsive" : "attractive");
            std::filesystem::create_directories(run_root);
            gpe::DrivenSimulationParameters p;
            p.F = force;
            p.G = coupling;
            p.t_final = FINAL_TIME;
            p.config_save_every = SAVE_INTERVAL;
            p.velocity_save_every = SAVE_INTERVAL;
            p.output_directory = run_root / "full";
            const auto full = gpe::GPE_solver_driven(p);
            p.t_final = CHECKPOINT_TIME;
            p.output_directory = run_root / "half";
            const auto half = gpe::GPE_solver_driven(p);
            p.t_initial = half.final_time;
            p.t_final = FINAL_TIME;
            p.initial_velocity = half.final_velocity;
            p.output_directory = run_root / "resumed";
            const auto resumed = gpe::GPE_solver_driven_contd(
                p, half.configurations_archive, "phi_10.dat");
            require((resumed.final_field - full.final_field).norm() < TOLERANCE,
                    "Driven compressed restart changes final field");
            require(std::abs(resumed.final_velocity - full.final_velocity) < TOLERANCE,
                    "Driven restart resets force clock or changes final velocity");
            p.initial_field = half.final_field;
            p.output_directory = run_root / "memory";
            const auto memory = gpe::GPE_solver_driven(p);
            require((memory.final_field - full.final_field).norm() < TOLERANCE
                    && std::abs(memory.final_velocity - full.final_velocity) < TOLERANCE,
                    "Driven in-memory restart changes evolution");
        }
        gpe::DrivenSimulationParameters free;
        free.F = force;
        free.G = ZERO_COUPLING;
        free.output_directory = force_root / "free";
        const auto result = gpe::GPE_solver_driven(free);
        require(std::abs(result.final_velocity - free.mass_ratio * force * free.t_final)
                    < TOLERANCE, "Free driven impurity acceleration is incorrect");
    }
}
