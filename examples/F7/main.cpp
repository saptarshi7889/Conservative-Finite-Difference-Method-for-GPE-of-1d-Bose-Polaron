#include <gpe/solver.hpp>
#include <iostream>

namespace {
constexpr double FORCE = 7.0;
constexpr double GAMMA = 1.0/400.0;
constexpr double DT = 0.01;
constexpr double LENGTH = 1500.0;
constexpr double DX = 1.0/2000.0;
constexpr double FINAL_TIME = 10.0;
constexpr double COUPLING = 15.0;
constexpr double MASS_RATIO = 1.0/3.0;
constexpr double TOLERANCE = 1e-8;
constexpr int ITERATIONS = 12;
constexpr int VELOCITY_INTERVAL = 10;
constexpr int CONFIG_INTERVAL = 50;
}
int main(int argc, char** argv) {
    try {
        gpe::SimulationParameters p;
        p.force=FORCE; p.gamma=GAMMA; p.dt=DT; p.length=LENGTH; p.dx=DX;
        p.t_final=FINAL_TIME; p.G=COUPLING; p.mass_ratio=MASS_RATIO; p.V_0=0;
        p.accept_unconverged_steps=true;
        p.picard_tolerance=TOLERANCE; p.max_picard_iterations=ITERATIONS;
        p.velocity_save_every=VELOCITY_INTERVAL; p.config_save_every=CONFIG_INTERVAL;
        p.output_directory=argc>1 ? argv[1] : "results_F7_t10_picard12";
        const auto result=gpe::GPE_solver(p);
        std::cout << "Completed t=" << result.final_time << "\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
