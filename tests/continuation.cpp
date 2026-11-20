#include <gpe/GPE_solver_contd.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

namespace {
constexpr double TEST_FORCE = 0.5;
constexpr double CHECKPOINT_TIME = 0.01;
constexpr int SAVE_INTERVAL = 10;
constexpr int OUTPUT_PRECISION = std::numeric_limits<double>::max_digits10;
}

int main(int argc, char** argv) {
    if (argc < 3) return 1;
    try {
        const std::filesystem::path root = argv[1];
        gpe::SimulationParameters p;
        p.force = TEST_FORCE;
        p.V_0 = std::stod(argv[2]);
        p.config_save_every = SAVE_INTERVAL;
        p.velocity_save_every = SAVE_INTERVAL;
        if (argc == 3) {
            p.output_directory = root / "full";
            gpe::GPE_solver(p);
            p.output_directory = root / "half";
            p.t_final = CHECKPOINT_TIME;
            const auto half = gpe::GPE_solver(p);
            std::ofstream velocity(root / "checkpoint_velocity.txt");
            velocity << std::setprecision(OUTPUT_PRECISION)
                     << half.final_velocity.real() << ' ' << half.final_velocity.imag();
        } else {
            if (argc != 6 && argc != 7) return 1;
            p.t_initial = CHECKPOINT_TIME;
            double re, im;
            std::ifstream velocity(root / "checkpoint_velocity.txt");
            if (!(velocity >> re >> im)) return 1;
            p.initial_velocity = std::complex<double>(re, im);
            p.output_directory = root / argv[4];
            const std::string member = std::string(argv[5]) == "-" ? "" : argv[5];
            if (argc == 7) {
                const std::string invalid = argv[6];
                if (invalid == "zero-time") p.t_initial = 0;
                if (invalid == "missing-velocity") p.initial_velocity.reset();
                if (invalid == "nonfinite-velocity") p.initial_velocity = std::numeric_limits<double>::infinity();
                if (invalid == "invalid-grid") p.dx = 0;
                if (invalid == "initial-field") p.initial_field = Eigen::VectorXcd::Ones(3);
            }
            gpe::GPE_solver_contd(p, argv[3], member);
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
