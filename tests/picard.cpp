#include <gpe/GPE_solver_contd.hpp>
#include <gpe/GPE_solver_driven_contd.hpp>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {
constexpr double CHECKPOINT_TIME = 0.01;
constexpr double FINAL_TIME = 0.02;
constexpr double FAILURE_TOLERANCE = 1e-30;
constexpr double INITIAL_VELOCITY = 0.0;
constexpr int FAILURE_ITERATIONS = 1;
constexpr double ONE_STEP = 0.001;

template<class Run>
void expect_failure(Run run) {
    try { run(); }
    catch (const std::runtime_error& error) {
        const std::string message = error.what();
        if (message.find("Nonlinear iteration failed at time ") == 0
            && message.find(" after 1 iterations; residual ") != std::string::npos
            && message.find(", tolerance ") != std::string::npos) return;
        throw;
    }
    throw std::runtime_error("Unconverged Picard step was accepted");
}

template<class Parameters, class Run, class Continue>
void check(const std::filesystem::path& root, Parameters p, Run run, Continue resume) {
    std::filesystem::create_directory(root);
    p.output_directory = root / "checkpoint";
    p.t_final = CHECKPOINT_TIME;
    const auto checkpoint = run(p);
    p.max_picard_iterations = FAILURE_ITERATIONS;
    p.picard_tolerance = FAILURE_TOLERANCE;
    p.t_final = FINAL_TIME;
    p.output_directory = root / "fresh_failed";
    expect_failure([&] { run(p); });
    p.t_initial = CHECKPOINT_TIME;
    p.initial_velocity = checkpoint.final_velocity;
    p.output_directory = root / "restart_failed";
    expect_failure([&] { resume(p, checkpoint.configurations_archive, "phi_10.dat"); });

    // Measure the iterations needed for a single converged step, then require
    // success when the cap is exactly that count (not one iteration more).
    Parameters single;
    single.t_final = ONE_STEP;
    single.output_directory = root / "single";
    const auto baseline = run(single);
    std::ifstream convergence(single.output_directory / "convergence.dat");
    std::string header;
    std::getline(convergence, header);
    double time, error;
    int iterations;
    if (!(convergence >> time >> iterations >> error))
        throw std::runtime_error("Missing converged-step diagnostic");
    single.max_picard_iterations = iterations;
    single.output_directory = root / "last_iteration";
    const auto capped = run(single);
    if (capped.final_velocity != baseline.final_velocity
        || (capped.final_field - baseline.final_field).norm() != 0)
        throw std::runtime_error("Convergence at iteration cap changed the solution");
}
}

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    const std::filesystem::path root = argv[1];
    gpe::SimulationParameters quenched;
    quenched.V_0 = INITIAL_VELOCITY;
    check(root / "quenched", quenched, gpe::GPE_solver, gpe::GPE_solver_contd);
    check(root / "driven", gpe::DrivenSimulationParameters{},
          gpe::GPE_solver_driven, gpe::GPE_solver_driven_contd);
}
