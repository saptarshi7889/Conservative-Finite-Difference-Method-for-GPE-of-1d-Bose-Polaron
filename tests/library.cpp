#include <gpe/solver.hpp>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {
constexpr double RESTART_TOLERANCE = 2e-13;
constexpr double CHECKPOINT_TIME = 0.01;
constexpr int CONFIG_INTERVAL = 10;
constexpr int VELOCITY_INTERVAL = 5;
constexpr double LIVE_END_TIME = 1000.0;
constexpr int LIVE_CONFIG_INTERVAL = 1000000;
constexpr double VELOCITIES[] = {2.0, -2.0};
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}
int main(int argc, char** argv) {
    if (argc==3 && std::string(argv[2])=="live") {
        gpe::SimulationParameters p;
        p.output_directory=std::filesystem::path(argv[1])/"live";
        p.t_final=LIVE_END_TIME;
        p.config_save_every=LIVE_CONFIG_INTERVAL;
        p.velocity_save_every=VELOCITY_INTERVAL;
        p.V_0=0; p.G=0;
        gpe::GPE_solver(p);
        return 0;
    }
    if (argc!=2) return 1;
    const std::filesystem::path root=argv[1];
    for (double velocity : VELOCITIES) {
        const std::string name=velocity>0 ? "positive" : "negative";
        gpe::SimulationParameters p;
        p.V_0=velocity;
        p.config_save_every=CONFIG_INTERVAL;
        p.velocity_save_every=VELOCITY_INTERVAL;
        p.output_directory=root/name;
        const auto full=gpe::GPE_solver(p);
        require(full.final_time==p.t_final, "Incorrect final time");
        require(!std::filesystem::exists(p.output_directory/"phi_0.dat"),
                "Library should write compressed snapshots only");
        p.t_final=CHECKPOINT_TIME;
        p.output_directory=root/(name+"_half");
        const auto half=gpe::GPE_solver(p);
        p.t_initial=half.final_time;
        p.t_final=full.final_time;
        p.initial_field=half.final_field;
        p.initial_velocity=half.final_velocity;
        p.output_directory=root/(name+"_restart");
        const auto resumed=gpe::GPE_solver(p);
        require((resumed.final_field-full.final_field).norm()<RESTART_TOLERANCE,
                "Library restart changes final field");
        require(std::abs(resumed.final_velocity-full.final_velocity)<RESTART_TOLERANCE,
                "Library restart changes final velocity");
    }
    gpe::SimulationParameters bad;
    bad.output_directory=root/"bad";
    bad.initial_field=Eigen::VectorXcd::Ones(3);
    bool rejected=false;
    try { gpe::GPE_solver(bad); } catch (const std::exception&) { rejected=true; }
    require(rejected, "Grid mismatch must throw");
    require(!std::filesystem::exists(bad.output_directory), "Invalid input created output");
}
