#include <gpe/GPE_solver_driven_contd.hpp>
#include <gpe/GPE_solver_contd.hpp>

namespace {
gpe::SimulationParameters settings_for(const gpe::DrivenSimulationParameters& parameters) {
    gpe::SimulationParameters settings;
    settings.V_0 = 0.0;
    settings.force = parameters.F;
    settings.t_final = parameters.t_final;
    settings.length = parameters.length;
    settings.mass_ratio = parameters.mass_ratio;
    settings.G = parameters.G;
    settings.gamma = parameters.gamma;
    settings.dx = parameters.dx;
    settings.dt = parameters.dt;
    settings.picard_tolerance = parameters.picard_tolerance;
    settings.max_picard_iterations = parameters.max_picard_iterations;
    settings.config_save_every = parameters.config_save_every;
    settings.velocity_save_every = parameters.velocity_save_every;
    settings.output_directory = parameters.output_directory;
    settings.initial_field = parameters.initial_field;
    settings.t_initial = parameters.t_initial;
    settings.initial_velocity = parameters.initial_velocity;
    settings.accept_unconverged_steps = parameters.accept_unconverged_steps;
    settings.write_plain_snapshots = parameters.write_plain_snapshots;
    return settings;
}
}

gpe::SimulationResult gpe::GPE_solver_driven(const DrivenSimulationParameters& parameters) {
    return GPE_solver(settings_for(parameters));
}

gpe::SimulationResult gpe::GPE_solver_driven_contd(
    const DrivenSimulationParameters& parameters,
    const std::filesystem::path& checkpoint_file,
    const std::string& checkpoint_member) {
    return GPE_solver_contd(settings_for(parameters), checkpoint_file, checkpoint_member);
}
