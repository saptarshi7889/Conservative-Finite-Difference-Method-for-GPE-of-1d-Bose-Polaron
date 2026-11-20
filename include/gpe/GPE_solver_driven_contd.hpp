#pragma once
#include <gpe/GPE_solver_driven.hpp>
#include <string>

namespace gpe {
// Compressed continuation: retain F and physical/grid settings; supply the
// absolute checkpoint time and its recorded initial_velocity. Leave
// initial_field empty. File formats/member rules match GPE_solver_contd.
SimulationResult GPE_solver_driven_contd(
    const DrivenSimulationParameters& parameters,
    const std::filesystem::path& checkpoint_file,
    const std::string& checkpoint_member = {});
}
