#pragma once
#include <gpe/solver.hpp>
#include <string>

namespace gpe {
// Continue from a compressed local checkpoint (x, Re(phi), Im(phi)).
// Supply t_initial > 0 and initial_velocity at that checkpoint, retaining V_0.
// initial_field must be empty: the field is loaded from checkpoint_file.
// For .dat.gz/.dat.xz, leave checkpoint_member empty. For a compressed tar
// archive, specify the exact member name, e.g. "phi_1000.dat".
// Decompression is streamed into the field; no expanded file is created.
// Uses the same evolution and output conventions as GPE_solver.
SimulationResult GPE_solver_contd(
    const SimulationParameters& parameters,
    const std::filesystem::path& checkpoint_file,
    const std::string& checkpoint_member = {});
}
