#include <utility>
#include <gpe/GPE_solver_contd.hpp>
#include <gpe/GPE_solver_driven_contd.hpp>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;

namespace {
template<class Parameters>
py::class_<Parameters> parameter_class(py::module_& module, const char* name) {
    return py::class_<Parameters>(module, name)
        .def(py::init<>())
        .def_readwrite("t_final", &Parameters::t_final)
        .def_readwrite("length", &Parameters::length)
        .def_readwrite("mass_ratio", &Parameters::mass_ratio)
        .def_readwrite("G", &Parameters::G)
        .def_readwrite("gamma", &Parameters::gamma)
        .def_readwrite("dx", &Parameters::dx)
        .def_readwrite("dt", &Parameters::dt)
        .def_readwrite("picard_tolerance", &Parameters::picard_tolerance)
        .def_readwrite("max_picard_iterations", &Parameters::max_picard_iterations)
        .def_readwrite("config_save_every", &Parameters::config_save_every)
        .def_readwrite("velocity_save_every", &Parameters::velocity_save_every)
        .def_readwrite("output_directory", &Parameters::output_directory)
        .def_property("initial_field",
            [](const Parameters& parameters) { return Eigen::VectorXcd(parameters.initial_field); },
            [](Parameters& parameters, const py::object& field) {
                if (field.is_none()) parameters.initial_field.resize(0);
                else parameters.initial_field = field.cast<Eigen::VectorXcd>();
            }, "Initial field as a NumPy complex vector. Assignment copies; None clears it.")
        .def_readwrite("t_initial", &Parameters::t_initial)
        .def_readwrite("initial_velocity", &Parameters::initial_velocity)
        .def_readwrite("write_plain_snapshots", &Parameters::write_plain_snapshots);
}

// Copy settings before releasing the GIL: another Python thread may then edit
// the caller's object without racing with the running C++ solver.
template<class Parameters, class Run>
gpe::SimulationResult run(Parameters parameters, Run solver) {
    py::gil_scoped_release release;
    return solver(parameters);
}

template<class Parameters, class Continue>
gpe::SimulationResult resume(Parameters parameters, const std::filesystem::path& checkpoint,
                             const std::string& member, Continue solver) {
    py::gil_scoped_release release;
    return solver(parameters, checkpoint, member);
}
}

PYBIND11_MODULE(_core, module) {
    module.doc() = "C++17 Bose-polaron GPE solver bindings.";
    module.attr("__version__") = GPE_VERSION;
    parameter_class<gpe::SimulationParameters>(module, "SimulationParameters")
        .def_readwrite("V_0", &gpe::SimulationParameters::V_0)
        .def_readwrite("force", &gpe::SimulationParameters::force);
    parameter_class<gpe::DrivenSimulationParameters>(module, "DrivenSimulationParameters")
        .def_readwrite("F", &gpe::DrivenSimulationParameters::F);
    py::class_<gpe::SimulationResult>(module, "SimulationResult")
        .def_readonly("configurations_archive", &gpe::SimulationResult::configurations_archive)
        .def_readonly("velocity_file", &gpe::SimulationResult::velocity_file)
        .def_property_readonly("final_field", [](const gpe::SimulationResult& result)
                -> const Eigen::VectorXcd& { return result.final_field; },
            py::return_value_policy::reference_internal,
            "Read-only NumPy view; keeps the result alive. Use .copy() for a writable array.")
        .def_readonly("final_velocity", &gpe::SimulationResult::final_velocity)
        .def_readonly("final_time", &gpe::SimulationResult::final_time);
    module.def("GPE_solver", [](gpe::SimulationParameters parameters) {
        return run(std::move(parameters), gpe::GPE_solver);
    }, py::arg("parameters"), "Run a fresh simulation or continue from a field in memory.");
    module.def("GPE_solver_driven", [](gpe::DrivenSimulationParameters parameters) {
        return run(std::move(parameters), gpe::GPE_solver_driven);
    }, py::arg("parameters"), "Run constant-force evolution with zero initial drive velocity.");
    module.def("GPE_solver_contd", [](gpe::SimulationParameters parameters,
            const std::filesystem::path& checkpoint, const std::string& member) {
        return resume(std::move(parameters), checkpoint, member, gpe::GPE_solver_contd);
    }, py::arg("parameters"), py::arg("checkpoint_file"), py::arg("checkpoint_member") = "",
       "Continue from a compressed checkpoint at its absolute time and recorded velocity.");
    module.def("GPE_solver_driven_contd", [](gpe::DrivenSimulationParameters parameters,
            const std::filesystem::path& checkpoint, const std::string& member) {
        return resume(std::move(parameters), checkpoint, member, gpe::GPE_solver_driven_contd);
    }, py::arg("parameters"), py::arg("checkpoint_file"), py::arg("checkpoint_member") = "",
       "Continue driven evolution from a compressed checkpoint without resetting the force clock.");
}
