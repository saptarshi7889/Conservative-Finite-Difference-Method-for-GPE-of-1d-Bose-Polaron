"""Small quenched/driven runs followed by compressed continuation."""
from pathlib import Path
import bose_polaron_gpe as gpe

OUTPUT = Path("python-example-results")
CONTINUATION_END_TIME = 0.04
CHECKPOINT_MEMBER = "phi_20.dat"
FORCE = 0.5


def main():
    OUTPUT.mkdir()
    for name, parameter_type, solve, resume in (
        ("quenched", gpe.SimulationParameters, gpe.GPE_solver, gpe.GPE_solver_contd),
        ("driven", gpe.DrivenSimulationParameters, gpe.GPE_solver_driven,
         gpe.GPE_solver_driven_contd),
    ):
        parameters = parameter_type()
        if name == "driven":
            parameters.F = FORCE
        parameters.output_directory = OUTPUT / name
        result = solve(parameters)
        parameters.t_initial = result.final_time
        parameters.t_final = CONTINUATION_END_TIME
        parameters.initial_velocity = result.final_velocity
        parameters.output_directory = OUTPUT / f"{name}-continued"
        continued = resume(parameters, result.configurations_archive, CHECKPOINT_MEMBER)
        print(name, continued.final_time, continued.final_velocity,
              continued.final_field.shape, continued.velocity_file)


if __name__ == "__main__":
    main()
