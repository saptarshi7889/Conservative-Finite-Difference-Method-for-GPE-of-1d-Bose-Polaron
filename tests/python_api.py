"""Installed Python package regression tests; NumPy and standard library only."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gc
import gzip
import importlib.metadata
import io
import lzma
from pathlib import Path
import subprocess
import tarfile
import tempfile
import time
import unittest

import numpy as np
import bose_polaron_gpe as gpe

CHECKPOINT_TIME = 0.01
FINAL_TIME = 0.02
ONE_STEP = 0.001
FIELD_TOLERANCE = 2e-13
FAILURE_TOLERANCE = 1e-30
FAILURE_ITERATIONS = 1
CHECKPOINT_MEMBER = "phi_10.dat"
FORCE = 0.5
INITIAL_VELOCITY = 2.0
THREAD_LENGTH = 512.0
THREAD_END_TIME = 1.0
THREAD_TIMEOUT = 10.0
POLL_INTERVAL = 0.001
ZERO_COUPLING = 0.0
CASES = (
    ("quenched", gpe.SimulationParameters, gpe.GPE_solver, gpe.GPE_solver_contd),
    ("driven", gpe.DrivenSimulationParameters, gpe.GPE_solver_driven,
     gpe.GPE_solver_driven_contd),
)
CLI = None


def checkpoint_bytes(path):
    with tarfile.open(path) as archive:
        return archive.extractfile(CHECKPOINT_MEMBER).read()


def numeric_rows(path):
    return np.loadtxt(path, ndmin=2)


class PythonAPI(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def assert_same(self, actual, expected):
        np.testing.assert_allclose(actual.final_field, expected.final_field,
                                   rtol=0, atol=FIELD_TOLERANCE)
        self.assertLessEqual(abs(actual.final_velocity - expected.final_velocity),
                             FIELD_TOLERANCE)
        self.assertEqual(actual.final_time, expected.final_time)

    def test_version_and_parameter_types(self):
        self.assertEqual(gpe.__version__, importlib.metadata.version("bose-polaron-gpe"))
        for _, parameter_type, _, _ in CASES:
            p = parameter_type()
            self.assertIsNone(p.initial_velocity)
            p.initial_velocity = complex(INITIAL_VELOCITY, FORCE)
            self.assertEqual(p.initial_velocity, complex(INITIAL_VELOCITY, FORCE))
            p.initial_velocity = None
            p.output_directory = self.root
            self.assertEqual(p.output_directory, self.root)
            p.output_directory = str(self.root)
            self.assertEqual(p.output_directory, self.root)
            p.initial_field = None
            self.assertEqual(p.initial_field.size, 0)
        self.assertFalse(hasattr(gpe.DrivenSimulationParameters(), "V_0"))
        self.assertFalse(hasattr(gpe.DrivenSimulationParameters(), "force"))

    def test_array_ownership_and_input_copy(self):
        for name, parameter_type, solve, _ in CASES:
            p = parameter_type()
            points = round(p.length / p.dx)
            source = np.ones(points, dtype=np.complex128)
            p.initial_field = source
            source[:] = 0
            np.testing.assert_array_equal(p.initial_field, np.ones(points))
            initial_copy = p.initial_field
            initial_copy[:] = 0
            np.testing.assert_array_equal(p.initial_field, np.ones(points))
            p.output_directory = self.root / name
            result = solve(p)
            field = result.final_field
            self.assertEqual(field.dtype, np.dtype("complex128"))
            self.assertEqual(field.shape, (points,))
            self.assertFalse(field.flags.writeable)
            with self.assertRaises(ValueError):
                field[0] = 0
            expected = field.copy()
            self.assertIsInstance(result.velocity_file, Path)
            self.assertIsInstance(result.configurations_archive, Path)
            del result
            gc.collect()
            np.testing.assert_array_equal(field, expected)
            independent = field.copy()
            independent[:] = 0
            np.testing.assert_array_equal(field, expected)
            # Also accept noncontiguous complex input by copying it.
            p.initial_field = np.ones(points * 2, dtype=np.complex128)[::2]
            self.assertEqual(p.initial_field.size, points)
            p.initial_field = None
            self.assertEqual(p.initial_field.size, 0)

    def test_restarts_and_compression_formats(self):
        for name, parameter_type, solve, resume in CASES:
            for sign in (-1, 1):
                run_root = self.root / f"{name}_{sign}"
                run_root.mkdir()
                p = parameter_type()
                if name == "driven":
                    p.F = sign * FORCE
                else:
                    p.V_0 = sign * INITIAL_VELOCITY
                p.output_directory = run_root / "full"
                full = solve(p)
                p.t_final = CHECKPOINT_TIME
                p.output_directory = run_root / "half"
                half = solve(p)
                raw = checkpoint_bytes(half.configurations_archive)
                (run_root / "field.dat.gz").write_bytes(gzip.compress(raw))
                (run_root / "field.dat.xz").write_bytes(lzma.compress(raw))
                with tarfile.open(run_root / "fields.tar.xz", "w:xz") as archive:
                    info = tarfile.TarInfo(CHECKPOINT_MEMBER)
                    info.size = len(raw)
                    archive.addfile(info, io.BytesIO(raw))
                p.t_initial = CHECKPOINT_TIME
                p.t_final = FINAL_TIME
                p.initial_velocity = half.final_velocity
                fixtures = [
                    (half.configurations_archive, CHECKPOINT_MEMBER),
                    (run_root / "fields.tar.xz", CHECKPOINT_MEMBER),
                    (run_root / "field.dat.gz", None),
                    (run_root / "field.dat.xz", None),
                ]
                for index, (source, member) in enumerate(fixtures):
                    before = source.read_bytes()
                    p.output_directory = run_root / f"restart_{index}"
                    if member is None:
                        result = resume(parameters=p, checkpoint_file=source)
                    else:
                        result = resume(parameters=p, checkpoint_file=source,
                                        checkpoint_member=member)
                    self.assert_same(result, full)
                    self.assertEqual(source.read_bytes(), before)
                p.initial_field = half.final_field
                p.output_directory = run_root / "memory"
                self.assert_same(solve(p), full)

    def test_nonconvergence_all_four_interfaces(self):
        for name, parameter_type, solve, resume in CASES:
            run_root = self.root / name
            run_root.mkdir()
            p = parameter_type()
            p.t_final = CHECKPOINT_TIME
            p.output_directory = run_root / "checkpoint"
            half = solve(p)
            p.t_final = FINAL_TIME
            p.max_picard_iterations = FAILURE_ITERATIONS
            p.picard_tolerance = FAILURE_TOLERANCE
            for continued in (False, True):
                p.output_directory = run_root / ("restart_failed" if continued else "fresh_failed")
                p.t_initial = CHECKPOINT_TIME if continued else 0.0
                p.initial_velocity = half.final_velocity if continued else None
                with self.assertRaisesRegex(RuntimeError, "Nonlinear iteration failed at time"):
                    if continued:
                        resume(p, half.configurations_archive, CHECKPOINT_MEMBER)
                    else:
                        solve(p)
                rows = numeric_rows(p.output_directory / "velocity.dat")
                self.assertEqual(rows.shape[0], 1)
                self.assertEqual(rows[0, 0], p.t_initial)
                with tarfile.open(p.output_directory / "configurations.tar.gz") as archive:
                    self.assertEqual(archive.getnames(), ["phi_0.dat"])

    def test_convergence_on_final_iteration(self):
        for name, parameter_type, solve, _ in CASES:
            p = parameter_type()
            p.t_final = ONE_STEP
            p.output_directory = self.root / f"{name}_baseline"
            baseline = solve(p)
            diagnostic = numeric_rows(p.output_directory / "convergence.dat")
            p.max_picard_iterations = int(diagnostic[0, 1])
            p.output_directory = self.root / f"{name}_capped"
            self.assert_same(solve(p), baseline)

    def test_validation_and_missing_checkpoint(self):
        for name, parameter_type, solve, resume in CASES:
            p = parameter_type()
            p.output_directory = self.root / name
            p.initial_field = np.ones(3, dtype=np.complex128)
            with self.assertRaisesRegex(RuntimeError, "Initial field"):
                solve(p)
            self.assertFalse(p.output_directory.exists())
            p.initial_field = None
            p.t_initial = CHECKPOINT_TIME
            p.initial_velocity = 0j
            with self.assertRaises(RuntimeError):
                resume(p, self.root / "missing.dat.gz")
            self.assertFalse(p.output_directory.exists())

    def test_gil_release_and_parameter_snapshot(self):
        p = gpe.SimulationParameters()
        p.G = ZERO_COUPLING
        p.V_0 = ZERO_COUPLING
        p.length = THREAD_LENGTH
        p.t_final = THREAD_END_TIME
        p.output_directory = self.root / "threaded"
        velocity_file = p.output_directory / "velocity.dat"
        with ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(gpe.GPE_solver, p)
            deadline = time.monotonic() + THREAD_TIMEOUT
            while not velocity_file.exists() and not future.done():
                self.assertLess(time.monotonic(), deadline)
                time.sleep(POLL_INTERVAL)
            self.assertFalse(future.done(), "Python did not run during C++ evolution")
            p.t_final = ONE_STEP
            p.initial_field = None
            p.output_directory = self.root / "edited_during_run"
            result = future.result(timeout=THREAD_TIMEOUT)
        self.assertEqual(result.final_time, THREAD_END_TIME)
        self.assertEqual(result.velocity_file, velocity_file)
        self.assertFalse(p.output_directory.exists())

    def test_cpp_cli_equivalence(self):
        if CLI is None:
            self.skipTest("Pass --cli to compare with a compiled C++ executable")
        for name, parameter_type, solve, _ in CASES:
            p = parameter_type()
            p.output_directory = self.root / f"{name}_python"
            result = solve(p)
            output = self.root / f"{name}_cpp"
            velocity = p.V_0 if name == "quenched" else 0.0
            force = p.force if name == "quenched" else p.F
            subprocess.run([CLI, "--output", str(output), "--velocity", str(velocity),
                            "--force", str(force)], check=True)
            np.testing.assert_allclose(numeric_rows(result.velocity_file),
                                       numeric_rows(output / "velocity.dat"),
                                       rtol=0, atol=FIELD_TOLERANCE)
            field = numeric_rows(output / "phi_20.dat")
            np.testing.assert_allclose(result.final_field, field[:, 1] + 1j * field[:, 2],
                                       rtol=0, atol=FIELD_TOLERANCE)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", type=lambda value: str(Path(value).resolve()))
    args, remaining = parser.parse_known_args()
    CLI = args.cli
    unittest.main(argv=[__file__, *remaining])
