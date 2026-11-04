#!/usr/bin/env python3
"""Offline checks of the gate itself; these do not execute GPU kernels."""

import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import signal
import shutil
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
import run
import build_deps

spec = importlib.util.spec_from_file_location("snapshot_export", Path(__file__).with_name("export.py"))
snapshot_export = importlib.util.module_from_spec(spec)
spec.loader.exec_module(snapshot_export)


def receipt(rank, ranks, backend="CUDA", transport="direct"):
    return {"backend": backend, "rank": rank, "ranks": ranks, "local_rank": rank,
            "local_ranks": ranks, "device": rank, "device_id": "0000:{:02x}:00.0".format(rank),
            "hostname": "fixture", "kokkos_version": 40401, "mpi_device_buffer": transport == "direct",
            "mpi_library": "fixture MPI", "device_name": "fixture GPU", "runtime_version": 12000,
            "driver_version": 12000}


def output(ranks, backend="CUDA", transport="direct", tests=3):
    lines = []
    for rank in range(ranks):
        lines.append(run.RECEIPT + json.dumps(receipt(rank, ranks, backend, transport)))
        lines.append("[==========] Running {} tests from 1 test suite.".format(tests))
        lines.append(run.TEST_RECEIPT + json.dumps({"rank": rank, "ranks": ranks, "cases": tests,
                                                   "passed": tests, "skipped": 0}))
    return "\n".join(lines)


class KitTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def junit(self, text=None, expected=None):
        expected = run.SUITES if expected is None else expected
        root = ET.Element("testsuite", tests=str(len(expected)))
        for name, ranks in expected.items():
            case = ET.SubElement(root, "testcase", name=name, status="run")
            base = name.split("_mpi")[0]
            ET.SubElement(case, "system-out").text = output(ranks, tests=run.GTEST_CASES.get(base, 3)) if text is None else text
        path = self.root / "results.xml"
        ET.ElementTree(root).write(path)
        return path

    def arguments(self, backend="cuda", transport="direct"):
        return SimpleNamespace(source=self.root / "snapshot/source", work=self.root / "run",
            backend=backend, transport=transport, cxx=self.root / "c++", gpu_compiler=self.root / "gpucc",
            kokkos_dir=self.root / "kokkos/lib/cmake/Kokkos", gtest_source=self.root / "gtest",
            launcher=self.root / "srun", rank_flag="-n", launcher_arg=["--export=ALL", "--cpu-bind=cores"],
            jobs=4, cmake="cmake", ctest="ctest", arch="90" if backend == "cuda" else "gfx942", plan_only=False)

    def test_exact_matrix_and_gpu_receipts_pass(self):
        result = run.audit_junit(self.junit(), run.SUITES, "cuda", "direct", True)
        self.assertEqual(len(result), 11)
        self.assertEqual(sum(map(len, result.values())), 23)
        self.assertEqual(sum(run.GTEST_CASES[name.split("_mpi")[0]] for name in run.SUITES), 159)
        with self.assertRaises(ValueError):
            run.audit_junit(self.junit(output(1, tests=12), {"mesh_ghost_layer_test": 1}),
                            {"mesh_ghost_layer_test": 1}, "cuda", "direct", True)

    def test_cuda_and_hip_staged_are_distinct(self):
        for backend in ("cuda", "hip"):
            for transport in ("direct", "staged"):
                run.audit_receipts(output(4, backend.upper(), transport), backend, 4, transport)
                with self.assertRaises(ValueError):
                    run.audit_receipts(output(4, backend.upper(), transport), backend, 4,
                                       "direct" if transport == "staged" else "staged")

    def test_missing_duplicate_cpu_wrong_version_and_shared_gpu_rejected(self):
        baseline = [receipt(rank, 4) for rank in range(4)]
        variants = [baseline[:-1], baseline + [baseline[0]]]
        for key, value in (("backend", "Serial"), ("kokkos_version", 40300), ("ranks", 2),
                           ("local_rank", 0), ("device", -1), ("device_id", baseline[0]["device_id"])):
            candidate = copy.deepcopy(baseline)
            candidate[1][key] = value
            variants.append(candidate)
        for records in variants:
            with self.assertRaises(ValueError):
                run.audit_receipts("\n".join(run.RECEIPT + json.dumps(record) for record in records), "cuda", 4, "direct")

    def test_cross_node_receipts_pass(self):
        records = [receipt(rank, 4) for rank in range(4)]
        for rank, record in enumerate(records):
            record.update(hostname="node{}".format(rank // 2), local_rank=rank % 2, local_ranks=2,
                          device=rank % 2, device_id="0000:{:02x}:00.0".format(rank % 2))
        run.audit_receipts("\n".join(run.RECEIPT + json.dumps(record) for record in records), "cuda", 4, "direct")

    def test_zero_missing_skipped_failed_and_disabled_tests_rejected(self):
        for cases in (0, 1, 3, 5):
            with self.assertRaises(ValueError):
                run.audit_junit(self.junit(output(1, tests=cases), {"morton_test": 1}),
                                {"morton_test": 1}, "cuda", "direct", True)
        text = output(1, tests=4).replace('"skipped": 0', '"skipped": 1')
        with self.assertRaises(ValueError):
            run.audit_junit(self.junit(text, {"morton_test": 1}), {"morton_test": 1}, "cuda", "direct", True)
        with self.assertRaises(ValueError):
            run.audit_junit(self.junit(expected={}), run.SUITES, "cuda", "direct", True)
        for tag in ("skipped", "failure", "error"):
            path = self.junit(expected={"morton_test": 1})
            tree = ET.parse(path)
            ET.SubElement(tree.getroot().find("testcase"), tag)
            tree.write(path)
            with self.assertRaises(ValueError):
                run.audit_junit(path, {"morton_test": 1}, "cuda", "direct", True)

    def test_source_mutation_and_extra_files_rejected(self):
        source = self.root / "source"
        source.mkdir()
        file = source / "CMakeLists.txt"
        file.write_text("project(fixture)\n")
        manifest = self.root / "source-manifest.json"
        manifest.write_text(json.dumps({"files": {"CMakeLists.txt": run.digest(file)}}))
        run.verify_snapshot(source)
        file.write_text("changed\n")
        with self.assertRaises(ValueError):
            run.verify_snapshot(source)
        manifest.write_text(json.dumps({"files": {"CMakeLists.txt": run.digest(file)}}))
        (source / "unexpected.txt").write_text("extra")
        with self.assertRaises(ValueError):
            run.verify_snapshot(source)

    def test_compiler_and_launcher_symlink_spelling_preserved(self):
        binary = self.root / "generic-driver"
        binary.write_text('#!/bin/sh\nprintf "%s\\n" "${0##*/}"\n')
        binary.chmod(0o755)
        compiler = self.root / "clang++"
        compiler.symlink_to(binary)
        launcher = self.root / "mpiexec"
        launcher.symlink_to(binary)
        args = self.arguments()
        args.cxx, args.launcher = compiler, launcher
        configure = dict(run.plan(args))["configure"]
        self.assertIn("-DCMAKE_CXX_COMPILER=" + str(compiler), configure)
        self.assertIn("-DMPIEXEC_EXECUTABLE=" + str(launcher), configure)
        self.assertIn("-DCMAKE_CUDA_HOST_COMPILER=" + str(compiler), dict(build_deps.commands(args))["configure"])
        version = run.tool_version(compiler, self.root / "version.log")
        self.assertEqual(version["version_banner"].strip(), "clang++")
        self.assertEqual(version["resolved_binary"], str(binary.resolve()))

    def test_binary_hash_targets_executable_after_launcher(self):
        for name, executable in (("morton_test", "morton_test"), ("structured_install", "mars_structured_consumer")):
            path = self.root / executable
            test = {"name": name, "command": ["/usr/bin/srun", "-n", "1", "--export=ALL", str(path)]}
            self.assertEqual(run.test_binary(test), path)
        with self.assertRaises(ValueError):
            run.test_binary({"name": "morton_test", "command": ["/usr/bin/srun", "-n", "1"]})

    def test_export_scope_excludes_other_backends_and_artifacts(self):
        for name in ("internal-notes/brief.md", "SYNC.md", "AGENTS.md", "docs/private.md", "mesh.bin",
                     "backend/serial/source.hpp", "backend/amr/source.hpp", "backend/distributed/unstructured/domain.hpp",
                     "backend/distributed/testsuite/unstructured/mars_domain_test.cu",
                     ".local-worktrees/file.hpp", "backend/distributed/mesh_manager/pump.hpp"):
            self.assertFalse(snapshot_export.eligible(name), name)
        for name in ("CMakeLists.txt", "cmake/MarsOptions.cmake", "core/mars.hpp",
                     "backend/distributed/testsuite/common/mars_structured_test.hpp",
                     "tests/structured/mars_structured_runtime.hpp", "tests/structured/cluster/run.py"):
            self.assertTrue(snapshot_export.eligible(name), name)

    def test_native_build_plan_and_offline_dependency_options(self):
        for backend in ("cuda", "hip"):
            for transport in ("direct", "staged"):
                args = self.arguments(backend, transport)
                steps = dict(run.plan(args))
                self.assertIn("-DMARS_ENABLE_{}=ON".format(backend.upper()), steps["configure"])
                self.assertIn("-DMARS_ENABLE_CUDAUVM=OFF", steps["configure"])
                self.assertIn("-DMARS_ENABLE_RDMA=" + ("ON" if transport == "direct" else "OFF"), steps["configure"])
                self.assertIn("-DFETCHCONTENT_FULLY_DISCONNECTED=ON", steps["configure"])
                self.assertIn("-DMARS_ENABLE_UNSTRUCTURED=OFF", steps["configure"])
                self.assertIn("-DMPIEXEC_PREFLAGS=--export=ALL;--cpu-bind=cores", steps["configure"])
                self.assertIn("--no-tests=error", steps["tests"])
                self.assertIn("--test-output-size-passed", steps["tests"])
                self.assertIn("--test-output-size-passed", steps["consumer-tests"])
                self.assertEqual(steps["tests"][steps["tests"].index("--parallel") + 1], "1")
                self.assertEqual(steps["consumer-tests"][steps["consumer-tests"].index("--parallel") + 1], "1")
                deps = dict(build_deps.commands(args))
                self.assertIn("-DKokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE=ON", deps["configure"])
                self.assertIn("-DCMAKE_INSTALL_LIBDIR=lib", deps["configure"])

    def test_plan_only_never_claims_pass(self):
        args = self.arguments()
        args.plan_only = True
        with patch.object(run, "verify_snapshot", return_value="fixture"), \
             patch.object(run, "verify_dependencies", return_value={}), \
             patch.object(run.platform, "platform", return_value="fixture"), \
             patch.object(run.subprocess, "run") as execute:
            result = run.run(args)
        self.assertEqual(execute.call_count, 0)
        self.assertEqual(result["status"], "PLANNED_NOT_VALIDATED")
        self.assertEqual(json.loads((args.work / "result.json").read_text())["status"], "PLANNED_NOT_VALIDATED")

    def test_subprocess_failure_keeps_nonpassing_receipt(self):
        args = self.arguments()
        with patch.object(run, "verify_snapshot", return_value="fixture"), \
             patch.object(run, "verify_dependencies", return_value={}), \
             patch.object(run.platform, "platform", return_value="fixture"), \
             patch.object(run, "tool_version", return_value={}), \
             patch.object(run.subprocess, "run", return_value=subprocess.CompletedProcess([], 7)):
            with self.assertRaises(ValueError):
                run.run(args)
        status = json.loads((args.work / "result.json").read_text())
        self.assertEqual(status["status"], "FAILED")
        self.assertIn("exit 7", status["error"])

    def test_srun_step_cannot_launch_nested_job(self):
        env = {"PATH": "/bin:/usr/bin", "SLURM_JOB_ID": "1", "SLURM_STEP_ID": "0"}
        command = ["bash", str(Path(__file__).with_name("job.sh")), "snapshot", "results", "cuda", "90", "c++", "nvcc", "direct"]
        process = subprocess.run(command, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(process.returncode, 2)
        self.assertIn("resource-holding srun step", process.stderr)


@unittest.skipUnless(hasattr(signal, "pthread_sigmask"), "POSIX signal masks required")
class ChildSignals(unittest.TestCase):
    def test_unblocks_only_child_sigchld(self):
        previous = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGCHLD, signal.SIGUSR1})
        try:
            command = [sys.executable, "-c", "import signal,json; print(json.dumps(sorted(int(s) for s in signal.pthread_sigmask(signal.SIG_BLOCK,set()))))"]
            result = run.run_child(command, stdout=subprocess.PIPE, text=True, timeout=10)
            self.assertEqual(result.returncode, 0)
            mask = json.loads(result.stdout)
            self.assertNotIn(int(signal.SIGCHLD), mask)
            self.assertIn(int(signal.SIGUSR1), mask)
            self.assertTrue(run.child_signal_policy()["parent_sigchld_blocked"])
        finally:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous)

    @unittest.skipUnless(shutil.which("cmake"), "CMake required")
    def test_cmake_configures_with_blocked_parent_sigchld(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source"
            source.mkdir()
            (source / "CMakeLists.txt").write_text("cmake_minimum_required(VERSION 3.22)\nproject(SignalCheck LANGUAGES CXX)\n")
            previous = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGCHLD})
            try:
                result = run.run_child(["cmake", "-S", str(source), "-B", str(Path(directory) / "build")],
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stdout)
                self.assertTrue(run.child_signal_policy()["parent_sigchld_blocked"])
            finally:
                signal.pthread_sigmask(signal.SIG_SETMASK, previous)


if __name__ == "__main__":
    unittest.main()
