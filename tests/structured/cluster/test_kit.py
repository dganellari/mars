#!/usr/bin/env python3
"""Offline checks of the gate itself; these do not execute GPU kernels."""

import copy
import importlib.util
import json
import os
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
            "driver_version": 12000, "transport": transport, "concurrency": 32}


def output(ranks, backend="CUDA", transport="direct", tests=3):
    lines = []
    for rank in range(ranks):
        lines.append(run.RECEIPT + json.dumps(receipt(rank, ranks, backend, transport)))
        lines.append("[==========] Running {} tests from 1 test suite.".format(tests))
        lines.append(run.TEST_RECEIPT + json.dumps({"rank": rank, "ranks": ranks, "cases": tests,
                                                   "passed": tests, "skipped": 0}))
    return "\n".join(lines)


def halo_output():
    return ("[ RUN      ] MarsNodeHaloTopologyTest.SingleRankAccessor\n" + "diagnostic\n" * 200 +
            "[       OK ] MarsNodeHaloTopologyTest.SingleRankAccessor\n[  PASSED  ] 1 test.")


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
        self.assertEqual(sum(run.GTEST_CASES[name.split("_mpi")[0]] for name in run.SUITES), 171)
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

    def test_cuda_dependency_uses_mpi_compatible_allocator(self):
        for backend in ("cuda", "hip", "serial", "openmp"):
            configure = dict(build_deps.commands(self.arguments(backend)))["configure"]
            option = "-DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF"
            self.assertEqual(option in configure, backend == "cuda")

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

    def test_cpu_receipts_threads_transport_and_nodes(self):
        for backend in ("serial", "openmp"):
            for transport in ("direct", "staged"):
                records = [receipt(rank, 4, backend.upper(), transport) for rank in range(4)]
                for rank, record in enumerate(records):
                    record.update(device=-1, device_id="cpu", device_name="cpu", runtime_version=0,
                                  driver_version=0, mpi_device_buffer=False, concurrency=2 if backend == "openmp" else 1,
                                  hostname="host" + str(rank // 2), local_rank=rank % 2, local_ranks=2)
                text = "\n".join(run.RECEIPT + json.dumps(x) for x in records)
                run.audit_receipts(text, backend, 4, transport, 2, 2 if backend == "openmp" else 1)
                with self.assertRaises(ValueError):
                    run.audit_receipts(text, backend, 4, transport, 3)
                with self.assertRaises(ValueError):
                    run.audit_receipts(text, backend, 4, transport, 2, 3)
                with self.assertRaises(ValueError):
                    run.audit_receipts(text, "cuda", 4, transport)
                args = self.arguments(backend, transport)
                args.gpu_compiler = None
                config = dict(run.plan(args))["configure"]
                self.assertIn("-DMARS_ENABLE_CUDA=OFF", config)
                self.assertIn("-DMARS_ENABLE_HIP=OFF", config)
                self.assertIn("-DMARS_STRUCTURED_FORCE_MPI=ON", config)
                self.assertFalse(any("None" in x for x in config))
                deps = dict(build_deps.commands(args))["configure"]
                self.assertIn("-DKokkos_ENABLE_OPENMP=" + ("ON" if backend == "openmp" else "OFF"), deps)

    def test_actual_scheduler_command_places_ranks(self):
        import mpi_launch
        for ranks in (1, 2, 4):
            args = SimpleNamespace(launcher="/usr/bin/srun", ranks=ranks,
                launcher_arg=["--uenv=fixture", "--cpu-bind=cores"], binary=["/build/test"])
            for layout in ("single", "multi"):
                argv = mpi_launch.command(args, layout)
                nodes = min(ranks, 2) if layout == "multi" else 1
                self.assertIn("--nodes=" + str(nodes), argv)
                self.assertIn("--ntasks-per-node=" + str(ranks // nodes), argv)
                self.assertIn("--export=ALL", argv)
                self.assertEqual(argv[-1], "/build/test")

    def test_performance_gate_rejects_incomplete_or_slow_results(self):
        import performance
        samples = []
        for ranks in (1, 2, 4):
            cases = {name: {metric: 1.0 for metric in performance.METRICS} for name in performance.CASES}
            for repeat in range(5):
                for variant in ("reference", "candidate"):
                    samples.append(dict(layout="single", ranks=ranks, round=repeat, variant=variant,
                                        audit={"cases": copy.deepcopy(cases)}))
        self.assertTrue(all(x["pass"] for x in performance.compare(samples)))
        with self.assertRaises(ValueError):
            performance.compare(samples[:-1])
        for s in samples:
            if s["variant"] == "candidate":
                s["audit"]["cases"]["hilbert-hex-full"]["halo_seconds"] = 2.0
        self.assertFalse(all(x["pass"] for x in performance.compare(samples)))
        rows = []
        for name, (nx, ny, nz) in performance.CASES.items():
            row = dict(case=name, nx=nx, ny=ny, nz=nz, cells=nx*ny*(nz or 1),
                       owned_dofs=(nx+1)*(ny+1)*(nz+1 if nz else 1), ranks=1, iterations=20, errors=0)
            row.update({metric: 1.0 for metric in performance.METRICS})
            rows.append(row)
        def checked(rows):
            return performance.audit(output(1) + "\n" + "\n".join(performance.PREFIX + json.dumps(x) for x in rows),
                                     1, "cuda", "direct", 1, 1)
        checked(rows)
        for key, value in (("errors",1),("owned_dofs",0),("halo_seconds",float("nan")),("nx",0),("iterations",0)):
            broken = copy.deepcopy(rows)
            broken[0][key] = value
            with self.assertRaises(ValueError): checked(broken)
        with self.assertRaises(ValueError): checked(rows[:-1])

    def test_release_regression_matrix_rejects_missing_and_skipped(self):
        import regression
        self.assertEqual(len(regression.NAMES), 50)
        root = ET.Element("testsuite")
        for name in regression.NAMES:
            case = ET.SubElement(root, "testcase", name=name, status="run")
            if name == "release_node_halo_single_rank":
                ET.SubElement(case, "system-out").text = halo_output()
        path = self.root / "release.xml"
        ET.ElementTree(root).write(path)
        self.assertEqual(len(regression.audit_junit(path)), 50)
        halo = root.find("testcase[@name='release_node_halo_single_rank']/system-out")
        for text in ("", "[  PASSED  ] 0 tests.", "[  SKIPPED ] 1 test. [ PASSED ] 1 test."):
            halo.text = text
            ET.ElementTree(root).write(path)
            with self.assertRaises(ValueError): regression.audit_junit(path)
        halo.text = halo_output()
        ET.SubElement(root.find("testcase"), "skipped")
        ET.ElementTree(root).write(path)
        with self.assertRaises(ValueError): regression.audit_junit(path)
        root.remove(root.find("testcase"))
        ET.ElementTree(root).write(path)
        with self.assertRaises(ValueError): regression.audit_junit(path)

    def test_release_truncated_halo_uses_full_case_log(self):
        import regression
        root = ET.Element("testsuite")
        for name in regression.NAMES:
            case = ET.SubElement(root, "testcase", name=name, status="run")
            if name == "release_node_halo_single_rank":
                ET.SubElement(case, "system-out").text = halo_output()[:1024]
        path, log = self.root / "release.xml", self.root / "release.log"
        ET.ElementTree(root).write(path)
        text = "Start 52: release_node_halo_single_rank\n" + "\n".join(
            "52: " + line for line in halo_output().splitlines())
        log.write_text(text)
        with self.assertRaises(ValueError): regression.audit_junit(path)
        self.assertEqual(len(regression.audit_junit(path, log)), 50)
        for broken in (text.replace("Start 52:", "Start 53:"), text + text,
                       text.replace("SingleRankAccessor", "DifferentCase"),
                       text.replace("[       OK ]", "[  SKIPPED ]"),
                       text.replace("[  PASSED  ] 1 test.", "[  PASSED  ] 0 tests.")):
            log.write_text(broken)
            with self.assertRaises(ValueError): regression.audit_junit(path, log)

    def test_performance_retains_failure_and_runs_remaining_pairs(self):
        import performance
        source, reference = self.root / 'source', self.root / 'reference'
        target, reuse = self.root / 'candidate', self.root / 'reuse'
        target.mkdir(); reuse.mkdir()
        checked = {'status': 'PASS', 'source_manifest_sha256': 'source', 'artifacts': {},
                   'backend': 'cuda', 'transport': 'direct', 'environment': {}}
        (target / 'result.json').write_text(json.dumps(checked))
        (target / 'commands.json').write_text('{}')
        (target / 'consumer-list.log').write_text(json.dumps({'tests': [
            {'name': 'structured_install_mpi' + str(n), 'command': ['launcher', '-n', str(n), '/mars_structured_consumer']}
            for n in (1, 2, 4)]}))
        binaries = {name: self.root / (name + '-binary') for name in ('reference', 'candidate')}
        for path in binaries.values(): path.write_text('binary fixture')
        (reuse / 'result.json').write_text(json.dumps({'status': 'BUILT_NOT_RUN',
            'source_manifest_sha256': 'source', 'reference_manifest_sha256': 'reference',
            'candidate_result_sha256': run.digest(target / 'result.json'),
            'binary_paths': {k: str(v) for k, v in binaries.items()},
            'binaries': {str(p): run.digest(p) for p in binaries.values()}}))
        calls = []
        def execute(command, log, env):
            calls.append(log.name)
            log.write_text('fixture')
            if log.name == 'multi-np2-0-reference.log': raise ValueError('retained MPI failure')
        args = SimpleNamespace(source=source, reference=reference, candidate=target, reuse=reuse,
            work=self.root / 'performance', threads=1, multi_node=True, prepare_only=False)
        with patch.object(run, 'verify_snapshot', side_effect=lambda p: p.name), \
             patch.object(performance, 'execute', side_effect=execute), \
             patch.object(performance, 'audit', return_value={}):
            result = performance.validate(args)
        self.assertEqual(len(calls), 50)
        self.assertIn('multi-np2-0-candidate.log', calls)
        self.assertEqual(len(result['samples']), 49)
        self.assertEqual(len(result['failures']), 1)
        self.assertEqual(result['status'], 'FAILED')
        self.assertNotIn('comparisons', result)

    def test_release_export_excludes_private_material(self):
        import export_release
        for name in ("docs/pump.md", "internal-notes/note.md", "tests/data/private/mesh.e", "pump.stl", "mesh.bin"):
            self.assertFalse(export_release.eligible(name), name)
        for name in ("tests/data/poiseuille/poiseuille_hex_14k_elem.e", "tests/release/ReleaseTests.cmake",
                     "backend/distributed/unstructured/domain.cu"):
            self.assertTrue(export_release.eligible(name), name)

    def test_execution_environment_preserves_live_allocation(self):
        recorded = {"SLURM_JOB_ID": "old", "SLURM_JOB_NODELIST": "oldnode", "SLURM_NTASKS": "4",
                    "CUDA_VISIBLE_DEVICES": "0", "HIP_VISIBLE_DEVICES": "0", "ROCR_VISIBLE_DEVICES": "0",
                    "MPICH_GPU_SUPPORT_ENABLED": "1", "MPICH_GPU_IPC_ENABLED": "0", "OMP_NUM_THREADS": "4",
                    "MARS_STRUCTURED_EXPECT_BACKEND": "CUDA", "PATH": "/old/path", "OLD_ONLY": "old"}
        current = dict(os.environ, SLURM_JOB_ID="new", SLURM_JOB_NODELIST="newnodes", SLURM_NTASKS="8",
                       CUDA_VISIBLE_DEVICES="2,3", HIP_VISIBLE_DEVICES="2,3", ROCR_VISIBLE_DEVICES="2,3",
                       GTEST_FILTER="wrong", MARS_STRUCTURED_DIAGNOSTIC="wrong")
        with patch.dict(os.environ, current, clear=True):
            env = run.execution_environment(recorded)
        for key in ("SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_NTASKS", "CUDA_VISIBLE_DEVICES",
                    "HIP_VISIBLE_DEVICES", "ROCR_VISIBLE_DEVICES", "PATH"):
            self.assertEqual(env[key], current[key])
        for key in run.TEST_ENVIRONMENT:
            self.assertEqual(env[key], recorded[key])
        self.assertNotIn("OLD_ONLY", env)
        self.assertNotIn("GTEST_FILTER", env)
        self.assertNotIn("MARS_STRUCTURED_DIAGNOSTIC", env)
        with patch.dict(os.environ, {}, clear=True):
            env = run.execution_environment(recorded)
        self.assertNotIn("SLURM_JOB_ID", env)
        self.assertNotIn("CUDA_VISIBLE_DEVICES", env)

    def test_validation_call_sites_launch_with_current_allocation(self):
        import release
        import performance
        class StopAfterProbe(Exception):
            pass
        live = dict(os.environ, SLURM_JOB_ID="222", SLURM_JOB_NODELIST="live1,live2", SLURM_NTASKS="4",
                    SLURM_STEP_ID="batch", CUDA_VISIBLE_DEVICES="2,3", HIP_VISIBLE_DEVICES="2,3")
        old = {"SLURM_JOB_ID": "111", "SLURM_JOB_NODELIST": "expired", "SLURM_NTASKS": "1",
               "CUDA_VISIBLE_DEVICES": "0", "HIP_VISIBLE_DEVICES": "0", "OMP_NUM_THREADS": "1",
               "MARS_STRUCTURED_EXPECT_BACKEND": "CUDA", "MPICH_GPU_SUPPORT_ENABLED": "1"}
        def probe(command, log, env):
            process = run.run_child([sys.executable, "-c",
                "import os,json; print(json.dumps({k:os.environ.get(k) for k in ['SLURM_JOB_ID','SLURM_JOB_NODELIST','SLURM_NTASKS','CUDA_VISIBLE_DEVICES','HIP_VISIBLE_DEVICES']}))"],
                env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=10)
            self.assertEqual(process.returncode, 0, process.stderr)
            child = json.loads(process.stdout)
            for key in child:
                self.assertEqual(child[key], live[key], key)
            raise StopAfterProbe()
        source, reference = self.root / "source", self.root / "reference"
        prep = self.root / "prepared"
        target = prep / "cuda-direct"
        target.mkdir(parents=True)
        (prep / "result.json").write_text(json.dumps({"status":"PREPARED_NOT_RELEASE_VALIDATED",
            "backend":"cuda", "source_manifest_sha256":"source", "reference_manifest_sha256":"reference"}))
        (target / "result.json").write_text(json.dumps({"status":"PASS", "source_manifest_sha256":"source",
            "backend":"cuda", "transport":"direct", "artifacts":{}, "dependencies":{"kokkos_files":{}}, "environment":old}))
        (target / "commands.json").write_text(json.dumps({"tests":["ctest","--output-junit","old.xml"]}))
        (target / "consumer-list.log").write_text(json.dumps({"tests":[{"name":"structured_install",
            "command":["srun","-n","1","/binary/mars_structured_consumer"]}]}))
        binaries = {}
        for name in ("candidate", "reference"):
            path=self.root / (name+"-binary"); path.write_text("fixture")
            binaries[name]=str(path)
        reuse = self.root / "performance-build"
        reuse.mkdir()
        (reuse / "result.json").write_text(json.dumps({"status":"BUILT_NOT_RUN",
            "source_manifest_sha256":"source", "reference_manifest_sha256":"reference",
            "candidate_result_sha256":run.digest(target/"result.json"), "binary_paths":binaries,
            "binaries":{path:run.digest(path) for path in binaries.values()}}))
        helper = run.execution_environment
        for legacy in (False, True):
            environment = (lambda recorded: dict(os.environ, **recorded)) if legacy else helper
            expected_error = AssertionError if legacy else StopAfterProbe
            args = SimpleNamespace(source=source, reference=reference, prepared=prep, backend="cuda",
                                   threads=4, work=self.root / ("distributed-old" if legacy else "distributed-new"))
            args.work.mkdir()
            record={"source_manifest_sha256":"source", "reference_manifest_sha256":"reference", "configurations":[]}
            with patch.dict(os.environ, live, clear=True), patch.object(run, "execution_environment", side_effect=environment), \
                 patch.object(release, "execute", side_effect=probe):
                with self.assertRaises(expected_error): release.validate(args, record)
            args=SimpleNamespace(source=source, reference=reference, candidate=target, reuse=reuse, threads=1,
                multi_node=True, prepare_only=False, work=self.root / ("performance-old" if legacy else "performance-new"))
            with patch.dict(os.environ, live, clear=True), patch.object(run, "execution_environment", side_effect=environment), \
                 patch.object(run, "verify_snapshot", side_effect=lambda path: "source" if path==source else "reference"), \
                 patch.object(performance, "execute", side_effect=probe):
                with self.assertRaises(expected_error): performance.validate(args)

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
