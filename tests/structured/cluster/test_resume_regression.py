#!/usr/bin/env python3
"""Exercise continuation with real CTest and synthetic executables, without a GPU."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
import regression
import resume_regression as resume
import run


class ResumeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.prepared = self.root / "prepared"
        self.build = self.prepared / "build"
        self.build.mkdir(parents=True)
        self.binary = self.build / "fixture"
        self.binary.write_text("#!/usr/bin/env python3\n" + '''import os,sys
from pathlib import Path
assert os.environ['MARS_GATE_LAYOUT']=='multi'
for arg in sys.argv[1:]:
    if arg.startswith('--rows='):
        p=Path(arg.split('=',1)[1]); p.parent.mkdir(parents=True,exist_ok=True); p.write_text('new result')
print('[ RUN      ] MarsNodeHaloTopologyTest.SingleRankAccessor')
print('diagnostic\\n'*200)
print('[       OK ] MarsNodeHaloTopologyTest.SingleRankAccessor')
print('[  PASSED  ] 1 test.')
''')
        self.binary.chmod(0o755)
        tests = []
        for name in sorted(regression.NAMES):
            tests.append({"name": name, "command": [str(self.binary), "--rows=" + str(self.build / "rows" / name)],
                          "properties": [{"name": "LABELS", "value": ["release"]},
                                         {"name": "TIMEOUT", "value": 30.0},
                                         {"name": "WORKING_DIRECTORY", "value": str(self.build / "work")} ]})
        tests[0]["properties"].append({"name": "FIXTURES_SETUP", "value": ["fixture"]})
        tests[1]["properties"].append({"name": "FIXTURES_REQUIRED", "value": ["fixture"]})
        tests[1]["properties"].append({"name": "DEPENDS", "value": [tests[0]["name"]]})
        (self.prepared / "inventory.json").write_text(json.dumps({"tests": tests}))
        (self.prepared / "result.json").write_text(json.dumps({"backend": "cuda", "source_manifest_sha256": "fixture",
            "error": "Single-rank halo regression did not execute its test"}))
        root = ET.Element("testsuite")
        for name in regression.NAMES:
            case = ET.SubElement(root, "testcase", name=name, status="run")
            ET.SubElement(case, "system-out").text = "truncated before summary"
        ET.ElementTree(root).write(self.prepared / "release-single.xml")
        (self.prepared / "release-single.log").write_text(
            "Start 52: release_node_halo_single_rank\n52: [ RUN ] MarsNodeHaloTopologyTest.SingleRankAccessor\n"
            "52: [ OK ] MarsNodeHaloTopologyTest.SingleRankAccessor\n52: [ PASSED ] 1 test.\n")
        files = {str(p): run.digest(p) for p in self.prepared.rglob("*") if p.is_file()}
        self.expected = {"backend": "cuda", "prepared": str(self.prepared), "source": str(self.root / "source"),
                         "source_manifest_sha256": "fixture", "dependencies": {}, "files": files,
                         "binaries": {str(self.binary): run.digest(self.binary)}}

    def tearDown(self):
        self.temp.cleanup()

    def test_continuation_preserves_prepared_files_and_rejects_binary_mutation(self):
        with patch.object(run, "verify_snapshot", return_value="fixture"):
            self.assertEqual(resume.resume(self.expected, self.root / "continuation"), 0)
            resume.verify_files(self.expected["files"])
            self.assertFalse((self.build / "rows").exists())
            self.assertTrue((self.root / "continuation/ctest/rows/release_node_halo_single_rank").is_file())
            result = json.loads((self.root / "continuation/result.json").read_text())
            self.assertEqual(len(result["single_node"]), 50)
            self.assertEqual(len(result["multi_node"]), 50)
            self.binary.write_text(self.binary.read_text() + "# mutation\n")
            self.assertEqual(resume.resume(self.expected, self.root / "rejected"), 1)
            rejected = json.loads((self.root / "rejected/result.json").read_text())
            self.assertIn("Prepared file changed", rejected["error"])
            self.assertNotIn("command", rejected)

    def test_real_ctest_truncation_replays_from_full_log(self):
        target = self.root / "truncation"
        target.mkdir()
        inventory = json.loads((self.prepared / "inventory.json").read_text())
        resume.mirror_inventory(inventory, self.build, target, self.expected["binaries"])
        xml, log = target / "result.xml", target / "result.log"
        with log.open("w") as stream:
            subprocess.run(["ctest", "--test-dir", str(target), "-L", "release", "--parallel", "1", "-V",
                            "--test-output-size-passed", "1024", "--output-junit", str(xml)],
                           stdout=stream, stderr=subprocess.STDOUT, env=dict(os.environ, MARS_GATE_LAYOUT="multi"), check=True)
        with self.assertRaises(ValueError):
            regression.audit_junit(xml)
        self.assertEqual(len(regression.audit_junit(xml, log)), 50)
        resume.verify_files(self.expected["files"])


if __name__ == "__main__":
    unittest.main()
