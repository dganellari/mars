"""Public profile/provenance rejection checks; numerical checks use the C++ gate."""
import json
from pathlib import Path
import tempfile
import unittest

from prepare_openaccel_geometry import DECK_SHA256, MESH_SHA256, read_frames, verify_run


class GeometryCaptureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.exports = self.root / 'exports'
        self.exports.mkdir()
        self.hashes = {'pressure.jsonl':'captured-hash'}
        self.record = dict(returncode=0,status='update_capture_completed',fixture='public_channel',
                           bundle={'sha256':{'input.i':DECK_SHA256,'channel.exo':MESH_SHA256}},
                           export_sha256=self.hashes,binary_sha256='binary-hash')

    def verify(self):
        (self.root/'run.json').write_text(json.dumps(self.record))
        return verify_run(self.exports,self.hashes,'export_sha256')

    def test_completed_public_profile(self):
        self.assertEqual(self.verify(),'binary-hash')

    def test_wrong_deck_or_mesh(self):
        for name in ('input.i','channel.exo'):
            original=self.record['bundle']['sha256'][name]
            self.record['bundle']['sha256'][name]='different'
            with self.assertRaisesRegex(ValueError,'profile'): self.verify()
            self.record['bundle']['sha256'][name]=original

    def test_failed_or_smoke_only_run(self):
        self.record['status']='runtime_smoke_completed'
        with self.assertRaisesRegex(ValueError,'did not complete'): self.verify()
        self.record['status']='update_capture_completed'; self.record['returncode']=1
        with self.assertRaisesRegex(ValueError,'did not complete'): self.verify()

    def test_export_changes(self):
        self.record['export_sha256']={'pressure.jsonl':'changed'}
        with self.assertRaisesRegex(ValueError,'hashes'): self.verify()

    def test_distributed_capture_rejected(self):
        (self.exports/'capture.jsonl').write_text(json.dumps(dict(kind='header',ranks=2,rank=0))+'\n')
        with self.assertRaisesRegex(ValueError,'one reference rank'): read_frames(self.exports)

    def test_missing_state_rejected(self):
        (self.exports/'capture.jsonl').write_text(json.dumps(
            dict(kind='header',ranks=1,rank=0,stage='pressure.interior',call=1))+'\n')
        with self.assertRaisesRegex(ValueError,'both captured fields'): read_frames(self.exports)


if __name__ == '__main__': unittest.main()
