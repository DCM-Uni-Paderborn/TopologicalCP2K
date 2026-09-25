"""Input and evidence-integrity checks for the independent material controls."""

import hashlib
import io
import json
import os
from pathlib import Path
import tarfile
import tempfile
import unittest

import numpy as np

from replay_stanene_controls import verify
from archive_stanene_controls import retain
from run_stanene_controls import cases, inputs


class ControlInputs(unittest.TestCase):
    def test_mesh_and_fixed_buckling(self):
        root = Path(os.environ["CP2K_ROOT"])
        controls = cases()
        self.assertEqual(len(controls), len({c["label"] for c in controls}))
        for case in controls:
            with self.subTest(case=case["label"]):
                files = inputs(case, root)
                inp = files["input.inp"]
                atoms = [line.split() for line in inp.splitlines() if line.strip().startswith("Sn ")]
                self.assertAlmostEqual((float(atoms[0][-1]) - float(atoms[1][-1])) * case["height"], .852, places=12)
                self.assertIn(f"BASIS_SET {case['basis']}-MOLOPT-GGA-GTH-q4", inp)
                self.assertIn(f"CUTOFF {case['cutoff']}\n", inp)
                self.assertIn(f"REL_CUTOFF {case['rel_cutoff']}\n", inp)
                self.assertIn(f"MONKHORST-PACK {case['scf_mesh']} {case['scf_mesh']} 1", inp)
                if case["mode"] != "bloch":
                    continue
                lines = files["mesh.nnkp"].splitlines()
                real = np.array([line.split() for line in lines[1:4]], float)
                reciprocal = np.array([line.split() for line in lines[6:9]], float)
                np.testing.assert_allclose(real @ reciprocal.T, 2 * np.pi * np.eye(3), atol=1e-14)
                self.assertEqual(real[2, 2], case["height"])
                start = lines.index("begin kpoints")
                nk = int(lines[start + 1])
                self.assertEqual(nk, case["size"]**2)
                points = np.array([line.split() for line in lines[start+2:start+2+nk]], float)
                start = lines.index("begin nnkpts")
                links = np.array([line.split() for line in lines[start+2:lines.index("end nnkpts")]], int)
                self.assertEqual(links.shape, (2 * nk, 5))
                for a, b, x, y, z in links:
                    step = points[b-1] + [x, y, z] - points[a-1]
                    self.assertAlmostEqual(sum(step), 1 / case["size"], places=14)
                    self.assertTrue(min(step) >= -1e-14)

    def test_changed_archive_rejected(self):
        with tempfile.TemporaryDirectory(prefix="stanene-integrity-") as tmp:
            path = Path(tmp) / "record.tar.gz"
            payload = b"complete evidence\n"
            manifest = dict(files={"result.txt": dict(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())})
            with tarfile.open(path, "w:gz") as archive:
                for name, data in (("result.txt", payload), ("manifest.json", json.dumps(manifest).encode())):
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    archive.addfile(info, io.BytesIO(data))
            expected = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(verify(path, expected), manifest)
            with self.assertRaises(AssertionError):
                verify(path, "0" * 64)
            path.write_bytes(path.read_bytes() + b"changed")
            with self.assertRaises(AssertionError):
                verify(path, expected)

    def test_archive_resume_checks_sources(self):
        with tempfile.TemporaryDirectory(prefix="stanene-resume-") as tmp:
            base = Path(tmp)
            source, archive = base / "input.txt", base / "case.tar.gz"
            source.write_text("verified source\n")
            files, metadata = {"input.txt": source}, {"case": "control"}
            record = retain(archive, files, metadata)
            self.assertEqual(retain(archive, files, metadata), record)
            source.write_text("changed source\n")
            with self.assertRaises(AssertionError):
                retain(archive, files, metadata)


if __name__ == "__main__":
    unittest.main()
