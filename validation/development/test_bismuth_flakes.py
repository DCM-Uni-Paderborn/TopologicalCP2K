"""Independent checks of flake geometry, metric integration and full-band gauges."""

from contextlib import redirect_stdout
import io
from itertools import product
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import numpy as np
from scipy.linalg import block_diag, eigh, eigvalsh

from analyze_bismuth_flakes import moments, native_comparison, operators, queries
from run_bismuth_flakes import geometry, inputs, nnkp

sys.path.insert(0, str(Path(os.environ["CP2K_ROOT"]) / "build-serial"))
from gaussian_states import Shell, primitive_integrals
from check_bloch_localizer import aii_skew, skew_sign


def s_basis(shift=0.):
    centers = np.array([[.1, .5, -.2], [1.1, -.2, 1.3]]) + shift
    alpha = .7
    norm = (2 * alpha / np.pi)**.75
    shells = [Shell(c, np.array([i]), np.array([alpha]), np.array([15.]),
                    np.zeros((1, 3), int), np.array([[[norm]]]), 15.)
              for i, c in enumerate(centers)]
    return SimpleNamespace(positions=centers, atom_sizes=[1, 1], periodic=np.zeros(3, int),
                           kpoints=np.zeros((1, 3)), shells=shells)


class FlakeChecks(unittest.TestCase):
    def test_native_comparison(self):
        text = """SPECTRAL_LOCALIZER| Position [bohr]: 0 1 2
SPECTRAL_LOCALIZER| Energy [hartree]: -0.1
SPECTRAL_LOCALIZER| Kappa [hartree/bohr]: 0.01
SPECTRAL_LOCALIZER| Gap bracket [hartree]: 0.02999999 0.03000001
SPECTRAL_LOCALIZER| Gap [hartree]: 0.02999999
SPECTRAL_LOCALIZER| Z2 index: 1
"""
        query = dict(position_bohr=[0., 1., 2.], energy_hartree=-.1,
                     kappa_hartree_bohr=.01, gap_hartree=.03, z2=1)
        self.assertTrue(native_comparison(text, [query])["accepted"])
        for key, value in (("gap_hartree", .031), ("z2", 0)):
            self.assertFalse(native_comparison(text, [dict(query, **{key: value})])["accepted"])
        with self.assertRaises(ValueError):
            native_comparison(text, [dict(query, energy_hartree=-.2)])
        with self.assertRaises(ValueError):
            native_comparison(text, [])
        dense = "\n".join(line for line in text.splitlines() if "Gap bracket" not in line)
        self.assertTrue(native_comparison(dense, [query])["accepted"])

    def test_geometry(self):
        for size in (0, 1, 3, 5):
            cell, atoms = geometry(size, 20.)
            atoms = np.array(atoms)
            self.assertEqual(len(atoms), 2 * max(1, size**2))
            self.assertAlmostEqual(np.ptp(atoms[:, 2]), 1.74, places=12)
            if size:
                np.testing.assert_allclose(np.mean(atoms, axis=0), np.diag(cell) / 2, atol=1e-14)
                np.testing.assert_allclose(np.min(atoms, axis=0), 10., atol=1e-14)
            else:
                self.assertAlmostEqual(np.linalg.norm(cell[1]), 4.33)
            nearest = np.linalg.norm(atoms[0] - atoms[1])
            self.assertAlmostEqual(nearest, np.sqrt(4.33**2/3 + 1.74**2), places=12)

    def test_input_boundaries(self):
        case = dict(size=0, vacuum=20., mode="wilson", basis="TZVP", scf_mesh=8,
                    temperature=300, cutoff=400)
        text = inputs(SimpleNamespace(**case))["input.inp"]
        bands = next(line.split()[1:] for line in text.splitlines() if line.strip().startswith("EXCLUDE_BANDS"))
        self.assertEqual(list(map(int, bands)), list(range(11, 69)))
        case.update(mode="spectrum", size=3)
        files = inputs(SimpleNamespace(**case))
        self.assertEqual(files["input.inp"].count("PERIODIC NONE"), 2)
        self.assertNotIn("&KPOINTS", files["input.inp"])
        self.assertIn("1 1 0 0 0", files["gamma.nnkp"])
        case.update(mode="localizer", export_spectrum=True, energy=[-.1], kappa=[.003],
                    offset=[0.], solver="TACHO", checkpoint=True, restart=Path("start.wfn"))
        files = inputs(SimpleNamespace(**case))
        for expected in ("&SPECTRAL_LOCALIZER", "&WANNIER90", "SCF_GUESS RESTART", "QS_SCF 5"):
            self.assertIn(expected, files["input.inp"])
        self.assertIn("gamma.nnkp", files)
        for key, value in (("vacuum", float("nan")), ("kappa", [-.1])):
            with self.assertRaises(ValueError):
                inputs(SimpleNamespace(**dict(case, **{key: value})))
        lines = nnkp(geometry(0, 20)[0]).splitlines()
        real = np.array([line.split() for line in lines[1:4]], float)
        reciprocal = np.array([line.split() for line in lines[6:9]], float)
        np.testing.assert_allclose(real @ reciprocal.T, 2*np.pi*np.eye(3), atol=1e-14)

    def test_analytic_s_moments_and_origin(self):
        snapshot = s_basis()
        integrated = moments(snapshot, primitive_integrals)
        distance = np.linalg.norm(snapshot.positions[0] - snapshot.positions[1])
        overlap = np.array([[1., np.exp(-.35*distance**2)], [np.exp(-.35*distance**2), 1.]])
        np.testing.assert_allclose(integrated[0], overlap, atol=1e-14)
        for axis in range(3):
            centers = snapshot.positions[:, axis]
            np.testing.assert_allclose(integrated[axis+1], overlap*(centers[:, None]+centers[None, :])/2,
                                       atol=1e-14)
        moved = moments(s_basis(3.), primitive_integrals)
        np.testing.assert_allclose(moved[0], integrated[0], atol=1e-14)
        np.testing.assert_allclose(moved[1:], integrated[1:]+3*integrated[0], atol=1e-14)

    def test_p_d_moments_against_quadrature(self):
        powers = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [2, 0, 0], [1, 1, 0]])
        a, b = np.array([.2, -.7, .8]), np.array([1.3, .2, -.3])
        alpha, beta = .6, .9
        total = alpha + beta
        center = (alpha*a+beta*b)/total
        nodes, weights = np.polynomial.hermite.hermgauss(4)
        grid = np.array(list(product(range(4), repeat=3)))
        points = center + nodes[grid]/np.sqrt(total)
        weights = np.prod(weights[grid], axis=1) * np.exp(-alpha*beta/total*np.sum((a-b)**2))/total**1.5
        left = np.prod((points[None, :, :]-a)**powers[:, None, :], axis=2)
        right = np.prod((points[None, :, :]-b)**powers[:, None, :], axis=2)
        overlap = primitive_integrals(a, b, alpha, beta, powers, powers, np.zeros(3))
        np.testing.assert_allclose(overlap, (left*weights)@right.T, atol=1e-14)
        for axis in range(3):
            raised = powers.copy()
            raised[:, axis] += 1
            analytic = primitive_integrals(a, b, alpha, beta, raised, powers, np.zeros(3))+a[axis]*overlap
            np.testing.assert_allclose(analytic, (left*(weights*points[:, axis]))@right.T, atol=1e-14)

    def test_degenerate_gauge_and_translation(self):
        snapshot = s_basis()
        integrated = moments(snapshot, primitive_integrals)
        w, v = eigh(integrated[0])
        rinv = (v / np.sqrt(w)) @ v.T.conj()
        energy = np.array([-.2, -.2, .4, .4])
        frame = np.eye(4)[:, [0, 2, 1, 3]].astype(complex)
        snapshot.bands = np.arange(4)
        snapshot.energies = energy[None, :]
        snapshot.coefficients = (block_diag(rinv, rinv) @ frame).reshape(1, 2, 2, 4)
        h, xyz, _ = operators(snapshot, integrated)
        rng = np.random.default_rng(1927)
        rotations = [np.linalg.qr(rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2)))[0] for _ in range(2)]
        snapshot.coefficients = (snapshot.coefficients.reshape(4, 4) @ block_diag(*rotations)).reshape(1, 2, 2, 4)
        changed, _, _ = operators(snapshot, integrated)
        np.testing.assert_allclose(h, changed, atol=1e-14)
        with redirect_stdout(io.StringIO()):
            original = queries(h, xyz, [.1], [.02], [np.zeros(3)], aii_skew, skew_sign)[0]
            translated = queries(h, [x+2*np.eye(4) for x in xyz], [.1], [.02],
                                 [2*np.ones(3)], aii_skew, skew_sign)[0]
        self.assertEqual(original["z2"], translated["z2"])
        self.assertAlmostEqual(original["gap_hartree"], translated["gap_hartree"], places=14)

    def test_covariant_metric_gap(self):
        integrated = moments(s_basis(), primitive_integrals)
        w, v = eigh(integrated[0])
        root = block_diag(*[(v*np.sqrt(w))@v.T.conj()]*2)
        inverse = np.linalg.inv(root)
        a = np.array([[.2, .03+.04j], [.03-.04j, -.1]])
        b = np.array([[0, .01-.02j], [-.01+.02j, 0]])
        h = np.block([[a, b], [-b.conj(), a.conj()]])
        overlap = block_diag(integrated[0], integrated[0])
        x, y = [block_diag(r, r) for r in integrated[1:3]]
        mass = root@h@root - .06*overlap
        cross = .02*(x-.3*overlap-1j*(y+.4*overlap))
        covariant = np.block([[mass, cross], [cross.conj().T, -mass]])
        metric = block_diag(overlap, overlap)
        inverse = block_diag(inverse, inverse)
        orthonormal = inverse@covariant@inverse
        np.testing.assert_allclose(eigvalsh(covariant, metric), eigvalsh(orthonormal), atol=1e-14)
        skew, _ = aii_skew(orthonormal)
        self.assertIn(skew_sign(skew), (-1, 1))


if __name__ == "__main__":
    unittest.main()
