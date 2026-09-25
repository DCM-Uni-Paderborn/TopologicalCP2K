"""Cross-basis metrics must not confuse a basis or gauge change with physics."""

import os
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import numpy as np
from scipy.linalg import solve

from compare_bismuth_basis import cross_moments, subspace_overlap

sys.path.insert(0, str(Path(os.environ["CP2K_ROOT"]) / "build-serial"))
from gaussian_states import Shell, primitive_integrals


class BasisChecks(unittest.TestCase):
    def test_gauge_and_metric_covariance(self):
        rng = np.random.default_rng(30192)
        a = np.linalg.qr(rng.normal(size=(12, 4)) + 1j*rng.normal(size=(12, 4)))[0]
        b = np.linalg.qr(rng.normal(size=(12, 4)) + 1j*rng.normal(size=(12, 4)))[0]
        reference = subspace_overlap(a, b, np.eye(12), np.eye(12), np.eye(12))
        bases = [rng.normal(size=(12, 12)) + 1j*rng.normal(size=(12, 12)) + 10*np.eye(12)
                 for _ in range(2)]
        rotations = [np.linalg.qr(rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4)))[0]
                     for _ in range(2)]
        left, right = [solve(t, c @ u) for t, c, u in zip(bases, (a, b), rotations)]
        x, y = bases
        changed = subspace_overlap(left, right, x.conj().T @ x, y.conj().T @ y, x.conj().T @ y)
        for key in ("overlap_singular_values", "projector_frobenius_squared", "projector_spectral_distance"):
            np.testing.assert_allclose(changed[key], reference[key], atol=1e-13)

    def test_rotated_degenerate_pair(self):
        a = np.eye(4, 2, dtype=complex)
        rotation = np.array([[1., 1j], [1j, 1.]]) / np.sqrt(2.)
        result = subspace_overlap(a, a @ rotation, np.eye(4), np.eye(4), np.eye(4))
        self.assertLess(abs(result["projector_frobenius_squared"]), 1e-14)
        np.testing.assert_allclose(result["overlap_singular_values"], 1., atol=1e-14)

    def test_unequal_ranks_and_unphysical_overlap(self):
        result = subspace_overlap(np.eye(3, 1), np.eye(3, 2), np.eye(3), np.eye(3), np.eye(3))
        self.assertEqual(result["projector_frobenius_squared"], 1.)
        self.assertEqual(result["projector_spectral_distance"], 1.)
        for cross in (np.eye(3)*2, np.full((3, 3), np.nan)):
            with self.assertRaises(ValueError):
                subspace_overlap(np.eye(3), np.eye(3), np.eye(3), np.eye(3), cross)
        with self.assertRaises(ValueError):
            subspace_overlap(np.eye(3)*2, np.eye(3), np.eye(3), np.eye(3), np.eye(3))

    def test_rectangular_gaussian_integrals(self):
        def snapshot(exponents):
            shells = [Shell(np.zeros(3), np.array([i]), np.array([a]), np.array([20.]),
                np.zeros((1, 3), int), np.array([[[(2*a/np.pi)**.75]]]), 20.)
                for i, a in enumerate(exponents)]
            return SimpleNamespace(positions=np.zeros((1, 3)), cell=np.eye(3)*20,
                kinds=np.array([1]), periodic=np.zeros(3, int), kpoints=np.zeros((1, 3)),
                atom_sizes=[len(exponents)], shells=shells)
        left, right = snapshot([.7]), snapshot([.7, 1.3])
        ml, mr, cross = cross_moments(left, right, primitive_integrals)
        self.assertEqual(cross.shape, (4, 1, 2))
        np.testing.assert_allclose(cross[0], [[1., (2*np.sqrt(.7*1.3)/2.)**1.5]], atol=1e-14)
        embedding = solve(mr[0], cross[0].conj().T, assume_a="pos")
        np.testing.assert_allclose(embedding.conj().T @ mr[0] @ embedding, ml[0], atol=1e-14)
        np.testing.assert_allclose(cross[1:], 0., atol=1e-14)
        right.positions[0, 0] += 1
        with self.assertRaises(ValueError):
            cross_moments(left, right, primitive_integrals)


if __name__ == "__main__":
    unittest.main()
