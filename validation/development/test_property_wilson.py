"""Independent sensitivity controls for property-frame comparisons."""

from pathlib import Path
import tempfile
import unittest

import numpy as np

from compare_property_wilson import circular_distance, frame_errors, overlaps


class PropertyComparisonTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(582)
        x = self.rng.normal(size=(6, 6)) + 1j * self.rng.normal(size=(6, 6))
        self.metric = x.conj().T @ x + np.eye(6)
        vectors = self.rng.normal(size=(12, 4)) + 1j * self.rng.normal(size=(12, 4))
        q, _ = np.linalg.qr(vectors)
        factor = np.linalg.cholesky(self.metric)
        self.frame = np.array([np.linalg.solve(factor.conj().T, v) for v in q.reshape(2, 6, 4)])
        u = self.rng.normal(size=(4, 4)) + 1j * self.rng.normal(size=(4, 4))
        self.gauge, _ = np.linalg.qr(u)

    def test_degenerate_spinor_gauge(self):
        sewing, errors = frame_errors(self.metric, self.frame, self.frame @ self.gauge)
        np.testing.assert_allclose(sewing, self.gauge, atol=3e-14)
        self.assertLess(max(errors.values()), 3e-14)

    def test_wrong_subspace(self):
        wrong = self.frame.copy()
        wrong[0, 0, 0] += 0.03j
        _, errors = frame_errors(self.metric, self.frame, wrong)
        self.assertGreater(errors["subspace_residual"], 0.01)

    def test_wrong_metric(self):
        _, errors = frame_errors(np.eye(6), self.frame, self.frame @ self.gauge)
        self.assertGreater(errors["frame_norm_error"], 0.1)

    def test_nonfinite(self):
        wrong = self.frame.copy()
        wrong[0, 0, 0] = np.nan
        with self.assertRaises(ValueError):
            frame_errors(self.metric, self.frame, wrong)

    def test_directed_link_covariance(self):
        matrix = self.rng.normal(size=(4, 4)) + 1j * self.rng.normal(size=(4, 4))
        left = self.gauge
        right = np.roll(self.gauge, 1, axis=1)
        transformed = left.conj().T @ matrix @ right
        np.testing.assert_allclose(left @ transformed @ right.conj().T, matrix, atol=2e-14)
        self.assertGreater(np.max(np.abs(transformed - left.conj().T @ matrix @ right.conj())), 0.1)

    def test_branch_cut_and_permutation(self):
        a = np.array([0.0, 0.2, 0.5, 0.5])
        b = np.array([0.5, 1.0 - 1e-12, 0.5, 0.2])
        self.assertLess(circular_distance(a, b), 2e-12)
        self.assertGreater(circular_distance(a, b + 0.07), 0.05)

    def test_overlap_parser(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / "states.mmn"
            p.write_text("test\n2 1 1\n1 1 1 0 0\n1 2\n3 4\n5 6\n7 8\n")
            matrix = overlaps(p)[1, 1, 1, 0, 0]
            np.testing.assert_array_equal(matrix, [[1 + 2j, 5 + 6j], [3 + 4j, 7 + 8j]])
            p.write_text(p.read_text().replace("1 2\n", "nan 2\n"))
            with self.assertRaises(ValueError):
                overlaps(p)


if __name__ == "__main__":
    unittest.main()
