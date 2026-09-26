"""Independent matrix, subdivision and rejection checks for localizer boxes."""

from itertools import product
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.linalg import eigvalsh, svdvals

from localizer_region import cover_region, variation_bound
from archive_bismuth_regions import case_archive, compare


class RegionChecks(unittest.TestCase):
    def test_joint_operator_bound(self):
        rng = np.random.default_rng(4427)
        n = 5
        xyz = rng.normal(size=(2, n, n)) + 1j * rng.normal(size=(2, n, n))
        xyz += xyz.conj().swapaxes(1, 2)
        origin = np.array([.3, -.8])
        eye = np.eye(n)
        norm = svdvals(xyz[0] - origin[0]*eye - 1j*(xyz[1] - origin[1]*eye))[0]

        def matrix(p):
            e, k, x, y = p
            b = k * (xyz[0] - x*eye - 1j*(xyz[1] - y*eye))
            return np.block([[-e*eye, b], [b.conj().T, e*eye]])

        for _ in range(12):
            center = np.r_[rng.normal(), rng.uniform(.2, .5), rng.normal(size=2)]
            half = rng.uniform(.01, .1, 4)
            lower, upper = center-half, center+half
            bound, _ = variation_bound(lower, upper, norm, origin)
            base = matrix(center)
            for signs in product((-1, 1), repeat=4):
                point = center+half*np.array(signs)
                actual = max(abs(eigvalsh(matrix(point)-base)))
                self.assertLessEqual(actual, bound*(1+1e-12))

    def test_exact_scalar_changes(self):
        # Position-only and energy-only bounds attain equality.
        self.assertAlmostEqual(variation_bound([-.2, 2, 0, 0], [.2, 2, 0, 0], 7, [0, 0])[0], .2)
        self.assertAlmostEqual(variation_bound([0, 2, -.3, -.4], [0, 2, .3, .4], 7, [0, 0])[0], 1.)

    def test_connected_cover(self):
        lo, hi = [-.1, .1, -.2, -.1], [.1, .3, .2, .1]
        result = cover_region(lo, hi, 2, [0, 0], lambda _: .2)
        self.assertTrue(result['resolved'])
        self.assertGreater(len(result['covered']), 1)
        volume = sum(np.prod(np.subtract(r['upper'], r['lower'])) for r in result['covered'])
        self.assertAlmostEqual(volume, np.prod(np.subtract(hi, lo)))
        self.assertTrue(all(r['lower_bound'] > 0 for r in result['covered']))

    def test_complete_partition_including_failed_leaves(self):
        lo, hi = np.array([-.1, .1, -.2, -.1]), np.array([.1, .3, .2, .1])
        result = cover_region(lo, hi, 2, [0, 0], lambda p: abs(p[0]), max_nodes=11)
        leaves = result['covered'] + result['unresolved']
        volume = sum(np.prod(np.subtract(row['upper'], row['lower'])) for row in leaves)
        self.assertAlmostEqual(volume, np.prod(hi-lo))
        self.assertFalse(result['resolved'])
        for i, a in enumerate(leaves):
            for b in leaves[i+1:]:
                intersection = np.minimum(a['upper'], b['upper']) - np.maximum(a['lower'], b['lower'])
                self.assertFalse(np.all(intersection > 0))

    def test_gapped_matrix_cover(self):
        # A commuting two-level Hamiltonian has exact localizer eigenvalues
        # +/-sqrt((h-E)^2+kappa^2*((X-x)^2+(Y-y)^2)).
        h, x, y = np.array([-.7, .9]), np.array([-1., 2.]), np.array([.5, -.4])
        rho = max(np.hypot(x, y))
        def analytic_gap(p):
            e, k, px, py = p
            return float(min(np.sqrt((h-e)**2+k**2*((x-px)**2+(y-py)**2))))
        result = cover_region([-.2, .1, -.5, -.5], [.2, .4, .5, .5], rho, [0, 0], analytic_gap)
        self.assertTrue(result['resolved'])
        for row in result['covered']:
            for point in product(*zip(row['lower'], row['upper'])):
                self.assertGreaterEqual(analytic_gap(point), row['lower_bound'])

    def test_closing_and_budget(self):
        for budget in (1, 63):
            result = cover_region([-.1, .1, 0, 0], [.1, .1, 0, 0], 1, [0, 0],
                                  lambda p: abs(p[0]), max_nodes=budget)
            self.assertFalse(result['resolved'])
            self.assertTrue(result['unresolved'])
            self.assertLessEqual(len(result['evaluations']), budget)
        point = cover_region([0, .1, 0, 0], [0, .1, 0, 0], 1, [0, 0], lambda _: 0)
        self.assertFalse(point['resolved'])

    def test_origin_translation(self):
        lo, hi = np.array([0, .1, -1, -2]), np.array([.1, .2, 2, 1])
        translation = np.array([0, 0, 31.3, -21.9])
        a = variation_bound(lo, hi, 3.4, [.3, -.2])[0]
        b = variation_bound(lo+translation, hi+translation, 3.4, translation[2:]+[.3, -.2])[0]
        self.assertAlmostEqual(a, b)

    def test_rejections(self):
        args = ([0, .1, 0, 0], [.1, .2, 1, 1], 1, [0, 0])
        for bad in (-1, float('nan')):
            with self.assertRaises(ValueError):
                variation_bound(args[0], args[1], bad, args[3])
        with self.assertRaises(ValueError):
            cover_region(*args, lambda _: float('nan'))
        with self.assertRaises(ValueError):
            cover_region(*args, lambda _: 1, max_nodes=0)
        with self.assertRaises(ValueError):
            variation_bound(args[1], args[0], 1, [0, 0])

    def test_replay_discrete_fields_are_exact(self):
        expected = dict(resolved=True, index=1, gap=0.3, samples=[1, 2], elapsed_seconds=10.)
        for changed in (dict(resolved=False), dict(index=0), dict(index=1.), dict(samples=[1]),
                        dict(gap=float('nan')), dict(gap=.31)):
            with self.assertRaises(ValueError):
                compare(expected, {**expected, **changed})
        self.assertLess(compare(expected, {**expected, 'gap': .3 + 1e-12,
                                            'elapsed_seconds': 42.}), 2e-12)

    def test_replay_catches_partition_change(self):
        expected = dict(covered=[dict(lower=[0., 1.], upper=[.5, 2.])], unresolved=[])
        with self.assertRaises(ValueError):
            compare(expected, dict(covered=[dict(lower=[.1, 1.], upper=[.5, 2.])], unresolved=[]))
        with self.assertRaises(ValueError):
            compare(expected, dict(covered=expected['covered'], unresolved=[[0., .1]]))

    def test_case_archive_stays_in_evidence_tree(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            self.assertEqual(case_archive(root / 'new-cases', 'flake5', root),
                             root / 'new-cases/flake5.tar.gz')
            for name in ('', '.', '..', '../flake5', '/flake5'):
                with self.assertRaises(ValueError):
                    case_archive(root / 'new-cases', name, root)
            with self.assertRaises(ValueError):
                case_archive(root.parent, 'flake5', root)
            outside = root / 'outside'
            outside.symlink_to(root.parent, target_is_directory=True)
            with self.assertRaises(ValueError):
                case_archive(outside, 'flake5', root)


if __name__ == '__main__':
    unittest.main()
