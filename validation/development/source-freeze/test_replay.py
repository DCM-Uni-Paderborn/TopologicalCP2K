"""Small negative controls for the portable MMN replay reader."""

import io
import unittest
import numpy as np

from replay import circle_error, read_loops


def sample():
    return 'example\n2 2 1\n1 2 0 0 0\n1 0\n0 2\n0 -3\n4 0\n2 1 1 0 0\n1 0\n0 0\n0 0\n1 0\n'


class ReplayTests(unittest.TestCase):
    def test_fortran_order(self):
        matrix = list(read_loops(io.StringIO(sample()), 1, 2))[0][0]
        np.testing.assert_array_equal(matrix, [[1, -3j], [2j, 4]])

    def test_bad_closure(self):
        with self.assertRaises(ValueError):
            list(read_loops(io.StringIO(sample().replace('2 1 1 0 0', '2 1 0 0 0')), 1, 2))

    def test_extra_data(self):
        with self.assertRaises(ValueError):
            list(read_loops(io.StringIO(sample() + 'trailing\n'), 1, 2))

    def test_truncated_matrix(self):
        with self.assertRaises(ValueError):
            list(read_loops(io.StringIO(sample()[:-5]), 1, 2))

    def test_circle_cut(self):
        self.assertLess(circle_error(np.array([0., .25]), np.array([.25, 1.])), 1e-14)
        self.assertGreater(circle_error(np.array([0., .25]), np.array([.25, .5])), 1.)


if __name__ == '__main__':
    unittest.main()
