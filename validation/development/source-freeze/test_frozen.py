"""Negative controls for numerical frozen-run comparison."""

from pathlib import Path
import tempfile
import unittest

from verify_frozen import differences, total_energy


class FrozenTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='wilson-check-test-')
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.old, self.new = self.root / 'old', self.root / 'new'
        for directory in (self.old, self.new):
            directory.mkdir()
            (directory / 'input.inp').write_text('same input\n')
            (directory / 'run.log').write_text('ENERGY| Total FORCE_EVAL ( QS ) energy [hartree] -1.5\n')
            (directory / 'neon.eig').write_text('1 1 -0.3\n2 1 0.2\n')
            (directory / 'neon.wilson').write_text('1 0.99 0 0 .2 .2 .4 .4 .6 .6\n')

    def test_unchanged(self):
        self.assertTrue(all(value == 0 for value in differences(self.old, self.new, 'neon').values()))

    def test_changed_eigenvalue(self):
        (self.new / 'neon.eig').write_text('1 1 -0.2\n2 1 0.2\n')
        self.assertAlmostEqual(differences(self.old, self.new, 'neon')['eigenvalues_eV'], 0.1)

    def test_changed_input(self):
        (self.new / 'input.inp').write_text('different input\n')
        with self.assertRaises(ValueError):
            differences(self.old, self.new, 'neon')

    def test_changed_indices(self):
        (self.new / 'neon.eig').write_text('1 2 -0.3\n2 2 0.2\n')
        with self.assertRaises(ValueError):
            differences(self.old, self.new, 'neon')

    def test_nonfinite_spectrum(self):
        (self.new / 'neon.eig').write_text('1 1 nan\n2 1 0.2\n')
        with self.assertRaises(ValueError):
            differences(self.old, self.new, 'neon')

    def test_duplicate_energy(self):
        log = self.new / 'run.log'
        log.write_text(log.read_text() * 2)
        with self.assertRaises(ValueError):
            total_energy(log)


if __name__ == '__main__':
    unittest.main()
