"""Independent checks of flake geometry, metric integration and full-band gauges."""

from contextlib import redirect_stdout
import io
from itertools import product
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
from scipy.linalg import block_diag, eigh, eigvalsh

from analyze_bismuth_flakes import moments, native_comparison, operators, queries
from archive_bismuth_validation import digest, retain_case
from diagnose_bismuth_operators import matrix, soc_components
from certify_bismuth_scale_window import cover_gap
from run_bismuth_flakes import geometry, inputs, nnkp, process_tree_rss
from replay_bismuth_neutral_window import compare_interval
from scan_bismuth_neutral_window import frontier_localization, window
from verify_soc_print_export import headerless

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
    def test_interval_replay(self):
        expected = dict(case="test", interval=[.1, .2], inputs={}, methods={}, script_sha256="method",
                        reference_scan_sha256="scan", anchor_kappa=.15, anchor_z2=1,
                        lipschitz_bohr=2., minimum_gap_lower_bound_hartree=.1,
                        bound=cover_gap(.1, .2, 2., lambda x: .3))
        actual = json.loads(json.dumps(expected))
        self.assertEqual(compare_interval(expected, actual)["maximum_absolute_difference"], 0.)
        actual["anchor_z2"] = 0
        with self.assertRaises(ValueError):
            compare_interval(expected, actual)
        actual = json.loads(json.dumps(expected))
        actual["bound"]["covered"][0]["gap"] += 1e-3
        with self.assertRaises(ValueError):
            compare_interval(expected, actual)
        actual = json.loads(json.dumps(expected))
        actual["bound"]["resolved"] = False
        with self.assertRaises(ValueError):
            compare_interval(expected, actual)
        actual = json.loads(json.dumps(expected))
        actual["bound"]["covered"][0]["lower_bound"] = float("nan")
        with self.assertRaises(ValueError):
            compare_interval(expected, actual)

    def test_lipschitz_gap_cover(self):
        covered = cover_gap(.1, 1., 2., lambda x: .1)
        self.assertTrue(covered["resolved"])
        intervals = covered["covered"]
        self.assertEqual(intervals[0]["lower"], .1)
        self.assertEqual(intervals[-1]["upper"], 1.)
        for first, second in zip(intervals, intervals[1:]):
            self.assertEqual(first["upper"], second["lower"])
        self.assertTrue(all(row["lower_bound"] > 0 for row in intervals))
        closing = cover_gap(.1, .9, 1., lambda x: abs(x-.5), max_nodes=127)
        self.assertFalse(closing["resolved"])
        self.assertTrue(closing["unresolved"])
        with self.assertRaises(ValueError):
            cover_gap(.1, .9, float("nan"), lambda x: 1.)

    def test_trimmed_termination(self):
        for size in (2, 3, 4):
            cell, original = geometry(size, 20.)
            trimmed_cell, trimmed = geometry(size, 20., trim=True)
            self.assertEqual(cell, trimmed_cell)
            self.assertEqual(len(trimmed), len(original) - 2)
            self.assertTrue(all(p in original for p in trimmed))
            np.testing.assert_allclose(np.mean(original, axis=0), np.mean(trimmed, axis=0), atol=1e-14)
            distances = np.linalg.norm(np.asarray(trimmed)[:, None] - trimmed, axis=2)
            coordination = np.sum((distances > 1e-8) & (distances < 3.3), axis=1)
            self.assertGreaterEqual(min(coordination), 2)
            self.assertLessEqual(max(coordination), 3)
        for size in (0, 1):
            with self.assertRaises(ValueError):
                geometry(size, 20., trim=True)

    def test_frontier_subspace_weights(self):
        positions = np.array([[0., 0., 0.], [3., 0., 0.],
                              [-1.5, 1.5*np.sqrt(3), 0.], [-1.5, -1.5*np.sqrt(3), 0.]]) / .52917720859
        rng = np.random.default_rng(1934)
        frame = np.linalg.qr(rng.normal(size=(24, 24)) + 1j*rng.normal(size=(24, 24)))[0]
        snapshot = SimpleNamespace(positions=positions, atom_sizes=[3]*4,
            coefficients=frame.reshape(1, 2, 12, 24), energies=np.repeat(np.arange(12.), 2)[None, :])
        first = frontier_localization(snapshot, np.eye(12))
        self.assertEqual(first["coordination"], [3, 1, 1, 1])
        self.assertEqual(first["edge_atoms_one_based"], [2, 3, 4])
        self.assertEqual(first["frontier"][0]["bands_one_based"], [19, 20])
        original_frame = frame.copy()
        for start in (18, 20):
            rotation = np.linalg.qr(rng.normal(size=(2, 2)) + 1j*rng.normal(size=(2, 2)))[0]
            frame[:, start:start+2] = frame[:, start:start+2] @ rotation
        self.assertGreater(np.max(abs(frame - original_frame)), .1)
        snapshot.coefficients = frame.reshape(1, 2, 12, 24)
        changed = frontier_localization(snapshot, np.eye(12))
        for original, rotated in zip(first["frontier"], changed["frontier"]):
            np.testing.assert_allclose(original["atom_weights"], rotated["atom_weights"], atol=1e-14)
            self.assertAlmostEqual(sum(rotated["atom_weights"]), 1.)
        trial = rng.normal(size=(12, 12))
        metric = trial @ trial.T + np.eye(12)
        eigenvalues, vectors = eigh(metric)
        inverse_root = (vectors / np.sqrt(eigenvalues)) @ vectors.T
        snapshot.coefficients = (block_diag(inverse_root, inverse_root) @ frame).reshape(1, 2, 12, 24)
        covariant = frontier_localization(snapshot, metric)
        for original, restored in zip(first["frontier"], covariant["frontier"]):
            np.testing.assert_allclose(original["atom_weights"], restored["atom_weights"], atol=1e-14)

    def test_process_tree_memory(self):
        snapshot = "100 1 10\n101 100 20\n102 101 30\n103 100 40\n999 1 5000\n"
        self.assertEqual(process_tree_rss(snapshot, 100), (100, 40, 4))
        self.assertEqual(process_tree_rss(snapshot, 101), (50, 30, 2))
        self.assertEqual(process_tree_rss(snapshot, 200), (0, 0, 0))
        self.assertEqual(process_tree_rss("", 100), (0, 0, 0))

    def test_neutral_window(self):
        snapshot = SimpleNamespace(energies=np.arange(12.)[None, :],
                                   positions=np.array([[0., 1., 2.], [4., 3., 6.]]))
        energies, points, record = window(snapshot, [-.25, 0., .25], [0., 1.5])
        np.testing.assert_allclose(energies, [9.25, 9.5, 9.75])
        np.testing.assert_allclose(points, [[2., 2., 4.], [5., 2., 4.]])
        self.assertEqual(record["gap_hartree"], 1.)
        snapshot.positions += 3
        snapshot.energies += 2
        shifted_e, shifted_p, _ = window(snapshot, [-.25, 0., .25], [0., 1.5])
        np.testing.assert_allclose(shifted_e, energies + 2)
        np.testing.assert_allclose(shifted_p, np.asarray(points) + 3)
        for fractions in ([], [.5], [float("nan")]):
            with self.assertRaises(ValueError):
                window(snapshot, fractions, [0.])
        snapshot.energies[0, 10] = snapshot.energies[0, 9]
        with self.assertRaises(ValueError):
            window(snapshot, [0.], [0.])

    def test_headerless_matrix(self):
        text = "MATRIX\n1.0 2.0\n4.0 5.0\n7.0 8.0\n\n3.0\n6.0\n9.0\n"
        np.testing.assert_allclose(headerless(text, "MATRIX", 3), np.arange(1., 10.).reshape(3, 3))
        with self.assertRaises(ValueError):
            headerless(text.replace("9.0", ""), "MATRIX", 3)

    def test_ao_matrix_reader(self):
        text = "\nOVERLAP MATRIX\n\n1 2\n1 1 Bi 2s 1.0 0.2\n2 1 Bi 3s 0.2 1.0\n"
        np.testing.assert_allclose(matrix(text, "OVERLAP MATRIX", 2), [[1., .2], [.2, 1.]])
        with self.assertRaises(ValueError):
            matrix(text.replace("2 1 Bi 3s 0.2 1.0", ""), "OVERLAP MATRIX", 2)
        with self.assertRaises(ValueError):
            matrix(text, "KOHN-SHAM MATRIX", 2)
        skew = "\nSOC\n\n1 2\n1 1 Bi 2s 0.0 0.2\n2 1 Bi 3s -0.2 0.0\n"
        np.testing.assert_allclose(matrix(skew, "SOC", 2, antisymmetric=True), [[0., .2], [-.2, 0.]])
        with self.assertRaises(ValueError):
            matrix(skew, "SOC", 2)

    def test_soc_spinor_components(self):
        rng = np.random.default_rng(503)
        scalar = rng.normal(size=(5, 5))
        scalar += scalar.T
        v = rng.normal(size=(3, 5, 5))
        v -= v.transpose(0, 2, 1)
        h = np.block([[scalar+1j*v[2], 1j*v[0]-v[1]],
                      [1j*v[0]+v[1], scalar-1j*v[2]]])
        actual_h, actual_v = soc_components(h)
        np.testing.assert_allclose(actual_h, scalar)
        np.testing.assert_allclose(actual_v, v)
        with self.assertRaises(ValueError):
            soc_components(h+np.diag([1.]*5+[-1.]*5))

    def test_retained_case_after_pruning(self):
        with tempfile.TemporaryDirectory(prefix="bismuth-archive-test-") as temporary:
            root = Path(temporary)
            work = root / "periodic-test"
            work.mkdir()
            (work / "runner.py").write_text("# retained test driver\n")
            (work / "input.inp").write_text("test input\n")
            (work / "output.out").write_text(
                "*** SCF run converged\nWilson surface sampling converged.\n"
                "Converged Z2 invariant: 1\nSampled indirect gap [eV]: 0.4\nPROGRAM ENDED AT\n")
            run = dict(completed=True, scf_converged=True, returncode=0,
                       options=dict(mode="wilson"),
                       provenance=dict(runner_sha256=digest(work / "runner.py")),
                       files={p.name: digest(p) for p in work.iterdir()})
            (work / "run.json").write_text(json.dumps(run))
            archive = root / "periodic-test.tar.gz"
            first = retain_case(work, archive)
            (work / "input.inp").unlink()
            self.assertEqual(retain_case(work, archive), first)
            self.assertFalse((work / "input.inp").exists())
            (work / "runner.py").write_text("# modified driver\n")
            with self.assertRaises(ValueError):
                retain_case(work, archive)

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
        case.update(mode="localizer", export_spectrum=True, ao_matrices=True, energy=[-.1], kappa=[.003],
                    offset=[0.], solver="TACHO", checkpoint=True, restart=Path("start.wfn"))
        files = inputs(SimpleNamespace(**case))
        for expected in ("&SPECTRAL_LOCALIZER", "&WANNIER90", "SCF_GUESS RESTART", "QS_SCF 5"):
            self.assertIn(expected, files["input.inp"])
        self.assertIn("gamma.nnkp", files)
        self.assertIn("        SOC T\n        NDIGITS 16", files["input.inp"])
        stack, matrix_paths = [], []
        for line in files["input.inp"].splitlines():
            token = line.strip().split()
            if not token:
                continue
            if token[0] == "&END":
                stack.pop()
            elif token[0].startswith("&"):
                stack.append(token[0][1:])
                if stack[-1] == "AO_MATRICES":
                    matrix_paths.append(tuple(stack))
        self.assertEqual(matrix_paths, [("FORCE_EVAL", "DFT", "PRINT", "AO_MATRICES")])
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
