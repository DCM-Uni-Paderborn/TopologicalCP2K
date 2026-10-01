"""Compare the retained Wilson evidence with independently rebuilt CP2K runs."""

import argparse
import hashlib
import json
import logging
from pathlib import Path
import re
import tempfile

import numpy as np

from replay import circle_error, native_check, require, unpack_checked


TOLERANCES = {'total_energy_hartree': 1e-8, 'eigenvalues_eV': 1e-6,
              'wilson_eigenvalue_chord': 1e-8, 'link_singular_value': 1e-8}


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def total_energy(path):
    matches = re.findall(r'ENERGY\| Total FORCE_EVAL .*?energy \[hartree\]\s+([-+0-9.Ee]+)',
                         path.read_text())
    require(len(matches) == 1, 'Expected exactly one total energy')
    energy = float(matches[0])
    require(np.isfinite(energy), 'Nonfinite energy')
    return energy


def differences(reference, fresh, name):
    require((reference / 'input.inp').read_bytes() == (fresh / 'input.inp').read_bytes(),
            'Input differs from retained calculation')
    old_e = np.loadtxt(reference / f'{name}.eig', ndmin=2)
    new_e = np.loadtxt(fresh / f'{name}.eig', ndmin=2)
    require(old_e.shape == new_e.shape and old_e.shape[1] == 3 and
            np.isfinite(old_e).all() and np.isfinite(new_e).all() and
            np.array_equal(old_e[:, :2], new_e[:, :2]), 'Invalid eigenvalue sampling')
    old_w = np.loadtxt(reference / f'{name}.wilson', ndmin=2)
    new_w = np.loadtxt(fresh / f'{name}.wilson', ndmin=2)
    require(old_w.shape == new_w.shape and old_w.shape[1] == 10 and
            np.isfinite(old_w).all() and np.isfinite(new_w).all() and
            np.array_equal(old_w[:, 0], new_w[:, 0]), 'Invalid Wilson sampling')
    return {
        'total_energy_hartree': abs(total_energy(fresh / 'run.log') -
                                    total_energy(reference / 'run.log')),
        'eigenvalues_eV': float(np.max(abs(new_e[:, 2] - old_e[:, 2]))),
        'wilson_eigenvalue_chord': max(circle_error(a, b) for a, b in
            zip(old_w[:, 2:], new_w[:, 2:], strict=True)),
        'link_singular_value': float(np.max(abs(new_w[:, 1] - old_w[:, 1]))),
    }


def check_provenance(root, manifest):
    record = json.loads((root / 'run-provenance.json').read_text())
    source = json.loads((root / 'source-reproduction/source-manifest.json').read_text())
    require(record['source_and_binary_stable'] and record['build']['exit_code'] == 0 and
            record['unit_test']['exit_code'] == 0, 'Incomplete source/build validation')
    require(source['target_commit'] == record['source']['commit'] and
            source['target_tree'] == source['reconstructed_tree'] == record['source']['tree'] and
            source['exact_tree_reconstruction'], 'Inconsistent source identity')
    require(digest(root / 'source-reproduction/source.patch') == source['patch_sha256'] and
            digest(root / 'source-reproduction/CP2K-LICENSE') == source['license_sha256'],
            'Source patch or license differs')
    require('Wilson-loop and Z2 unit tests passed.' in (root / 'wilson-unit-test.log').read_text(),
            'Missing native unit-test completion')
    require(manifest['files']['CMakeCache.txt'] == record['cmake_cache'], 'Build configuration differs')
    for name in ('neon', 'stanene'):
        run = record['runs'][name]
        require(run['exit_code'] == 0, 'Failed native run')
        expected = dict(run['outputs'], **{'input.inp': run['input']})
        for filename, fingerprint in expected.items():
            require(manifest['files'][f'native/{name}/{filename}'] == fingerprint,
                    'Native output differs from recorded provenance')
        log = (root / f'native/{name}/run.log').read_text()
        threads = run['runtime_settings']['OMP_NUM_THREADS']
        require(re.search(r'Number of threads for this process\s+' + threads + r'\s', log)
                is not None, 'Thread count differs from recorded run settings')
        require(re.search(r'CP2K\| source code revision number:\s+' +
                          re.escape(source['target_commit'][:7]) + r'\s', log) is not None,
                'Native log revision differs from frozen source')
    return record['source']


def control_checks(reference, fresh):
    """Replay the failed initial configuration and the scalar-space controls too."""
    report = {'cases': {}}
    controls = json.loads((fresh / 'diagnostics/thread-controls.json').read_text())
    with tempfile.TemporaryDirectory(prefix='wilson-control-view-') as temporary:
        view = Path(temporary)
        (view / 'native').mkdir()
        link = view / 'native/neon'

        def replay_case(name, directory):
            link.symlink_to(directory, target_is_directory=True)
            report['cases'][name] = native_check(view, 'neon', 0)
            link.unlink()

        for name, record in controls['cases'].items():
            directory = fresh / 'diagnostics' / name
            delta = differences(reference / 'native/neon', directory, 'neon')
            for key, expected in record['reference_differences'].items():
                require(np.isclose(delta[key], expected, rtol=1e-5, atol=1e-14),
                        'Thread control disagrees with retained diagnostics')
            replay_case(name, directory)
            report['cases'][name]['historical_differences'] = delta
        replay_case('initial-omp2', fresh / 'diagnostics/initial-neon-omp2')
        for group, label in (('neon-scalar-space-controls', 'soc-full'),
                             ('neon-triplet-space-controls', 'soc-triplet')):
            root = fresh / 'diagnostics' / group
            record = json.loads((root / 'comparison.json').read_text())
            base = root / f'{label}-omp1-repeat1'
            for name in record['full_space_differences']:
                delta = differences(base, root / name, 'neon')
                require(all(delta[key] <= tol for key, tol in TOLERANCES.items()),
                        'Complete scalar-space thread control is not stable')
                replay_case(f'{group}/{name}', root / name)
                report['cases'][f'{group}/{name}']['same_space_differences'] = delta
    scalar_path = fresh / 'diagnostics/neon-scalar-space-controls/scalar-full/neon.eig'
    spectrum = np.loadtxt(scalar_path)
    require(spectrum.shape == (40 * 13, 3) and np.isfinite(spectrum).all(), 'Invalid scalar control')
    require(np.array_equal(spectrum[:, 0], np.tile(np.arange(1, 14), 40)) and
            np.array_equal(spectrum[:, 1], np.repeat(np.arange(1, 41), 13)), 'Scalar point indices differ')
    energy = spectrum[:, 2].reshape(40, 13)
    report['minimum_scalar_cutoff_gaps_eV'] = {
        '5': float(np.min(energy[:, 5] - energy[:, 4])),
        '7': float(np.min(energy[:, 7] - energy[:, 6]))}
    recorded = json.loads((fresh / 'diagnostics/scalar-space-diagnosis.json').read_text())
    for key, value in report['minimum_scalar_cutoff_gaps_eV'].items():
        require(np.isclose(value, recorded['minimum_scalar_cutoff_gaps_eV'][key], rtol=1e-10, atol=1e-14),
                'Scalar cutoff-gap diagnosis differs')
    report['scope'] = ('Historical five-state thread sensitivity retained, not hidden or retoleranced. '
                       'Seven- and thirteen-state controls are separate changed-input diagnostics.')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference_archive', type=Path)
    parser.add_argument('fresh_archive', type=Path)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    logging.getLogger('z2pack').setLevel(logging.ERROR)
    result = {'reference_sha256': digest(args.reference_archive),
              'fresh_sha256': digest(args.fresh_archive), 'tolerances': TOLERANCES,
              'materials': {}, 'new_dft_calculations_during_replay': 0}
    with tempfile.TemporaryDirectory(prefix='wilson-frozen-check-') as temporary:
        reference, fresh = Path(temporary) / 'reference', Path(temporary) / 'fresh'
        previous_manifest = unpack_checked(args.reference_archive, reference)
        manifest = unpack_checked(args.fresh_archive, fresh)
        result['source'] = check_provenance(fresh, manifest)
        result['files_verified'] = len(previous_manifest['files']) + len(manifest['files'])
        for name, z2 in (('neon', 0), ('stanene', 1)):
            for folder in (reference, fresh):
                log = (folder / f'native/{name}/run.log').read_text()
                expected_threads = '1' if name == 'neon' else '2'
                require(re.search(r'Number of threads for this process\s+' + expected_threads + r'\s', log)
                        is not None, 'Thread count differs from retained material settings')
            delta = differences(reference / 'native' / name, fresh / 'native' / name, name)
            require(all(delta[key] <= tol for key, tol in TOLERANCES.items()),
                    f'{name}: fresh calculation differs from retained evidence: {delta}')
            result['materials'][name] = {'maximum_differences': delta,
                'reference_replay': native_check(reference, name, z2),
                'fresh_replay': native_check(fresh, name, z2)}
            print(name, json.dumps(result['materials'][name]), flush=True)
        result['diagnostics'] = control_checks(reference, fresh)
    require(result['reference_sha256'] == digest(args.reference_archive) and
            result['fresh_sha256'] == digest(args.fresh_archive), 'Archive changed during replay')
    result['all_passed'] = True
    result['scope'] = ('Archive-only replay of fresh pinned-build runs versus retained finite-basis inputs. '
                       'No material-convergence claim and no recovery of the historical executable. '
                       'Build provenance is retained evidence, not a new build performed by this verifier.')
    args.report.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
