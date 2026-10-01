"""Verify retained cold-build evidence without CP2K or third-party Python modules."""

import argparse
import hashlib
import json
from pathlib import PurePosixPath
import re
import tarfile


def verify(path):
    payload = {}
    with tarfile.open(path) as archive:
        for entry in archive:
            name = PurePosixPath(entry.name)
            assert entry.isfile() and not name.is_absolute() and '..' not in name.parts
            assert entry.name not in payload
            payload[entry.name] = archive.extractfile(entry).read()
    manifest = json.loads(payload.pop('manifest.json'))
    assert set(payload) == set(manifest['files'])
    for name, data in payload.items():
        record = manifest['files'][name]
        assert len(data) == record['bytes']
        assert hashlib.sha256(data).hexdigest() == record['sha256']

    def read(name):
        return json.loads(payload[name])

    old = read('original/source-provenance.json')
    new = read('corrected/source-provenance.json')
    assert set(old['files']) == set(new['files']) and len(old['files']) == 9746
    changed = sorted(name for name in old['files'] if old['files'][name] != new['files'][name])
    expected = ['src/qs_rho_atom_methods.F', 'src/quadratic_pseudospectrum.F',
                'src/quadratic_pseudospectrum_unittest.F']
    assert changed == expected == sorted(new['changed_files'])
    assert hashlib.sha256(payload['corrected/combined-fix.patch']).hexdigest() == new['patch_sha256']
    for name in expected:
        assert hashlib.sha256(payload['corrected/files/' + name]).hexdigest() == new['files'][name]['sha256']
    config = read('corrected/build-provenance.json')
    assert '-fcheck=bounds' in config['options']['CMAKE_Fortran_FLAGS'].split()
    assert config['options']['CP2K_USE_MPI'] == config['options']['CP2K_USE_MUMPS'] == 'ON'
    assert config['options']['CP2K_USE_TACHO'] == 'ON'
    assert config['source_provenance_sha256'] == hashlib.sha256(payload['corrected/source-provenance.json']).hexdigest()

    baseline, corrected = read('original/validation-result.json'), read('corrected/validation-result.json')
    assert not baseline['passed'] and corrected['passed']
    assert len(corrected['unit_runs']) == 12 and len(corrected['regression_runs']) == 2
    assert len([row for row in baseline['unit_runs'] if not row['passed']]) == 2
    units = ['spectral_localizer_unittest', 'spectral_localizer_flatten_unittest',
             'spectral_localizer_sparse_unittest', 'quadratic_pseudospectrum_unittest',
             'quadratic_dbcsr_unittest']
    expected_units = {(unit + '.psmp', ranks, 1) for unit in units for ranks in (1, 2)}
    expected_units.update({('spectral_localizer_sparse_unittest.psmp', 4, 1),
                           ('spectral_localizer_sparse_unittest.psmp', 2, 2)})
    observed_units = set()
    for row in corrected['unit_runs']:
        command = row['command']
        observed_units.add((PurePosixPath(command[command.index('-np') + 2]).name,
                            int(command[command.index('-np') + 1]),
                            int(row['environment']['OMP_NUM_THREADS'])))
    assert observed_units == expected_units
    counts = []
    for row in corrected['unit_runs'] + corrected['regression_runs'] + [corrected['version']]:
        assert row['passed'] and row['returncode'] == 0
        data = payload['corrected/' + row['log']]
        assert hashlib.sha256(data).hexdigest() == row['log_sha256']
        text = data.decode()
        assert not re.search(r'Fortran runtime error|ERROR STOP|runtime error:|RUNTIME FAIL|Status: FAILED', text)
    assert [row['log'] for row in corrected['regression_runs']] == ['regtests-2r1t.log', 'regtests-1r2t.log']
    for row, (ranks, threads) in zip(corrected['regression_runs'], [(2, 1), (1, 2)]):
        command = row['command']
        assert command[command.index('--mpiranks') + 1] == str(ranks)
        assert command[command.index('--ompthreads') + 1] == str(threads)
        assert row['environment']['OMP_NUM_THREADS'] == str(threads)
        text = payload['corrected/' + row['log']].decode()
        matched = re.search(r'Summary: correct: (\d+) / (\d+);', text)
        assert matched and matched[1] == matched[2] and int(matched[1]) > 100
        assert 'Launched 5 test directories' in text and 'Status: OK' in text
        counts.append(int(matched[1]))
    for row in baseline['regression_runs']:
        assert not row['passed'] and row['returncode'] != 0
        data = payload['original/' + row['log']]
        assert hashlib.sha256(data).hexdigest() == row['log_sha256']
        text = data.decode()
        assert 'Number of FAILED  tests 19' in text and 'Number of WRONG   tests 0' in text
    controls = read('corrected/control-result.json')
    assert controls['passed'] and len(controls['runs']) == 8
    assert hashlib.sha256(payload['corrected/fix.patch']).hexdigest() == controls['fix_patch_sha256']
    for row in controls['runs']:
        assert row['passed']
        assert (row['returncode'] == 0) == (row['version'] == 'fixed')
        name = f"corrected/controls/{row['version']}-{row['optimization']}/regression-{row['threads']}t.log"
        assert hashlib.sha256(payload[name]).hexdigest() == row['log_sha256']
        if row['version'] == 'original':
            assert b'Partial residual block did not fill the available space' in payload[name]
    pointer = read('pointer/result.json')
    assert len(pointer['runs']) == 16 and pointer['all_alias_variants_passed']
    observed_pointer = set()
    for row in pointer['runs']:
        assert (row['returncode'] == 0) == row['use_alias']
        observed_pointer.add((row['optimization'], row['openmp'], row['use_alias'], row['threads']))
        name = (f"pointer/{row['optimization']}-omp{int(row['openmp'])}"
                f"-alias{int(row['use_alias'])}-{row['threads']}t/run.log")
        assert hashlib.sha256(payload[name]).hexdigest() == row['log_sha256']
        assert row['environment']['OMP_NUM_THREADS'] == str(row['threads'])
    expected_pointer = {(opt, omp, alias, threads) for opt in ('O0', 'O3')
                        for omp in (False, True) for alias in (False, True)
                        for threads in ((1, 2, 4) if omp else (1,))}
    assert observed_pointer == expected_pointer
    assert hashlib.sha256(payload['tools/gapw_pointer_probe.F']).hexdigest() == pointer['source_sha256']
    return {'passed': True, 'files': len(payload), 'source_paths_verified': len(old['files']),
            'changed_source_files': changed, 'corrected_unit_runs': 12,
            'corrected_regression_checks': counts,
            'scope': 'Integrity and consistency of retained execution evidence, not an independent rerun or scaling claim.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('archive')
    args = parser.parse_args()
    print(json.dumps(verify(args.archive), indent=2))
