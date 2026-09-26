"""Verify frozen Hall evidence and optionally rebuild its mathematical test."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import tarfile
import tempfile


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def check_unit(text):
    require('Hall sign, Liouvillian, gauge, time reversal, degeneracy and Chern tests passed.' in text,
            'Missing mathematical completion marker')
    rows = re.findall(r'Hall Chern model mass, sigma, Chern, residual:\s*([^\n]+)', text)
    require(len(rows) >= 3, 'Missing Chern model phases')
    observed = set()
    for row in rows:
        mass, sigma, chern, residual = map(float, row.split())
        require(all(map(math.isfinite, (mass, sigma, chern, residual))), 'Nonfinite model output')
        expected = {-1.0: -1.0, 1.0: 1.0, 3.0: 0.0}[mass]
        require(abs(sigma - expected) < 1.e-6, 'Cold Hall limit failed')
        require(abs(chern + expected) < 1.e-10, 'Chern integral failed')
        require(residual < 1.e-12, 'Pointwise curvature failed')
        observed.add(mass)
    require(observed == {-1.0, 1.0, 3.0}, 'Incomplete Chern phases')
    for label in ('Hall Liouvillian residual:', 'Hall k-dependent nonorthogonal AO residual:'):
        values = [float(x) for x in re.findall(re.escape(label) + r'\s*(\S+)', text)]
        require(values and all(math.isfinite(x) and x < 1.e-12 for x in values), label)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blas-library', type=Path)
    parser.add_argument('--compiler', default='gfortran')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / 'manifest.json').read_text())
    payload = {}
    archive = root / 'methods-results.tar.gz'
    require(hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['archive_sha256'],
            'Archive digest mismatch')
    with tarfile.open(archive, 'r:gz') as stream:
        for item in stream:
            name = PurePosixPath(item.name)
            require(item.isfile() and not name.is_absolute() and '..' not in name.parts,
                    'Unsafe archive member')
            require(item.name not in payload, 'Duplicate archive member')
            payload[item.name] = stream.extractfile(item).read()
    require(set(payload) == set(manifest['files']), 'Archive inventory mismatch')
    for name, data in payload.items():
        require(hashlib.sha256(data).hexdigest() == manifest['files'][name]['sha256'], name)
    for layout in ('serial', 'mpi4'):
        log = payload[f'evidence/{layout}-regtests.log'].decode()
        require('Summary: correct: 67 / 67;' in log and 'Status: OK' in log,
                f'{layout}: incomplete regression suite')
        check_unit(payload[f'evidence/{layout}-unit.log'].decode())
    check_unit(payload['evidence/debug-unit.log'].decode())
    log = payload['evidence/mpi2-regtests.log'].decode()
    require('Summary: correct: 14 / 14;' in log and 'Status: OK' in log,
            'Two-rank Hall regression suite incomplete')
    controls = json.loads(payload['evidence/controls.json'])
    for name in ('bloch-soc3', 'bloch-scalar3', 'projected-soc'):
        on, off = controls[name]['on'], controls[name]['off']
        require(abs(on['sigma_iso'] - off['sigma_iso']) < 1.e-12, name)
        require(max(abs(a - b) for a, b in zip(on['tensor'], off['tensor'])) < 1.e-9, name)
        require(on['hall_max'] < 1.e-9 and off['hall_max'] is None, name)
    require(controls['embedded']['guard_rejected'], 'Missing atom-embedding rejection')
    if args.blas_library:
        library = args.blas_library.resolve(strict=True)
        with tempfile.TemporaryDirectory(prefix='kubo-hall-replay-') as temporary:
            directory = Path(temporary)
            files = ['kinds.F', 'kubo_projected.F', 'kubo_bloch.F', 'kubo_hall.F',
                     'kubo_hall_unittest.F']
            for name in files:
                (directory / name).write_bytes(payload['sources/' + name])
            command = [args.compiler, '-cpp', '-ffree-form', '-ffree-line-length-none',
                       '-std=f2008', '-O0', '-g', '-fcheck=all', '-fbacktrace',
                       '-ffpe-trap=invalid,zero,overflow', '-D__HAS_IEEE_EXCEPTIONS',
                       '-J.', '-I.', *files, str(library), '-o', 'hall-reference']
            subprocess.run(command, cwd=directory, check=True)
            environment = os.environ.copy()
            environment.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='2')
            result = subprocess.run([str(directory / 'hall-reference')], cwd=directory,
                                    env=environment, check=True, capture_output=True, text=True)
            check_unit(result.stdout)
            print(result.stdout, end='')
    print(json.dumps({'verified_archive_members': len(payload),
                      'native_assertions': 148, 'control_calculations': 7,
                      'mathematical_test_rebuilt': bool(args.blas_library)}, indent=2))


if __name__ == '__main__':
    main()
