"""Replay archived CP2K/Wilson and adaptive Z2Pack evidence without running DFT."""

import argparse
import hashlib
import io
import json
import logging
from importlib.metadata import version
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import tempfile

import numpy as np
import z2pack


def require(condition, message):
    if not condition:
        raise ValueError(message)


def circle_error(left, right):
    left, right = np.sort(left), np.sort(right)
    require(left.shape == right.shape, 'Different Wilson spectrum sizes')
    return min(float(np.max(abs(np.exp(2j * np.pi * left) -
                               np.exp(2j * np.pi * np.roll(right, shift)))))
               for shift in range(len(left)))


def read_loops(stream, nloops, nbands=8):
    """Validate directed MMN connections, including each reciprocal closure."""
    require(bool(stream.readline()), 'Missing MMN comment')
    bands, nk, neighbours = map(int, stream.readline().split())
    require(bands == nbands and nk > 0 and nk % nloops == 0 and neighbours == 1,
            'Unexpected MMN dimensions')
    per_loop = nk // nloops
    require(per_loop >= 2, 'A loop needs at least two links')
    for loop in range(nloops):
        matrices = []
        for link in range(per_loop):
            source = loop * per_loop + link + 1
            last = link == per_loop - 1
            target = loop * per_loop + 1 if last else source + 1
            expected = (source, target, int(last), 0, 0)
            require(tuple(map(int, stream.readline().split())) == expected,
                    f'Invalid MMN connection at source {source}')
            values = []
            for _ in range(bands * bands):
                real, imag = map(float, stream.readline().replace('D', 'E').split())
                values.append(complex(real, imag))
            matrix = np.array(values).reshape((bands, bands), order='F')
            require(np.isfinite(matrix).all(), 'Nonfinite MMN matrix')
            matrices.append(matrix)
        yield np.array(matrices)
    require(not stream.read().strip(), 'Trailing MMN data')


def polar(matrices):
    u, s, vh = np.linalg.svd(matrices)
    require(np.min(s) > 1e-10, 'Singular overlap link')
    return u @ vh, float(np.min(s))


def surface(lines):
    result = z2pack.surface.SurfaceData()
    for t, matrices in lines:
        data = z2pack.line.OverlapLineData(matrices)
        result.add_line(t, z2pack.line.LineResult(data, [], []))
    return result


def native_check(root, name, expected_z2):
    directory = root / 'native' / name
    native = np.loadtxt(directory / f'{name}.wilson', ndmin=2)
    log = (directory / 'run.log').read_text()
    require(f'Converged Z2 invariant: {expected_z2}' in log and 'PROGRAM ENDED AT' in log,
            'Native calculation did not finish with the expected index')
    require(np.array_equal(native[:, 0], np.arange(1, len(native) + 1)), 'Wrong loop numbering')
    raw_lines, polar_lines = [], []
    error, singular_error, minimum_singular = 0., 0., 1.
    with (directory / f'{name}.mmn').open() as stream:
        for i, matrices in enumerate(read_loops(stream, len(native))):
            unitary, minimum = polar(matrices)
            t = i / (len(native) - 1)
            raw_lines.append((t, matrices))
            polar_lines.append((t, unitary))
            centers = z2pack.line.OverlapLineData(unitary).wcc
            error = max(error, circle_error(native[i, 2:], centers))
            singular_error = max(singular_error, abs(minimum - native[i, 1]))
            minimum_singular = min(minimum_singular, minimum)
    require(error < 1e-9 and singular_error < 1e-10, 'Native overlap replay mismatch')
    raw_z2 = z2pack.invariant.z2(surface(raw_lines))
    polar_z2 = z2pack.invariant.z2(surface(polar_lines))
    require(raw_z2 == polar_z2 == expected_z2, 'Native/Z2Pack index mismatch')
    return {'loops': len(native), 'points_per_loop': len(raw_lines[0][1]),
            'raw_z2': raw_z2, 'polar_z2': polar_z2,
            'maximum_wilson_eigenvalue_error': error,
            'maximum_singular_value_error': singular_error,
            'minimum_link_singular_value': minimum_singular}


class ArchivedSystem(z2pack.system.OverlapSystem):
    def __init__(self, root):
        self.root = root
        self.requests = json.loads((root / 'adaptive/requests.json').read_text())
        self.used = []

    def get_mmn(self, kpt):
        points = np.array(kpt)
        matches = [entry for entry in self.requests
                   if np.shape(entry['points']) == points.shape and
                   np.allclose(entry['points'], points, atol=1e-13, rtol=0)]
        require(len(matches) == 1, 'Requested adaptive loop is missing or ambiguous')
        entry = matches[0]
        require(entry['polar'] is True, 'Archived adapter did not use polar links')
        path = self.root / 'adaptive/lines' / entry['directory'] / 'loop.mmn'
        with path.open() as stream:
            loops = list(read_loops(stream, 1))
        require(len(loops[0]) + 1 == len(points), 'Adaptive point count mismatch')
        self.used.append(entry['directory'])
        return polar(loops[0])[0]


def adaptive_check(root):
    system = ArchivedSystem(root)
    original = z2pack.io.load(str(root / 'adaptive/surface.json'), serializer=json)
    result = z2pack.surface.run(
        system=system, surface=lambda s, t: [t, s / 2, 0], num_lines=11,
        min_neighbour_dist=1e-4, iterator=[16, 32, 64, 128, 256, 512], pos_tol=1e-3)
    for item in (original, result):
        require(z2pack.invariant.z2(item) == 1, 'Adaptive index mismatch')
        require(not item.convergence_report['line']['PosCheck']['FAILED'] and
                not item.convergence_report['line']['PosCheck']['MISSING'], 'Unconverged lines')
        for name in ('MoveCheck', 'GapCheck'):
            require(not item.convergence_report['surface'][name]['FAILED'], 'Unconverged surface')
    require(np.allclose(original.t, result.t, atol=1e-13, rtol=0), 'Adaptive surfaces differ')
    error = max(circle_error(a, b) for a, b in zip(original.wcc, result.wcc, strict=True))
    require(error < 1e-9, 'Adaptive checkpoint spectrum differs')
    require(len(system.used) == len(system.requests) == 42, 'Unexpected adaptive replay length')
    require(set(system.used) == {entry['directory'] for entry in system.requests},
            'Not all adaptive requests replayed')
    return {'z2': 1, 'final_lines': len(result.t), 'replayed_requests': len(system.used),
            'maximum_checkpoint_eigenvalue_error': error, 'all_checks_passed': True}


def unpack_checked(archive, destination):
    hashes = {}
    with tarfile.open(archive, 'r:gz') as bundle:
        for member in bundle:
            path = PurePosixPath(member.name)
            require(member.isfile() and not path.is_absolute() and '..' not in path.parts,
                    'Unsafe bundle member')
            require(member.name not in hashes, 'Duplicate bundle member')
            target = destination / member.name
            target.parent.mkdir(parents=True, exist_ok=True)
            with bundle.extractfile(member) as source, target.open('wb') as output:
                shutil.copyfileobj(source, output)
            hashes[member.name] = hashlib.sha256(target.read_bytes()).hexdigest()
    manifest = json.loads((destination / 'manifest.json').read_text())
    require(set(hashes) == set(manifest['files']) | {'manifest.json'}, 'Incomplete bundle manifest')
    for name, record in manifest['files'].items():
        require(hashes[name] == record['sha256'], f'Checksum failure: {name}')
        require((destination / name).stat().st_size == record['bytes'], f'Size failure: {name}')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    logging.getLogger('z2pack').setLevel(logging.ERROR)
    with tempfile.TemporaryDirectory(prefix='wilson-replay-') as temporary:
        root = Path(temporary)
        manifest = unpack_checked(args.archive, root)
        result = {'versions': {name: version(name) for name in ('numpy', 'scipy', 'z2pack')},
                  'archive_sha256': hashlib.sha256(args.archive.read_bytes()).hexdigest(),
                  'files_verified': len(manifest['files']), 'new_dft_calculations': 0}
        for name, expected in [('neon', 0), ('stanene', 1)]:
            result[name] = native_check(root, name, expected)
            print(name, json.dumps(result[name]), flush=True)
        result['adaptive'] = adaptive_check(root)
    result['scope'] = ('Replay of retained overlaps and convergence history. '
                       'No new CP2K calculation, material-convergence study, or exact source-build reconstruction.')
    args.report.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
