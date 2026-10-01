"""Verify retained hashes and numerical comparisons without a CP2K installation."""

import argparse
import hashlib
import io
import json
import math
from pathlib import Path, PurePosixPath
import re
import tarfile

REFERENCE_SHA256 = 'b8d380dd7b02fa1cc519f3d0afdfb16a7c8e37fe18afdecd919cc01c5f4ec897'
SOURCE_COMMIT = '851fe2178997f1fbb2a8c2141a1da5fff40e644d'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def unpack(data):
    files = {}
    with tarfile.open(fileobj=io.BytesIO(data)) as archive:
        for member in archive:
            path = PurePosixPath(member.name)
            require(member.isfile() and not path.is_absolute() and '..' not in path.parts,
                    'Unsafe archive member: ' + member.name)
            require(member.name not in files, 'Duplicate archive member')
            require(member.size < 64*1024**2, 'Unexpectedly large evidence member')
            files[member.name] = archive.extractfile(member).read()
    require('manifest.json' in files, 'Missing manifest')
    manifest = json.loads(files.pop('manifest.json'))
    require(set(files) == set(manifest['files']), 'Manifest member set mismatch')
    for name, record in manifest['files'].items():
        require(len(files[name]) == record['bytes'] and digest(files[name]) == record['sha256'],
                'Manifest mismatch: ' + name)
    return files


def verify(files):
    reference = files['reference/methods-results.tar.gz']
    require(digest(reference) == REFERENCE_SHA256, 'Historical reference was replaced')
    oldfiles = unpack(reference)
    oldindex = json.loads(files['reference/index.json'])
    provenance = json.loads(files['native/source-provenance.json'])
    require(oldindex['source_commit'] == provenance['source_commit'] == SOURCE_COMMIT,
            'Source revision mismatch')
    for name, expected in oldindex['source_sha256'].items():
        require(digest(oldfiles[name]) == expected == provenance['source_sha256'][name],
                'Source hash mismatch: ' + name)
    source = json.loads(files['source-reproduction/source-manifest.json'])
    require(source['exact_tree_reconstruction'] and source['target_tree'] == source['reconstructed_tree'],
            'Source tree was not reconstructed')
    require(source['target_commit'] == SOURCE_COMMIT and
            source['patch_sha256'] == digest(files['source-reproduction/source.patch']),
            'Source patch mismatch')
    build = json.loads(files['native/build-result.json'])
    require(build['runtime_check_symbols'] and any('ubsan_handle' in name for name in build['runtime_check_symbols']),
            'Missing UBSan instrumentation evidence')
    closure = json.loads(files['native/dependency-closure.json'])
    require(closure['library_sha256'] == build['sha256'], 'Dependency closure belongs to another library')
    commands = json.loads(files['native/compile_commands.json'])
    require(len(commands) == 1 and '-fsanitize=undefined' in commands[0]['command'] and
            '-fno-sanitize-recover=all' in commands[0]['command'], 'Unexpected compiler options')
    require(build['compile_commands_sha256'] == digest(files['native/compile_commands.json']),
            'Compile-command hash mismatch')
    commands_checked = 0
    for name, data in files.items():
        if not name.startswith('native/') or not name.endswith('.command.json'):
            continue
        command = json.loads(data)
        log = files['native/' + command['log']]
        require(command['returncode'] == 0 and command['log_sha256'] == digest(log),
                'Failed command or changed log: ' + name)
        require(not re.search(rb'runtime error:|Assertion.*failed', log), 'Runtime failure in ' + name)
        commands_checked += 1
    require(commands_checked == 15, 'Expected configure/build, three synthetic and ten material commands')
    synthetic = []
    for label, total, resolved, singular in [('independent', 600, 468, 132), ('stress', 636, 336, 300)]:
        actual = json.loads(files['native/' + label + '.json'])
        old = json.loads(oldfiles[f'build-mpi/delayed-native-verified-{label}.json'])
        require(len(actual['results']) == len(old['results']) == total, 'Synthetic query count')
        for new, previous in zip(actual['results'], old['results'], strict=True):
            for key in ('case', 'scale', 'expected', 'singular'):
                require(new[key] == previous[key], 'Synthetic identity mismatch: ' + label)
            if new['singular']:
                require(new['sign'] is None, 'Singular query was accepted')
            else:
                require(new['sign'] == new['expected'] == previous['sign'], 'Wrong Pfaffian sign')
                require(math.isfinite(new['residual']) and 0 <= new['residual'] <= 1e-10,
                        'Unresolved factor residual')
        require(sum(not row['singular'] for row in actual['results']) == resolved and
                sum(row['singular'] for row in actual['results']) == singular, 'Wrong singular counts')
        synthetic.append({'queries': total, 'resolved': resolved, 'singular': singular})
    pattern = json.loads(files['native/roundoff-chain.log'])
    oldpattern = json.loads(oldfiles['build-mpi/delayed-native-pattern-chain-after.json'])
    require(len(pattern) == len(oldpattern) == 930, 'Roundoff pattern count')
    for new, previous in zip(pattern, oldpattern, strict=True):
        require((new['i'], new['j']) == (previous['i'], previous['j']), 'Roundoff coordinate mismatch')
        require(new['status'] == 0 and new['sign'] == 1 and math.isfinite(new['residual'])
                and 0 <= new['residual'] <= 1e-10, 'Roundoff pattern failed')
    api = json.loads(files['native/api-boundaries.json'])
    require(api['library_sha256'] == build['sha256'] and api['passed'] and len(api['results']) == 63,
            'Incomplete ABI test record')
    for row in api['results']:
        require(row['status'] == row['expected_status'] and row['sign'] == row['expected_sign'],
                'ABI expectation failed')
        require(math.isfinite(row['residual']) and row['residual'] >= 0, 'Invalid ABI residual')
    case_index = json.loads(files['reference/stanene-convergence-index.json'])
    require(len(case_index['cases']) == 10, 'Material case count')
    differences, factor_residuals, eigen_residuals, dense_errors = [], [], [], []
    for case in case_index['cases']:
        label = case['label']
        old = json.loads(oldfiles[f'build-mpi/delayed-native-materials-verified/{label}.json'])
        actual = json.loads(files[f'native/{label}.json'])
        require(actual['library_sha256'] == build['sha256'], 'Material library mismatch')
        require(actual['inputs'] == old['inputs'] and actual['script_sha256'] == old['script_sha256'],
                'Material inputs or analysis changed')
        require(len(actual['queries']) == len(old['queries']), 'Material query count')
        for new, previous in zip(actual['queries'], old['queries'], strict=True):
            for key in ('energy', 'eta', 'flatten_scale', 'nu'):
                require(new[key] == previous[key], 'Material mismatch: ' + label + ' ' + key)
            require(new['nu'] in (0, 1), 'Unresolved material index')
            gap_error = abs(new['gap']-previous['gap'])
            require(math.isfinite(gap_error) and gap_error < 1e-10, 'Changed material gap')
            require(math.isfinite(new['factor_residual']) and 0 <= new['factor_residual'] <= 1e-10,
                    'Material factor residual')
            require(math.isfinite(new['eigen_residual']) and 0 <= new['eigen_residual'] < 1e-8,
                    'Material eigen residual')
            attempts = new['factor_attempts']
            require(len(attempts) == 1 and attempts[0]['ordering'] == 'original' and attempts[0]['resolved'],
                    'Material needed extra orderings')
            differences.append(gap_error)
            factor_residuals.append(new['factor_residual'])
            eigen_residuals.append(new['eigen_residual'])
            if label == 'mesh3':
                require(new['dense_nu'] == new['nu'] and math.isfinite(new['dense_gap_error']) and
                        0 <= new['dense_gap_error'] < 1e-10, 'Independent dense check failed')
                dense_errors.append(new['dense_gap_error'])
    require(len(differences) == 33 and len(dense_errors) == 3, 'Incomplete material comparison')
    return {'passed': True, 'synthetic': synthetic, 'roundoff_patterns': 930, 'api_queries': 63,
            'material_cases': 10, 'material_queries': 33, 'max_gap_error': max(differences),
            'max_factor_residual': max(factor_residuals), 'max_eigen_residual': max(eigen_residuals),
            'independent_dense_queries': 3, 'max_dense_gap_error': max(dense_errors),
            'scope': 'Archive integrity and independent re-evaluation of retained numerical comparisons. '
                     'No native executable, SCF, factorization or MPI process is rerun by this checker.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    data = args.archive.read_bytes()
    files = unpack(data)
    result = verify(files)
    result.update(archive_sha256=digest(data), verified_files=len(files))
    text = json.dumps(result, indent=2) + '\n'
    if args.report:
        with args.report.open('x') as handle:
            handle.write(text)
    print(text, end='')


if __name__ == '__main__':
    main()
