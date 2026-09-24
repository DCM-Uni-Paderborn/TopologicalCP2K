"""Create numeric native-test fixtures from pinned independent affine reference data.

Usage: cache-dir output-fixture provenance-jsonl [space-group ...]
Only arithmetic AST nodes and x/y/z are accepted; source strings are not evaluated.
No third-party code is imported. Generated files contain derived reference data.
"""
import ast
import hashlib
import json
import re
import sys
import time
import urllib.request
from fractions import Fraction as F
from pathlib import Path

REVISION = "9cfb644a1a4cc1c7baec457b321885459bb98bd1"
BASE = f"https://raw.githubusercontent.com/thchr/Crystalline.jl/{REVISION}/"
cache = Path(sys.argv[1])
cache.mkdir(parents=True, exist_ok=True)
groups = list(map(int, sys.argv[4:])) or list(range(1, 231))
provenance = []


def download(path, name):
    target = cache / name
    url = BASE + path
    if not target.exists():
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=30) as response:
                    data = response.read()
                target.write_bytes(data)
                break
            except OSError:
                if attempt == 2:
                    raise
                time.sleep(2)
    data = target.read_bytes()
    provenance.append(dict(url=url, revision=REVISION, sha256=hashlib.sha256(data).hexdigest()))
    return data.decode()


def expression(text):
    text = re.sub(r"(?<=[0-9])(?=[xyz])", "*", text)
    def parse(node):
        if isinstance(node, ast.Constant) and type(node.value) in (int, float):
            return [F(str(node.value)), F(0), F(0), F(0)]
        if isinstance(node, ast.Name) and node.id in ("x", "y", "z"):
            return [F(0)] + [F(node.id == var) for var in ("x", "y", "z")]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            sign = -1 if isinstance(node.op, ast.USub) else 1
            return [sign * v for v in parse(node.operand)]
        if isinstance(node, ast.BinOp):
            a, b = parse(node.left), parse(node.right)
            if isinstance(node.op, (ast.Add, ast.Sub)):
                sign = -1 if isinstance(node.op, ast.Sub) else 1
                return [x + sign*y for x, y in zip(a, b)]
            if isinstance(node.op, ast.Mult):
                if not any(a[1:]):
                    return [a[0]*v for v in b]
                if not any(b[1:]):
                    return [b[0]*v for v in a]
            if isinstance(node.op, ast.Div) and not any(b[1:]) and b[0]:
                return [v/b[0] for v in a]
        raise ValueError(f"Not an affine arithmetic expression: {text!r}")
    return parse(ast.parse(text.strip(), mode="eval").body)


def affine(text):
    rows = [expression(part) for part in text.split(",")]
    assert len(rows) == 3
    linear = tuple(v for row in rows for v in row[1:])
    assert all(v.denominator == 1 for v in linear)
    return tuple(map(int, linear)), tuple(row[0] for row in rows)


def close_group(generators):
    identity = (1, 0, 0, 0, 1, 0, 0, 0, 1)
    # All supplied crystallographic translations have denominators dividing 24.
    gens = []
    for r, t in generators:
        assert all((24*v).denominator == 1 for v in t)
        gens.append((r, tuple(int(24*v) % 24 for v in t)))
    work = [(identity, (0, 0, 0))]
    seen = set(work)
    for r, t in work:
        for a, b in gens:
            ar = tuple(sum(a[3*i+k]*r[3*k+j] for k in range(3)) for i in range(3) for j in range(3))
            at = tuple((sum(a[3*i+k]*t[k] for k in range(3)) + b[i]) % 24 for i in range(3))
            item = (ar, at)
            if item not in seen:
                seen.add(item)
                work.append(item)
                assert len(work) <= 192
    return work


def column_major(values):
    return [values[3*i+j] for j in range(3) for i in range(3)]


lines = [str(len(groups))]
for sg in groups:
    assert 1 <= sg <= 230
    generators = download(f"test/data/xyzt/generators/sgs/3d/{sg}.csv", f"generators-{sg}.csv")
    group = close_group([affine(line) for line in generators.splitlines() if line.strip()])
    rows = download(f"data/wyckpos/3d/{sg}.csv", f"wyckoff-{sg}.csv").splitlines()
    lines.append(f"{sg} {len(group)} {len(rows)}")
    for r, t in group:
        lines.append(" ".join(map(str, column_major(r) + [float(F(v, 24)) for v in t])))
    for row in rows:
        mult, label, text = row.split("|")
        directions, position = affine(text)
        lines.append(" ".join(map(str, [label, int(mult)] + list(map(float, position)) + column_major(directions))))
    print(f"Prepared group {sg}: {len(group)} operations, {len(rows)} families", flush=True)
Path(sys.argv[2]).write_text("\n".join(lines) + "\n")
Path(sys.argv[3]).write_text("".join(json.dumps(row) + "\n" for row in provenance))
