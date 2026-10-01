"""Verify the exact CP2K source tree using a separate temporary Git index."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("bundle", type=Path)
parser.add_argument("repository", type=Path)
args = parser.parse_args()
manifest = json.loads((args.bundle / "source-manifest.json").read_text())
patch = (args.bundle / "source.patch").read_bytes()
assert len(patch) == manifest["patch_bytes"]
assert hashlib.sha256(patch).hexdigest() == manifest["patch_sha256"]
assert hashlib.sha256((args.bundle / "CP2K-LICENSE").read_bytes()).hexdigest() == manifest["license_sha256"]


def git(*command, data=None, env=None):
    return subprocess.check_output(["git", *command], cwd=args.repository, input=data, env=env)


before = {"head": git("rev-parse", "HEAD"), "status": git("status", "--porcelain=v1", "-z"),
          "index": git("ls-files", "--stage", "-z")}
with tempfile.TemporaryDirectory(prefix="topology-source-replay-") as directory:
    env = os.environ | {"GIT_INDEX_FILE": str(Path(directory) / "index")}
    git("read-tree", manifest["base_commit"], env=env)
    git("apply", "--cached", "--binary", "-", data=patch, env=env)
    actual = git("write-tree", env=env).decode().strip()
assert actual == manifest["target_tree"] == manifest["reconstructed_tree"]
assert before == {"head": git("rev-parse", "HEAD"), "status": git("status", "--porcelain=v1", "-z"),
                  "index": git("ls-files", "--stage", "-z")}
print(json.dumps({"accepted": True, "tree": actual, "checkout_and_index_unchanged": True}))
