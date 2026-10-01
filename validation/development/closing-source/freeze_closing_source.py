"""Preserve and verify the complete closing CP2K source tree without changing it."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile


root = Path("/Users/tkuehne/cp2k-spectral-localizer")
output = Path(__file__).resolve().parent / "closing-source"
base = "92574dc7b41b039d88f32e23bbd59a71b0eb8836"


def git(*args, data=None, env=None):
    return subprocess.check_output(["git", *args], cwd=root, input=data, env=env)


assert not git("status", "--porcelain"), "Do not freeze a dirty source checkout"
target = git("rev-parse", "HEAD").decode().strip()
tree = git("rev-parse", "HEAD^{tree}").decode().strip()
patch = git("diff", "--binary", "--full-index", base, target)
with tempfile.TemporaryDirectory(prefix="topology-source-verification-") as temporary:
    env = os.environ | {"GIT_INDEX_FILE": str(Path(temporary) / "index")}
    git("read-tree", base, env=env)
    git("apply", "--cached", "--binary", "-", data=patch, env=env)
    reconstructed = git("write-tree", env=env).decode().strip()
assert reconstructed == tree
assert not git("status", "--porcelain")
assert git("rev-parse", "HEAD").decode().strip() == target
output.mkdir(exist_ok=False)
(output / "source.patch").write_bytes(patch)
license_text = git("show", target + ":LICENSE")
(output / "CP2K-LICENSE").write_bytes(license_text)
manifest = {
    "base_commit": base, "base_url": "https://github.com/cp2k/cp2k/commit/" + base,
    "target_commit": target, "target_tree": tree, "reconstructed_tree": reconstructed,
    "exact_tree_reconstruction": True, "source_checkout_unchanged": True,
    "patch_sha256": hashlib.sha256(patch).hexdigest(), "patch_bytes": len(patch),
    "license_sha256": hashlib.sha256(license_text).hexdigest(),
    "changed_files": git("diff", "--name-only", base, target).decode().splitlines(),
    "scope": "Complete tracked CP2K source tree, not an upstream release, dependency bundle or claim of full-suite validation",
    "dependency_notices": "../source-freeze/topology-source-notices.tar.gz",
    "historical_runs": "Retain their own earlier source identities; this snapshot does not retrospectively reassign them",
}
(output / "source-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps({"exact_tree_reconstruction": True, "changed_files": len(manifest["changed_files"]),
                  "patch_bytes": len(patch), "target": target}))
