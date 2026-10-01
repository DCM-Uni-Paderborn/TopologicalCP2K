Consolidated closing CP2K source snapshot
=======================================

This bundle reconstructs the complete tracked source tree at the closing
implementation checkpoint, not only selected Fortran modules. It does not
contain build artifacts, binaries, Git metadata or untracked files.

The public upstream base and exact target/tree identities are recorded in
source-manifest.json. The full-index, binary-safe patch changes 296 files.
CP2K-LICENSE is copied verbatim from the target tree. Dependency notices
remain in ../source-freeze/topology-source-notices.tar.gz. This bundle does
not assign new terms to manuscript text, analysis scripts or dependencies.

Verification requires a CP2K Git repository containing the public base
commit in the manifest. The target commit is not required.

  python verify_closing_source.py BUNDLE_DIRECTORY CP2K_REPOSITORY

The verifier checks patch and license fingerprints, reads the public base
into a separate temporary Git index, applies the patch there, and compares
the resulting tree identity to the recorded target. It confirms that HEAD,
the working-tree status and the normal index are unchanged. Git may add
reconstructed blobs/trees to the repository's object database.

For an actual build, apply source.patch to a separate clean checkout of
the stated base, then use the retained build configurations and dependency
records appropriate to the desired solver. The verification command itself
does not apply changes to that checkout or run a build.

This is not a released CP2K version or a claim that the complete regression
suite was repeated. The manuscript distinguishes selected native checks,
archive-only reanalysis and material convergence. Historical calculations
continue to refer to their own earlier pinned source and runtime records;
this snapshot does not retroactively change their provenance.

Verified reconstruction, 1 October 2026:
  Target tree: 58ed6fe33f3b2533cd4c48b948e90620fea15419
  Checkout and normal index unchanged: true
