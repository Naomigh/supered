# Third-party source provenance

The upstream files below are vendored as ordinary files, so a normal clone can
build supered without submodules or access to the original workspace.

## Edlib

- Upstream: https://github.com/Martinsos/edlib
- Commit: `0ddc23ea06abd95ab6227442acb4c27bd7607b68`
- Included: original implementation, public header, README, and MIT license.
- License: [edlib_origin/LICENSE](edlib_origin/LICENSE).
- `generated/original_nw_preallocated.hpp` is derived from this Edlib source;
  `scripts/generate_original.py` regenerates it with stack block storage.

## Corrected LEAP

- Upstream: https://github.com/YourongZhou/LEAP
- Branch: `codex/leap-correctness-fixes`
- Commit: `6a711282f3494b6b9886236e58c7dfba6efa885a`
- Included: the seven implementation files, their headers, and upstream README
  needed for the backend used by supered. Unrelated Bowtie2 trees and data are
  omitted. The upstream README also describes targets from the full LEAP tree.
- Original source notices are retained. No replacement license is assigned to
  the LEAP sources by this packaging step.

Each vendor directory contains `UPSTREAM.json` with the pinned commit and
SHA256 hashes of its unmodified upstream files.
