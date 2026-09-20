#!/usr/bin/env python3
"""Freeze the power-of-two wrapped-EC parent certificate candidate.

This is a thin integration driver for the expander-code project's existing
certificate generator.  It fixes the parent parameters needed by the
shorten-and-puncture wrapper in EXPANDER_CODE_ALTERNATIVE_AUDIT.md.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


WORKSTREAM = Path(__file__).resolve().parent
DEFAULT_EXPANDER_ROOT = Path(
    "C:/Users/peter/.codex/worktrees/30ff/permute_conv/expander_codes"
)
DEFAULT_OUTPUT = WORKSTREAM / "expander_ec_d10_m15_k20_parent_candidate.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expander-root", type=Path, default=DEFAULT_EXPANDER_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    scripts = args.expander_root.resolve() / "scripts"
    if not scripts.is_dir():
        parser.error(f"missing expander-code scripts directory: {scripts}")
    sys.path.insert(0, str(scripts))
    from binary_biregular_ec_certificate import generate_certificate

    certificate = generate_certificate(
        k=1_048_585,
        left_degree=10,
        right_degree=5,
        cutoff=228_607,
        memory=15,
        exact_limit=32,
        intermediate_relative_width=0.2,
        dense_relative_width=0.02,
        target_bits=40,
        precision_bits=192,
    )
    args.output.write_text(json.dumps(certificate, indent=2) + "\n", encoding="utf-8")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
