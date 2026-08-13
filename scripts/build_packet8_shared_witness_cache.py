#!/usr/bin/env python3
"""Materialize an upgraded cache plus retuned reports as one witness cache."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from probe_packet8_shared_witness_io import load_shared_witness_arrays


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upgraded-cache", type=Path, required=True)
    parser.add_argument(
        "--extra-shared-report", type=Path, action="append", default=[]
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    names, constants, charges = load_shared_witness_arrays(
        args.upgraded_cache, args.extra_shared_report
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        names=np.asarray(names),
        constants=constants,
        charges=charges,
    )
    print(f"materialized_witnesses={len(constants)} output={args.output}")


if __name__ == "__main__":
    main()
