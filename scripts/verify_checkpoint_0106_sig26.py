#!/usr/bin/env python3
"""Run the stronger sigma=26 finite-n dense+dense checkpoint certificate."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def run(args: list[str]) -> int:
    print("==> " + " ".join(Path(a).name if a.endswith(".py") else a for a in args), flush=True)
    proc = subprocess.run([sys.executable, *args], check=False)
    print(flush=True)
    return proc.returncode


def main() -> int:
    here = Path(__file__).resolve().parent
    prefix = here / "verify_checkpoint_prefix.py"
    residual = here / "verify_checkpoint_residual.py"
    summary = here / "global_episode_cover_suffix_be_exactsurv30_k1048576_sig26_d0106_h600_r16_coarse_summary.csv"
    tailcert = here / "global_episode_cover_suffix_be_exactsurv30_k1048576_sig26_d0106_h600_r16_coarse_tailcert.csv"

    checks = [
        [
            str(prefix),
            "--summary",
            str(summary),
            "--tailcert",
            str(tailcert),
            "--sigma",
            "26",
            "--h-max",
            "600",
            "--r-cap",
            "16",
            "--prefix-r-max",
            "16",
            "--block-ratio",
            "1.02",
            "--suffix-block-ratio",
            "1.1",
            "--max-total-log2",
            "-0.9",
            "--max-tail-log2",
            "-300",
            "--max-ratio-log2",
            "-0.49",
            "--max-prefix-large-r-low-union-log2",
            "-50",
            "--max-prefix-large-r-high-union-log2",
            "-1000",
        ],
        [
            str(residual),
            "--sigma",
            "26",
            "--h-min",
            "600",
            "--r-cap",
            "16",
            "--min-ratio-slack",
            "0.002",
            "--max-low-slot-ledger",
            "-130",
            "--max-union-bound",
            "-50",
        ],
    ]

    ok = True
    for args in checks:
        ok = run(args) == 0 and ok
    print(f"SIGMA=26 CHECKPOINT STATUS: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
