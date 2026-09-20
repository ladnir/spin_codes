#!/usr/bin/env python3
"""Run the sigma=32 finite-n checkpoint with a terminated outer parity tail."""

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
    summary = here / "global_episode_cover_tail32_suffix_be_exactsurv30_k1048576_sig32_d0106_h500_r12_coarse_summary.csv"
    tailcert = here / "global_episode_cover_tail32_suffix_be_exactsurv30_k1048576_sig32_d0106_h500_r12_coarse_tailcert.csv"

    checks = [
        [
            str(prefix),
            "--summary",
            str(summary),
            "--tailcert",
            str(tailcert),
            "--n",
            str(2 * 1048576 + 32),
            "--d",
            "222301",
            "--parity-extra",
            "32",
            "--sigma",
            "32",
            "--h-max",
            "500",
            "--r-cap",
            "12",
            "--prefix-r-max",
            "12",
            "--block-ratio",
            "1.02",
            "--suffix-block-ratio",
            "1.1",
            "--peak-h",
            "8",
            "--max-total-log2",
            "-6.2",
            "--max-tail-log2",
            "-300",
            "--max-ratio-log2",
            "-0.68",
            "--max-prefix-large-r-low-union-log2",
            "-100",
            "--max-prefix-large-r-high-union-log2",
            "-1000",
        ],
        [
            str(residual),
            "--sigma",
            "32",
            "--parity-extra",
            "32",
            "--h-min",
            "500",
            "--r-cap",
            "12",
            "--min-ratio-slack",
            "0.01",
            "--max-low-slot-ledger",
            "-200",
            "--max-union-bound",
            "-100",
        ],
    ]

    ok = True
    for args in checks:
        ok = run(args) == 0 and ok
    print(f"SIGMA=32 TERMINATED-TAIL CHECKPOINT STATUS: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
