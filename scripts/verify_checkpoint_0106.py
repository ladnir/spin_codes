#!/usr/bin/env python3
"""Run the complete 0.106 dense+dense checkpoint certificate.

This is the small top-level verifier for the current finite-n checkpoint:

    k = 2^20, N = 2^21, sigma = 25, delta = 0.106.

It checks the finite prefix summary and the analytic residual constants.  The
two component checks are kept separate so their proof obligations stay visible.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def run(script: Path) -> int:
    print(f"==> {script.name}", flush=True)
    proc = subprocess.run([sys.executable, str(script)], check=False)
    print(flush=True)
    return proc.returncode


def main() -> int:
    here = Path(__file__).resolve().parent
    scripts = [
        here / "verify_checkpoint_prefix.py",
        here / "verify_checkpoint_residual.py",
    ]
    ok = True
    for script in scripts:
        rc = run(script)
        ok = ok and rc == 0
    print(f"CHECKPOINT STATUS: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
