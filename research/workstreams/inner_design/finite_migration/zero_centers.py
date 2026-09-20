"""Diagnose unresolved centers with the integer all-one split backend."""
import argparse
from pathlib import Path
from types import SimpleNamespace

import budget_probe
import zero_constant_dense


def run(seed,output,count):
    original = budget_probe.budget.short_dense
    budget_probe.budget.short_dense = SimpleNamespace(
        Checker=zero_constant_dense.Checker,mixed_dense=zero_constant_dense.parent.short_dense.mixed_dense)
    try:
        budget_probe.run(seed,output,count)
    finally:
        budget_probe.budget.short_dense = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--count',type=int,default=5)
    a = p.parse_args()
    run(a.seed.resolve(),a.output.resolve(),a.count)
