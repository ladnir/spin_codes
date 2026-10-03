"""Cache floating local operators from one actual 24-bit map census.

These small, source-pinned receipts support weight-tilt searches without
repeating the 2^24-state census. They are proposals, not outward endpoints.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'larger_state'))
import screen24
import numpy as np


def source_pins():
    pins = screen24.sources()
    path = Path(__file__).resolve()
    pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


def main(args):
    tilts = tuple(Fraction(x) for x in args.tilts)
    if (not tilts or min(tilts) <= 0 or len(set(tilts)) != len(tilts)
            or args.output_dir.exists()):
        raise ValueError('distinct positive tilts and a fresh output directory required')
    pins = source_pins()
    started = monotonic()
    data = screen24.maps24.make_maps()
    counts = screen24.maps24.census(data, include_character=True)
    census_seconds = monotonic()-started
    args.output_dir.mkdir(parents=True)
    print(f'actual24 census: {census_seconds:.2f}s', flush=True)
    for tilt in tilts:
        local_started = monotonic()
        local, diagnostic = screen24.local_operators(data, counts,
            np.exp(-float(tilt)), progress=True)
        if source_pins() != pins:
            raise ArithmeticError('source changed during cache construction')
        result = dict(schema='packet8-actual24-local-cache-proposal-1',
            tilt=str(tilt), map_record=data['record'],
            local_operator_matrices=local.tolist(), local_diagnostics=diagnostic,
            source_sha256=pins, source_pins_verified_at_finish=True,
            proposal_only=True, whole_code_certificate=False,
            has_outward_endpoints=False, actual_larger_state_maps=True,
            arithmetic='float64 finite formulas, explicit diagnostic clipping',
            census_seconds=census_seconds, local_seconds=monotonic()-local_started,
            elapsed_seconds=monotonic()-started)
        path = args.output_dir/f'theta_{tilt.numerator}_{tilt.denominator}.json'
        path.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(f'saved theta={tilt}: local={result["local_seconds"]:.2f}s', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts', nargs='+', default=['.1','.2','.3','.4','.5','.6','.7','.8'])
    parser.add_argument('--output-dir', type=Path, required=True)
    main(parser.parse_args())
