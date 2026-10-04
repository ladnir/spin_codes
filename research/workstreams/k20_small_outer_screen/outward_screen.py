"""Outward selected-occupancy screen, not a whole-code certificate."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import sys
from time import monotonic

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'k18_margin60'))
import replay_fractional as replay


def scalar(regional, q, tilt, alpha, k):
    groups = k // 128
    beta = F(1 << 256, 65535**8)
    power = regional ** 64
    moment = sum((power[0, j] for j in range(power.ncols())), arb(0))
    prefactor = (replay.aq(alpha) * (
        q * replay.aq(beta).log() + replay.aq(tilt) * ((2*k)//10))).exp()
    return replay.upper(arb(comb(groups, q)) * prefactor * moment)


def run(local, output, alphas, qs, *, k=1 << 20):
    groups, remainder = divmod(k, 128)
    if (remainder or groups % 32 or output.exists() or not qs
            or min(qs) < 1 or max(qs) > groups):
        raise ValueError('natural geometry, fresh output, and valid occupancies required')
    ctx.prec = 256
    physical, tilt, source = replay.load_local(local)
    pins = dict(source['source_sha256'])
    for path in (Path(local), Path(__file__), Path(replay.__file__), Path(replay.placement.__file__)):
        pins[str(path.resolve())] = replay.digest(path)
    trials = []
    result = dict(schema='k20-small-outer-selected-outward-1',
        whole_code_certificate=False, all_occupancies_covered=False,
        message_bits=k, output_bits=2*k, cutoff=(2*k)//10,
        outer_groups=groups, outer_dimension=128, outer_length=256,
        regions=64, slots_per_region=groups, steps_per_macro=2,
        local_receipt=str(local), map_record=source['map_record'],
        tilt=str(tilt), precision=256, selected_occupancies=list(qs),
        source_sha256=pins, trials=trials)
    for alpha in alphas:
        alpha = F(alpha)
        start = monotonic()
        blocks = replay.fractional_blocks(physical, alpha)
        print(f'G2 operators ready: alpha={alpha}', flush=True)
        regionals = replay.placement.placement_power(blocks,
            epochs=groups//32, windows=32, maximum_groups=max(qs))
        endpoints = {str(q): replay.endpoint(scalar(regionals[q], q, tilt, alpha, k)) for q in qs}
        trial = dict(alpha=str(alpha), endpoints=endpoints,
            margins={q: replay.bits(replay.fraction(value)) for q, value in endpoints.items()},
            seconds=monotonic()-start)
        trials.append(trial)
        replay.authenticate(pins)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(json.dumps({key: trial[key] for key in ('alpha', 'margins', 'seconds')}), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--alphas', nargs='+', default=['.8', '.9', '1'])
    parser.add_argument('--qs', type=int, nargs='+', default=[64,128,192,256,320,356,400,512])
    parser.add_argument('--k', type=int, default=1 << 20)
    args = parser.parse_args()
    run(args.local, args.output, args.alphas, args.qs, k=args.k)
