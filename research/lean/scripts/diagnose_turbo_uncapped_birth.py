"""Bounded proof-witness diagnostic; not a Lean or whole-code certificate.

Reuse retained outward live/zero endpoints, and replace selected newborn
rows by a coarse arbitrary-live mass g(z)^j. No state census, search, raw
array export, construction change, or encoder benchmark is performed.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
import sys
import subprocess
from time import monotonic

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'research/workstreams/k18_margin60'))
import replay_fractional as fractional
from flint import arb, arb_mat, ctx


def run(args):
    if args.output.exists():
        raise FileExistsError('fresh diagnostic output required')
    if args.groups not in (2048, 8192) or not all(3 <= q <= args.groups for q in args.q):
        raise ValueError('K18/K20 groups and feasible positive tail occupancies required')
    if not 0 < args.alpha <= 1:
        raise ValueError('fractional exponent in (0,1] required')
    ctx.prec = 256
    physical, tilt, source = fractional.load_local(args.physical)
    source_hash = hashlib.sha256(args.physical.read_bytes()).hexdigest()
    started = monotonic()
    records = []
    for mode in args.mode:
        if monotonic() - started >= args.max_seconds:
            break
        local = [arb_mat(m) for m in physical]
        z = (-fractional.aq(tilt)).exp()
        g = ((1 + z) ** 4 - 1) / 15
        changed = []
        for j in range(1, 17):
            if (mode == 'all' or (mode == 'capped' and source['birth_row_choices'][j] > 0)
                    or (mode.startswith('from') and j >= int(mode[4:]))):
                changed.append(j)
                for k in range(1, 11):
                    local[j][0, k] = 0
                local[j][0, 1] = fractional.upper(g ** j)
        blocks = fractional.fractional_blocks(local, args.alpha)
        regional = fractional.placement.placement_power(blocks,
            epochs=args.groups // 32, windows=32, maximum_groups=max(args.q))
        endpoints = {}
        beta = Q(1 << 256, 65535 ** 8)
        for q in args.q:
            power = regional[q] ** 64
            moment = sum((power[0, k] for k in range(11)), arb(0))
            prefactor = (fractional.aq(args.alpha) *
                (q * fractional.aq(beta).log() + fractional.aq(tilt) * args.cutoff)).exp()
            endpoint = fractional.endpoint(arb(comb(args.groups, q)) * prefactor * moment)
            endpoints[str(q)] = dict(upper=endpoint,
                margin_bits_display=fractional.bits(fractional.fraction(endpoint)))
        records.append(dict(mode=mode, changed_birth_occupancies=changed,
            endpoints=endpoints, elapsed_seconds=monotonic() - started))
        print(json.dumps(records[-1]), flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        partial = dict(schema='turbo-coarse-birth-diagnostic-1', status='partial',
            scope='Outward diagnostic only; not a complete certificate.',
            source=str(args.physical.resolve()), source_sha256=source_hash,
            groups=args.groups, cutoff=args.cutoff, alpha=str(args.alpha), tilt=str(tilt),
            requested_modes=args.mode, completed_modes=len(records), results=records,
            elapsed_seconds=monotonic() - started)
        args.output.write_text(json.dumps(partial, indent=2) + '\n', encoding='utf-8')
    result = dict(schema='turbo-coarse-birth-diagnostic-1',
        status='complete' if len(records) == len(args.mode) else 'budget_stopped',
        scope='Outward Python proof-witness diagnostic using retained live/zero entries; '
              'not Lean-authenticated physical rows or a complete occupation certificate.',
        source=str(args.physical.resolve()), source_sha256=source_hash,
        groups=args.groups, cutoff=args.cutoff, alpha=str(args.alpha), tilt=str(tilt),
        requested_modes=args.mode, completed_modes=len(records),
        results=records, elapsed_seconds=monotonic() - started)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', encoding='utf-8') as out:
        json.dump(result, out, indent=2)
        out.write('\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--physical', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--groups', type=int, default=2048)
    parser.add_argument('--cutoff', type=int, default=52428)
    parser.add_argument('--alpha', type=Q, default=Q(4, 5))
    parser.add_argument('--q', type=int, nargs='+', required=True)
    parser.add_argument('--mode', choices=('retained', 'capped', 'all', 'from2', 'from3', 'from4', 'from5', 'from6'), nargs='+',
                        default=['capped', 'all'])
    parser.add_argument('--max-seconds', type=float, default=120)
    parser.add_argument('--_worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.max_seconds <= 0:
        parser.error('--max-seconds must be positive')
    if args._worker:
        run(args)
    else:
        try:
            result = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                *sys.argv[1:], '--_worker'], timeout=args.max_seconds)
            raise SystemExit(result.returncode)
        except subprocess.TimeoutExpired:
            partial = (json.loads(args.output.read_text(encoding='utf-8'))
                       if args.output.exists() else {'results': []})
            partial.update(schema='turbo-coarse-birth-diagnostic-1', status='timed_out',
                scope='Incomplete outward diagnostic; no whole-code claim.',
                timeout_seconds=args.max_seconds)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(partial, indent=2) + '\n', encoding='utf-8')
            print('Numeric diagnostic worker exceeded its wall-clock limit', file=sys.stderr)
            raise SystemExit(2)
