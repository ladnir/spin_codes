"""Bounded, pinned all-profile FLOAT gates from preserved GF256 witnesses.

This is a research driver. It does not certify gaps between selected weights.
It restores the complete old marker bank, not a thinned subset.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import monotonic

from hd_field128_lift_alias_screen import AliasBankScreen, source_pins
from hd_field128_lift_screen import Budget


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--centers', type=int, nargs='+', required=True)
    parser.add_argument('--distance', default='7/20')
    parser.add_argument('--max-seconds', type=float, default=50.)
    parser.add_argument('--method', choices=('runs', 'marginal'), default='runs')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists')
    directory = Path(__file__).resolve().parent
    manifest_path = directory/'hd_field128_lift_tail_seed_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    seeds = {row['central_weight']: row for row in manifest['seeds']}
    if any(r not in seeds or r == 128 for r in args.centers):
        parser.error('This driver handles the saved JSON banks: 320,1024,2048,3584,3936')
    budget = Budget(args.max_seconds)
    rows, failure = [], None
    pins = source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    pins[manifest_path.name] = sha256(manifest_path.read_bytes()).hexdigest()
    for r in args.centers:
        try:
            budget.check()
            seed = seeds[r]
            source = directory/seed['source_path']
            digest = sha256(source.read_bytes()).hexdigest()
            if digest != seed['source_sha256']:
                raise RuntimeError('Preserved source pin mismatch: '+source.name)
            old = json.loads(source.read_text())['fixed_witness']
            z = F(seed['source_exact_z'])**16
            alpha = F(seed['source_exact_alpha'])/16
            ratio = log2(F(seed['source_exact_potential']))
            remaining = args.max_seconds-(monotonic()-budget.started)
            screen = AliasBankScreen(float(z), float(alpha), old['marker_log_odds'], ratio,
                                     distance=F(args.distance), max_seconds=remaining)
            budget.check()
            row = screen.evaluate(r, method=args.method, optimize_potential=True,
                                  max_seconds=args.max_seconds-(monotonic()-budget.started))
            row.update(witness=screen.witness(), source_witness_path=source.name,
                       source_witness_sha256=digest, transformed_exact_z=str(z),
                       transformed_exact_alpha=str(alpha))
            rows.append(row)
            print(f"r={r} log2_bound={row['log2_upper_screen']:.9f} "
                  f"diagnostics={row['regional_alias_diagnostics']}", flush=True)
        except (TimeoutError, ArithmeticError) as error:
            failure = str(error)
            break
    result = dict(requested_weights=args.centers, results=rows,
                  completed_requested_weights=len(rows)==len(args.centers), failure=failure,
                  source_sha256=pins, pin_scope='listed sources, not all transitive imports',
                  elapsed_seconds=monotonic()-budget.started, max_seconds=args.max_seconds,
                  target_distance=args.distance, field_bits=128, directed=False,
                  whole_code_certificate=False,
                  scope='all supports only at the explicitly completed selected weights')
    with args.output.open('x', encoding='utf-8') as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write('\n')


if __name__ == '__main__':
    main()
