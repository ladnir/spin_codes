"""Bounded neighboring-weight FLOAT screens from one exact positive witness."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import hd_field128_gamma_midpoint_enrich as enrich


provider, mode, HERE = enrich.provider, enrich.mode, enrich.HERE


def run(source_path, weights, *, max_seconds=90.):
    weights = tuple(weights)
    if (not 1 <= len(weights) <= 8 or len(set(weights)) != len(weights)
            or any(type(r) is not int or not 2 <= r < 4096 for r in weights)
            or not isfinite(max_seconds) or not 0 < max_seconds <= 90):
        raise ValueError('Expected1..8 distinct weights2..4095 and at most90 seconds')
    path = Path(source_path).resolve()
    raw = path.read_bytes()
    source = json.loads(raw)
    state, distance = source['state_dimension'], F(source['target_distance'])
    if source['construction'] != provider.gamma.construction(state):
        raise ValueError('Source construction mismatch')
    witness = source['witness']
    positives, marker_key = enrich.positive_markers(witness)
    z, alpha = F(witness['z']), F(witness['alpha'])
    odds = [tuple(log2(float(v)) for v in pair) for pair in positives]
    pins = enrich.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    for name, digest in source['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting source pin: '+name)
        pins[name] = digest
    def check():
        if path.read_bytes() != raw:
            raise RuntimeError('Source witness changed')
        for name, digest in pins.items():
            if sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Source closure changed: '+name)
    check()
    started, budget = monotonic(), provider.old.Budget(max_seconds)
    completed, status, pending, metadata = [], 'complete', None, None
    try:
        coeff = provider.AliasCoefficients(float(z), float(alpha), state=state, _budget=budget)
        bank = provider.AliasBankScreen(float(z), float(alpha), odds, 0., state=state, distance=distance,
                    grid=witness.get('grid', 128), floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
        metadata = bank.metadata()
        for weight in weights:
            budget.check()
            pending = weight
            regions = bank.regional(weight)
            row = mode.support_sum(4096, weight, regions, max_singletons=64, budget=budget)
            prefix = bank.constant+bank.alpha*bank.density*weight
            row.update(active_groups=weight, log2_outer_prefix=prefix,
                       log2_upper_screen=prefix+row['log2_support_sum_screen'],
                       regional_alias_diagnostics=dict(bank.regional_diagnostics),
                       regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                       complete_selected_weight_screen=True, elapsed_seconds=monotonic()-started)
            completed.append(row)
            pending = None
            print(json.dumps(dict(weight=weight, score=row['log2_upper_screen'], largest=row['largest_profile'],
                                  diagnostics=row['regional_alias_diagnostics'], elapsed_seconds=row['elapsed_seconds'])), flush=True)
    except TimeoutError:
        status = 'timeout'
    check()
    return dict(schema='gamma-fixed-positive-marker-neighbors-FLOAT-v1', status=status,
                construction=source['construction'], state_dimension=state, target_distance=str(distance),
                source_calibration_active_groups=source['active_groups'], requested_weights=list(weights),
                completed_weights=completed, incomplete_weight=pending,
                pending_weights=[r for r in weights if r not in {row['active_groups'] for row in completed}],
                source=dict(path=str(path), sha256=sha256(raw).hexdigest()), source_sha256=dict(sorted(pins.items())),
                source_marker_key=marker_key, marker_count=len(positives), bank_metadata=metadata,
                witness=dict(z=str(z), alpha=str(alpha), markers=[list(map(str, pair)) for pair in positives],
                             marker_log_odds_exact=[[str(F(v)) for v in pair] for pair in odds], grid=witness.get('grid', 128)),
                marker_conversion='Saved exact positive rationals retained; FLOAT log2(float(rational)) only',
                one_shared_float_bank=True, regional_normalizer_and_outer_prefix_fresh_per_weight=True,
                max_weights=8, max_seconds=max_seconds, elapsed_seconds=monotonic()-started,
                common_potential_used=False, old_numerical_bounds_reused=False, no_weight_interpolation=True,
                directed=False, selected_weight_certificate=False, whole_code_certificate=False,
                scope='Complete FLOAT sums only at explicitly completed weights; no intervening-weight or contiguous claim')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--weights', type=int, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=90.)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = run(args.source, args.weights, max_seconds=args.max_seconds)
    print(json.dumps(dict(status=result['status'], artifact=provider.search.write_new(args.output, result))), flush=True)


if __name__ == '__main__':
    main()
