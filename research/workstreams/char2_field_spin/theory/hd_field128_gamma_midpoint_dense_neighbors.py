"""FLOAT neighbors at original weights through8192; complement counts only."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import hd_field128_gamma_midpoint_neighbors as frozen


provider, mode, HERE = frozen.provider, frozen.mode, frozen.HERE


def run(source_path, weights, *, max_seconds=90.):
    weights = tuple(weights)
    if (not 1 <= len(weights) <= 8 or len(set(weights)) != len(weights)
            or any(type(r) is not int or not 2 <= r <= 8192 for r in weights)
            or not isfinite(max_seconds) or not 0 < max_seconds <= 90):
        raise ValueError('Expected1..8 distinct original weights2..8192, at most90 seconds')
    path = Path(source_path).resolve()
    raw = path.read_bytes()
    source = json.loads(raw)
    state, distance = source['state_dimension'], F(source['target_distance'])
    if source['construction'] != provider.gamma.construction(state):
        raise ValueError('Source construction mismatch')
    witness = source['witness']
    positives, marker_key = frozen.enrich.positive_markers(witness)
    calibration = source.get('active_groups')
    parent_pin = None
    if calibration is None:
        # An actual-bank metadata file may bind its literal proposal inline.
        proposal = source['witness_proposal']
        parent_positive, unused = frozen.enrich.positive_markers(proposal['witness'])
        if (parent_positive != positives or any(F(proposal['witness'][key]) != F(witness[key]) for key in ('z', 'alpha'))
                or F(proposal['new_target_distance']) != distance):
            raise ValueError('Embedded exact proposal differs from bank witness')
        calibration = proposal['source_calibration_weight']
        parent_pin = proposal['source']
    if type(calibration) is not int or not 1 <= calibration <= 8192:
        raise ValueError('Explicit source weight calibration required')
    z, alpha = F(witness['z']), F(witness['alpha'])
    odds = [tuple(log2(float(v)) for v in pair) for pair in positives]
    pins = frozen.enrich.source_pins()
    for filename in (Path(frozen.__file__), Path(__file__)):
        pins[filename.name] = sha256(filename.read_bytes()).hexdigest()
    for name, digest in source['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting source pin: '+name)
        pins[name] = digest
    def check():
        if path.read_bytes() != raw:
            raise RuntimeError('Source witness changed')
        if parent_pin and sha256(Path(parent_pin['path']).read_bytes()).hexdigest() != parent_pin['sha256']:
            raise RuntimeError('Embedded literal-proposal source changed')
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
            full = provider.alias.full_counts(4096, weight)
            regions = bank.regional(weight)
            if regions.shape != (len(full), 2, 2):
                raise ArithmeticError('Original-weight feasible regional range mismatch')
            counts = mode.MidpointCounts(4096, weight, max_singletons=64, budget=budget)
            row = mode.support_sum(4096, weight, regions, counts=counts, budget=budget)
            if row['profile_pairs_evaluated'] != len(full)**2:
                raise ArithmeticError('Incomplete original-weight profile sum')
            prefix = bank.constant+bank.alpha*bank.density*weight
            row.update(active_groups=weight, feasible_full_range=[int(full[0]), int(full[-1])],
                       canonical_count_weight=min(weight, 8192-weight), count_full_shift=max(0, weight-4096),
                       log2_outer_prefix=prefix, log2_upper_screen=prefix+row['log2_support_sum_screen'],
                       regional_alias_diagnostics=dict(bank.regional_diagnostics),
                       regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                       complete_selected_weight_screen=True, elapsed_seconds=monotonic()-started)
            completed.append(row)
            pending = None
            print(json.dumps(dict(weight=weight, score=row['log2_upper_screen'], full_range=row['feasible_full_range'],
                                  profiles=row['profile_pairs_evaluated'], largest=row['largest_profile'],
                                  diagnostics=row['regional_alias_diagnostics'], elapsed_seconds=row['elapsed_seconds'])), flush=True)
    except TimeoutError:
        status = 'timeout'
    check()
    return dict(schema='gamma-original-weight-dense-neighbors-FLOAT-v1', status=status,
                construction=source['construction'], state_dimension=state, target_distance=str(distance),
                source_calibration_active_groups=calibration, requested_weights=list(weights), completed_weights=completed,
                incomplete_weight=pending, pending_weights=[r for r in weights if r not in {row['active_groups'] for row in completed}],
                source=dict(path=str(path), sha256=sha256(raw).hexdigest()), embedded_literal_source=parent_pin,
                source_sha256=dict(sorted(pins.items())), source_marker_key=marker_key, marker_count=len(positives), bank_metadata=metadata,
                witness=dict(z=str(z), alpha=str(alpha), markers=[list(map(str, pair)) for pair in positives],
                             marker_log_odds_exact=[[str(F(v)) for v in pair] for pair in odds], grid=witness.get('grid', 128)),
                count_complement_only=True, original_weight_used_for_histograms_markers_matrices_and_prefix=True,
                one_shared_float_bank=True, max_weights=8, max_seconds=max_seconds, elapsed_seconds=monotonic()-started,
                common_potential_used=False, old_numerical_bounds_reused=False, no_weight_interpolation=True,
                directed=False, selected_weight_certificate=False, whole_code_certificate=False,
                scope='Complete FLOAT sums at listed original weights only; no complement transfer or intervening-weight claim')


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
