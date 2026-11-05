"""Rebuild a FLOAT bank and retain every selected regional-entry marker."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import monotonic

import numpy as np

import hd_field128_gamma_midpoint_enrich as enrich


provider, mode, HERE = enrich.provider, enrich.mode, enrich.HERE


def run(source_path, weight):
    if type(weight) is not int or not 2 <= weight < 4096:
        raise ValueError('Expected weight2..4095')
    path = Path(source_path).resolve()
    raw = path.read_bytes()
    source = json.loads(raw)
    if source.get('active_groups') != weight:
        raise ValueError('Source must explicitly match the selected weight')
    state, distance = source['state_dimension'], F(source['target_distance'])
    if source['construction'] != provider.gamma.construction(state):
        raise ValueError('Construction mismatch')
    pins = enrich.source_pins()
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    for name, digest in source['source_sha256'].items():
        if name in pins and pins[name] != digest:
            raise ValueError('Source conflict: '+name)
        pins[name] = digest
    def check():
        if path.read_bytes() != raw:
            raise RuntimeError('Source artifact changed')
        for name, digest in pins.items():
            if sha256((HERE/name).read_bytes()).hexdigest() != digest:
                raise RuntimeError('Source changed: '+name)
    check()
    witness = source['witness']
    positives, key = enrich.positive_markers(witness)
    odds = [tuple(log2(float(value)) for value in row) for row in positives]
    z, alpha = F(witness['z']), F(witness['alpha'])
    started, budget = monotonic(), provider.old.Budget(30.)
    coeff = provider.AliasCoefficients(float(z), float(alpha), state=state, _budget=budget)
    bank = provider.AliasBankScreen(float(z), float(alpha), odds, 0., state=state, distance=distance,
                grid=witness.get('grid', 128), floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
    regions = bank.regional(weight)
    full = np.arange(weight//2+1)
    single, empty = weight-2*full, 4096-weight+full
    factorial = provider.alias.LOG_FACTORIAL
    normalizer = factorial[4096]-factorial[empty]-factorial[single]-factorial[full]
    best, winners = np.full((4, len(full)), np.inf), np.full((4, len(full)), -1, dtype=int)
    reliability = np.zeros(best.shape, dtype=bool)
    def candidate(index):
        budget.check()
        marker, table = bank.odds[index], bank.tables[index]
        reliable = table.reliable[:, single % bank.grid, full % bank.grid]
        numerator = np.where(reliable, table.logs[:, single % bank.grid, full % bank.grid], table.positive_numerator[:, None])
        return numerator-normalizer-single*marker[0]-full*marker[1], reliable
    for index in range(len(bank.tables)):
        values, reliable = candidate(index)
        improved = values < best
        winners[improved], reliability[improved] = index, reliable[improved]
        np.minimum(best, values, out=best)
    if not np.array_equal(best.T.reshape(regions.shape), regions) or np.any(winners < 0):
        raise ArithmeticError('Winner reconstruction differs from frozen regional caps')
    retained = sorted(set(map(int, winners.ravel())))
    pruned = np.full(best.shape, np.inf)
    for index in retained:
        np.minimum(pruned, candidate(index)[0], out=pruned)
    pruned = pruned.T.reshape(regions.shape)
    if not np.array_equal(regions, pruned):
        raise ArithmeticError('Pruning changed a regional cap')
    counts = mode.MidpointCounts(4096, weight, max_singletons=64, budget=budget)
    original_result = mode.support_sum(4096, weight, regions, counts=counts, budget=budget)
    pruned_result = mode.support_sum(4096, weight, pruned, counts=counts, budget=budget)
    prefix = bank.constant+bank.alpha*bank.density*weight
    for result in (original_result, pruned_result):
        result.update(log2_outer_prefix=prefix, log2_upper_screen=prefix+result['log2_support_sum_screen'])
    if original_result['log2_upper_screen'] != pruned_result['log2_upper_screen']:
        raise ArithmeticError('Complete support score changed after pruning')
    source_indices = [odds.index(bank.odds[index]) for index in retained]
    check()
    return dict(schema='gamma-selected-entry-winner-pruned-proposal-FLOAT-v1', active_groups=weight,
                construction=source['construction'], state_dimension=state, target_distance=str(distance),
                source=dict(path=str(path), sha256=sha256(raw).hexdigest()), source_sha256=dict(sorted(pins.items())),
                source_marker_key=key, source_marker_count=len(positives), float_unique_marker_count=len(bank.odds),
                marker_count=len(retained), retained_bank_marker_indices=retained, retained_source_marker_indices=source_indices,
                selected_marker_indices_by_entry=winners.tolist(), selected_reliable_entries=int(reliability.sum()),
                selected_positive_fallback_entries=int((~reliability).sum()), all_regional_entries_retained=True,
                regional_caps_bitwise_equal=True, complete_scores_equal=True,
                regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                rebuilt_original=original_result, rebuilt_pruned=pruned_result,
                witness=dict(z=str(z), alpha=str(alpha), markers=[list(map(str, positives[i])) for i in source_indices],
                             grid=bank.grid, state_dimension=state,
                             marker_log_odds_exact=[[str(F(v)) for v in odds[i]] for i in source_indices]),
                source_log_witness_literals=witness.get('marker_log_odds_exact'),
                marker_selection='Union of first strict minimum winners over all4 entries and all feasible full counts; positive fallback included',
                old_numerical_bounds_reused=False, initial_bank_rebuilt=True, max_seconds=30,
                elapsed_seconds=monotonic()-started, directed=False, selected_weight_certificate=False,
                whole_code_certificate=False, proposal_only=True,
                scope='Exact positive marker subset proposal; cap equality is checked only in FLOAT at the listed weight')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = run(args.source, args.r)
    pin = provider.search.write_new(args.output, result)
    print(json.dumps(dict(artifact=pin, markers=result['marker_count'],
                          log2_upper_screen=result['rebuilt_pruned']['log2_upper_screen'])), flush=True)


if __name__ == '__main__':
    main()
