"""Bounded FLOAT witness refinement for the XOR2 (a1,a2)=(1,16) cohort.

Every group keeps its constituent row set across all32 columns. The source
matching geometry and pair-multiplicity bins are counted outside the RS
support polynomial. This is a proposal, not a directed certificate.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import minimize

from hd_xor_geometry_strata import MultiplicityCapsFloat, geometry_stratum_counts
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_gate import inherited_witnesses, ordered_region_tables, sparse_log_count
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION


BOUNDS = ((.99, .9998), (-100., -20.))
MAX_EVALUATIONS = 32
MU_BINS = ((1, 8), (9, 14), (15, 15), (16, 16))


@lru_cache(None)
def source_counts():
    return matching_joint_counts(REGION_LENGTH, 17)


@lru_cache(None)
def raw_tables(mark):
    return pair_regions_screen(float(mark), 17)


def geometry_count(low, high):
    if (low, high) not in MU_BINS:
        raise ValueError('Expected one of the four fixed multiplicity bins')
    counts = geometry_stratum_counts(16)
    return sum(counts[mu] for mu in range(low, high+1))


def evaluate_at(source, mu_range, mark, log_potential):
    source, mu_range = tuple(source), tuple(mu_range)
    z, logp = F(mark), F(log_potential)
    if source not in source_counts() or not 0 < z < 1 or logp > 0:
        raise ValueError('Valid source class and Chernoff/potential witness required')
    geometries = geometry_count(*mu_range)
    tables = ordered_region_tables(raw_tables(z), 17, *source)
    caps = MultiplicityCapsFloat(tables, float(logp), 17)(mu_range[1])
    value = (log2(source_counts()[source])+log2(17)+log2(8)+log2(geometries)
             -log2(Q-1)-LIMIT*log2(float(z))-float(logp)
             +sparse_log_count(caps, 1, 16))
    return dict(source_matching_counts=list(source), mu_range=list(mu_range),
        mark=str(z), log2_potential=str(logp), source_supports=source_counts()[source],
        cohort_colorings=17, one_row_choices=8, two_row_geometries=geometries,
        failure_log2_screen=value)


def inherited_literal(source, mu_range):
    inherited, _ = inherited_witnesses()
    low, high = mu_range
    candidates = [(z, p) for lo, hi, z, p in inherited[tuple(source)]
                  if lo <= low and hi >= high]
    if len(candidates) != 1:
        raise ValueError('Expected one inherited literal for the requested bin')
    return candidates[0]


def fit_cell(source=(8, 0), mu_range=(9, 14), progress=None):
    initial = inherited_literal(source, mu_range)
    proposals, cache = [], {}

    def objective(point):
        key = tuple(float(value) for value in point)
        if key in cache:
            return cache[key]
        if len(proposals) >= MAX_EVALUATIONS:
            raise ArithmeticError('The bounded witness search exceeded its budget')
        row = evaluate_at(source, mu_range, *map(F, key))
        row['evaluation'] = len(proposals)+1
        proposals.append(row)
        cache[key] = row['failure_log2_screen']
        if progress:
            progress(dict(stage='bounded-cell-fit', **row))
        return cache[key]

    start = np.array(list(map(float, initial)))
    fit = minimize(objective, start, method='Nelder-Mead', bounds=BOUNDS,
        options=dict(maxfev=MAX_EVALUATIONS, maxiter=MAX_EVALUATIONS,
            initial_simplex=np.array([start, start+[-.001, 0], start+[0, 8.]]),
            xatol=1e-10, fatol=.01, adaptive=False))
    best = min(proposals, key=lambda row: row['failure_log2_screen'])
    point = [float(F(best[key])) for key in ('mark', 'log2_potential')]
    boundary = [name for name, value, (lo, hi) in zip(('mark', 'log2_potential'), point, BOUNDS)
                if min(abs(value-lo), abs(value-hi)) < 1e-10]
    return dict(best=best, proposals=proposals, evaluations=len(proposals),
        maximum_evaluations=MAX_EVALUATIONS, optimizer_message=str(fit.message),
        boundary_coordinates=boundary)


def screen(bank, *, all_sources=False, progress=None):
    """Test a fixed witness bank without any additional fitting."""
    sources = sorted(source_counts()) if all_sources else [(0, 0), (0, 8), (8, 0), (8, 8)]
    candidates = tuple(sorted({(F(z), F(p)) for z, p in bank}))
    entries = []
    for source in sources:
        for mu_range in MU_BINS:
            literals = set(candidates) | {inherited_literal(source, mu_range)}
            rows = [evaluate_at(source, mu_range, *literal) for literal in sorted(literals)]
            selected = min(rows, key=lambda row: row['failure_log2_screen'])
            selected['witnesses_evaluated'] = len(rows)
            entries.append(selected)
            if progress:
                progress(dict(stage='fixed-bank-screen', completed=len(entries),
                    total=len(sources)*len(MU_BINS), **selected))
    return dict(entries=entries, source_classes=len(sources), all_source_classes=len(source_counts()),
        complete_cohort=all_sources,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries])))


def metadata():
    return dict(construction=CONSTRUCTION, active_groups=17, cohort_counts=[1, 16, 0],
        field_order=Q, forbidden_weight=LIMIT, sparse_active_rows=33,
        sparse_support_floor=561, sparse_support_normalizer_exponent=528,
        pair_swaps=True, fixed_geometry_across_columns=True,
        projective_normalizations=1, endpoint_potential_factors=1,
        directed=False, proposal_requires_fresh_directed_replay=True,
        whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bank-from', type=Path)
    parser.add_argument('--all-sources', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite frozen proposal artifacts')
    started = perf_counter()
    progress = lambda row: print(json.dumps(row, allow_nan=False), flush=True)
    result = metadata()
    if args.bank_from:
        prior = json.loads(args.bank_from.read_text(encoding='utf-8'))
        if any(prior[key] != result[key] for key in ('construction', 'cohort_counts', 'active_groups')):
            parser.error('Witness source has a different construction or cohort')
        fit = prior['fit']
        result['witness_source_sha256'] = sha256(args.bank_from.read_bytes()).hexdigest()
    else:
        fit = fit_cell(progress=progress)
    result['fit'] = fit
    best = fit['best']
    result['screen'] = screen([(best['mark'], best['log2_potential'])],
        all_sources=args.all_sources, progress=progress)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(output=str(args.output.resolve()),
        bound=result['screen']['failure_log2_screen'], elapsed_seconds=result['elapsed_seconds'])), flush=True)
