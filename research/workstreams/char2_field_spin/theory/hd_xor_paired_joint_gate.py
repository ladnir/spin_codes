"""Fixed-r17 XOR2 proof gate retaining both matching counts and pair geometry.

Unlike hd_xor_paired_gate's max(a,b) rectangles, each source stratum keeps
the exact ordered pair (a,b). The regional output mask is a subset of this
fixed source support, so both lower and upper full-pair bounds apply.
The saved schedule uses literal Chernoff witnesses per constituent stratum.
It does not change the construction or choose a witness after observing
setup randomness. No mixed-message or whole-code claim is made.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_group_geometry import matching_joint_counts, subset_full_pair_range
from hd_xor_geometry_strata import (geometry_stratum_counts, MultiplicityCapsFloat,
                                     MultiplicityCapsArb)
from hd_xor_pair_route import pair_regions_screen, arb_pair_regions
from hd_xor_joint_columns import joint_log_count, joint_count_arb, arb_upper_max
from hd_xor_paired_gate import CONSTRUCTION
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT


GROUPS = 17
SCHEDULE = (
    (1, 8, '996248217/1000000000', '-127/2'),
    (9, 14, '996248217/1000000000', '-127/2'),
    (15, 15, '996248217/1000000000', '-127/2'),
    (16, 16, '996548217/1000000000', '-127/2'),
    (17, 17, '996598217/1000000000', '-127/2'),
)


def ordered_tables(raw, a, b, directed=False):
    """Dominate every regional subset of a fixed17-group source support."""
    if any(type(x) is not int or not 0 <= x <= GROUPS//2 for x in (a, b)):
        raise ValueError('Invalid r17 matching counts')
    tables = []
    for full in (a, b):
        family = []
        for m in range(GROUPS+1):
            low, high = subset_full_pair_range(GROUPS, full, m)
            if directed:
                entries = list(raw[m, low])
                for f in range(low+1, high+1):
                    entries = [arb_upper_max(x, y) for x, y in zip(entries, raw[m, f])]
                family.append(tuple(entries))
            else:
                family.append(np.maximum.reduce([raw[m, f] for f in range(low, high+1)]))
        tables.append(family)
    return [tables[i % 2] for i in range(8)] if directed else np.array([tables[i % 2] for i in range(8)])


def checked_schedule(schedule):
    covered = set()
    result = []
    for low, high, mark, potential in schedule:
        if type(low) is not int or type(high) is not int or not 1 <= low <= high <= GROUPS:
            raise ValueError('Invalid constituent-pair stratum')
        values = set(range(low, high+1))
        if values & covered:
            raise ValueError('Constituent-pair strata overlap')
        covered |= values
        mark, potential = F(mark), F(potential)
        if not 0 < mark < 1 or potential > 0:
            raise ValueError('Invalid literal proof witness')
        result.append((low, high, mark, potential))
    return result, covered == set(range(1, GROUPS+1))


def literal_schedules(assigned, schedule=SCHEDULE, witness_rows=None):
    """Validate fixed witnesses; saved numerical bounds/counts are never trusted."""
    if witness_rows is None:
        checked, complete = checked_schedule(schedule)
        return {source: checked for source in assigned}, complete
    schedules = {}
    for row in witness_rows:
        source = tuple(row['source_matching_counts'])
        if len(source) != 2 or any(type(x) is not int or not 0 <= x <= 8 for x in source):
            raise ValueError('Invalid saved source stratum')
        low, high = row['mu_range']
        schedules.setdefault(source, []).append((low, high, row['mark'], row['log2_potential']))
    if not set(assigned) <= schedules.keys():
        raise ValueError('Saved witnesses omit an assigned source stratum')
    checked = {source: checked_schedule(rows) for source, rows in schedules.items()}
    return {source: checked[source][0] for source in assigned}, all(checked[source][1] for source in assigned)


def run(*, precision=None, source_pairs=None, schedule=SCHEDULE, witness_rows=None, progress=None):
    if precision is not None and (type(precision) is not int or precision < 128):
        raise ValueError('Directed replay requires at least128 bits')
    counts = matching_joint_counts(REGION_LENGTH, GROUPS)
    pair_counts = geometry_stratum_counts(GROUPS)
    assigned = set(counts) if source_pairs is None else set(map(tuple, source_pairs))
    if not assigned or not assigned <= set(counts):
        raise ValueError('Select a nonempty subset of feasible source classes')
    schedules, complete_mu = literal_schedules(assigned, schedule, witness_rows)
    complete = complete_mu and assigned == set(counts)
    started = perf_counter()
    raw_tables, rows = {}, []
    marks = {z for rows in schedules.values() for _, _, z, _ in rows}
    common = dict(construction=CONSTRUCTION, active_groups=GROUPS, active_rows_per_group=2,
        source_dimension=2**20, output_length=2**21, forbidden_weight=LIMIT,
        pair_swaps=True, exact_ordered_source_matching_counts=True,
        constituent_geometry_fixed_across_columns=True,
        complete_geometry_partition=complete, whole_code_certificate=False,
        support_scope='fixed-r17 homogeneous two-constituent class',
        projective_normalizations=1, individual_support_floors_relaxed=True)
    if precision is None:
        for z in marks:
            raw_tables[z] = pair_regions_screen(float(z), GROUPS)
        for a, b in sorted(assigned):
            # Discard each source class's resource cache after use. Its bounds
            # can differ from every other class; no max/average swap is used.
            evaluators = {}
            for low, high, z, logp in schedules[a, b]:
                key = z, logp
                if key not in evaluators:
                    evaluators[key] = MultiplicityCapsFloat(
                        ordered_tables(raw_tables[z], a, b), float(logp), GROUPS)
                caps = evaluators[key](high)
                geometries = sum(pair_counts[mu] for mu in range(low, high+1))
                value = (log2(counts[a, b])+log2(geometries)-log2(Q-1)
                         -LIMIT*log2(float(z))-float(logp)+joint_log_count(caps, GROUPS))
                row = dict(source_matching_counts=[a, b], mu_range=[low, high],
                           source_supports=counts[a, b], constituent_pairs=geometries,
                           mark=str(z), log2_potential=str(logp), failure_log2_screen=value)
                rows.append(row)
                if progress:
                    progress(row)
        return dict(common, directed=False, entries=rows,
                    failure_log2_screen=float(np.logaddexp2.reduce([x['failure_log2_screen'] for x in rows])),
                    elapsed_seconds=perf_counter()-started)

    from flint import arb, ctx
    from hd_shifted_coverage import rational_ceiling
    with ctx.workprec(precision):
        for z in marks:
            raw_tables[z] = arb_pair_regions(z, GROUPS)
        total = arb(0)
        for a, b in sorted(assigned):
            evaluators = {}
            for low, high, z, logp in schedules[a, b]:
                key = z, logp
                potential = arb(2)**(arb(logp.numerator)/logp.denominator)
                if key not in evaluators:
                    evaluators[key] = MultiplicityCapsArb(
                        ordered_tables(raw_tables[z], a, b, True), potential, GROUPS, labelled_regions=True)
                caps = evaluators[key](high)
                geometries = sum(pair_counts[mu] for mu in range(low, high+1))
                value = (counts[a, b]*geometries*joint_count_arb(caps, GROUPS)
                         *(arb(z.numerator)/z.denominator)**(-LIMIT)/(Q-1)/potential)
                if not value > 0 or not value.is_finite():
                    raise ArithmeticError('Invalid geometry contribution')
                logarithm = value.log()/arb(2).log()
                row = dict(source_matching_counts=[a, b], mu_range=[low, high],
                           source_supports=counts[a, b], constituent_pairs=geometries,
                           mark=str(z), log2_potential=str(logp),
                           failure_log2_upper=str(rational_ceiling(logarithm)))
                rows.append(row)
                total += value
                if progress:
                    progress(row)
        logarithm = total.log()/arb(2).log()
        return dict(common, directed=True, entries=rows, precision_bits=precision,
                    failure_log2_upper=str(rational_ceiling(logarithm)),
                    homogeneous_class_reaches_40_bits=bool(complete and logarithm < -40),
                    elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int)
    parser.add_argument('--source-pairs', nargs=2, type=int, action='append')
    parser.add_argument('--witnesses', type=Path,
                        help='Read literal class witnesses from entries; recompute every bound')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Frozen gate artifacts are never overwritten')
    source_hash = sha256(Path(__file__).read_bytes()).hexdigest()
    witness_bytes = args.witnesses.read_bytes() if args.witnesses else None
    witness_rows = json.loads(witness_bytes)['entries'] if witness_bytes else None
    result = run(precision=args.precision, source_pairs=args.source_pairs, witness_rows=witness_rows,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = source_hash
    if witness_bytes:
        result['witness_sha256'] = sha256(witness_bytes).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, allow_nan=False, indent=2)
        target.write('\n')
    print(json.dumps(dict(path=str(args.output.resolve()), directed=result['directed'],
                          bound=result.get('failure_log2_upper', result.get('failure_log2_screen')),
                          elapsed_seconds=result['elapsed_seconds'])), flush=True)
