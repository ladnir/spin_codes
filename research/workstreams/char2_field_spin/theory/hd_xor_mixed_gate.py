"""Bounded FLOAT mixed-cohort XOR2 paired-route proof proposals.

Cohorts are (n1,n2,nH): one active constituent, two active constituents,
and between three and eight active constituents per group. Dense groups
use the proved product domination, summed over row sets ONCE per group.
Its iid nonzero probability is eta=p*255/256. Average those occupancies
inside each labelled regional operator before applying sparse column caps.
The sparse RS support polynomial retains J>=17*(n1+2*n2), and the fixed
two-row geometry multiplicity is counted once outside its 32-column power.

This is an alpha=1 FLOAT gate, not a directed or whole-code certificate.
Only inherited literal marks/potentials are tried; no optimizer is called.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, factorial, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_column_screen import group_density_bits
from hd_xor_geometry_strata import MultiplicityCapsFloat, geometry_stratum_counts
from hd_xor_group_geometry import matching_joint_counts, subset_full_pair_range
from hd_xor_pair_route import pair_regions_screen
from hd_xor_paired_gate import CONSTRUCTION
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT


SOURCE = Path(__file__).with_name('hd_xor_paired_joint_gate17_potential_screen.json')
DENSE_P = '9815016183/10000000000'
PROFILES = {'one-plus-pairs': (1, 16, 0), 'dense-plus-pairs': (0, 16, 1),
            'dense-one': (0, 0, 1), 'dense-seventeen': (0, 0, 17)}


def profile_counts(profile):
    if len(profile) != 3 or any(type(n) is not int or n < 0 for n in profile) or not 1 <= sum(profile) <= 17:
        raise ValueError('Three nonnegative cohort counts with total1..17 required')
    return tuple(profile)


def ordered_region_tables(raw, total_groups, a, b):
    """Universal subset bounds conditional on the fixed whole-source geometry."""
    if any(type(x) is not int or not 0 <= x <= total_groups//2 for x in (a, b)):
        raise ValueError('Invalid ordered source matching counts')
    tables = []
    for full in (a, b):
        family = []
        for m in range(total_groups+1):
            low, high = subset_full_pair_range(total_groups, full, m)
            family.append(np.maximum.reduce([raw[m, f] for f in range(low, high+1)]))
        tables.append(family)
    return np.array([tables[position % 2] for position in range(8)])


def dense_density_log2(probability):
    return float(np.logaddexp2.reduce([group_density_bits(active, probability) for active in range(3, 9)]))


def integrate_dense_regions(tables, sparse_groups, dense_groups, probability):
    """Average dense surrogate occupancy before taking any column potential cap.

    The source coloring need not be sampled: the whole-source subset table
    dominates every coloring and every sparse output subset of that coloring.
    """
    tables = np.asarray(tables, dtype=float)
    if (tables.shape != (8, sparse_groups+dense_groups+1, 2, 2)
            or not 0 < probability <= 1 or min(sparse_groups, dense_groups) < 0):
        raise ValueError('Invalid regional table or dense cohort')
    eta = probability*(Q-1)/Q
    result = np.full((8, sparse_groups+1, 2, 2), -np.inf)
    for occupied in range(dense_groups+1):
        weight = (log2(comb(dense_groups, occupied))+occupied*log2(eta)
                  +(dense_groups-occupied)*log2(1-eta))
        result = np.logaddexp2(result, tables[:, occupied:occupied+sparse_groups+1]+weight)
    return result


def mixed_column_log_coefficients(caps, n1, n2):
    """Exact binary source-support multiplicities for fixed sparse row sets."""
    polynomial = np.full(n1+2*n2+1, -np.inf)
    for one in range(n1+1):
        for single in range(n2+1):
            for double in range(n2-single+1):
                degree = one+single+2*double
                multiplicity = comb(n1, one)*comb(n2, single)*comb(n2-single, double)*2**single
                polynomial[degree] = np.logaddexp2(polynomial[degree],
                    log2(multiplicity)+caps[one+single, double])
    return polynomial


def sparse_log_count(caps, n1, n2, n=32, k=16, q=256):
    if not 0 < k <= n or q < 2 or min(n1, n2) < 0:
        raise ValueError('Invalid dimensions or cohort counts')
    active_rows = n1+2*n2
    polynomial = mixed_column_log_coefficients(caps, n1, n2)
    polynomial += np.arange(active_rows+1)*log2(q-1)
    expanded = np.array([0.])
    for _ in range(n):
        following = np.full(len(expanded)+active_rows, -np.inf)
        for degree, coefficient in enumerate(polynomial):
            following[degree:degree+len(expanded)] = np.logaddexp2(
                following[degree:degree+len(expanded)], expanded+coefficient)
        expanded = following
    floor = (n-k+1)*active_rows
    return float(np.logaddexp2.reduce(expanded[floor:])-(n-k)*active_rows*log2(q-1))


def source_coloring_count(profile):
    profile = profile_counts(profile)
    return factorial(sum(profile))//np.prod([factorial(n) for n in profile], dtype=object)


def inherited_witnesses(path=SOURCE):
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    if data['construction'] != CONSTRUCTION or data['active_groups'] != 17 or data['directed']:
        raise ValueError('Expected the fixed r17 FLOAT proposal source')
    rows = {}
    for row in data['entries']:
        source = tuple(row['source_matching_counts'])
        low, high = row['mu_range']
        z, potential = F(row['mark']), F(row['log2_potential'])
        if not 1 <= low <= high <= 17 or not 0 < z < 1 or potential > 0:
            raise ValueError('Invalid inherited literal witness')
        rows.setdefault(source, []).append((low, high, z, potential))
    if set(rows) != set(matching_joint_counts(REGION_LENGTH, 17)):
        raise ValueError('The inherited source classes are incomplete')
    for group in rows.values():
        if sorted(mu for low, high, _, _ in group for mu in range(low, high+1)) != list(range(1, 18)):
            raise ValueError('Inherited multiplicity strata must partition1..17')
    return rows, sha256(Path(path).read_bytes()).hexdigest()


def run(profile, *, all_sources=False, source_path=SOURCE, progress=None):
    n1, n2, nH = profile_counts(profile)
    r, sparse = sum(profile), n1+n2
    p = float(F(DENSE_P))
    density = dense_density_log2(p) if nH else 0.
    counts = matching_joint_counts(REGION_LENGTH, r)
    selected = sorted(counts) if all_sources else sorted({(0, 0), (r//2, 0), (0, r//2), (r//2, r//2)} & set(counts))
    witness_rows, source_hash = inherited_witnesses(source_path)
    pair_counts = geometry_stratum_counts(n2) if n2 else {0: 1}
    raw_tables, entries = {}, []
    started = perf_counter()
    for a, b in selected:
        inherited = witness_rows[a, b]
        if n2:
            scheduled = [(low, min(high, n2), z, potential) for low, high, z, potential in inherited if low <= n2]
        else:
            # With no two-row cohort, all inherited literals are candidate
            # witnesses for this one class, NOT distinct union contributions.
            scheduled = [(0, 0, z, potential) for z, potential in sorted({(z, p) for _, _, z, p in inherited})]
        candidates, evaluators = [], {}
        for low, high, z, potential in scheduled:
            if z not in raw_tables:
                raw_tables[z] = pair_regions_screen(float(z), r)
            key = z, potential
            if key not in evaluators:
                original = ordered_region_tables(raw_tables[z], r, a, b)
                tables = integrate_dense_regions(original, sparse, nH, p)
                evaluators[key] = MultiplicityCapsFloat(tables, float(potential), sparse) if sparse else tables
            if sparse:
                caps = evaluators[key](max(1, high))
            else:
                vector = np.array([0., float(potential)])
                for table in reversed(evaluators[key]):
                    vector = np.logaddexp2(table[0, :, 0]+vector[0], table[0, :, 1]+vector[1])
                caps = {(0, 0): float(np.max(vector-[0., float(potential)]))}
            geometries = sum(pair_counts[mu] for mu in range(low, high+1))
            value = (log2(counts[a, b])+log2(source_coloring_count(profile))+3*n1+log2(geometries)
                     +nH*density-log2(Q-1)-LIMIT*log2(float(z))-float(potential)
                     +sparse_log_count(caps, n1, n2))
            row = dict(source_matching_counts=[a, b], mu_range=[low, high],
                       source_supports=counts[a, b], cohort_colorings=int(source_coloring_count(profile)),
                       one_row_choices=8**n1, two_row_geometries=geometries,
                       mark=str(z), log2_potential=str(potential), failure_log2_screen=value)
            candidates.append(row)
            if n2 and progress:
                progress(row)
        chosen = candidates if n2 else [min(candidates, key=lambda row: row['failure_log2_screen'])]
        entries.extend(chosen)
        if not n2 and progress:
            progress(chosen[0])
    return dict(construction=CONSTRUCTION, cohort_counts=[n1, n2, nH], cohort_order=['a1', 'a2', 'a3-through-a8'],
        active_groups=r, dense_probability=DENSE_P, dense_group_density_log2=density,
        sparse_active_rows=n1+2*n2, sparse_support_floor=17*(n1+2*n2), alpha=1,
        witness_source_sha256=source_hash, source_classes=len(selected), all_source_classes=len(counts),
        complete_class=bool(set(selected) == set(counts)), exact_source_color_count=True,
        fixed_geometry_across_columns=True, projective_normalizations=1,
        directed=False, proposal_requires_fresh_directed_replay=True, whole_code_certificate=False,
        failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in entries])),
        entries=entries, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', choices=PROFILES, required=True)
    parser.add_argument('--all-sources', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite frozen diagnostic artifacts')
    result = run(PROFILES[args.profile], all_sources=args.all_sources,
                 progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps(dict(output=str(args.output.resolve()), bound=result['failure_log2_screen'],
                          elapsed_seconds=result['elapsed_seconds'], complete_class=result['complete_class'])), flush=True)
