"""Bounded FLOAT gate for sixteen RS[16,8] rows and an affine two-sum mixer.

The group has 128 input and 256 output GF256 coordinates, as before. The
new labels are 0x06+GF16 inside AES GF256. Each singleton reaches all16
outputs; a fixed pair gives full support239/255 and each one-hole mask1/255.
Row sets remain fixed through16 columns. The old regional operator applies
pointwise to each fixed regional subset, but no old outer receipt transfers.
Only r17 sparse cohorts are screened, with one inherited global witness.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from gf8_packet_maps import MUL, INV, matmul, xor_sum
from hd_two_sum_affine_outer import column_action_float
from hd_xor_group_geometry import matching_joint_counts, subset_full_pair_range
from hd_xor_joint_columns import categories
from hd_xor_mixed_gate import sparse_log_count
from hd_xor_multigroup import LIMIT
from hd_xor_pair_route import pair_regions_screen


BASIS = (1, 0x0c, 0x50, 0xb0)
LABELS = tuple(6 ^ xor_sum(b for j, b in enumerate(BASIS) if i >> j & 1)
               for i in range(16))
CONSTRUCTION = 'two-sum-gf16-coset-g16-rs16x8-h4-swapped-pairs-k20-rate-half'
MARK, LOG_POTENTIAL = F('996248217/1000000000'), F('-127/2')
GROUPS, REGION_LENGTH = 17, 8192


@lru_cache(None)
def matrix():
    return tuple(tuple(int(i == j) ^ 1 ^ MUL[a][b] for j, b in enumerate(LABELS))
                 for i, a in enumerate(LABELS))


def direct(values):
    if len(values) != 16 or any(type(x) is not int or not 0 <= x < 256 for x in values):
        raise ValueError('Sixteen AES-field bytes required')
    s = xor_sum(values)
    t = xor_sum(MUL[u][x] for u, x in zip(LABELS, values))
    return tuple(x ^ s ^ MUL[u][t] for u, x in zip(LABELS, values))


def circuit(values):
    """Scalar oracle for the shared-sum eight-fixed-product circuit."""
    if len(values) != 16 or any(type(x) is not int or not 0 <= x < 256 for x in values):
        raise ValueError('Sixteen AES-field bytes required')
    pairs = tuple(values[2*i] ^ values[2*i+1] for i in range(8))
    quarters = tuple(pairs[2*i] ^ pairs[2*i+1] for i in range(4))
    halves = (quarters[0] ^ quarters[1], quarters[2] ^ quarters[3])
    s = halves[0] ^ halves[1]
    a0 = xor_sum(values[1::2])
    a1 = xor_sum(pairs[1::2])
    a2 = quarters[1] ^ quarters[3]
    a3 = halves[1]
    t = MUL[6][s] ^ a0 ^ MUL[0x0c][a1] ^ MUL[0x50][a2] ^ MUL[0xb0][a3]
    basis_t = (t, MUL[0x0c][t], MUL[0x50][t], MUL[0xb0][t])
    corrections = [s ^ MUL[6][t]]
    for scaled in basis_t:
        corrections += [x ^ scaled for x in corrections]
    return tuple(x ^ correction for x, correction in zip(values, corrections))


def sparse_law(arity):
    if type(arity) is not int or not 0 <= arity <= 2:
        raise ValueError('Only sparse arities0,1,2 are covered')
    if arity == 0:
        return {0: F(1)}
    if arity == 1:
        return {65535: F(1)}
    return {65535: F(239, 255), **{65535 ^ (1 << i): F(1, 255) for i in range(16)}}


def exact_local_checks():
    h = matrix()
    pairs = list(combinations(range(16), 2))
    if set(LABELS).intersection(INV[x] for x in LABELS):
        raise ArithmeticError('Coset meets its inverse')
    if matmul(h, h) != tuple(tuple(int(i == j) for j in range(16)) for i in range(16)):
        raise ArithmeticError('Mixer is not an involution')
    if not all(x for row in h for x in row):
        raise ArithmeticError('Singleton support is not full')
    for i, j in pairs:
        for k, ell in pairs:
            if not MUL[h[k][i]][h[ell][j]] ^ MUL[h[k][j]][h[ell][i]]:
                raise ArithmeticError('A two-column minor is zero')
        counts = Counter(sum(bool(h[k][i] ^ MUL[h[k][j]][ratio]) << k for k in range(16))
                         for ratio in range(1, 256))
        if counts != {mask: int(p*255) for mask, p in sparse_law(2).items()}:
            raise ArithmeticError('Projective pair histogram disagrees')
    return dict(nonzero_entries=256, nonzero_two_column_minors=14400,
                projective_pairs_checked=30600, involution=True)


def ordered_regions(raw, total, a, b):
    """Sixteen chronological regions, alternating the same two matchings."""
    if any(type(x) is not int or not 0 <= x <= total//2 for x in (a, b)):
        raise ValueError('Invalid source matching counts')
    families = []
    for full in (a, b):
        family = []
        for m in range(total+1):
            low, high = subset_full_pair_range(total, full, m)
            family.append(np.maximum.reduce([raw[m, f] for f in range(low, high+1)]))
        families.append(family)
    return np.array([families[i % 2] for i in range(16)])


def screen(*, all_sources=False, budget_seconds=120, progress=None):
    if not 0 < budget_seconds <= 120:
        raise ValueError('Use a positive total budget at most120 seconds')
    started = perf_counter()
    counts = matching_joint_counts(REGION_LENGTH, GROUPS)
    if len(counts) != 81 or sum(counts.values()) != comb(REGION_LENGTH, GROUPS):
        raise ArithmeticError('Source classes do not partition source supports')
    selected = sorted(counts) if all_sources else [(0, 0), (0, 8), (8, 0), (8, 8)]
    raw = pair_regions_screen(float(MARK), GROUPS)
    logp = float(LOG_POTENTIAL)
    common = -log2(255)-LIMIT*log2(float(MARK))-logp
    sources, exhausted = [], False
    for source in selected:
        if perf_counter()-started > budget_seconds:
            exhausted = True
            break
        regions = ordered_regions(raw, GROUPS, *source)
        caps = {}
        for single, double in categories(GROUPS):
            value = column_action_float(regions, logp, single, double)
            caps[single, double] = float(np.max(value-[0., logp]))
        entries = []
        for n1 in range(18):
            n2 = GROUPS-n1
            value = (log2(counts[source])+log2(comb(GROUPS, n1))+4*n1+n2*log2(120)
                     +common+sparse_log_count(caps, n1, n2, n=16, k=8))
            entries.append(dict(cohort_counts=[n1, n2, 0], failure_log2_screen=value,
                                source_supports=counts[source], source_colorings=comb(GROUPS, n1),
                                singleton_row_choices=16**n1, pair_row_choices=120**n2,
                                sparse_support_floor=9*(n1+2*n2),
                                sparse_normalizer_exponent=8*(n1+2*n2)))
        row = dict(source_matching_counts=list(source), source_supports=counts[source], entries=entries,
                   failure_log2_screen=float(np.logaddexp2.reduce([e['failure_log2_screen'] for e in entries])))
        sources.append(row)
        if progress:
            progress({k: v for k, v in row.items() if k != 'entries'})
    complete_selected = len(sources) == len(selected) and not exhausted
    return dict(construction=CONSTRUCTION, field_order=256, field_polynomial='0x11b', labels=LABELS,
                active_groups=GROUPS, total_columns=16, constituent_dimension=8,
                forbidden_weight=LIMIT, mark=str(MARK), log2_potential=str(LOG_POTENTIAL),
                witness_banks=1, fitted_new_witnesses=0, raw_regional_banks=1,
                source_classes=len(sources), expected_source_classes=81, requested_source_classes=len(selected),
                complete_source_coverage=complete_selected and all_sources,
                complete_sparse_cohort_coverage=complete_selected and all_sources,
                source_count_partition=comb(REGION_LENGTH, GROUPS),
                cohort_order=['a1', 'a2', 'a3-through-a16'], selected_cells=18*len(sources),
                fixed_rows_across_columns=True, fixed_source_geometry=True,
                local_mask_mixture_exact=True, regional_matching_count_envelopes_used=True,
                column_occupancy_resource_maximum_used=False, pair_swaps=True,
                potential_quotient_after_hole_mixture=True, projective_normalizations=1,
                endpoint_potential_factors=1, directed=False, whole_code_certificate=False,
                old_certificates_transferred=False, budget_exhausted=exhausted,
                budget_seconds=budget_seconds, elapsed_seconds=perf_counter()-started,
                failure_log2_screen=float(np.logaddexp2.reduce([s['failure_log2_screen'] for s in sources])) if sources else None,
                sources=sources)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all-sources', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a research artifact')
    checks = exact_local_checks()
    result = screen(all_sources=args.all_sources, progress=lambda row: print(json.dumps(row), flush=True))
    result['exact_local_checks'] = checks
    result['source_sha256'] = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in (Path(__file__).name, 'hd_two_sum_affine_outer.py', 'hd_xor_mixed_gate.py',
                     'hd_xor_pair_route.py', 'hd_xor_group_geometry.py')}
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'sources'}, allow_nan=False))
