"""Exact density moments modulo a proved subgroup of Hlambda symmetries.

Labels are vectors in F2^4. Let L permute their four bits, optionally
followed by tau(i)=i+parity(i)*(1,1,1,1). Both maps are orthogonal.
Consequently Hlambda[Li,Lj]=Hlambda[i,j], including its diagonal identity.
For a mapped complement O, the code ker Hlambda[O,:] is therefore carried
to ker Hlambda[LO,:] by the SAME coordinate permutation L. Independent
nonzero output diagonals preserve this equivalence. Every density moment
agrees on each orbit. No larger affine symmetry is assumed. This only
reuses local moments; it does not quotient scheduled whole-word chronology
or remove the original labelled rowset multiplicities from a later bound.

This module keeps an explicit representative-to-complement permutation.
It computes the representative's exact support counts once, then computes
several exact rational moments. Complement sizes0..4 mean arities12..16.
There are no regional, source, RS, or whole-code certificate factors.
"""
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

import hd_lambda_renyi_support as support

DIRECTORY = Path(__file__).resolve().parent
PINS = {
    'hd_lambda_renyi_support.py': 'c4662bb016f248faa7ffcc1972e01349bc66306c40762e12df9245fc04a3c03e',
    'hd_binary_outer_relations.py': 'ec1c13cd4a685760983b33ebd8e722efc3735f6335ed55bfe444aea01b1a20f2',
}
MAX_MOMENT = 16


def _binary_entry(i, j):
    return int(i == j) ^ ((i & j).bit_count() & 1)


@lru_cache(maxsize=1)
def orthogonal_maps():
    """The explicit48-map subgroup; each tuple maps old labels to new ones."""
    maps = []
    for bit_order in permutations(range(4)):
        moved = tuple(sum(((i >> b) & 1) << bit_order[b] for b in range(4)) for i in range(16))
        for use_tau in (False, True):
            image = tuple(j ^ (15 if use_tau and j.bit_count() & 1 else 0) for j in moved)
            if sorted(image) != list(range(16)):
                raise ArithmeticError('Claimed label map is not a permutation')
            # Equality of lambda plus the binary entry is symbolic in lambda.
            if any(_binary_entry(image[i], image[j]) != _binary_entry(i, j)
                   for i in range(16) for j in range(16)):
                raise ArithmeticError('Claimed map does not preserve full Hlambda')
            maps.append(image)
    result = tuple(sorted(set(maps)))
    if len(result) != 48:
        raise ArithmeticError('Expected the48 explicit orthogonal maps')
    return result


def complement_orbits(complement_size):
    """Partition all labelled complements; members contain actual witnesses.

    Each returned record has representative and members. A member is the
    pair (complement, permutation), where permutation maps representative
    to complement. All label tuples are sorted; the coordinate map is not.
    """
    if type(complement_size) is not int or not 0 <= complement_size <= 4:
        raise ValueError('Complement size must be0..4')
    remaining = set(combinations(range(16), complement_size))
    result = []
    maps = orthogonal_maps()
    while remaining:
        representative = min(remaining)
        members = {}
        for image in maps:
            target = tuple(sorted(image[i] for i in representative))
            members.setdefault(target, image)
        if not set(members) <= remaining:
            raise ArithmeticError('Computed subgroup orbits overlap')
        remaining.difference_update(members)
        result.append(dict(representative=representative,
                           members=tuple(sorted(members.items()))))
    if sum(len(row['members']) for row in result) != comb(16, complement_size):
        raise ArithmeticError('Labelled complement partition is incomplete')
    return tuple(result)


def _orders(orders):
    orders = tuple(orders)
    if (not 1 <= len(orders) <= 8 or len(set(orders)) != len(orders)
            or any(type(p) is not int or not 1 <= p <= MAX_MOMENT for p in orders)):
        raise ValueError('Use one to eight distinct integer moment orders1..16')
    return tuple(sorted(orders))


def exact_moments(counts, codimension, *, orders=(3, 4), q=256):
    """Return p->Fraction for normalized density, using exact integer sums.

    For n coordinates and codimension b, the formula is
    I_p=q^(p*b-n)*sum_S n_C(S)^p/(q-1)^((p-1)*|S|).
    It is evaluated by support weight. No logarithms or FLOAT powers enter.
    """
    orders = _orders(orders)
    counts = tuple(counts)
    if (type(q) is not int or q not in (4, 256) or not counts
            or len(counts) & (len(counts)-1) or not 2 <= len(counts) <= 65536
            or any(type(n) is not int or n < 0 for n in counts)):
        raise ValueError('GF4/GF256 labelled exact counts, width1..16, are required')
    width = len(counts).bit_length()-1
    if (type(codimension) is not int or not 0 <= codimension <= width
            or counts[0] != 1 or sum(counts) != q**(width-codimension)):
        raise ValueError('Wrong code mass or codimension')
    weight_counts = [[] for _ in range(width+1)]
    for mask, count in enumerate(counts):
        if count:
            weight_counts[mask.bit_count()].append(count)
    moments = {}
    for p in orders:
        value = sum((Fraction(sum(n**p for n in bucket), (q-1)**((p-1)*weight))
                     for weight, bucket in enumerate(weight_counts)), Fraction())
        exponent = p*codimension-width
        value *= Fraction(q**exponent) if exponent >= 0 else Fraction(1, q**(-exponent))
        if value < 1 or (p == 1 and value != 1):
            raise ArithmeticError('Exact normalized density check failed')
        moments[p] = value
    return moments


def representative_moments(complement, *, orders=(3, 4), q=256):
    """One exact support table, reused for every requested moment order."""
    complement = tuple(complement)
    orders = _orders(orders)
    data = support.support_counts(complement, q=q)
    return dict(complement=complement, arity=16-len(complement),
                moments=exact_moments(data['counts'], len(complement), orders=orders, q=q),
                raw_code_mass=data['mass'])


def checked_dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
    return dict(PINS, **{Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()})


def inventory(*, complement_sizes=(0, 1, 2, 3, 4), orders=(3, 4),
              max_seconds=360, progress=None):
    """Bounded serial representative inventory; partial output is explicit.

    A support calculation is indivisible. The time budget is checked
    before each representative; it is not an interrupting timeout.
    Every output row records the coordinate permutation witnessing reuse.
    """
    sizes = tuple(complement_sizes)
    orders = _orders(orders)
    if (not sizes or len(set(sizes)) != len(sizes)
            or any(type(b) is not int or not 0 <= b <= 4 for b in sizes)):
        raise ValueError('Use distinct complement sizes from0..4')
    if (type(max_seconds) not in (int, float) or not isfinite(max_seconds)
            or not 0 < max_seconds <= 600):
        raise ValueError('Use a positive budget at most600 seconds')
    before = checked_dependencies()
    plan = tuple((b, orbit) for b in sorted(sizes) for orbit in complement_orbits(b))
    started, rows, representatives = perf_counter(), [], []
    maxima = {b: {p: Fraction(1) for p in orders} for b in sizes}
    stop = None
    for b, orbit in plan:
        if perf_counter()-started >= max_seconds:
            stop = 'Budget exhausted: partial orbit inventory is not complete'
            break
        data = representative_moments(orbit['representative'], orders=orders)
        values = {str(p): str(data['moments'][p]) for p in orders}
        index = len(representatives)
        representatives.append(dict(orbit_index=index, complement=list(orbit['representative']),
            arity=16-b, multiplicity=len(orbit['members']), moments=values))
        for complement, image in orbit['members']:
            rows.append(dict(complement=list(complement), arity=16-b,
                representative=list(orbit['representative']), orbit_index=index,
                representative_to_complement_permutation=list(image), moments=values))
        for p in orders:
            maxima[b][p] = max(maxima[b][p], data['moments'][p])
        if progress:
            progress(dict(completed_representatives=index+1, total_representatives=len(plan),
                          completed_labelled_complements=len(rows), elapsed_seconds=perf_counter()-started))
    if checked_dependencies() != before:
        raise ValueError('Source changed during exact orbit inventory')
    complete = len(representatives) == len(plan)
    return dict(field_order=256, field_polynomial='0x11b', lambda_aes_byte=0xbc,
        complement_sizes=sorted(sizes), arities=sorted(16-b for b in sizes), moment_orders=list(orders),
        exact_rational_arithmetic=True, complete=complete, stop_reason=stop,
        full_parity_includes_identity=True, code='ker(Hlambda[complement,:])',
        output_diagonal='independent uniform nonzero field multipliers',
        normalized_likelihood='256^|complement|*n_C(S)/255^|S|',
        subgroup_size=48, symmetry='bit permutations and i->i+parity(i)*1111',
        no_larger_symmetry_claim=True, all_maps_check_full256_entries=True,
        maps=[list(image) for image in orthogonal_maps()],
        representatives=representatives, rows=sorted(rows, key=lambda row: (len(row['complement']), row['complement'])),
        maximum_moments_by_arity={str(16-b): {str(p): str(value) for p, value in values.items()}
                                 for b, values in maxima.items()},
        expected_labelled_complements=sum(comb(16, b) for b in sizes),
        expected_representatives=len(plan), actual_support_tables=len(representatives),
        reviewed_source_sha256=before, source_pins_unchanged=True,
        regional_provider_runs=0, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--complement-sizes', default='0,1,2,3,4')
    parser.add_argument('--orders', default='3,4')
    parser.add_argument('--max-seconds', type=float, default=360)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = inventory(complement_sizes=tuple(map(int, args.complement_sizes.split(','))),
        orders=tuple(map(int, args.orders.split(','))), max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'representatives', 'maps')}))
