"""Exact Hlambda support moments for mapped complements of sizes5..7.

The full parity matrix is H[O,:]=B+lambda*J, where
B[o,i]=1[o=i]+dot_F2(o,i). Its identity entries are retained.
Every minor of a rank-one perturbation is affine in lambda over F2.
For lambda outside F2 its rank therefore equals max(rank B,rank(B+J)).

Two packed binary-basis subset DPs compute those ranks. A basis uses at
most7 slots of7 bits; it is not an index into a byte-sized state table.
Once either binary rank is full, every larger support also has full rank.
Exact integer Mobius inversion then gives the labelled support counts.

The frozen48-map orthogonal subgroup transports local density moments.
Every labelled complement retains its actual coordinate-map witness.
This does not quotient scheduled words or remove rowset multiplicities.
No regional provider, RS bound, or whole-code certificate is computed.
"""
from array import array
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

import hd_lambda_renyi_orbits as small

DIRECTORY = Path(__file__).resolve().parent
MAX_WIDTH, MIN_COMPLEMENT, MAX_COMPLEMENT = 16, 5, 7
PINS = dict(small.PINS, **{
    'hd_lambda_renyi_orbits.py': '08e828b3d227fca874bce4622266a994550686f9df7997933890baaaae44799a',
})


def subset_ranks(binary_rows, width):
    """Exact ranks of (B+lambda*J)[:,S], with0..7 full binary rows."""
    rows = tuple(binary_rows)
    if (type(width) is not int or not 1 <= width <= MAX_WIDTH or len(rows) > MAX_COMPLEMENT
            or any(type(row) is not int or not 0 <= row < 1 << width for row in rows)):
        raise ValueError('At most7 binary rows and width1..16 are required')
    b, size = len(rows), 1 << width
    ranks = bytearray(size)
    if not b:
        return bytes(ranks)
    columns = tuple(sum(((row >> j) & 1) << i for i, row in enumerate(rows)) for j in range(width))
    ones = (1 << b)-1
    first, second = array('Q', [0])*size, array('Q', [0])*size
    first_rank, second_rank = bytearray(size), bytearray(size)
    for mask in range(1, size):
        bit = mask & -mask
        previous = mask ^ bit
        if ranks[previous] == b:
            ranks[mask] = b
            continue
        column = columns[bit.bit_length()-1]
        basis, rank, vector = first[previous], first_rank[previous], column
        while vector:
            shift = (vector.bit_length()-1)*b
            pivot = (basis >> shift) & ones
            if pivot:
                vector ^= pivot
            else:
                basis |= vector << shift
                rank += 1
                break
        if rank == b:
            ranks[mask] = b
            continue
        first[mask], first_rank[mask] = basis, rank
        basis, other_rank, vector = second[previous], second_rank[previous], column ^ ones
        while vector:
            shift = (vector.bit_length()-1)*b
            pivot = (basis >> shift) & ones
            if pivot:
                vector ^= pivot
            else:
                basis |= vector << shift
                other_rank += 1
                break
        second[mask], second_rank[mask] = basis, other_rank
        ranks[mask] = max(rank, other_rank)
    return bytes(ranks)


def support_counts_from_binary(binary_rows, width, *, q=256):
    """Labelled exact counts for a full B+lambda*J parity matrix."""
    if type(q) is not int or q not in (4, 256):
        raise ValueError('Only GF4 toy and AES-GF256 are supported')
    ranks = subset_ranks(binary_rows, width)
    powers = tuple(q**i for i in range(width+1))
    counts = [powers[mask.bit_count()-rank] for mask, rank in enumerate(ranks)]
    for j in range(width):
        step = 1 << j
        for start in range(0, 1 << width, step*2):
            for offset in range(step):
                counts[start+step+offset] -= counts[start+offset]
    codimension = ranks[-1]
    if counts[0] != 1 or min(counts) < 0 or sum(counts) != powers[width-codimension]:
        raise ArithmeticError('Exact Mobius inversion failed its mass/positivity check')
    return dict(width=width, q=q, codimension=codimension, subset_ranks=ranks,
                counts=tuple(counts), mass=powers[width-codimension])


def parity_binary_rows(complement):
    complement = tuple(complement)
    if (not MIN_COMPLEMENT <= len(complement) <= MAX_COMPLEMENT
            or len(set(complement)) != len(complement)
            or any(type(i) is not int or not 0 <= i < 16 for i in complement)):
        raise ValueError('Use5..7 distinct mapped complement labels')
    return tuple(sum((int(o == i) ^ ((o & i).bit_count() & 1)) << i for i in range(16))
                 for o in complement)


def support_counts(complement, *, q=256):
    """Actual code ker H[O,:], for mapped complement O BEFORE physical Q."""
    complement = tuple(complement)
    data = support_counts_from_binary(parity_binary_rows(complement), 16, q=q)
    if data['codimension'] != len(complement):
        raise ArithmeticError('Rows of involutive full Hlambda lost independence')
    return dict(data, complement=complement, arity=16-len(complement),
                lambda_aes_byte=0xbc, includes_diagonal_identity=True)


def representative_moments(complement, *, orders=(2, 3, 4), q=256):
    """Compute support counts once, then exact normalized density moments."""
    complement = tuple(complement)
    orders = small._orders(orders)
    data = support_counts(complement, q=q)
    return dict(complement=complement, arity=data['arity'],
                moments=small.exact_moments(data['counts'], len(complement), orders=orders, q=q),
                raw_code_mass=data['mass'])


def complement_orbits(complement_size):
    """Exact labelled partition under the frozen48-map subgroup only."""
    if type(complement_size) is not int or not MIN_COMPLEMENT <= complement_size <= MAX_COMPLEMENT:
        raise ValueError('Complement size must be5..7')
    remaining = set(combinations(range(16), complement_size))
    result = []
    for representative in combinations(range(16), complement_size):
        if representative not in remaining:
            continue
        members = {}
        for image in small.orthogonal_maps():
            target = tuple(sorted(image[i] for i in representative))
            members.setdefault(target, image)
        if not set(members) <= remaining:
            raise ArithmeticError('Computed subgroup orbits overlap')
        remaining.difference_update(members)
        result.append(dict(representative=representative,
                           members=tuple(sorted(members.items()))))
    if remaining or sum(len(row['members']) for row in result) != comb(16, complement_size):
        raise ArithmeticError('Labelled complement partition is incomplete')
    return tuple(result)


def checked_dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
    return dict(PINS, **{Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()})


def inventory(*, complement_sizes=(5,), orders=(2, 3, 4), max_complement_size=5,
              max_seconds=180, progress=None):
    """Bounded serial inventory; the default admits only complement size5.

    Sizes6/7 require explicitly increasing max_complement_size. The budget
    is checked before each indivisible representative support calculation.
    A partial return retains its rows but is not a completed inventory.
    """
    sizes, orders = tuple(complement_sizes), small._orders(orders)
    if (type(max_complement_size) is not int or not 5 <= max_complement_size <= 7
            or not sizes or len(set(sizes)) != len(sizes)
            or any(type(b) is not int or not 5 <= b <= max_complement_size for b in sizes)):
        raise ValueError('Use distinct complement sizes5..cap, with explicit cap5..7')
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
    return dict(field_order=256, field_polynomial='0x11b', lambda_aes_byte=0xbc,
        complement_sizes=sorted(sizes), arities=sorted(16-b for b in sizes), moment_orders=list(orders),
        max_complement_size=max_complement_size, exact_rational_arithmetic=True,
        complete=len(representatives) == len(plan), stop_reason=stop,
        full_parity_includes_identity=True, code='ker(Hlambda[complement,:])',
        output_diagonal='independent uniform nonzero field multipliers',
        normalized_likelihood='256^|complement|*n_C(S)/255^|S|',
        subgroup_size=48, symmetry='bit permutations and i->i+parity(i)*1111',
        no_larger_symmetry_claim=True, all_maps_check_full256_entries=True,
        maps=[list(image) for image in small.orthogonal_maps()], representatives=representatives,
        rows=sorted(rows, key=lambda row: (len(row['complement']), row['complement'])),
        maximum_moments_by_arity={str(16-b): {str(p): str(value) for p, value in values.items()}
                                 for b, values in maxima.items()},
        expected_labelled_complements=sum(comb(16, b) for b in sizes), expected_representatives=len(plan),
        actual_support_tables=len(representatives), rank_engine='two packed binary-basis subset DPs',
        reviewed_source_sha256=before, source_pins_unchanged=True,
        regional_provider_runs=0, whole_code_certificate=False, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--complement-sizes', default='5')
    parser.add_argument('--max-complement-size', type=int, default=5, choices=(5, 6, 7))
    parser.add_argument('--orders', default='2,3,4')
    parser.add_argument('--max-seconds', type=float, default=180)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = inventory(complement_sizes=tuple(map(int, args.complement_sizes.split(','))),
        orders=tuple(map(int, args.orders.split(','))), max_complement_size=args.max_complement_size,
        max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'representatives', 'maps')}))
