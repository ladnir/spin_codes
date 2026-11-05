"""Isolated Hlambda variant with one binary quiet transvection.

Let v be the indicator of a four-label set V whose label XOR is zero.
For H=I+lambda*J+sum_k b_k*b_k^T, these conditions give Hv=v and v^T v=0.
Thus T=I+v*v^T commutes with H, T^2=I, and

    H' = T H = H+v*v^T

is a symmetric involution. Every entry remains lambda or lambda+1, hence
nonzero when lambda is outside F2. This is a NEW construction, not a
performance-preserving implementation of the measured Hlambda matrix.

The naive circuit applies H then a four-output parity update: seven extra
XORs. For default V=(0,1,2,3), the existing H circuit already computes
q0=x0+x1+x2+x3. Since v^T Hx=v^T x=q0, four output XORs suffice. The fused
scalar circuit has61 XORs and one fixed product, excluding all diagonals
and routing. No native cost or timing is inferred from these counts.

For a mapped complement O, the parity matrix is H'[O,:]=B'+lambda*J,
where B'=I+binary-dot+v*v^T includes its identity entries. The reviewed
packed subset-rank and integer Mobius engine therefore applies unchanged.
Only maps in the old48-map group that stabilize V are used for local
moment reuse. For V0123 this subgroup has four maps, not48. Every labelled
complement remains in an inventory; no scheduled word is quotiented.

If O is disjoint from V, its parity matrix is exactly unchanged. There
are C(12,b) such complements of size b, including792 at b=5. Local gains
therefore do not imply a better worst-case moment. Old numerical receipts,
distance claims, and the measured3.442ms timing do not transfer.
In particular, the old support-distance table needs a new proof; the
affine-image premise behind that table changes after adding v*v^T.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

import hd_lambda_outer as original
import hd_lambda_renyi_middle as middle
import hd_lambda_renyi_orbits as symmetry

DIRECTORY = Path(__file__).resolve().parent
DEFAULT_SUPPORT = (0, 1, 2, 3)
DEFAULT_LAMBDA = 0xbc
CONSTRUCTION = 'g16-rs16x8-hlambda-bc-plus-quiet-transvection-v000f-research'
PINS = {
    'hd_lambda_outer.py': 'f0c063843e295a160a1186c6c084a07dd1c2a3745005ed68057fde59913c99e7',
    'hd_outer_pair_layers.py': '017c5ba20747ae8e25acf02393e9702a03b086038e1814ac2424cdfe49d66f60',
    'hd_lambda_renyi_middle.py': '1e43e08877c2f23195605c293d810f40063caee10073abf8383a7388f63fb115',
    'hd_lambda_renyi_orbits.py': '08e828b3d227fca874bce4622266a994550686f9df7997933890baaaae44799a',
}


def _support(support):
    support = tuple(support)
    if (len(support) != 4 or any(type(i) is not int or not 0 <= i < 16 for i in support)
            or len(set(support)) != 4):
        raise ValueError('A binary weight-four support on labels0..15 is required')
    value = 0
    for i in support:
        value ^= i
    if value:
        raise ValueError('The four labels must XOR to zero so that Hv=v')
    return tuple(sorted(support))


def _complement(complement):
    complement = tuple(complement)
    if (len(complement) > 7 or any(type(i) is not int or not 0 <= i < 16 for i in complement)
            or len(set(complement)) != len(complement)):
        raise ValueError('Use at most seven distinct mapped complement labels')
    return tuple(sorted(complement))


def construction(support=DEFAULT_SUPPORT):
    support = _support(support)
    return f'g16-rs16x8-hlambda-bc-plus-quiet-transvection-v{sum(1 << i for i in support):04x}-research'


@lru_cache(maxsize=32)
def _matrix(lambda_value, q, support):
    old = original.matrix(lambda_value, q)
    active = set(support)
    return tuple(tuple(old[i][j] ^ int(i in active and j in active) for j in range(16))
                 for i in range(16))


def matrix(lambda_value=DEFAULT_LAMBDA, q=256, *, support=DEFAULT_SUPPORT):
    """Actual new matrix; GF4 toys use lambda_value=2 or3 explicitly."""
    return _matrix(lambda_value, q, _support(support))


def _values(values, lambda_value, q):
    original._parameters(lambda_value, q)
    values = tuple(values)
    if len(values) != 16 or any(type(x) is not int or not 0 <= x < q for x in values):
        raise ValueError('Sixteen values in the declared field are required')
    return values


def direct(values, lambda_value=DEFAULT_LAMBDA, q=256, *, support=DEFAULT_SUPPORT):
    values = _values(values, lambda_value, q)
    output = []
    for row in matrix(lambda_value, q, support=support):
        value = 0
        for coefficient, x in zip(row, values):
            value ^= original._multiply(coefficient, x, q)
        output.append(value)
    return tuple(output)


def _naive_circuit(x, add, multiply, lambda_value, support):
    output = list(original._circuit(x, add, multiply, lambda_value))
    parity = output[support[0]]
    for i in support[1:]:
        parity = add(parity, output[i])
    for i in support:
        output[i] = add(output[i], parity)
    return tuple(output)


def _fused_default_circuit(x, add, multiply, lambda_value):
    # Preserve the old fixed16 circuit and reuse its existing q[0].
    p = tuple(add(x[2*i], x[2*i+1]) for i in range(8))
    a0 = add(add(add(x[1], x[3]), add(x[5], x[7])), add(add(x[9], x[11]), add(x[13], x[15])))
    a1 = add(add(p[1], p[3]), add(p[5], p[7]))
    q = tuple(add(p[2*i], p[2*i+1]) for i in range(4))
    a2, a3 = add(q[1], q[3]), add(q[2], q[3])
    total = multiply(lambda_value, add(add(q[0], q[1]), a3))
    syndromes = (a0, a1, a2, a3)
    correction = [total]*16
    for i in range(1, 16):
        low = i & -i
        correction[i] = add(correction[i ^ low], syndromes[low.bit_length()-1])
    output = [add(x[i], correction[i]) for i in range(16)]
    for i in range(4):
        output[i] = add(output[i], q[0])
    return tuple(output)


def circuit(values, lambda_value=DEFAULT_LAMBDA, q=256, *, support=DEFAULT_SUPPORT,
            fuse_default=True):
    """Use the fused default circuit, or the literal post-H transvection."""
    support = _support(support)
    values = _values(values, lambda_value, q)
    if type(fuse_default) is not bool:
        raise ValueError('fuse_default must be Boolean')
    add = lambda a, b: a ^ b
    multiply = lambda a, b: original._multiply(a, b, q)
    if fuse_default and support == DEFAULT_SUPPORT:
        return _fused_default_circuit(values, add, multiply, lambda_value)
    return _naive_circuit(values, add, multiply, lambda_value, support)


def primitive_counts(*, support=DEFAULT_SUPPORT, fuse_default=True):
    """Scalar operation counts only; no SIMD, key-load, or timing claim."""
    support = _support(support)
    if type(fuse_default) is not bool:
        raise ValueError('fuse_default must be Boolean')
    counts = Counter()
    def add(a, b):
        counts['binary_xors'] += 1
        return a ^ b
    def multiply(a, b):
        counts['fixed_products'] += 1
        return original._multiply(a, b, 256)
    if fuse_default and support == DEFAULT_SUPPORT:
        _fused_default_circuit((0,)*16, add, multiply, DEFAULT_LAMBDA)
    else:
        _naive_circuit((0,)*16, add, multiply, DEFAULT_LAMBDA, support)
    return dict(construction=construction(support), support=support,
        binary_xors=counts['binary_xors'], fixed_products=counts['fixed_products'], variable_products=0,
        extra_xors_over_Hlambda=counts['binary_xors']-original.primitive_counts()['binary_xors'],
        fused_existing_q0=fuse_default and support == DEFAULT_SUPPORT,
        includes_diagonals=False, includes_route=False, native_implementation_exists=False,
        native_performance_claim=False, measured_current_construction_unchanged=True)


def parity_binary_rows(complement, *, support=DEFAULT_SUPPORT):
    """Full B' parity rows before physical Q, including diagonal identity."""
    complement, support = _complement(complement), _support(support)
    active = set(support)
    return tuple(sum((int(o == i) ^ ((o & i).bit_count() & 1)
                       ^ int(o in active and i in active)) << i for i in range(16))
                 for o in complement)


def support_counts(complement, *, q=256, support=DEFAULT_SUPPORT):
    """Exact labelled counts of ker H'[O,:], for mapped O of size0..7."""
    complement, support = _complement(complement), _support(support)
    data = middle.support_counts_from_binary(parity_binary_rows(complement, support=support), 16, q=q)
    if data['codimension'] != len(complement):
        raise ArithmeticError('Rows of the new involution lost independence')
    return dict(data, complement=complement, arity=16-len(complement),
                construction=construction(support), transvection_support=support,
                lambda_aes_byte=DEFAULT_LAMBDA, includes_diagonal_identity=True)


def support_moments(complement, *, orders=(2, 3, 4), q=256, support=DEFAULT_SUPPORT):
    """One exact support table, reused for each normalized density moment."""
    complement, support = _complement(complement), _support(support)
    orders = symmetry._orders(orders)
    data = support_counts(complement, q=q, support=support)
    return dict(complement=complement, arity=16-len(complement),
        construction=construction(support), transvection_support=support, field_order=q,
        moments=symmetry.exact_moments(data['counts'], len(complement), orders=orders, q=q),
        zero_moment_contributions={p: F(q)**(p*len(complement)-16) for p in orders},
        raw_code_mass=data['mass'], parity_unchanged_from_Hlambda=not set(complement).intersection(support),
        old_receipts_transferred=False, regional_provider_runs=0, whole_code_certificate=False)


@lru_cache(maxsize=16)
def _stabilizer(support):
    active = set(support)
    maps = tuple(image for image in symmetry.orthogonal_maps()
                 if {image[i] for i in support} == active)
    for image in maps:
        for i in range(16):
            for j in range(16):
                old = int(i == j) ^ ((i & j).bit_count() & 1) ^ int(i in active and j in active)
                moved = (int(image[i] == image[j]) ^ ((image[i] & image[j]).bit_count() & 1)
                         ^ int(image[i] in active and image[j] in active))
                if old != moved:
                    raise ArithmeticError('A claimed stabilizer map changes the new matrix')
    if support == DEFAULT_SUPPORT and len(maps) != 4:
        raise ArithmeticError('The default stabilizer must have exactly four maps')
    return maps


def stabilizer_maps(support=DEFAULT_SUPPORT):
    """Only the stabilizer of V inside the reviewed48-map group is claimed."""
    return _stabilizer(_support(support))


def complement_orbits(size, *, support=DEFAULT_SUPPORT):
    if type(size) is not int or not 0 <= size <= 7:
        raise ValueError('Complement size0..7 required')
    support = _support(support)
    maps = stabilizer_maps(support)
    remaining = set(combinations(range(16), size))
    result = []
    for representative in combinations(range(16), size):
        if representative not in remaining:
            continue
        members = {}
        for image in maps:
            target = tuple(sorted(image[i] for i in representative))
            members.setdefault(target, image)
        if not set(members) <= remaining:
            raise ArithmeticError('Stabilizer orbits overlap')
        remaining.difference_update(members)
        result.append(dict(representative=representative, members=tuple(sorted(members.items()))))
    if remaining or sum(len(row['members']) for row in result) != comb(16, size):
        raise ArithmeticError('Labelled complement partition is incomplete')
    return tuple(result)


def checked_dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
    return {**middle.checked_dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def inventory(*, complement_sizes=(0, 1, 2, 3, 4, 5), orders=(2, 3, 4),
              support=DEFAULT_SUPPORT, max_seconds=180, progress=None):
    """Optional bounded local inventory; no regional bank is ever invoked.

    Each representative is indivisible. Budget exhaustion returns partial
    evidence without a completeness claim. No previously measured moment
    or timing is substituted for this variant's actual calculation.
    """
    support, sizes, orders = _support(support), tuple(complement_sizes), symmetry._orders(orders)
    if (not sizes or len(set(sizes)) != len(sizes)
            or any(type(b) is not int or not 0 <= b <= 5 for b in sizes)):
        raise ValueError('The bounded inventory allows distinct complement sizes0..5 only')
    if (type(max_seconds) not in (int, float) or not isfinite(max_seconds)
            or not 0 < max_seconds <= 600):
        raise ValueError('Use a positive bounded budget at most600 seconds')
    pins = checked_dependencies()
    started = perf_counter()
    plan = tuple((b, orbit) for b in sorted(sizes) for orbit in complement_orbits(b, support=support))
    rows, representatives, stop = [], [], None
    maxima = {b: {p: F(1) for p in orders} for b in sizes}
    for b, orbit in plan:
        if perf_counter()-started >= max_seconds:
            stop = 'Time budget: partial local inventory is not complete'
            break
        data = support_moments(orbit['representative'], orders=orders, support=support)
        index = len(representatives)
        values = {str(p): str(data['moments'][p]) for p in orders}
        representatives.append(dict(orbit_index=index, complement=list(orbit['representative']),
            arity=16-b, multiplicity=len(orbit['members']), moments=values))
        for complement, image in orbit['members']:
            rows.append(dict(complement=list(complement), arity=16-b, orbit_index=index,
                representative=list(orbit['representative']),
                representative_to_complement_permutation=list(image), moments=values))
        for p in orders:
            maxima[b][p] = max(maxima[b][p], data['moments'][p])
        if progress:
            progress(dict(completed_representatives=index+1, total_representatives=len(plan),
                completed_labelled_complements=len(rows), elapsed_seconds=perf_counter()-started))
    if checked_dependencies() != pins:
        raise ValueError('Frozen dependencies changed during variant inventory')
    return dict(construction=construction(support), transvection_support=list(support),
        field_order=256, field_polynomial='0x11b', lambda_aes_byte=DEFAULT_LAMBDA,
        complement_sizes=sorted(sizes), arities=sorted(16-b for b in sizes), moment_orders=list(orders),
        exact_rational_arithmetic=True, complete=len(representatives) == len(plan), stop_reason=stop,
        full_parity_includes_identity=True, code='ker((Hlambda+v*v^T)[complement,:])',
        output_diagonal='independent uniform nonzero field multipliers',
        normalized_likelihood='256^|complement|*n_C(S)/255^|S|',
        subgroup_size=len(stabilizer_maps(support)), old_subgroup_size=48,
        symmetry='stabilizer of transvection support inside proved orthogonal subgroup',
        no_larger_symmetry_claim=True, all_maps_check_full256_entries=True,
        maps=[list(image) for image in stabilizer_maps(support)],
        representatives=representatives, rows=sorted(rows, key=lambda row: (len(row['complement']), row['complement'])),
        maximum_moments_by_arity={str(16-b): {str(p): str(v) for p, v in values.items()} for b, values in maxima.items()},
        unchanged_disjoint_complements_by_size={str(b): comb(12, b) for b in sizes},
        expected_labelled_complements=sum(comb(16, b) for b in sizes), expected_representatives=len(plan),
        actual_support_tables=len(representatives), original_rowsets_not_resampled=True,
        reviewed_source_sha256=pins, source_pins_unchanged=True,
        circuit_counts=primitive_counts(support=support), old_receipts_transferred=False,
        regional_provider_runs=0, native_benchmarks=0, native_performance_claim=False,
        whole_code_certificate=False, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--complement-sizes', default='0,1,2,3,4,5')
    parser.add_argument('--orders', default='2,3,4')
    parser.add_argument('--support', default='0,1,2,3')
    parser.add_argument('--max-seconds', type=float, default=180)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = inventory(complement_sizes=tuple(map(int, args.complement_sizes.split(','))),
        orders=tuple(map(int, args.orders.split(','))), support=tuple(map(int, args.support.split(','))),
        max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'representatives', 'maps')}))
