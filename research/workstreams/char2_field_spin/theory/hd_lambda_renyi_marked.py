"""Exact input-weight polynomials for the column likelihood of Hlambda.

Fix a mapped source set A and allow every input in GF(q)^A, including zero.
Let n(w,S) count inputs of weight w with output support S before the final
independent uniform nonzero diagonal. Put K_S(X)=sum_w n(w,S) X^w.
The raw counting density relative to uniform GF(q)^16 has p-th moment

    J_p(X) = q^(16*(p-1)) sum_S K_S(X)^p/(q-1)^((p-1)*|S|).

The power is taken after summing input weights for the SAME labelled S.
Coefficients are exact Fractions, ordered by degree zero through a*p.
In particular J_p(1)=q^(a*p) I_p for the normalized column density.

For a fixed original A, a caller may multiply the sixteen scheduled J_p
polynomials and retain degrees at least 9*a*p. This relaxes the row floors
in all p moment copies. The resulting p-th root still needs 255^(-8*a).
No extra q factor belongs outside this raw-moment formula.

Only local mapped sets are identified by the proved 48-map subgroup.
The inventory retains every labelled A and a representative-to-A witness.
It does not compose columns, impose RS floors, average source sets, or
provide a whole-code certificate. Existing proof helpers stay unchanged.
"""
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

import hd_lambda_arity5_relations as relations
import hd_lambda_renyi_orbits as orbits

DIRECTORY = Path(__file__).resolve().parent
MAX_ARITY = 5
MOMENT_ORDERS = (2, 3, 4)
LAMBDA_AES_BYTE = 0xbc
PINS = {
    'hd_lambda_arity5_relations.py': 'f01f060d861c976a10a05dd22e3a9e79e674e1031ae752462ee980d5cb0cc781',
    'hd_quaternary_outer_relations.py': 'd8be4e7a1ab77b98c76695f13c1bc9a9cf867df2aae6408ac3e08c6fc12346c0',
    'hd_lambda_outer.py': 'f0c063843e295a160a1186c6c084a07dd1c2a3745005ed68057fde59913c99e7',
    'hd_pair_scheduled_outer.py': '64dc50bf18171d371ae16f4484087e9d42d89cc7c9b55a7576be55b1251e4741',
    'hd_outer_pair_layers.py': '017c5ba20747ae8e25acf02393e9702a03b086038e1814ac2424cdfe49d66f60',
    'hd_binary_syndrome_outer.py': 'b1dc2dcb6e212a28526f7c7dc7f7a55999be7e6f5f53370f78768d705f2bda96',
    'hd_lambda_renyi_orbits.py': '08e828b3d227fca874bce4622266a994550686f9df7997933890baaaae44799a',
    'hd_lambda_renyi_support.py': 'c4662bb016f248faa7ffcc1972e01349bc66306c40762e12df9245fc04a3c03e',
    'hd_binary_outer_relations.py': 'ec1c13cd4a685760983b33ebd8e722efc3735f6335ed55bfe444aea01b1a20f2',
}


def _orders(orders):
    orders = tuple(orders)
    if (not orders or len(set(orders)) != len(orders)
            or any(type(p) is not int or p not in MOMENT_ORDERS for p in orders)):
        raise ValueError('Use distinct integer moment orders from 2,3,4')
    return tuple(sorted(orders))


def _arity(arity):
    if type(arity) is not int or not 1 <= arity <= MAX_ARITY:
        raise ValueError('This interface supports arities 1..5 only')


def _field(q):
    if type(q) is not int or q not in (4, 256):
        raise ValueError('Use GF4 toys or the GF256 construction')


def _source(source):
    source = tuple(source)
    _arity(len(source))
    if (len(set(source)) != len(source)
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('Distinct mapped source labels in 0..15 are required')
    # Only TOTAL input weight is marked, so the source tuple order is immaterial.
    return tuple(sorted(source))


def _integer_convolution(left, right):
    result = [0]*(len(left)+len(right)-1)
    for i, x in enumerate(left):
        if x:
            for j, y in enumerate(right):
                if y:
                    result[i+j] += x*y
    return result


def marked_polynomials_from_joint(joint, arity, *, width=16, q=256,
                                  orders=MOMENT_ORDERS):
    """Return p -> immutable exact coefficients from raw (weight,mask) counts.

    The caller supplies an injective map's exact allowed-input counts.
    Validation checks mass, ranges, and its unique zero-input/output branch;
    these checks alone do not establish that arbitrary supplied counts came
    from a linear map. Tiny tests may use widths below sixteen.
    """
    _arity(arity)
    _field(q)
    orders = _orders(orders)
    if type(width) is not int or not arity <= width <= 16:
        raise ValueError('Use arity <= width <= 16')
    if not hasattr(joint, 'items'):
        raise ValueError('Supply a mapping of (input weight,output mask) to count')
    by_mask = {}
    mass = 0
    for key, count in joint.items():
        if (not isinstance(key, tuple) or len(key) != 2
                or type(key[0]) is not int or not 0 <= key[0] <= arity
                or type(key[1]) is not int or not 0 <= key[1] < 1 << width
                or type(count) is not int or count <= 0):
            raise ValueError('Counts must be positive integers with valid weight/mask keys')
        weight, mask = key
        if (weight == 0 or mask == 0) and (key != (0, 0) or count != 1):
            raise ValueError('An injective map has exactly one zero branch')
        by_mask.setdefault(mask, [0]*(arity+1))[weight] = count
        mass += count
    if joint.get((0, 0)) != 1 or mass != q**arity:
        raise ValueError('Wrong raw allowed-input mass or missing zero branch')

    # Integer buckets avoid a Fraction allocation for each output support.
    # Masks with equal support WEIGHT are combined only AFTER K_S is powered.
    buckets = {p: [[0]*(arity*p+1) for _ in range(width+1)] for p in orders}
    for mask, polynomial in by_mask.items():
        terms = tuple((i, n) for i, n in enumerate(polynomial) if n)
        weight = mask.bit_count()
        if len(terms) == 1:
            degree, count = terms[0]
            for p in orders:
                buckets[p][weight][p*degree] += count**p
            continue
        power = polynomial
        for p in range(2, max(orders)+1):
            power = _integer_convolution(power, polynomial)
            if p in buckets:
                target = buckets[p][weight]
                for degree, count in enumerate(power):
                    target[degree] += count
    result = {}
    for p in orders:
        scale = q**(width*(p-1))
        coefficients = [Fraction(0) for _ in range(arity*p+1)]
        for weight, integer_coefficients in enumerate(buckets[p]):
            denominator = (q-1)**((p-1)*weight)
            for degree, count in enumerate(integer_coefficients):
                if count:
                    coefficients[degree] += Fraction(scale*count, denominator)
        if coefficients[0] != scale:
            raise ArithmeticError('The raw zero-column moment must equal q^(width*(p-1))')
        result[p] = tuple(coefficients)
    return result


@lru_cache(maxsize=7000)
def _representative(source):
    maps = orbits.orthogonal_maps()
    representative = min(tuple(sorted(image[i] for i in source)) for image in maps)
    image = next(image for image in maps
                 if tuple(sorted(image[i] for i in representative)) == source)
    return representative, image


@lru_cache(maxsize=512)
def _representative_polynomials(representative, q, orders):
    joint = relations.joint_counts(representative, 'static', 0,
                                   lambda_value=LAMBDA_AES_BYTE, q=q)
    polynomials = marked_polynomials_from_joint(joint, len(representative), q=q, orders=orders)
    return (tuple(polynomials.items()), len(joint), len({mask for _, mask in joint}))


def mapped_source_polynomials(mapped_source, *, q=256, orders=MOMENT_ORDERS):
    """Exact local J_p, reusing only a proved simultaneous label permutation.

    Q_j merely permutes output masks and therefore leaves these moments
    unchanged. A scheduled caller still supplies its actual mapped P_j(A)
    separately in all sixteen columns, for the SAME original A.
    """
    source = _source(mapped_source)
    _field(q)
    orders = _orders(orders)
    representative, image = _representative(source)
    pairs, joint_entries, output_masks = _representative_polynomials(representative, q, orders)
    return dict(mapped_source=source, arity=len(source), width=16, field_order=q,
        lambda_aes_byte=LAMBDA_AES_BYTE, representative=representative,
        representative_to_source_permutation=image, polynomials=dict(pairs),
        raw_allowed_input_mass=q**len(source), joint_entries=joint_entries,
        output_masks=output_masks, rs_normalization_included=False,
        whole_code_certificate=False)


def source_orbits(arity):
    """Partition labelled mapped source sets; do not quotient word chronology."""
    _arity(arity)
    remaining = set(combinations(range(16), arity))
    result = []
    while remaining:
        representative = min(remaining)
        members = {}
        for image in orbits.orthogonal_maps():
            target = tuple(sorted(image[i] for i in representative))
            members.setdefault(target, image)
        if not set(members) <= remaining:
            raise ArithmeticError('Proved subgroup orbits overlap')
        remaining.difference_update(members)
        result.append(dict(representative=representative,
                           members=tuple(sorted(members.items()))))
    if sum(len(row['members']) for row in result) != comb(16, arity):
        raise ArithmeticError('Labelled source partition is incomplete')
    return tuple(result)


def checked_dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
    return dict(PINS, **{Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()})


def inventory(*, arities=(1, 2, 3, 4, 5), orders=MOMENT_ORDERS,
              max_seconds=360, progress=None):
    """Bounded exact representative inventory, with explicit partial output.

    Each representative calculation is indivisible. The budget is checked
    before it starts. Coefficients are stored once per representative;
    every labelled source retains its index and an explicit map witness.
    """
    arities = tuple(arities)
    orders = _orders(orders)
    if not arities or len(set(arities)) != len(arities):
        raise ValueError('Use distinct arities from 1..5')
    for arity in arities:
        _arity(arity)
    if (type(max_seconds) not in (int, float) or not isfinite(max_seconds)
            or not 0 < max_seconds <= 600):
        raise ValueError('Use a positive budget at most 600 seconds')
    before = checked_dependencies()
    started = perf_counter()
    plan = tuple((a, orbit) for a in sorted(arities) for orbit in source_orbits(a))
    rows, representatives, stop = [], [], None
    for arity, orbit in plan:
        if perf_counter()-started >= max_seconds:
            stop = 'Budget exhausted: partial local inventory is not complete'
            break
        data = mapped_source_polynomials(orbit['representative'], orders=orders)
        index = len(representatives)
        representatives.append(dict(orbit_index=index, arity=arity,
            representative=list(orbit['representative']), multiplicity=len(orbit['members']),
            coefficients={str(p): [str(c) for c in data['polynomials'][p]] for p in orders},
            raw_allowed_input_mass=data['raw_allowed_input_mass'],
            joint_entries=data['joint_entries'], output_masks=data['output_masks']))
        for source, image in orbit['members']:
            rows.append(dict(mapped_source=list(source), arity=arity, orbit_index=index,
                representative_to_source_permutation=list(image)))
        if progress:
            progress(dict(completed_representatives=index+1, total_representatives=len(plan),
                          completed_labelled_sources=len(rows), elapsed_seconds=perf_counter()-started))
    if checked_dependencies() != before:
        raise ValueError('Source changed during exact marked inventory')
    expected = sum(comb(16, a) for a in arities)
    keys = {(row['arity'], tuple(row['mapped_source'])) for row in rows}
    if len(keys) != len(rows):
        raise ArithmeticError('Duplicate labelled source in inventory')
    complete = len(representatives) == len(plan) and len(rows) == expected
    return dict(field_order=256, field_polynomial='0x11b', lambda_aes_byte=LAMBDA_AES_BYTE,
        width=16, arities=sorted(arities), moment_orders=list(orders),
        coefficients_order='ascending original input degree 0..arity*moment',
        formula='J_p(X)=256^(16*(p-1))*sum_S K_S(X)^p/255^((p-1)*popcount(S))',
        raw_likelihood=True, normalized_column_identity='J_p(1)=256^(arity*p)*I_p',
        exact_rational_arithmetic=True, complete=complete, stop_reason=stop,
        subgroup_size=48, symmetry='bit permutations and i->i+parity(i)*1111',
        no_larger_symmetry_claim=True, whole_word_chronology_quotiented=False,
        original_source_multiplicity_removed=False, rs_normalization_included=False,
        scheduled_columns_composed=False, input_floor_applied=False,
        maps=[list(image) for image in orbits.orthogonal_maps()],
        representatives=representatives, rows=sorted(rows, key=lambda row: (row['arity'], row['mapped_source'])),
        expected_labelled_sources=expected, expected_representatives=len(plan),
        completed_representatives=len(representatives), completed_labelled_sources=len(rows),
        reviewed_source_sha256=before, source_pins_unchanged=True,
        regional_provider_runs=0, whole_code_certificate=False,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arities', default='1,2,3,4,5')
    parser.add_argument('--orders', default='2,3,4')
    parser.add_argument('--max-seconds', type=float, default=360)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = inventory(arities=tuple(map(int, args.arities.split(','))),
        orders=tuple(map(int, args.orders.split(','))), max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'representatives', 'maps')}))
