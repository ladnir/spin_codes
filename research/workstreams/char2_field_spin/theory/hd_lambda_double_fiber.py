"""Exact labelled fibre factors for the double-transvection matrix.

Fix a mapped source set A and O=A^c. C=Hnew[O,A] has exact image and kernel
support enumerators C_A and Z_A. Symmetric involution implies that Hnew maps
ker C bijectively to itself. For output atoms in[0,1], Fourier positivity
of product weights gives K_A<=C_A Z_A. This comparison is unmarked.

The zero-outside fibre need not be pointwise quiet. Put N=Hnew-I and mark
i in A as protected exactly when adding row N[i,A] to C preserves rank.
On ker C, input and output agree at each protected coordinate. For a
nonzero zero-outside output mask T, inflating every unprotected coordinate
gives input-degree upper bound |T intersect I_protected|+|A minus I_protected|.
The zero output is exactly zero input and retains degree0. Nonzero outside
outputs are assigned degree|A| before the unmarked coset comparison.

All exact counts allow zero input coordinates inside A. Original A is kept
fixed across16 scheduled columns; both factors are transported through Q.
No RS normalization, regional provider, or old outer certificate is used.
"""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path

import hd_lambda_double_transvection as outer
from hd_lambda_fiber_outer import _support_enumerators, _labels, _permute_mask
from hd_pair_scheduled_outer import column_permutations


DIRECTORY = Path(__file__).resolve().parent
FULL = 65535
SCHEDULE = 'input-inverse-output-xor'
PINS = {
    'hd_lambda_double_transvection.py': '973175e25051621470bd19f8c6932806e4a830b908a8a1eb04186e4b616430bd',
    'hd_lambda_fiber_outer.py': 'ae4d54b67e8bf42990e481b36d1a085d1dd56aa34da472ee11aacbd40499edf6',
    'hd_pair_scheduled_outer.py': '64dc50bf18171d371ae16f4484087e9d42d89cc7c9b55a7576be55b1251e4741',
}
POINTS = tuple(1 | (i << 1) | (int(i in outer.SUPPORTS[0]) << 5)
               | (int(i in outer.SUPPORTS[1]) << 6) for i in range(16))


def _mask(mask):
    if type(mask) is not int or not 0 <= mask <= FULL:
        raise ValueError('A16-coordinate mask is required')
    return mask


def _source(source):
    source = tuple(source)
    if (not source or len(source) > 16 or len(set(source)) != len(source)
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('Distinct ordered source labels from0..15 required')
    return source


@lru_cache(maxsize=65536)
def syndrome_span(mask):
    """Canonical binary row basis for the seven-component augmented labels."""
    _mask(mask)
    basis = [0]*7
    for i in _labels(mask):
        value = POINTS[i]
        while value:
            pivot = value.bit_length()-1
            if basis[pivot]:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    for low in range(7):
        if basis[low]:
            for high in range(low+1, 7):
                if basis[high] & (1 << low):
                    basis[high] ^= basis[low]
    return tuple(value for value in reversed(basis) if value)


@lru_cache(maxsize=65536)
def _rank_spans(left, right):
    # On canonical spans, the lambda term is u*v^T, not necessarily J.
    # Every minor is affine in lambda, so max(rank B,rank(B+u*v^T)) is exact.
    binary = tuple(sum((((u >> 1) & (v >> 1)).bit_count() & 1) << j
                       for j, v in enumerate(right)) for u in left)
    rank_one = sum((v & 1) << j for j, v in enumerate(right))
    shifted = tuple(row ^ (rank_one if u & 1 else 0) for row, u in zip(binary, left))
    return max(outer.core.binary_rank(binary), outer.core.binary_rank(shifted))


def rank_update_masks(rows, columns):
    """rank((Hnew-I)[rows,columns]); it is Hnew's rank only for disjoint sets."""
    left, right = syndrome_span(_mask(rows)), syndrome_span(_mask(columns))
    if not left or not right:
        return 0
    return _rank_spans(min(left, right), max(left, right))


@lru_cache(maxsize=512)
def _fiber_cached(amask, q):
    source, outside = tuple(_labels(amask)), tuple(_labels(FULL ^ amask))
    omask = FULL ^ amask
    rank = rank_update_masks(omask, amask)
    outside_counts, zero_counts = _support_enumerators(source, outside, rank, rank_update_masks, q)
    protected = sum(1 << i for i in source if rank_update_masks(omask | (1 << i), amask) == rank)
    return rank, protected, outside_counts, zero_counts


def fiber_counts(mapped_source, *, q=256, max_subsets=2048):
    source = _source(mapped_source)
    if type(q) is not int or q not in (4, 256):
        raise ValueError('Only exact GF4 toy or GF256 counts are in scope')
    if (type(max_subsets) is not int or not 1 <= max_subsets <= 65536
            or max(1 << len(source), 1 << (16-len(source))) > max_subsets):
        raise ValueError('Both factor subset tables must fit the explicit budget')
    amask = sum(1 << i for i in source)
    rank, protected, outside, zero = _fiber_cached(amask, q)
    return dict(mapped_source=source, mapped_source_mask=amask, mapped_outside_mask=FULL ^ amask,
        arity=len(source), cross_rank=rank, syndrome_rank=len(syndrome_span(amask)),
        protected_mask=protected, all_zero_fiber_coordinates_protected=protected == amask,
        outside_counts=dict(outside), zero_fiber_counts=dict(zero),
        outside_mass=q**rank, zero_fiber_mass=q**(len(source)-rank), raw_column_mass=q**len(source),
        field_order=q, construction=outer.CONSTRUCTION,
        pointwise_atom_domain='[0,1]', old_outer_receipts_transferred=False)


def column_fiber_counts(original_source, column, *, q=256, max_subsets=2048):
    source = _source(original_source)
    before, after = column_permutations(SCHEDULE, column)
    data = fiber_counts(tuple(before[i] for i in source), q=q, max_subsets=max_subsets)
    data.update(original_source=source, column=column, input_permutation=before, output_permutation=after,
        physical_inside_mask=_permute_mask(data['mapped_source_mask'], after),
        physical_outside_mask=_permute_mask(data['mapped_outside_mask'], after),
        physical_protected_mask=_permute_mask(data['protected_mask'], after))
    for key in ('outside_counts', 'zero_fiber_counts'):
        data[key] = {_permute_mask(mask, after): count for mask, count in data[key].items()}
    return data


def protected_degree(mask, data):
    """Upper bound on original input weight in the actual zero-outside fibre."""
    if type(mask) is not int or mask < 0 or mask & ~data['physical_inside_mask']:
        raise ValueError('Physical zero-fibre mask must lie in the physical inside set')
    if not mask:
        return 0
    protected = data['physical_protected_mask']
    return (mask & protected).bit_count()+data['arity']-protected.bit_count()


def checked_dependencies():
    for name, expected in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed dependency changed: '+name)
    return {**outer.checked_dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}
