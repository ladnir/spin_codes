"""Exact labelled fibre factors for the scheduled Hlambda outer.

Fix a mapped input set A and its complement O before the output permutation.
The symmetric involution Hlambda has cross map C=Hlambda[O,A]. Its raw
allowed-input column enumerator is K_A(t)=sum_x prod_i t_i**[Hx_i != 0].
This sum allows zero coordinates in A; it is not a forced-support law.

Write Z_A for the support enumerator of ker C on A, and C_A for the
support enumerator of im C on O. For 0 <= t_i <= 1,

                         K_A(t) <= C_A(t_O) Z_A(t_A).

Indeed, Hlambda maps ker C bijectively to itself since Hlambda**2=I.
Every inside output fibre is consequently a translate of ker C. The
additive Fourier factors of the product weight are 1+(q-1)t_i and 1-t_i.
They are nonnegative, so the zero translate maximizes the weighted sum.
This proves a pointwise inequality, not a coefficientwise inequality.

The exact allowed-support counts are q**(c-rank C[O\\S,A]) outside
and q**(|S|-rank C[O,S]) inside, where c=rank C. Subset Mobius inversion
produces labelled support counts. The factor masses are q**c and
q**(|A|-c); their product is q**|A|, without an extra fibre multiplier.

Rank computations use binary matrices. For binary B and theta outside
{0,1}, rank(B+theta*J)=max(rank_F2(B),rank_F2(B+J)). Every minor is an
affine polynomial in theta with binary coefficients, which proves this
identity. For disjoint cross coordinates, affine hulls reduce B to at
most 5x5. No affine symmetry of the actual labelled Hlambda is assumed.

The scheduled API retains the original A, applies its input permutation
in each column, and transports BOTH factors through the output permutation.
It drops original-input/activity floors. No input-weight marker, RS factor,
regional provider, numerical gate, native change, or certificate is supplied.
"""
from functools import lru_cache
from math import prod

from hd_bit_syndrome_flat_geometry import flat_hull_mask
from hd_pair_scheduled_outer import GF16_POLYNOMIAL, column_permutations


ROWS = 16
FULL = (1 << ROWS)-1
SCHEDULE = 'input-inverse-output-xor'
LAMBDA_AES = 0xbc
CONSTRUCTION = 'G16-RS16x8-scheduled-GF4-Hlambda-no-outer-cross-research'


def _mask(mask):
    if type(mask) is not int or not 0 <= mask <= FULL:
        raise ValueError('A sixteen-coordinate integer mask is required')
    return mask


def _source(source):
    source = tuple(source)
    if (len(source) > ROWS or len(set(source)) != len(source)
            or any(type(i) is not int or not 0 <= i < ROWS for i in source)):
        raise ValueError('At most sixteen distinct ordered source labels are required')
    return source


def _parameters(q, max_subsets):
    if type(q) is not int or q not in (4, 256):
        raise ValueError('Use GF4 toy counts or AES-GF256 counts')
    if type(max_subsets) is not int or not 1 <= max_subsets <= 65536:
        raise ValueError('The explicit per-factor subset budget must lie in 1..65536')


def _labels(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask ^= bit


def _binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            lead = row.bit_length()-1
            if lead in pivots:
                row ^= pivots[lead]
            else:
                pivots[lead] = row
                break
    return len(pivots)


def binary_rank_one_rank(binary_rows, columns):
    """Exact rank of B+theta*J for any theta outside the binary subfield.

    B is supplied as integer-packed rows with `columns` bits. This generic
    helper does not assume Hlambda structure, disjoint sets, or square B.
    """
    rows = tuple(binary_rows)
    if (type(columns) is not int or not 0 <= columns <= ROWS
            or len(rows) > ROWS
            or any(type(row) is not int or not 0 <= row < 1 << columns for row in rows)):
        raise ValueError('At most sixteen binary rows and columns are required')
    ones = (1 << columns)-1
    return max(_binary_rank(rows), _binary_rank(row ^ ones for row in rows))


@lru_cache(maxsize=308)
def _affine_basis(flat):
    pivots, basis = {}, []
    for label in _labels(flat):
        value = ROWS | label
        while value:
            lead = value.bit_length()-1
            if lead in pivots:
                value ^= pivots[lead]
            else:
                pivots[lead] = value
                basis.append(label)
                break
    return tuple(basis)


@lru_cache(maxsize=47586)
def _rank_hulls(left, right):
    rows, columns = _affine_basis(left), _affine_basis(right)
    # Hlambda-I = (1+binary_dot) + (lambda+1)*J. Hull replacement
    # preserves the GF256 spans of the augmented binary label vectors.
    binary = tuple(sum((1 ^ ((i & j).bit_count() & 1)) << k
                       for k, j in enumerate(columns)) for i in rows)
    return binary_rank_one_rank(binary, len(columns))


def rank_cross_masks(row_mask, column_mask):
    """rank((Hlambda-I)[R,S]); equals rank(Hlambda[R,S]) if disjoint.

    Overlap is allowed ONLY with this bilinear meaning. The diagonal I is
    not included. The lazy symmetric cache is keyed by both affine hulls.
    """
    _mask(row_mask)
    _mask(column_mask)
    if not row_mask or not column_mask:
        return 0
    left, right = flat_hull_mask(row_mask), flat_hull_mask(column_mask)
    return _rank_hulls(min(left, right), max(left, right))


def _subset_masks(labels):
    result = [0]
    for label in labels:
        bit = 1 << label
        result += [mask | bit for mask in result]
    return result


def _exact_supports(labels, allowed_counts, mass):
    counts = list(allowed_counts)
    if len(counts) != 1 << len(labels):
        raise AssertionError('The allowed-support table has the wrong shape')
    for bit in range(len(labels)):
        step = 1 << bit
        for base in range(0, len(counts), 2*step):
            for offset in range(step):
                counts[base+step+offset] -= counts[base+offset]
    if min(counts) < 0 or counts[0] != 1 or sum(counts) != mass:
        raise AssertionError('Exact support inversion has invalid mass or coefficients')
    return tuple((mask, count) for mask, count in zip(_subset_masks(labels), counts) if count)


def _support_enumerators(source, outside, cross_rank, rank_cross, q):
    """Internal generic rank interface, also used by tiny involution tests."""
    amask, omask = sum(1 << i for i in source), sum(1 << i for i in outside)
    outside_subsets, inside_subsets = _subset_masks(outside), _subset_masks(source)
    outside_allowed = (q**(cross_rank-rank_cross(omask ^ subset, amask))
                       for subset in outside_subsets)
    inside_allowed = (q**(subset.bit_count()-rank_cross(omask, subset))
                      for subset in inside_subsets)
    return (_exact_supports(outside, outside_allowed, q**cross_rank),
            _exact_supports(source, inside_allowed, q**(len(source)-cross_rank)))


@lru_cache(maxsize=512)
def _fiber_cached(source_mask, q):
    source = tuple(_labels(source_mask))
    outside = tuple(_labels(FULL ^ source_mask))
    cross_rank = rank_cross_masks(FULL ^ source_mask, source_mask)
    outside_counts, zero_counts = _support_enumerators(
        source, outside, cross_rank, rank_cross_masks, q)
    return cross_rank, outside_counts, zero_counts


def fiber_counts(source, *, q=256, max_subsets=4096):
    """Exact pre-Q support factors for an allowed mapped source set.

    Return labelled mask-to-integer dictionaries outside_counts and
    zero_fiber_counts. Counts are cached only as immutable tuples; callers
    receive fresh dictionaries. Both factor lengths must meet max_subsets.
    """
    source = _source(source)
    _parameters(q, max_subsets)
    arity = len(source)
    if max(1 << arity, 1 << (ROWS-arity)) > max_subsets:
        raise ValueError('The source exceeds the explicit per-factor subset budget')
    amask = sum(1 << i for i in source)
    cross_rank, outside, zero = _fiber_cached(amask, q)
    return dict(source=source, source_mask=amask, outside_mask=FULL ^ amask,
        arity=arity, cross_rank=cross_rank, fiber_dimension=arity-cross_rank,
        outside_counts=dict(outside), zero_fiber_counts=dict(zero),
        outside_mass=q**cross_rank, zero_fiber_mass=q**(arity-cross_rank),
        raw_column_mass=q**arity, field_order=q, metadata=metadata())


def _permute_mask(mask, permutation):
    return sum(1 << permutation[i] for i in _labels(mask))


def column_fiber_counts(original_source, column, *, schedule=SCHEDULE,
                        q=256, max_subsets=4096):
    """Apply P_j before Hlambda and Q_j to BOTH output factor masks.

    All returned count masks use physical output labels. Original source
    ordering is retained as metadata, never averaged or replaced by an orbit.
    The default is the frozen inverse-input/output-xor schedule, without
    an outer pair layer. The other existing schedules are explicit controls.
    """
    original_source = _source(original_source)
    before, after = column_permutations(schedule, column)
    mapped = tuple(before[i] for i in original_source)
    data = fiber_counts(mapped, q=q, max_subsets=max_subsets)
    data.update(original_source=original_source,
                original_source_mask=sum(1 << i for i in original_source),
                mapped_source=mapped, mapped_source_mask=data['source_mask'],
                physical_inside_mask=_permute_mask(data['source_mask'], after),
                physical_outside_mask=_permute_mask(data['outside_mask'], after),
                input_permutation=before, output_permutation=after,
                schedule=schedule, column=column)
    for key in ('outside_counts', 'zero_fiber_counts'):
        data[key] = {_permute_mask(mask, after): count for mask, count in data[key].items()}
    data['metadata'].update(schedule=schedule, column=column,
        output_masks_are_physical=True, original_A_fixed_across_columns=True,
        source_order_preserved=True, label_field_polynomial=GF16_POLYNOMIAL)
    return data


def evaluate_fiber_bound(data, atoms):
    """Evaluate C_A(t)Z_A(t) for sixteen nonnegative atoms at most one.

    Fractions remain exact. FLOAT inputs give only a FLOAT evaluation.
    This evaluator contains no r-th moment or original-input marker.
    """
    atoms = tuple(atoms)
    if len(atoms) != ROWS or any(not 0 <= value <= 1 for value in atoms):
        raise ValueError('Sixteen atoms in [0,1] are required')
    def evaluate(counts):
        return sum(count*prod(atoms[i] for i in _labels(mask))
                   for mask, count in counts.items())
    return evaluate(data['outside_counts'])*evaluate(data['zero_fiber_counts'])


def metadata():
    return dict(construction=CONSTRUCTION, lambda_aes_byte=LAMBDA_AES,
        lambda_gf4_symbol=2, exact_labelled_counts=True,
        raw_allowed_input_counts=True, forced_input_support=False,
        pointwise_atom_domain='[0,1]', coefficientwise_domination_claimed=False,
        original_input_floor_included=False, input_weight_marker_supported=False,
        rs_normalization_included=False, source_count_included=False,
        crossed_pair_layer_included=False, regional_provider_runs=0,
        whole_code_certificate=False)
