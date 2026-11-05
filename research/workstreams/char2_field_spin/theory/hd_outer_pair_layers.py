"""Exact support laws for crossed random-ratio pair layers in one RS block.

This helper describes NEW ensembles, not a regional provider or distance
certificate.  Coordinates are (row,column) in {0,...,15}^2.  For each
matched pair, sample alpha independently from GF(q)* once at setup and
apply (a,b)->(a+alpha*b, a+c*alpha*b), where c is fixed outside {0,1}.
The existing independent nonzero output diagonal follows the last layer.

Conditioned on ANY prior values, one layer's output-mask law depends only
on its input mask.  Fresh alpha for each pair and layer therefore permits
exact mask composition without intermediate diagonal masks.  Output
amplitudes are not asserted independent until the final diagonal.
All inputs/messages are fixed before the setup probability is evaluated.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import comb

SIZE = 16
BASELINE_DELTAS = (2, 2, 2, 2)
VARIED_DELTAS = (2, 4, 8, 14)
COMPONENT_FIRST_PAIRS = ((0, 1), (2, 3))
COMPONENT_SECOND_PAIRS = ((0, 2), (1, 3))


def _field_order(q):
    if type(q) is not int or q < 4 or q & (q-1):
        raise ValueError('A characteristic-two field order q>=4 is required')


@lru_cache(maxsize=64)
def pair_mask_law(mask, q=256):
    """Exact probabilities in output-mask order00,01,10,11.

    Bit0 is the first coordinate.  For two active inputs, the two distinct
    cancellation keys are alpha=a/b and alpha=a/(c*b).  No input-amplitude
    distribution, final scaling, or extra orientation factor is used.
    """
    _field_order(q)
    if type(mask) is not int or not 0 <= mask < 4:
        raise ValueError('A two-coordinate support mask is required')
    if not mask:
        return (Fraction(1), Fraction(0), Fraction(0), Fraction(0))
    if mask != 3:
        return (Fraction(0), Fraction(0), Fraction(0), Fraction(1))
    return (Fraction(0), Fraction(1, q-1), Fraction(1, q-1), Fraction(q-3, q-1))


def layer_mask_law(mask, pairs, *, width=4, q=256):
    """Exact sparse mask distribution for one perfect matching."""
    _field_order(q)
    pairs = tuple(tuple(pair) for pair in pairs)
    if (type(width) is not int or width < 2 or width % 2 or type(mask) is not int
            or not 0 <= mask < 1 << width or any(len(pair) != 2 for pair in pairs)
            or sorted(i for pair in pairs for i in pair) != list(range(width))):
        raise ValueError('A labelled perfect matching and in-range mask are required')
    result = {0: Fraction(1)}
    for first, second in pairs:
        local = ((mask >> first) & 1) | (((mask >> second) & 1) << 1)
        following = Counter()
        for previous, probability in result.items():
            for output, mass in enumerate(pair_mask_law(local, q)):
                if mass:
                    placed = ((output & 1) << first) | (((output >> 1) & 1) << second)
                    following[previous | placed] += probability*mass
        result = dict(following)
    return result


@lru_cache(maxsize=64)
def component_mask_law(mask, q=256):
    """Two refreshed-ratio stages on one four-coordinate square.

    Returns16 rational coefficients of the multilinear OUTPUT support
    polynomial.  The two stages use four independent alpha keys total.
    No diagonal is needed between them.
    """
    answer = [Fraction(0) for _ in range(16)]
    for middle, first_mass in layer_mask_law(mask, COMPONENT_FIRST_PAIRS, q=q).items():
        for final, second_mass in layer_mask_law(middle, COMPONENT_SECOND_PAIRS, q=q).items():
            answer[final] += first_mass*second_mass
    if sum(answer) != 1 or (mask and answer[0]) or (not mask and answer[0] != 1):
        raise AssertionError('Invertibility or probability mass changed')
    return tuple(answer)


def weight_law(mask, q=256):
    result = [Fraction(0) for _ in range(5)]
    for output, mass in enumerate(component_mask_law(mask, q)):
        result[output.bit_count()] += mass
    return tuple(result)


def evaluate_component(mask, atoms, q=256):
    """Evaluate the exact labelled polynomial at nonnegative rational atoms."""
    atoms = tuple(Fraction(t) for t in atoms)
    if len(atoms) != 4 or min(atoms) < 0:
        raise ValueError('Four nonnegative rational output markers are required')
    answer = Fraction(0)
    for output, mass in enumerate(component_mask_law(mask, q)):
        term = mass
        for i, atom in enumerate(atoms):
            if output >> i & 1:
                term *= atom
        answer += term
    return answer


def backward_support_kernel(output_weights, q=256):
    """Exact expectation of a positive output function for each input mask."""
    weights = tuple(Fraction(value) for value in output_weights)
    if len(weights) != 16 or min(weights) < 0:
        raise ValueError('Sixteen nonnegative output-mask weights required')
    return tuple(sum(probability*weight for probability, weight in zip(component_mask_law(mask, q), weights))
                 for mask in range(16))


def backward_activity_kernel(input_weights, output_weights, *, q=256, column_labels=(0, 1, 2, 3)):
    """Positive transfer retaining ORIGINAL input-column activity flags.

    input_weights is caller-owned counting/probability mass on the actual
    pre-layer masks.  It is not a claim of independent H* coordinates.
    Each local coordinate belongs to one of four distinct columns.  The
    result's flag is an OR mask, not a final-output activity flag.  Combining
    components must OR their flags before imposing any column floor.
    """
    inputs = tuple(Fraction(value) for value in input_weights)
    labels = tuple(column_labels)
    if len(inputs) != 16 or min(inputs) < 0 or sorted(labels) != list(range(4)):
        raise ValueError('Sixteen positive input masses and four distinct column labels required')
    backward = backward_support_kernel(output_weights, q)
    result = [Fraction(0) for _ in range(16)]
    for mask, mass in enumerate(inputs):
        flag = sum(((mask >> i) & 1) << labels[i] for i in range(4))
        result[flag] += mass*backward[mask]
    return tuple(result)


def field_multiply(a, b, *, modulus=0x11b):
    """Tiny exact arithmetic for algebra tests, not a native implementation."""
    bits = modulus.bit_length()-1
    if bits < 2 or not 0 <= a < 1 << bits or not 0 <= b < 1 << bits:
        raise ValueError('Field elements must fit the specified binary modulus')
    answer = 0
    for _ in range(bits):
        if b & 1:
            answer ^= a
        b >>= 1
        a <<= 1
        if a >> bits:
            a ^= modulus
    return answer


def apply_pair(a, b, alpha, *, c=2, modulus=0x11b):
    q = 1 << (modulus.bit_length()-1)
    if not 0 < alpha < q or not 1 < c < q:
        raise ValueError('alpha must be nonzero and c outside{0,1}')
    scaled = field_multiply(alpha, b, modulus=modulus)
    return a ^ scaled, a ^ field_multiply(c, scaled, modulus=modulus)


def apply_component(values, alphas, *, c=2, modulus=0x11b):
    values, alphas = tuple(values), tuple(alphas)
    if len(values) != 4 or len(alphas) != 4:
        raise ValueError('Four input values and four independent pair keys required')
    out = list(values)
    for key, (i, j) in zip(alphas, COMPONENT_FIRST_PAIRS+COMPONENT_SECOND_PAIRS):
        out[i], out[j] = apply_pair(out[i], out[j], key, c=c, modulus=modulus)
    return tuple(out)


def _deltas(deltas):
    deltas = tuple(deltas)
    if len(deltas) != 4 or any(type(d) is not int or not 1 < d < 16 for d in deltas):
        raise ValueError('Four nonzero row deltas distinct from1 are required')
    return deltas


def one_layer_partner(coordinate):
    """Preferred single stage: row xor(8+floor(column/2)), column xor1."""
    if type(coordinate) is not int or not 0 <= coordinate < 256:
        raise ValueError('A16x16 coordinate is required')
    row, column = divmod(coordinate, 16)
    return 16*(row ^ (8+(column >> 1)))+(column ^ 1)


def one_layer_support_profile(row_mask, column_mask=65535):
    """Exact pair-mask counts for the explicit product support A×J."""
    if any(type(mask) is not int or not 0 <= mask < 1 << 16 for mask in (row_mask, column_mask)):
        raise ValueError('Two16-coordinate support masks are required')
    result = Counter()
    for point in range(256):
        other = one_layer_partner(point)
        if point < other:
            local = 0
            for i, coordinate in enumerate((point, other)):
                row, column = divmod(coordinate, 16)
                local |= bool((row_mask >> row & 1) and (column_mask >> column & 1)) << i
            result[local] += 1
    return dict(sorted(result.items()))


def one_layer_weight_law(row_mask, column_mask=65535, q=256):
    """Exact output-weight PGF for A×J, valid for arbitrary fixed NZ values."""
    _field_order(q)
    profile = one_layer_support_profile(row_mask, column_mask)
    singles, doubles = profile.get(1, 0)+profile.get(2, 0), profile.get(3, 0)
    result = [Fraction(0) for _ in range(257)]
    for full in range(doubles+1):
        result[2*singles+doubles+full] += (comb(doubles, full)
            * Fraction(2, q-1)**(doubles-full)*Fraction(q-3, q-1)**full)
    if sum(result) != 1:
        raise AssertionError('The pair-product weight law lost mass')
    return tuple(result)


def partner(coordinate, stage, second_row_deltas=VARIED_DELTAS):
    """Integer coordinate16*row+column; stage0=(xor1,xor1).

    Stage1 is (row xor d_b,column xor2), where b=column//4.  The optional
    varied schedule changes wiring only, and is a distinct construction.
    """
    deltas = _deltas(second_row_deltas)
    if type(coordinate) is not int or not 0 <= coordinate < 256 or stage not in (0, 1):
        raise ValueError('A16x16 coordinate and one of two stages are required')
    row, column = divmod(coordinate, 16)
    delta, shift = (1, 1) if stage == 0 else (deltas[column // 4], 2)
    return 16*(row ^ delta)+(column ^ shift)


def components(second_row_deltas=VARIED_DELTAS):
    """The64 disjoint squares, in local order0,first,second,both."""
    deltas = _deltas(second_row_deltas)
    seen, answer = set(), []
    for point in range(256):
        if point in seen:
            continue
        first, second = partner(point, 0, deltas), partner(point, 1, deltas)
        both = partner(first, 1, deltas)
        if both != partner(second, 0, deltas):
            raise AssertionError('The local matching square does not commute')
        component = (point, first, second, both)
        if len(set(component)) != 4 or seen.intersection(component):
            raise AssertionError('Components overlap or fail to have four positions')
        seen.update(component)
        answer.append(component)
    return tuple(answer)


def closed_rowset(mask, second_row_deltas=VARIED_DELTAS):
    """Whether A×all16columns is invariant under both pair matchings."""
    deltas = _deltas(second_row_deltas)
    if type(mask) is not int or not 0 <= mask < 1 << 16:
        raise ValueError('A16-row support mask is required')
    generators = {1, *deltas}
    return all(bool(mask >> row & 1) == bool(mask >> (row ^ d) & 1)
               for d in generators for row in range(16))


def support_profile(row_mask, second_row_deltas=VARIED_DELTAS):
    """Counts local input-mask types for A×all16columns; no RS counting."""
    if type(row_mask) is not int or not 0 <= row_mask < 1 << 16:
        raise ValueError('A16-row support mask is required')
    return dict(sorted(Counter(sum(bool(row_mask >> (point // 16) & 1) << i
                                   for i, point in enumerate(component))
                               for component in components(second_row_deltas)).items()))


def primitive_counts(layers=2, *, precompute_c_alpha=True):
    """Literal pair circuit counts, before any native fusion or gauge saving."""
    if type(layers) is not int or layers < 1:
        raise ValueError('A positive layer count is required')
    pairs = 128*layers
    return dict(pairs=pairs, fresh_nonzero_alpha_keys=pairs, binary_xors=2*pairs,
                variable_products=(2 if precompute_c_alpha else 1)*pairs,
                fixed_products=(0 if precompute_c_alpha else 1)*pairs,
                final_diagonal_products=256, final_diagonal_keys=256,
                product_total_including_final_diagonal=2*pairs+256,
                includes_Hstar=False, includes_D1=False, includes_RS=False,
                native_performance_claim=False)
