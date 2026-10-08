"""Exact small-column spectra for the dominant identity refresh branch.

The matrix G=F+A C is fixed over GF256 and invertible over GF65536.
For flat-state propagation, its zero-coset moment bounds every translated
column. These local facts are not a whole-code distance certificate.
"""

from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb, log2

from hd_blockwise_shared_maps import matrices, multiplication_table
from hd_blockwise_maps import matmul, power
from hd_randomized_outer_state_screen import add, total, subtract


def validate_field(field):
    if not isinstance(field, int) or field < 256 or field & (field-1) or (field.bit_length()-1) % 8:
        raise ValueError('Need a power of 256, the cardinality of a GF256 extension')


@lru_cache(None)
def base_inverses():
    return tuple(0 if not x else power(x, 254) for x in range(256))


@lru_cache(None)
def coupled_matrix():
    feedthrough, emission, feedback, unused, unused2 = matrices(64)
    perturbation = matmul(emission, feedback)
    return tuple(tuple(a ^ b for a, b in zip(row, change))
                 for row, change in zip(feedthrough, perturbation))


def pair_projective_weights(first, second, field=65536):
    """Projective weight histogram for two independent GF256 columns.

    Representatives are (1,t) for t in GF65536 and (0,1). All cancellations
    in a nonzero coordinate require t in GF256. The other q-256 slopes have
    full union support. No extension-field enumeration is needed.
    """
    validate_field(field)
    if len(first) != len(second):
        raise ValueError('Need equal column lengths')
    table = multiplication_table()
    inverses = base_inverses()
    zero_slopes = Counter()
    support = 0
    for a, b in zip(first, second):
        if a or b:
            support += 1
        if b:
            zero_slopes[table[a][inverses[b]]] += 1
    weights = Counter({support: field-256})
    for slope in range(256):
        weights[support-zero_slopes[slope]] += 1
    weights[sum(x != 0 for x in second)] += 1
    return +weights


@lru_cache(None, typed=True)
def aggregate_enumerator(potential, field=65536):
    """Integer weight counts summed over all supports of the declared size."""
    validate_field(field)
    if potential not in (0, 1, 2):
        raise ValueError('Only zero, one and two columns are implemented')
    matrix = coupled_matrix()
    if potential == 0:
        return (1,)+(0,)*64
    columns = tuple(zip(*matrix))
    coefficients = Counter({0: comb(64, potential)})
    if potential == 1:
        for column in columns:
            coefficients[sum(x != 0 for x in column)] += field-1
    else:
        for first, second in combinations(columns, 2):
            for weight, lines in pair_projective_weights(first, second, field).items():
                coefficients[weight] += (field-1)*lines
    return tuple(coefficients[weight] for weight in range(65))


def moment_log(potential, log_z, field=65536):
    counts = aggregate_enumerator(potential, field)
    return (total(log2(count)+weight*log_z for weight, count in enumerate(counts) if count)
            -log2(comb(64, potential))-potential*log2(field))


def refine_columns(local, log_z, field=65536):
    """Replace flat rows' epsilon gamma by a gamma_I+b gamma_L.

    For j<=2 retain the stronger of this exact I-column spectrum and the
    previous column cap. Other entries and every j>2 operator are unchanged.
    Input local is a three-component operator already averaged over subsets.
    """
    # Recover the old gamma contribution from its two mixture terms without
    # subtracting floating upper bounds: compute its explicit occupancy law.
    d = field**4-1
    a = (field**3-1)/d
    b = (field-1)/(field*d)
    c = 1-1/field
    log_q = log2(field)
    result = list(local)
    for j in (1, 2):
        column = moment_log(j, log_z, field)
        t = add(-log_q, log2(1-1/field)+log_z)
        if j == 1:
            occupancies = ((1, 1.0),)
        else:
            same = 8*comb(8, 2)/comb(64, 2)
            occupancies = ((1, same), (2, 1-same))
        gamma = total(log2(probability)+min(j*t, max(0, 64-8*u-3)*log_z)
                      for u, probability in occupancies)
        gamma0 = total(log2(probability)+min(subtract(j*t, -j*log_q),
                        max(0, 64-8*u-3)*log_z+log2(1-field**(-j)))
                       for u, probability in occupancies)
        # B->B = c V + epsilon gamma in the old operator. Instead of undoing
        # that sum, recover V by direct punctured-state/tile enumeration.
        from hd_randomized_outer_state_profile import punctured_state_moment
        from hd_randomized_outer_state_screen import mds_moment
        terms = []
        for u, probability in occupancies:
            tile = mds_moment(j, log_z, field) if u == 1 else 2*mds_moment(1, log_z, field)
            terms.append(log2(probability)+tile+punctured_state_moment(64-8*u, log_z, field))
        uniform = total(terms)
        identity_gamma = min(gamma, column)
        identity_gamma0 = min(gamma0, subtract(column, -j*log_q))
        new_bulk = total((log2(c)+uniform, log2(a)+identity_gamma, log2(b)+gamma))
        new_zero = total((log2(c)+uniform, log2(a)+identity_gamma0, log2(b)+gamma0))-log2(d)
        rows = [list(row) for row in local[j]]
        # A caller may have improved V separately. This whole-row alternative
        # cannot undo that improvement; take entrywise minima of valid bounds.
        rows[1][0] = min(rows[1][0], new_zero)
        rows[1][1] = min(rows[1][1], new_bulk)
        result[j] = tuple(tuple(row) for row in rows)
    return tuple(result)


if __name__ == '__main__':
    import json
    for potential in (1, 2):
        counts = aggregate_enumerator(potential)
        print(json.dumps(dict(potential=potential,
                              minimum_nonzero_weight=min(w for w in range(1, 65) if counts[w]),
                              counts_sum=sum(counts),
                              expected_sum=comb(64, potential)*65536**potential,
                              nonzero_weights=[w for w in range(1, 65) if counts[w]],
                              whole_code_certificate=False)))
