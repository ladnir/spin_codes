"""Exact fixed-message support kernels for prospective packed binary mixers.

These are local distribution lemmas, not a SPIN distance certificate. All
random maps below are independent per block/group and sampled once at setup.
CDF transfers require cumulative count caps, not arbitrary shell majorants.
Only the Python standard library is used.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, prod


def integer(value, low, high, name):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f'{name} must be an integer in [{low},{high}]')
    return value


def gaussian(n, r):
    """Number of r-dimensional subspaces of binary n-space."""
    if not 0 <= r <= n:
        return 0
    return prod((1 << n) - (1 << i) for i in range(r)) // prod(
        (1 << r) - (1 << i) for i in range(r))


def spanning_nonzero_columns(rank, columns):
    """Ordered nonzero columns in F2^rank spanning the whole space."""
    return sum((-1)**(rank-j) * (1 << ((rank-j)*(rank-j-1)//2))
               * gaussian(rank, j) * ((1 << j)-1)**columns
               for j in range(rank+1))


@lru_cache(None, typed=True)
def shared_gl(width, rank):
    """Union support after one uniform GL(width,2), shared across rows.

    The fixed input's row rank is `rank`. Its image subspace in the row
    alphabet is preserved, so this is not a rank-free packet randomizer.
    """
    integer(width, 1, 256, 'width')
    integer(rank, 0, width, 'rank')
    if rank == 0:
        return (Q(1),) + (Q(0),)*width
    denominator = prod((1 << width)-(1 << i) for i in range(rank))
    return tuple(Q(comb(width, w)*spanning_nonzero_columns(rank, w), denominator)
                 for w in range(width+1))


@lru_cache(None, typed=True)
def independent_rows(width, active_rows):
    """Independent uniform GL(width,2) on each nonzero row fragment.

    Independent uniform nonzero GF(2^width) scalars have the SAME law for
    each fixed input. Sharing a scalar across rows does not have this law.
    """
    integer(width, 1, 256, 'width')
    integer(active_rows, 0, 16, 'active_rows')
    if active_rows == 0:
        return (Q(1),) + (Q(0),)*width
    denominator = ((1 << width)-1)**active_rows
    return tuple(Q(comb(width, w)*sum(
        (-1)**j * comb(active_rows, j) * ((1 << (active_rows-j))-1)**w
        for j in range(active_rows+1)), denominator) for w in range(width+1))


@lru_cache(None, typed=True)
def full_block(width, rows=4):
    """Support of a nonzero block after uniform GL(rows*width,2).

    A uniform nonzero GF(2^(rows*width)) scalar has exactly this fixed-input
    law. Conditional on the complete packet support, all active packet
    labels are independent uniform nonzero rows-bit vectors. No extra
    packet randomizer is needed if this is the final local random map.
    """
    integer(width, 1, 256, 'width')
    integer(rows, 1, 16, 'rows')
    denominator = (1 << (rows*width))-1
    return (Q(0),) + tuple(Q(comb(width, w)*((1 << rows)-1)**w, denominator)
                           for w in range(1, width+1))


@lru_cache(None, typed=True)
def independent_pieces(width, piece_rows, active_pieces):
    """Independent uniform full maps on row-pieces, e.g. two GL16 row-pairs."""
    integer(width, 1, 256, 'width')
    integer(piece_rows, 1, 16, 'piece_rows')
    integer(active_pieces, 0, 16, 'active_pieces')
    if active_pieces == 0:
        return (Q(1),) + (Q(0),)*width
    denominator = ((1 << (piece_rows*width))-1)**active_pieces
    return tuple(Q(comb(width, w)*sum(
        (-1)**(w-j)*comb(w, j)*((1 << (piece_rows*j))-1)**active_pieces
        for j in range(w+1)), denominator) for w in range(width+1))


def checked_pmf(values):
    values = tuple(map(Q, values))
    if not values or min(values) < 0 or sum(values) != 1:
        raise ValueError('an exact normalized nonnegative PMF is required')
    return values


def convolve(left, right):
    result = [Q(0)]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    result[i+j] += a*b
    return tuple(result)


def repeated(kernel, copies):
    checked_pmf(kernel)
    integer(copies, 0, 256, 'copies')
    result = (Q(1),)
    for _ in range(copies):
        result = convolve(result, kernel)
    return result


def mean(kernel):
    return sum((i*p for i, p in enumerate(checked_pmf(kernel))), Q(0))


def cdf(kernel):
    total, result = Q(0), []
    for p in checked_pmf(kernel):
        total += p
        result.append(total)
    return tuple(result)


def moment(kernel, z):
    z = Q(z)
    if not 0 <= z <= 1:
        raise ValueError('support moment requires 0 <= z <= 1')
    return sum((p*z**i for i, p in enumerate(checked_pmf(kernel))), Q(0))


@lru_cache(None, typed=True)
def random_partition(total, width, support):
    """Exact number of hit width-column blocks after a fresh uniform partition.

    Equivalently, shuffle a fixed support uniformly BEFORE canonical blocks.
    The normal SPIN shuffle AFTER mixing cannot justify this distribution.
    """
    integer(total, 1, 4096, 'total')
    integer(width, 1, total, 'width')
    integer(support, 0, total, 'support')
    if total % width:
        raise ValueError('width must divide total')
    groups = total//width
    denominator = comb(total, support)
    values = []
    for h in range(groups+1):
        count = comb(groups, h)*sum(
            (-1)**(h-j)*comb(h, j)*comb(width*j, support)
            for j in range(h+1) if width*j >= support)
        values.append(Q(count, denominator))
    return checked_pmf(values)


def support_moment_kernel(total, width, local, z, random_partition_first=False):
    """Upper support moments from total union support alone.

    For full-block maps, canonical blocks have at least ceil(u/width)
    active blocks. For shared/independent-row maps, using shared_gl(width,1)
    supplies a conservative stochastic lower support law for every nonempty
    block. Random-partition mode requires a genuinely independent preceding
    permutation. Neither version assumes ranks/activity profiles from u.
    """
    integer(total, 1, 4096, 'total')
    integer(width, 1, total, 'width')
    if total % width or len(local) != width+1:
        raise ValueError('matching local width and integral block count required')
    if type(random_partition_first) is not bool:
        raise ValueError('random_partition_first must be bool')
    local = checked_pmf(local)
    if local[0]:
        raise ValueError('each nonempty input block must stay nonempty')
    f = moment(local, z)
    if random_partition_first:
        return tuple(sum((p*f**h for h, p in enumerate(random_partition(total, width, u))), Q(0))
                     for u in range(total+1))
    return tuple(f**((u+width-1)//width) for u in range(total+1))


def transfer_cdf_moment(caps, kernel):
    """Abel/summation-by-parts bound for decreasing nonnegative kernel.

    Input caps[u] bound cumulative NONZERO-message counts through support u.
    The final cap bounds their total mass. CDF differences are used only
    against this checked decreasing kernel; they are NOT shell bounds.
    The caller remains responsible for authenticating the input CDF caps.
    """
    caps, kernel = tuple(map(Q, caps)), tuple(map(Q, kernel))
    if len(caps) < 2 or len(caps) != len(kernel) or caps[0] != 0:
        raise ValueError('matching nonzero-message CDF and kernel required')
    if any(a > b or a < 0 for a, b in zip(caps, caps[1:])) or caps[-1] < 0:
        raise ValueError('nonnegative nondecreasing CDF caps required')
    if any(a < b or b < 0 for a, b in zip(kernel, kernel[1:])) or kernel[0] < 0:
        raise ValueError('nonnegative nonincreasing kernel required')
    return sum(((caps[u]-caps[u-1])*kernel[u] for u in range(1, len(caps))), Q(0))


def two_shear_row(q, first_active, second_active):
    """Exact row-activity kernel for A=x+a*y, B=y+b*A over GF(q).

    a,b are independent uniform nonzero scalars. Return probabilities in
    mask order 00,10,01,11, using bit0 for A and bit1 for B. A field of
    this cardinality must actually be supplied by the implementation.
    """
    integer(q, 2, 1 << 32, 'q')
    if q & (q-1):
        raise ValueError('this model requires a binary extension field')
    if type(first_active) is not bool or type(second_active) is not bool:
        raise ValueError('activity indicators must be bool')
    d = q-1
    if not first_active and not second_active:
        return (Q(1), Q(0), Q(0), Q(0))
    if first_active and not second_active:
        return (Q(0), Q(0), Q(0), Q(1))
    if not first_active:
        return (Q(0), Q(1, d), Q(0), Q(d-1, d))
    return (Q(0), Q(d-1, d*d), Q(1, d), Q((d-1)**2, d*d))


def two_shear_masks(q, first_mask, second_mask, rows=4):
    """Exact transition on the two row-activity masks; fresh scalars per row."""
    integer(rows, 1, 8, 'rows')
    integer(first_mask, 0, (1 << rows)-1, 'first_mask')
    integer(second_mask, 0, (1 << rows)-1, 'second_mask')
    result = {(0, 0): Q(1)}
    for i in range(rows):
        law = two_shear_row(q, bool(first_mask & (1 << i)), bool(second_mask & (1 << i)))
        next_result = {}
        for (a, b), mass in result.items():
            for outcome, probability in enumerate(law):
                if probability:
                    pair = (a | ((outcome & 1) << i), b | (((outcome >> 1) & 1) << i))
                    next_result[pair] = next_result.get(pair, Q(0)) + mass*probability
        result = next_result
    return result


if __name__ == '__main__':
    print('Exact local laws only; no whole-code certificate or timing claim.')
    for rank in range(1, 5):
        print('shared GL8 rank', rank, 'mean', float(mean(shared_gl(8, rank))))
        print('row-wise GL8/GF256 active rows', rank, 'mean', float(mean(independent_rows(8, rank))))
    print('full GL32/GF2^32 nonempty block mean', float(mean(full_block(8))))
    for u in (32, 64, 112, 128, 192, 256):
        random_hits = mean(random_partition(256, 8, u))
        print('support', u, 'canonical minimum hits', (u+7)//8,
              'random-partition mean hits', float(random_hits),
              'full-block mean canonical minimum', float(((u+7)//8)*mean(full_block(8))),
              'full-block mean with random partition', float(random_hits*mean(full_block(8))))
