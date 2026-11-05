"""Grouped MDS outer interface; exact toy checks, not a native benchmark.

See HD_GROUPED_OUTER.md for the conditional-information-set proof. The
envelope concerns expected counts of all messages under one sampled setup.
"""
from fractions import Fraction
from itertools import product


def atom_cap(field=65536, group=4, length=64, dimension=32):
    if field < 2 or group < 1 or not 0 < dimension < length:
        raise ValueError('Invalid grouped MDS dimensions')
    return Fraction(1, (field-1)**(group*(length-dimension)))


def envelope(field=65536, group=4, length=64, dimension=32):
    return field**(group*length)*atom_cap(field, group, length, dimension)


def gf4_mul(a, b):
    result = 0
    for _ in range(2):
        if b & 1: result ^= a
        b >>= 1
        a <<= 1
        if a & 4: a ^= 7
    return result


def toy_mixer_counts(two_diagonals=True, normalized=False):
    """All D2 H D1 maps, H interpolating nodes0,1 at2,3 over GF4.

    Vectors (x0,x1) are packed as x0+4*x1. With one diagonal omitted,
    this becomes a deliberate negative control, not an approved family.
    """
    counts = [[0]*16 for _ in range(16)]
    before = tuple(product(range(1,4),repeat=2)) if two_diagonals else ((1,1),)
    if normalized:
        if not two_diagonals:raise ValueError('Normalization is not diagonal removal')
        before=tuple((1,b) for b in range(1,4))
    for x in range(16):
        for a,b in before:
            u,v=gf4_mul(a,x&3),gf4_mul(b,x>>2)
            low=gf4_mul(3,u)^gf4_mul(2,v)
            high=gf4_mul(2,u)^gf4_mul(3,v)
            for c,d in product(range(1,4),repeat=2):
                y=gf4_mul(c,low)+4*gf4_mul(d,high)
                counts[x][y]+=1
    return counts, len(before)*9


def field_multiplier_cap(field, group, length, dimension):
    if field<2 or group<1 or not 0<dimension<length:raise ValueError('Invalid dimensions')
    return Fraction(1,(field**group-1)**(length-dimension))


def toy_group_counts(two_diagonals=True):
    """Exact expected atoms for two GF4 RS[4,2] rows, including correlations.

    Sum the 255 nonzero message images over independent column ensembles.
    Return integer numerators over the common denominator, not floats.
    """
    counts,den=toy_mixer_counts(two_diagonals)
    support=[tuple((y,c) for y,c in enumerate(row) if c) for row in counts]
    atoms=[0]*65536
    for a,b,c,d in product(range(4),repeat=4):
        if not(a or b or c or d):continue
        columns=[(a^gf4_mul(b,j))+4*(c^gf4_mul(d,j)) for j in range(4)]
        for y0,c0 in support[columns[0]]:
            for y1,c1 in support[columns[1]]:
                low=y0+16*y1; mass=c0*c1
                for y2,c2 in support[columns[2]]:
                    mid=low+256*y2; m=mass*c2
                    for y3,c3 in support[columns[3]]:
                        atoms[mid+4096*y3]+=m*c3
    return atoms,den**4
