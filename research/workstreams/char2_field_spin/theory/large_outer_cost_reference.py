"""Fixed-coefficient upper circuits for large GF256 systematic GRS outers.

Arithmetic/site counts only: no native instruction schedule or timing claim.
Checks every proposed XOR-circulant entry by independent interpolation.
"""

import argparse
from fractions import Fraction
import json

from cauchy_joint_return import MUL, INV


def interpolation(dimension, targets):
    weights = []
    for source in range(dimension):
        denominator = 1
        for other in range(dimension):
            if source != other:
                denominator = MUL[denominator][source ^ other]
        weights.append(INV[denominator])
    result = []
    for point in targets:
        numerator = 1
        for other in range(dimension):
            numerator = MUL[numerator][point ^ other]
        result.append([MUL[MUL[numerator][INV[point ^ source]]][weights[source]]
                       for source in range(dimension)])
    return result


def zeta(values):
    result = list(values)
    for index in range((len(result)-1).bit_length()):
        bit = 1 << index
        for mask in range(len(result)):
            if not mask & bit:
                result[mask] ^= result[mask | bit]
    return result


def circulant_apply(values, kernel):
    transformed, gamma = zeta(values), zeta(kernel)
    result = [0]*len(values)
    for mask in range(len(values)):
        subset = mask
        while True:
            result[mask] ^= MUL[gamma[subset]][transformed[mask ^ subset]]
            if not subset:
                break
            subset = (subset-1) & mask
    return zeta(result)


def coefficient_basis(needed):
    span = {0: 0, 1: 1}
    rank = 1
    for coefficient in sorted(set(needed)-{0}):
        if coefficient not in span:
            span.update({value ^ coefficient: mask | (1 << rank)
                         for value, mask in tuple(span.items())})
            rank += 1
    # A coefficient with t basis terms needs ceil((t-1)/2) ternary XORs.
    reconstruction = sum(span[c].bit_count()//2 for c in set(needed)-{0})
    return rank-1, reconstruction


def circulant_cost(dimension, block_size, check=False):
    matrix = interpolation(dimension, range(dimension, dimension+block_size))
    blocks = dimension//block_size
    kernels = [matrix[0][j*block_size:(j+1)*block_size] for j in range(blocks)]
    for row in range(block_size):
        for column in range(dimension):
            assert matrix[row][column] == kernels[column//block_size][row ^ (column % block_size)]
    if check:
        for kernel in kernels:
            for column in range(block_size):
                values = [int(j == column) for j in range(block_size)]
                assert circulant_apply(values, kernel) == [kernel[row ^ column]
                                                          for row in range(block_size)]
    gammas = [zeta(values) for values in kernels]
    counts, reconstruction = [], 0
    for input_mask in range(block_size):
        needed = {gamma[subset] for gamma in gammas for subset in range(block_size)
                  if subset & input_mask == 0}
        products, xors = coefficient_basis(needed)
        counts.append(products)
        reconstruction += xors
    accumulation = sum(sum(1 for subset in range(block_size)
                           if subset & mask == subset and gamma[subset])//2
                       for gamma in gammas for mask in range(block_size))
    transforms = (blocks+1)*(block_size//2)*(block_size.bit_length()-1)
    return dict(products=sum(counts), counts=counts,
                reconstruction=reconstruction, accumulation=accumulation,
                transforms=transforms, booleans=reconstruction+accumulation+transforms,
                entry_checks=dimension*block_size, gamma0=[gamma[0] for gamma in gammas])


def cost(length, check=False):
    n = length//4
    missing = n//4
    dimension = 3*n-missing
    base = circulant_cost(3*n, n, check)
    shortening = circulant_cost(dimension, missing, check)
    products = base['products']+shortening['products']
    booleans = base['booleans']+shortening['booleans']+3*n
    return dict(outer=[length, dimension], base=base, shortening=shortening,
                total_products_generic=products, total_booleans_generic=booleans,
                state4_volume=str(Fraction(products, 4*dimension)+Fraction(length, dimension)),
                state8_volume=str(Fraction(products, 4*dimension)+Fraction(91*length, 64*dimension)),
                arithmetic_only=True)


def further_shorten(old, dimension):
    length, old_dimension = old['outer']
    matrix = interpolation(dimension, range(dimension, old_dimension))
    counts = [coefficient_basis(row) for row in matrix]
    extra = sum(item[0] for item in counts)
    reconstruction = sum(item[1] for item in counts)
    # Add every correction to the already available old transpose output.
    accumulation = dimension*((old_dimension-dimension+1)//2)
    products = old['total_products_generic']+extra
    return dict(outer=[length, dimension], via=old['outer'], extra_products=extra,
                extra_ranks=[item[0] for item in counts],
                extra_reconstruction_booleans=reconstruction,
                extra_accumulation_booleans=accumulation,
                total_products=products,
                total_booleans=old['total_booleans_generic']+reconstruction+accumulation,
                state4_volume=str(Fraction(products, 4*dimension)+Fraction(length, dimension)),
                state8_volume=str(Fraction(products, 4*dimension)+Fraction(91*length, 64*dimension)),
                arithmetic_only=True)


def state8_width32_cost(check=False):
    from state_capacity_reference import LABELS, BETA, coefficients
    phases = []
    for group in ((0, 1), (2, 3)):
        labels = sum((LABELS[phase] for phase in group), ())
        kernels = [coefficients(labels[first] ^ BETA, True)
                   for first in range(0, 32, 8)]
        gammas = [zeta(kernel) for kernel in kernels]
        assert all(gamma[0] == 1 for gamma in gammas)
        if check:
            for kernel in kernels:
                for column in range(8):
                    values = [int(j == column) for j in range(8)]
                    assert circulant_apply(values, kernel) == [kernel[row ^ column]
                                                              for row in range(8)]
        counts = [coefficient_basis({gamma[subset] for gamma in gammas
                                     for subset in range(8) if not subset & mask})
                  for mask in range(8)]
        products = sum(item[0] for item in counts)
        reconstruction = sum(item[1] for item in counts)
        accumulation = sum(sum(1 for subset in range(8)
                               if subset & mask == subset and gamma[subset])//2
                           for gamma in gammas for mask in range(8))
        transforms = 5*12
        booleans = reconstruction+accumulation+transforms
        packet = Fraction(products, 2)+12
        # Conservative binary-XOR bound: expand each ternary XOR to two.
        # Transposition changes the binary addition count by outputs-inputs.
        phases.append(dict(old_phases=group, scalar_products_per_map=products,
                           product_counts=[item[0] for item in counts],
                           emission_booleans=booleans,
                           pair_binary_xor_upper=4*booleans+24,
                           packet_GFNI512_volume=str(packet),
                           full_256_176_volume=str(Fraction(503, 704)+packet/Fraction(22))))
    return dict(width=32, state=8, phases=phases,
                mean_packet_GFNI512_volume="121/4",
                mean_full_256_176_volume="1471/704", arithmetic_only=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    for length in (32, 64, 128, 256):
        result = cost(length, args.check)
        if args.check:
            expected = {32: (34, 8), 64: (82, 19), 128: (186, 46), 256: (400, 103)}
            if length in expected:
                assert (result['base']['products'], result['shortening']['products']) == expected[length]
        print(json.dumps(result), flush=True)
        if length == 256:
            for dimension in (171, 170):
                print(json.dumps(further_shorten(result, dimension)), flush=True)
    packet = state8_width32_cost(args.check)
    if args.check:
        assert [phase['scalar_products_per_map'] for phase in packet['phases']] == [36, 37]
    print(json.dumps(packet), flush=True)
