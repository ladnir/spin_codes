"""Directed complete-class replay for the Cauchy/diagonal-refresh variant.

The CHANGED inner has four fresh nonzero GF256 state scalars per packet,
superregular Cauchy-form emission, and dual injection. This is not a
whole-code certificate. It covers every inverse-input weight in [.11N,.24N);
the independent dense theorem covers [.24N,N], and lower weights remain open.
"""

import argparse
from fractions import Fraction
import json
import math

from diagonal_refresh_rank import rank_potential_bound
from return_cover_certificate import (Witness, WITNESSES as OLD_WITNESSES,
                                      interval, row_upper, terminal_power)
from subfield_certificate import ONE, ZERO, integer_margin


WITNESSES = tuple((Fraction(density), Witness(*(Fraction(value) for value in values)))
    for density, values in (
        ('11/100', ('32627/9127', '44/291', '54371/3437', '2722/4197')),
        ('3/25', ('40279/8997', '1080/8903', '164117/9570', '221/388')),
        ('13/100', ('30545/5532', '871/8874', '181840/9657', '2941/5883')),
        ('7/50', ('63199/9375', '557/6975', '83327/3988', '3676/8385')),
        ('3/20', ('67482/8233', '553/8508', '163349/6975', '3159/8230')),
        ('4/25', ('13579/1361', '499/9477', '218523/8237', '2837/8465')),
        ('17/100', ('85275/6982', '176/4165', '251913/8281', '2023/6943')),
        ('9/50', ('136405/9038', '319/9530', '109699/3100', '2248/8929')))) + OLD_WITNESSES


def transfer(witness):
    x, y, tau = (interval(value) for value in (witness.x, witness.y, witness.tau))
    if not witness.x > 0 or not 0 < witness.y <= 1 or not 1 <= witness.tau <= 255:
        raise ValueError('invalid witness')
    xp, yp = [x**i for i in range(17)], [y**i for i in range(17)]
    gp = [(ONE/tau)**i for i in range(21)]
    gamma = [interval(max(rank_potential_bound(r, z, witness.tau)
                          for r in range(1, 5))) for z in range(4)]
    matrix = [[ONE, ZERO], [ZERO, ZERO]]
    for w in range(1, 17):
        term = math.comb(16, w)*xp[w]*yp[w]
        matrix[0][1] += term
        if w >= 5:
            matrix[0][0] += term*gp[4]
    for z in range(4):
        for only_input in range(17-z):
            for both in range(17-z-only_input):
                only_output = 16-z-only_input-both
                w, v = only_input+both, only_output+both
                count = (math.comb(16, z)*math.comb(16-z, only_input)
                         *math.comb(16-z-only_input, both))
                term = count*xp[w]*yp[v]*gp[only_input]*gamma[z]
                matrix[1][1] += term
                if w and v and both and (v != 1 or w == 16):
                    matrix[1][0] += term*gp[min(both, 4)]
    return matrix


def certify(k_log, first=1100, last=2400, grid=10000):
    if not 0 < first < last <= grid:
        raise ValueError('invalid density range')
    rows = ((2**k_log+22*64-1)//(22*64))*64
    n, threshold = 32*rows, (32*rows+9)//10-1
    prefactors = {density: (interval(rows+1)**32
        * terminal_power(transfer(witness), n//16)
        / (interval(witness.tau)*interval(witness.y)**threshold))
        for density, witness in WITNESSES}
    total, next_weight, cells = ZERO, (first*n+grid-1)//grid, 0
    for index in range(first, last):
        lo = (index*n+grid-1)//grid
        hi = ((index+1)*n+grid-1)//grid-1
        if lo > hi:
            continue
        assert lo == next_weight
        next_weight = hi+1
        lower, upper = Fraction(lo, n), Fraction(hi, n)
        density, witness = min(WITNESSES, key=lambda item: abs(item[0]-(lower+upper)/2))
        base = witness.x*witness.tilt
        endpoint = lo if base >= 1 else hi
        bound = ((hi-lo+1)*prefactors[density]*row_upper(lower, upper, witness)**rows
                 / interval(base)**endpoint)
        if not 0 < bound.hi < 1:
            raise ArithmeticError(f'failed density cell {index}/{grid}: {bound.hi}')
        total += bound
        cells += 1
    assert next_weight == (last*n+grid-1)//grid
    margin = integer_margin(total.hi)
    if margin < 40:
        raise ArithmeticError(f'class margin below40: {margin}')
    return dict(construction='cauchy4-nonzero-diagonal-state-refresh-every16',
                outer=[32, 22], rate='11/16', k_log=k_log, k=22*rows, n=n, rows=rows,
                minimum_input_weight=(first*n+grid-1)//grid,
                maximum_input_weight=next_weight-1, output_threshold=threshold,
                density_min=str(Fraction(first, grid)), density_max_exclusive=str(Fraction(last, grid)),
                checked_cells=cells, failure_upper=str(total.hi),
                class_margin_bits_lower=margin, rigorous=True, whole_code=False,
                scope='complete inverse-input-weight class only; inputs below .11N remain open')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--k-log', nargs='+', type=int, default=[16, 18, 20])
    parser.add_argument('--first', type=int, default=1100)
    args = parser.parse_args()
    for k_log in args.k_log:
        print(json.dumps(certify(k_log, first=args.first)), flush=True)
