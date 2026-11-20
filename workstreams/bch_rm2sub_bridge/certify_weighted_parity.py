"""Exact degree/coefficient budgets for weighted all-region pair parity.

Polynomial variables are row averages of [0,1]-valued row-local functions.
Degree d marks at most d rows. No independence between regions is assumed.
"""
from fractions import Fraction as F
from pathlib import Path
import bridge as base
from certify_tail_parity_mixing import kraw_values
import math


def floor_negative_log2(value):
    assert 0 < value < 1
    bits = value.denominator.bit_length()-value.numerator.bit_length()
    while value > F(1, 1 << bits):
        bits -= 1
    while value <= F(1, 1 << (bits+1)):
        bits += 1
    return bits


def error_budget(q, degree, single, paired):
    assert 0 <= degree <= q and 0 < single < 1 and single*single <= paired < 1
    characters = (1 << 256)-2
    return characters*single**(q-degree)+F(characters**2, 4)*paired**(q-degree)


def run():
    path = base.HERE/'generated/weight80_family.json'
    family = base.read(path)
    for name, digest in family['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    single = max(F(abs(v), math.comb(256, j))
                 for j, v in enumerate(kraw_values(256, 80)) if 0 < j < 256)
    assert single == F(3, 8)
    far_tables = [kraw_values(256, w) for w in range(82, 161, 2)]
    far = max(F(abs(values[j]), math.comb(256, j))
              for values in far_tables for j in range(1, 256))
    assert far == F(23, 64)
    close = base.decode(family['close_pair_probability_upper'])
    paired = close+(1-close)*far
    rows = []
    for degree in [0, 32, 64, 128, 256, 512, 1024, 1536]:
        error = error_budget(2620, degree, single, paired)
        bits = floor_negative_log2(error)
        rows.append(dict(degree=degree, coefficient_l1_scaled_error_upper=base.encode(error),
                         unweighted_budget_bits=bits))
        print('degree', degree, 'coefficient-L1 budget error <=2^-'+str(bits), flush=True)
    target = error_budget(2620, 512, single, paired)*(1 << 2048)
    assert target < F(1, 1 << 500)
    result = dict(status='EXACT_POLYNOMIAL_WEIGHTED_PAIR_PARITY', occupation=2620,
                  row_family='actual fixed BCH weight80 shell',
                  single_character_bound=base.encode(single),
                  paired_character_bound=base.encode(paired), rows=rows,
                  uniform_over_occupied_supports=True,
                  polynomial_variables='averages of [0,1]-valued functions of individual row-pair variables',
                  degree512_coefficient_to_mean_ratio_cap=base.encode(F(1 << 2048)),
                  degree512_relative_error_upper=base.encode(target),
                  full_second_moment_certified=False,
                  local_sha256={str(p.relative_to(base.HERE)):base.sha(p)
                                for p in [Path(__file__), path, base.HERE/'certify_tail_parity_mixing.py']})
    output = base.HERE/'generated/weight80_weighted_parity.json'
    if output.exists():
        assert result == base.read(output)
    else:
        base.write_new(output, result)
    print('Degree<=512 and coefficient-L1/mean <=2^2048: relative parity error <2^-500', flush=True)


if __name__ == '__main__':
    run()
