"""Sharper central witness via a degree-six log-kernel polynomial and exact variance.

Numerical fitting only chooses rational coefficients. Every integer point
is checked with Arb. Exact additive moments use only one/two-column laws.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from numpy.polynomial import Chebyshev, Polynomial
from flint import arb, arb_poly, ctx
import bridge as base
import tightened_occupancy as tight
import fixed_slice_additive_moments as exact


def run(verify=False):
    output = base.HERE/'generated/zero_state_polynomial_witness.json'
    old = base.read(output) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    n, q, degree = 256, 2620, 6
    mu = F(q*5, 16)
    source = base.HERE/'generated/zero_state_central_witness.json'
    central = base.read(source)
    for name, digest in central['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    family_path = base.HERE/'generated/weight80_family.json'
    family = base.read(family_path)
    lo, hi = central['even_region_count_window']
    assert [lo, hi] == [740, 900]
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    logs = {j:(kernel[j]/math.comb(8192, j)).log() for j in range(lo, hi+1, 2)}
    if old:
        coefficients = [base.decode(v) for v in old['proposal_coefficients']]
    else:
        xs = np.array([float(F(j)-mu) for j in logs])
        fit = Chebyshev.fit(xs, [float(value) for value in logs.values()], degree).convert(kind=Polynomial)
        coefficients = [F.from_float(float(c)) for c in fit.coef]
    assert len(coefficients) == degree+1
    ball = lambda value:arb(value.numerator)/value.denominator
    from audit_bch_q1_full_arb import rational
    residual = min((value-sum(ball(c*(F(j)-mu)**r) for r, c in enumerate(coefficients))).lower()
                   for j, value in logs.items())
    intercept = (ball(coefficients[0])+residual).lower()
    random_coefficients = [F(0), F(0)]+coefficients[2:]
    mean, variance = exact.additive_mean_variance(n, 80, q, random_coefficients)
    assert mean > 0 and variance > 0
    norm = n*sum(abs(c)*(q+mu)**r for r, c in enumerate(random_coefficients))
    characters = (1 << 256)-2
    eps0 = characters*F(3, 8)**q
    epsd = characters*F(3, 8)**(q-degree)
    eps2d = characters*F(3, 8)**(q-2*degree)
    mean_parity = ((ball(mean-norm*epsd))/(1+ball(eps0))).lower()
    assert mean_parity > 0
    second_parity = ball(mean**2+variance+norm**2*eps2d).upper()
    variance_parity = (second_parity-mean_parity**2).upper()
    gamma = ball(base.decode(central['good_given_parity_probability_lower']))
    conditional_mean = (mean_parity-(variance_parity*(1-gamma)/gamma).sqrt()).lower()
    probability = rational((arb(2)**-255*gamma*(n*intercept+conditional_mean).exp()).lower())
    lower = math.comb(8192, q)*family['shell_size_lower']**q*probability
    previous = base.decode(central['first_moment_lower'])
    assert lower > previous
    if old:
        assert probability >= base.decode(old['restricted_zero_state_probability_lower'])
        assert lower >= base.decode(old['first_moment_lower'])
        print('512-bit degree-six central witness replay passed', flush=True)
        return
    local = [Path(__file__), Path(exact.__file__), source, family_path, Path(tight.__file__),
             base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']
    base.write_new(output, dict(status='ACTUAL_POLYNOMIAL_CENTRAL_COUNT_WITNESS_LOWER', configuration=tight.NAME,
        occupation=q, row_weight=80, even_region_count_window=[lo, hi], polynomial_degree=degree,
        proposal_coefficients=[base.encode(c) for c in coefficients],
        certified_constant_term_lower=base.encode(rational(intercept)),
        additive_residual_exact_mean=base.encode(mean), additive_residual_exact_variance=base.encode(variance),
        parity_polynomial_coefficient_norm_upper=base.encode(norm),
        conditional_residual_mean_lower=base.encode(rational(conditional_mean)),
        restricted_zero_state_probability_lower=base.encode(probability), first_moment_lower=base.encode(lower),
        same_encoder_and_setup_distribution=True, full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local}))
    log2 = lambda value:math.log2(value.numerator)-math.log2(value.denominator)
    print('Degree-six central witness log2 mean lower', log2(lower), 'gain', log2(lower/previous), 'bits;',
          'residual mean', float(mean), 'variance', float(variance), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    run(parser.parse_args().verify)
