"""A central-count bad-message witness with an outward first-moment lower.

Restrict the counted messages, not the encoder: every region count lies in
740..900. Polynomial parity mixing, exact count-energy variance, and a
quadratic envelope retain most of the Jensen denominator while removing
all outside-window region weights from this witness's second moment.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from flint import arb, arb_poly, ctx
import bridge as base
import tightened_occupancy as tight


def run(verify=False):
    path = base.HERE/'generated/zero_state_central_witness.json'
    old = base.read(path) if verify else None
    if old:
        for name, digest in old['local_sha256'].items():
            assert base.sha(base.HERE/name) == digest
    else:
        assert not path.exists()
    ctx.prec = 512 if verify else 256
    n, q, lo, hi = 256, 2620, 740, 900
    p, mu = F(5, 16), F(2620*5, 16)
    v = p*(1-p)
    source = base.HERE/'generated/zero_state_jensen_lower.json'
    jensen = base.read(source)
    for name, digest in jensen['local_sha256'].items():
        assert base.sha(base.HERE/name) == digest
    family_path = base.HERE/'generated/weight80_family.json'
    family = base.read(family_path)
    from audit_bch_q1_full_arb import rational
    ball = lambda value: arb(value.numerator)/value.denominator
    characters = (1 << 256)-2
    eps = {degree:characters*F(3, 8)**(q-degree) for degree in [0, 2, 4]}
    mean_energy = q*n*v
    variance_energy = 2*q*(q-1)*F(n*n, n-1)*v*v
    norm = n*(q+mu)**2
    mean_lower = ((ball(mean_energy)-ball(norm*eps[2]))/(1+ball(eps[0]))).lower()
    second_upper = ball(mean_energy**2+variance_energy+norm**2*eps[4]).upper()
    variance_upper = (second_upper-mean_lower**2).upper()
    mass = (arb(11)/16)**q
    outside = arb(0)
    for j in range(q+1):
        if j % 2 == 0 and not lo <= j <= hi:
            outside += mass
        if j < q:
            mass *= arb(q-j)*5/((j+1)*11)
    conditional_error = base.decode(jensen['conditional_parity_error_upper'])
    good_relative = (1-512*(1+ball(conditional_error))*outside).lower()
    good_given_parity = (good_relative/(1+ball(eps[0]))).lower()
    assert good_given_parity > arb(4)/5
    energy_good_lower = (mean_lower-(variance_upper*(1-good_given_parity)/good_given_parity).sqrt()).lower()
    assert energy_good_lower > 0
    spectrum = tight.kernel_spectrum()
    kernel = arb_poly([spectrum.get(j, 0) for j in range(129)])**64
    logs = {j:(kernel[j]/math.comb(8192, j)).log() for j in range(lo, hi+1, 2)}
    if old:
        slope, curvature = base.decode(old['linear_slope']), base.decode(old['quadratic_coefficient'])
    else:
        deltas = np.array([float(F(j)-mu) for j in logs])
        solution = linprog([-n, 0, -float(energy_good_lower)],
            A_ub=np.column_stack([np.ones(len(deltas)), deltas, deltas*deltas]),
            b_ub=[float(value) for value in logs.values()],
            bounds=[(None, None), (None, None), (0, None)], method='highs')
        assert solution.success
        slope, curvature = [F.from_float(float(x)) for x in solution.x[1:]]
    assert curvature > 0
    b, c = ball(slope), ball(curvature)
    intercept = min((value-b*ball(F(j)-mu)-c*ball((F(j)-mu)**2)).lower() for j, value in logs.items())
    probability = (arb(2)**-255*good_relative*(n*intercept+c*energy_good_lower).exp()).lower()
    lower = rational(probability)
    mean = math.comb(8192, q)*family['shell_size_lower']**q*lower
    assert mean > F(1 << 19000)
    if old:
        assert lower >= base.decode(old['restricted_zero_state_probability_lower'])
        assert mean >= base.decode(old['first_moment_lower'])
        print('512-bit central-count witness replay passed', flush=True)
        return
    local = [Path(__file__), source, family_path, Path(tight.__file__),
             base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']
    base.write_new(path, dict(status='ACTUAL_CENTRAL_COUNT_BAD_MESSAGE_WITNESS', configuration=tight.NAME,
        occupation=q, row_weight=80, even_region_count_window=[lo, hi],
        same_encoder_and_setup_distribution=True, restricted_message_count=True,
        unconditional_energy_mean=base.encode(mean_energy), unconditional_energy_variance=base.encode(variance_energy),
        good_given_parity_probability_lower=base.encode(rational(good_given_parity)),
        conditional_energy_lower=base.encode(rational(energy_good_lower)),
        linear_slope=base.encode(slope), quadratic_coefficient=base.encode(curvature),
        quadratic_intercept_lower=base.encode(rational(intercept)),
        restricted_zero_state_probability_lower=base.encode(lower), first_moment_lower=base.encode(mean),
        full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local}))
    log2 = lambda value:math.log2(value.numerator)-math.log2(value.denominator)
    print('Central witness: good-given-parity >=', float(good_given_parity),
          'mean lower log2', log2(mean), 'loss vs unrestricted Jensen bits',
          log2(base.decode(jensen['first_moment_lower'])/mean),
          'curvature', float(curvature), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    run(parser.parse_args().verify)
