"""Exact scalar certificate for the setup-only temporal selector statistic.

Inject quartet_selector_capacity.py definitions with its main block removed.
Its directed rational logarithm routines are the only dependencies.  This
helper never constructs a129-state matrix or a floating-point witness.

For iid H-slot probability p, threshold r and beta>1, put
q0=Pr[Bin(4,p)<r], q1=E[1_{j>=r} beta^(8(j-r))].  Choose z>=1 and
lambda=max(q0+q1*z,q0/z+q1*beta^4).  Log-convexity proves

  (q0+q1*z*beta^c)^4 <=lambda^4*z^c, c=0,...,4.

The positive coefficient identity then gives K(z^a)<=lambda^128 z^a
for the129-state setup operator.  Initial population0 and free endpoint
imply E beta^(E_r+8Z_r)<=lambda^(2^19).

Regional quota conditioning and the H-set union give the exact tail
bound C(8192,h)*pi_h(p)^(-256)*lambda^(2^19)*beta^(-score).
The score required by an actual profile is
Ls-512-f+8*((4-r)*(F-f)/4-R+128h).
Thresholds describe entire setup-geometry classes, not a pointwise
profile union.  This remains a partial-class result, not a full code
certificate.  The parent alone runs this source, serially and file-free.
"""

import argparse
from fractions import Fraction
from math import comb
from operator import index
import os
import shutil
import threading
import time


TSC_M = 8192
TSC_N = 1 << 21
TSC_J = TSC_N // 4
TSC_REGIONS = 256
TSC_WIDTH = 512
TSC_SCALE = 1 << 48
TSC_P = Fraction(5, 16)
TSC_BETA = Fraction(12, 11)
TSC_Z = Fraction(77, 64)
TSC_LAMBDA = Fraction(10613108988597669204677, 8760324352986371325952)
TSC_PROFILES = ((3584, 819200, 845414, 4, 1048576),
                (3072, 819200, 837222, 4, 1015808))


def _tsc_ceildiv(a, b):
    if b <= 0:
        raise ValueError('Positive divisor required')
    return -((-a) // b)


def _tsc_dependencies():
    for name in ('selector_log2_interval', 'selector_binomial_log_interval'):
        if name not in globals():
            raise RuntimeError('Inject quartet_selector_capacity.py directed log helpers')
    if globals().get('SC_SCALE') != TSC_SCALE:
        raise RuntimeError('The injected directed-log scale must equal2^48')


def temporal_scalar_parameters(p=TSC_P, r=2, beta=TSC_BETA, z=TSC_Z):
    """Five exact local coefficient comparisons, with no approximate inputs."""
    p, beta, z, r = Fraction(p), Fraction(beta), Fraction(z), index(r)
    if not 0 < p < 1 or not 1 <= r <= 3 or beta <= 1 or z < 1:
        raise ValueError('Require0<p<1,r in1..3,beta>1,z>=1')
    probabilities = [comb(4, j) * p ** j * (1 - p) ** (4 - j) for j in range(5)]
    q0 = sum(probabilities[:r], Fraction(0))
    q1 = sum((probabilities[j] * beta ** (8 * (j - r)) for j in range(r, 5)), Fraction(0))
    lam = max(q0 + q1 * z, q0 / z + q1 * beta ** 4)
    for c in range(5):
        if (q0 + q1 * z * beta ** c) ** 4 > lam ** 4 * z ** c:
            raise ArithmeticError('Exact local supersolution comparison failed')
    return dict(p=p, r=r, beta=beta, z=z, q0=q0, q1=q1, lam=lam,
                exact_local_comparisons=5)


def temporal_scalar_threshold(h, p=TSC_P, r=2, beta=TSC_BETA, z=TSC_Z,
                              tail_bits=92):
    """Certified threshold for the entire H-unioned setup tail at fixedh,r."""
    _tsc_dependencies()
    h, tail_bits = index(h), index(tail_bits)
    if not 1 <= h < TSC_M or tail_bits < 1:
        raise ValueError('Require1<=h<8192 and a positive tail budget')
    parameters = temporal_scalar_parameters(p, r, beta, z)
    p, beta, lam = parameters['p'], parameters['beta'], parameters['lam']
    lambda_upper = selector_log2_interval(lam.numerator, lam.denominator)[1]
    beta_lower, beta_upper = selector_log2_interval(beta.numerator, beta.denominator)
    p_lower = selector_log2_interval(p.numerator, p.denominator)[0]
    complement = 1 - p
    complement_lower = selector_log2_interval(complement.numerator, complement.denominator)[0]
    binomial_lower = selector_binomial_log_interval(h)[0]
    if beta_lower <= 0:
        raise ArithmeticError('Positive beta-log lower bound required')
    # log C(M,h)-256log pi_h(p) = -255log C(M,h)
    #                                  -256h log p-256(M-h)log(1-p).
    overhead_upper = (-(TSC_REGIONS - 1) * binomial_lower
                      - TSC_REGIONS * h * p_lower
                      - TSC_REGIONS * (TSC_M - h) * complement_lower)
    numerator = TSC_J * lambda_upper + overhead_upper
    threshold = _tsc_ceildiv(numerator + TSC_SCALE * tail_bits, beta_lower)
    if numerator - threshold * beta_lower > -tail_bits * TSC_SCALE:
        raise ArithmeticError('Directed threshold verification failed')
    return dict(h=h, r=parameters['r'],
                p=(p.numerator, p.denominator),
                beta=(beta.numerator, beta.denominator),
                z=(parameters['z'].numerator, parameters['z'].denominator),
                lambda_fraction=(lam.numerator, lam.denominator),
                exact_local_comparisons=parameters['exact_local_comparisons'],
                setup_moment_log2_integer_upper=_tsc_ceildiv(TSC_J * lambda_upper, TSC_SCALE),
                threshold=threshold, tail_bits=tail_bits,
                log_lambda_upper=lambda_upper, log_beta_lower=beta_lower,
                log_beta_upper=beta_upper, log_binomial_lower=binomial_lower,
                log_p_lower=p_lower, log_complement_lower=complement_lower,
                tail_numerator_upper=numerator,
                directed_log_scale=TSC_SCALE,
                exact_integer_certificate=True, no_state_history_union=True,
                no_profile_or_death_union=True, initial_empty_population=True,
                free_endpoint=True, full_code_certificate=False)


def temporal_scalar_profile(h, R, Ls, factual, F, **markers):
    """Evaluate a profile's required score in an already global setup class."""
    h, R, Ls, factual, F = map(index, (h, R, Ls, factual, F))
    if not 0 <= factual <= Ls <= F <= TSC_N or R < 128 * h:
        raise ValueError('Invalid actual-state profile')
    G = F - factual
    if G % 4:
        raise ValueError('Actual G=F-f must be divisible by four')
    row = temporal_scalar_threshold(h, **markers)
    selector = (4 - row['r']) * (G // 4) - R + 128 * h
    score = Ls - TSC_WIDTH - factual + 8 * selector
    beta_log = row['log_beta_lower'] if score >= 0 else row['log_beta_upper']
    exponent = _tsc_ceildiv(row['tail_numerator_upper'] - score * beta_log, TSC_SCALE)
    return dict(**row, R=R, Ls=Ls, factual=factual, F=F, G=G,
                selector_score=selector, temporal_score_required=score,
                exact_log2_tail_upper=exponent,
                qualifies_temporal_class=score >= row['threshold'])


def _tsc_polynomial_power(local, exponent):
    coefficients = [Fraction(1)]
    for _ in range(exponent):
        result = [Fraction(0)] * (len(coefficients) + len(local) - 1)
        for i, first in enumerate(coefficients):
            for j, second in enumerate(local):
                result[i + j] += first * second
        coefficients = result
    return coefficients


def temporal_scalar_checks():
    """Exact independent small-row and coefficient tests; no productionmatrix."""
    local_checks = row_checks = 0
    for p in (Fraction(1, 2), TSC_P):
        for r in (1, 2, 3):
            for beta in (Fraction(3, 2), TSC_BETA, Fraction(35, 32)):
                parameters = temporal_scalar_parameters(p, r, beta, TSC_Z)
                q0, q1, lam = (parameters[key] for key in ('q0', 'q1', 'lam'))
                local_checks += 5
                local = [comb(4, c) * (q0 + q1 * TSC_Z * beta ** c) for c in range(5)]
                for width in (4, 8):
                    m = width // 4
                    coefficients = _tsc_polynomial_power(local, m)
                    for a in range(m + 1):
                        u = 4 * a
                        direct = Fraction(0)
                        for b in range(m + 1):
                            v = 4 * b
                            overlap = sum(Fraction(comb(u, ell) * comb(width - u, v - ell), comb(width, v))
                                          * beta ** ell
                                          for ell in range(max(0, v - width + u), min(u, v) + 1))
                            weight = comb(m, b) * q0 ** (m - b) * q1 ** b
                            direct += weight * overlap * TSC_Z ** b
                        assert direct == coefficients[u] / comb(width, u)
                        assert direct <= lam ** m * TSC_Z ** a
                        row_checks += 1
    parameters = temporal_scalar_parameters()
    assert parameters['lam'] == TSC_LAMBDA
    assert TSC_Z >= 1 and TSC_J == 128 * 4096
    print(dict(exact_quartet_comparisons=local_checks,
               exact_small_scalar_rows=row_checks,
               fixed_lambda=(TSC_LAMBDA.numerator, TSC_LAMBDA.denominator),
               production_matrix_constructed=False, full_code_certificate=False), flush=True)


def temporal_scalar_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--certify', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, action='append',
                        metavar=('H', 'R', 'LS', 'F_BIRTH', 'F'))
    parser.add_argument('--p', type=Fraction, default=TSC_P)
    parser.add_argument('--r', type=int, default=2)
    parser.add_argument('--beta', type=Fraction, default=TSC_BETA)
    parser.add_argument('--z', type=Fraction, default=TSC_Z)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            temporal_scalar_checks()
        else:
            for profile in args.profile or TSC_PROFILES:
                print(temporal_scalar_profile(*profile, p=args.p, r=args.r,
                                              beta=args.beta, z=args.z), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    temporal_scalar_main()
