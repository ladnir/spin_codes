"""Exact factorial-moment refinements of the frozen V50 return cap.

Fix a potential support S, |S|=j, a good global gamma, X uniform on F^S,
R=CX, and Y=(I+gamma AC)X. Let K count zero coordinates of Y. All measures
below include the indicator R!=0; they are not conditioned probabilities.
For k>=0, summing V50's joint-zero bound over every k-set gives

    E[binom(K,k) 1_{R!=0}] <= H[k]
      = sum_{c+e=k} binom(j,c) binom(64-j,e) p[j,c,e].

H[0]=1-q^(-min(j,32)) is exact. If the actual zero set has size k, its
own joint-zero event has positive probability, so H[k]>0. Thus the last
positive H index M bounds K on R!=0. This uses exact rational zeros,
not a numerical cutoff or a claim that a positive bound is attainable.

For 0<=a<=b<=M, binom(K,b)/binom(M,b) <=
binom(K,a)/binom(M,a), including K<b. This clips each factorial moment.
Also Pr[R!=0,K>=k] <= min_{0<=l<=k} H[l]/binom(k,l). Finally,

    z^(64-K) = z^64 + (1-z) sum_{k=1}^K z^(64-k).

These two inequalities give a clipped positive Bernstein cap and a tail
cap. Their minimum with V50 is valid pointwise under the SAME global
principal event. No new gamma exclusion, independent-zero assumption,
fractional power, coefficient calculation, or production run is introduced.
Source-only consumers may import this module; it never writes artifacts.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb
from pathlib import Path

import hd_field128_translated_packet_v50 as provider

HERE = Path(__file__).resolve().parent
NAME = Path(__file__).name
PROVIDER_SHA256 = '8ac2a1346bf4ce75f668560f7f00192b91ece98d0cd8252e1e96232b8aefe985'
FIELD_ORDER, WIDTH, STATE = provider.FIELD_ORDER, provider.WIDTH, provider.STATE


def verify_frozen():
    path = Path(provider.__file__).resolve()
    if path.parent != HERE or sha256(path.read_bytes()).hexdigest() != PROVIDER_SHA256:
        raise RuntimeError('Frozen V50 provider changed')
    provider.verify_frozen()


def source_pins():
    verify_frozen()
    pins = provider.source_pins()
    if pins.get(Path(provider.__file__).name) != PROVIDER_SHA256:
        raise RuntimeError('V50 closure identity mismatch')
    digest = sha256(Path(__file__).read_bytes()).hexdigest()
    if NAME in pins and pins[NAME] != digest:
        raise RuntimeError('Conflicting source identity')
    pins[NAME] = digest
    verify_frozen()
    return dict(sorted(pins.items()))


def check_source_pins(expected):
    if source_pins() != expected:
        raise RuntimeError('Tail-moment source closure changed')


def _exact_nonnegative(value):
    if isinstance(value, (bool, float)):
        raise ValueError('Exact nonnegative rational moments are required')
    value = F(value)
    if value < 0:
        raise ValueError('Negative factorial-moment bound')
    return value


def factorial_envelope(bounds):
    """Pure generic finite-measure helper; the caller proves the input bounds.

    bounds must contain orders 0..n for a variable K in 0..n. A zero upper
    bound at the last positive index+1 forces the stated support cutoff.
    Exact all-zero measure is supported. No construction is inferred here.
    """
    bounds = tuple(_exact_nonnegative(x) for x in bounds)
    if not 1 <= len(bounds) <= 65:
        raise ValueError('Expected all factorial orders for width0..64')
    maximum = max((k for k, value in enumerate(bounds) if value), default=0)
    clipped = [F(0)]*len(bounds)
    running = bounds[0]
    for k in range(maximum+1):
        choose = comb(maximum, k)
        running = min(running, bounds[k]/choose)
        clipped[k] = choose*running
    # Every entry is a separate proved upper bound; selecting its minimum
    # is not a probabilistic independence or simultaneous-attainment claim.
    tails = [bounds[0]]
    for k in range(1, maximum+1):
        tails.append(min(clipped[l]/comb(k, l) for l in range(k+1)))
    return maximum, tuple(clipped), tuple(tails)


def evaluate_envelope(bounds, z):
    """Exact cap diagnostics for the generic finite measure above.

    z is an exact rational in (0,1). The measure need not have mass one.
    The original positive expansion is a control, not the chosen result.
    """
    z = provider._mark(z)
    bounds = tuple(_exact_nonnegative(x) for x in bounds)
    maximum, clipped, tails = factorial_envelope(bounds)
    width, w = len(bounds)-1, 1-z
    zp, wp = [F(1)], [F(1)]
    for unused in range(width):
        zp.append(zp[-1]*z)
        wp.append(wp[-1]*w)
    original = sum((h*zp[width-k]*wp[k] for k, h in enumerate(bounds)), F(0))
    bernstein = sum((h*zp[width-k]*wp[k] for k, h in enumerate(clipped)), F(0))
    tail = zp[width]*bounds[0]+w*sum(
        (zp[width-k]*tails[k] for k in range(1, maximum+1)), F(0))
    if bernstein > original or min(original, bernstein, tail) < 0:
        raise ArithmeticError('Exact factorial-envelope invariant failed')
    return dict(maximum_zeros=maximum, original_bernstein=original,
                clipped_bernstein=bernstein, tail=tail,
                refined=min(original, bernstein, tail))


@lru_cache(maxsize=65)
def _factorial_bounds(j):
    # V50 already groups its positive joint-zero probabilities by c+e.
    # Its omitted c=j and e>=32 terms are EXACT zeros. Reuse those exact
    # coefficients, not FLOAT logs or a newly reconstructed rank formula.
    bounds = provider._return_coefficients(j)
    expected = F(0) if j == 0 else 1-F(1, FIELD_ORDER**min(j, STATE))
    if len(bounds) != WIDTH+1 or bounds[0] != expected:
        raise ArithmeticError('Original q**j normalization / total mass changed')
    return bounds


def factorial_bounds(j):
    provider._geometry(j)
    verify_frozen()
    try:
        return _factorial_bounds(j)
    finally:
        verify_frozen()


def cap_components(j, z):
    """Return exact caps before division by q**32-1, including old basic caps."""
    provider._geometry(j)
    z = provider._mark(z)
    verify_frozen()
    try:
        old = provider.exact_moments(j, z)
        result = evaluate_envelope(_factorial_bounds(j), z)
        result.update(v50_return=old[2], selected=min(old[2], result['refined']))
        return result
    finally:
        verify_frozen()


def exact_moments(j, z):
    """Same interface as V50; quiet, birth and live are unchanged."""
    provider._geometry(j)
    z = provider._mark(z)
    verify_frozen()
    try:
        quiet, birth, ret, live = provider.exact_moments(j, z)
        cap = evaluate_envelope(_factorial_bounds(j), z)['refined']
        return quiet, birth, min(ret, cap), live
    finally:
        verify_frozen()


def return_cap(j, z):
    """The minimum of V50's full return cap and both moment refinements."""
    return exact_moments(j, z)[2]


def exact_packet(j, z):
    z = provider._mark(z)
    quiet, birth, ret, live = exact_moments(j, z)
    d = FIELD_ORDER**STATE-1
    density = F(d, (FIELD_ORDER-1)**STATE)
    return (quiet, density*birth, ret/d,
            min(density*live, z**(WIDTH-STATE+1)) if j == 0 else density*live)


def principal_bad_gamma_upper():
    return provider.principal_bad_gamma_upper()


def metadata():
    meta = provider.metadata()
    meta.update(schema='translated-state32-return-tail-moments-v57',
        source_sha256=source_pins(), mathematical_provider_sha256=PROVIDER_SHA256,
        return_refinement='min(V50, clipped factorial Bernstein, factorial tail)',
        exact_nonzero_syndrome_mass=True, forced_zero_cutoff_from_exact_bounds=True,
        factorial_moments_not_independent=True, construction_changed=False,
        numerical_improvement_measured=False)
    return meta
