"""Kernel-aware transfer and adaptive counting for arbitrary occupations.

All calculations here are nearest binary64 diagnostics. The formulas apply
through Q=L; callers must record the actual evaluated range separately.
"""
import math

import numpy as np

from activation_q1 import matrix_product


def epoch_logs(t, s, spectrum, kernel, lam, maximum):
    if (not 0 <= maximum <= t or not math.isfinite(lam) or lam <= 0
            or sum(spectrum.values()) != (1 << s)-1 or s < 2
            or any(not 0 < w <= t or n <= 0 for w,n in spectrum.items())
            or len(kernel) != t+1 or kernel[0] != 1
            or sum(kernel) != 1 << (t-s)
            or any(not 0 <= n <= math.comb(t,j) for j,n in enumerate(kernel))):
        raise ValueError('invalid map spectra or epoch parameters')
    den = (1 << s)-1
    log_den = math.log(den)
    log_kappa = math.log1p(1/(den-1))
    output = np.full((maximum+1,3,3), -np.inf)
    for j in range(maximum+1):
        total = math.comb(t,j)
        log_total = math.log(total)
        arbitrary = -math.inf
        uniform = -math.inf
        for w,n in spectrum.items():
            terms = [math.log(math.comb(w,v)*math.comb(t-w,j-v))-log_total-lam*(w+j-2*v)
                     for v in range(max(0,j-t+w), min(w,j)+1)]
            moment = float(np.logaddexp.reduce(terms))
            arbitrary = max(arbitrary, moment)
            uniform = float(np.logaddexp(uniform, math.log(n)-log_den+moment))
        # Use integer subtraction before taking logs; beta can be very near 1.
        nonkernel = math.log(total-kernel[j])-log_total if kernel[j] < total else -math.inf
        beta = math.log(kernel[j])-log_total if kernel[j] else -math.inf
        output[j,0,0] = beta-lam*j
        output[j,0,1] = nonkernel-lam*j
        output[j,1,0] = min(nonkernel,arbitrary)-log_den
        output[j,1,2] = arbitrary
        output[j,2,0] = log_kappa+min(nonkernel,uniform)-log_den
        output[j,2,2] = log_kappa+uniform
    return output


def polynomial_product(left, right, maximum):
    """Truncated log-domain positive matrix-polynomial convolution."""
    limit = min(maximum, len(left)+len(right)-2)
    result = np.full((limit+1,3,3), -np.inf)
    for degree in range(limit+1):
        lo, hi = max(0,degree-len(right)+1), min(degree,len(left)-1)
        products = matrix_product(left[lo:hi+1], right[degree-hi:degree-lo+1][::-1])
        result[degree] = np.logaddexp.reduce(products,axis=0)
    return result


def region_logs(epoch, t, length, maximum):
    if length < t or length % t or not 0 <= maximum <= length or len(epoch) != min(t,maximum)+1:
        raise ValueError('non-native or incomplete region input')
    power = epoch + np.array([math.log(math.comb(t,j)) for j in range(len(epoch))])[:,None,None]
    current = np.full((1,3,3), -np.inf)
    current[0,0,0] = current[0,1,1] = current[0,2,2] = 0.
    exponent = length//t
    while exponent:
        if exponent & 1:
            current = polynomial_product(current,power,maximum)
        exponent >>= 1
        if exponent:
            power = polynomial_product(power,power,maximum)
    current -= np.array([math.log(math.comb(length,j)) for j in range(maximum+1)])[:,None,None]
    return current


def terminal_log(matrix, length):
    """Log of e_zero M^length 1, retaining state across regions."""
    if length < 1:
        raise ValueError('positive constituent length required')
    current = np.full((1,3,3), -np.inf)
    current[0,0,0] = current[0,1,1] = current[0,2,2] = 0.
    power = matrix[None,:,:]
    while length:
        if length & 1:
            current = matrix_product(current,power)
        length >>= 1
        if length:
            power = matrix_product(power,power)
    return float(np.logaddexp.reduce(current[0,0]))


def bands_for(counts, block, band_count=8):
    """Deterministic contiguous rank bands, retaining singleton endpoints."""
    if band_count < 1 or not counts or any(not 0 < w <= block or n <= 0 for w,n in counts.items()):
        raise ValueError('nonzero integer spectrum/caps required')
    weights = sorted(w for w in counts if w != block)
    bands = [list(map(int, group)) for group in np.array_split(weights,min(band_count,len(weights))) if len(group)] if weights else []
    if block in counts:
        bands.append([block])
    return bands


def density_roots(counts, block, bands, shift=0.):
    """Return log(rho_g), log(p_g), log(1-p_g) for fixed counting measures."""
    if sorted(w for band in bands for w in band) != sorted(counts):
        raise ValueError('bands must partition the supported spectrum')
    roots, active, inactive = [], [], []
    for band in bands:
        if band == [block]:
            # An exact all-one support is already a Bernoulli(1) atom.
            roots.append(math.log(counts[block])/block)
            active.append(0.)
            inactive.append(-math.inf)
            continue
        midpoint = (min(band)+max(band))/(2*block)
        logit = math.log(midpoint)-math.log1p(-midpoint)+shift
        log_p = -float(np.logaddexp(0.,-logit))
        log_not_p = -float(np.logaddexp(0.,logit))
        cost = max(math.log(counts[w])-math.log(math.comb(block,w))-w*log_p-(block-w)*log_not_p for w in band)
        roots.append(cost/block)
        active.append(log_p)
        inactive.append(log_not_p)
    return np.array(roots), np.array(active), np.array(inactive)


def adaptive_logs(regions, roots, active, inactive):
    """Yield a bound for every integer Q in the supplied coefficient range."""
    current = regions.copy()
    for occupation in range(1,len(regions)):
        updated = np.full((len(current)-1,3,3), -np.inf)
        for rho,p,not_p in zip(roots,active,inactive):
            np.maximum(updated,rho+np.logaddexp(not_p+current[:-1],p+current[1:]),out=updated)
        current = updated
        yield occupation,current[0]


def occupation_bounds(regions, counts, block, length, cutoff, lam, shift=0., band_count=8):
    """Return log first-moment bounds for Q=1..len(regions)-1.

    Counts are deterministic exact multiplicities or simultaneous caps for
    one reused constituent. Expected spectra are deliberately not accepted
    as a separate model here; use random_spectrum_caps for that ensemble.
    """
    if any(not isinstance(n,int) or n <= 0 for n in counts.values()):
        raise ValueError('deterministic integer multiplicities/caps required')
    bands = bands_for(counts,block,band_count)
    roots,active,inactive = density_roots(counts,block,bands,shift)
    result = []
    for q,matrix in adaptive_logs(regions,roots,active,inactive):
        bound = math.log(math.comb(length,q))+q*math.log(len(bands))+cutoff*lam+terminal_log(matrix,block)
        # Counting all messages is also a valid first-moment bound.
        trivial = math.log(math.comb(length,q))+q*math.log(sum(counts.values()))
        result.append(min(bound,trivial))
    return np.array(result)


def random_spectrum_caps(block, dimension, failure_bits=60):
    """Simultaneous Markov caps for one uniform full-rank binary subspace.

    With probability >=1-2^-failure_bits every shell obeys these integer
    caps. The setup-failure charge belongs once in the final union, never
    once per occupation. Zero caps exclude a shell on the good setup event.
    """
    if not 0 < dimension < block or failure_bits < 0:
        raise ValueError('invalid random-subspace parameters')
    mass = (1 << dimension)-1
    denominator = (1 << block)-1
    multiplier = block*(1 << failure_bits)
    return {w: cap for w in range(1,block+1)
            if (cap := min(mass,math.comb(block,w),
                           math.comb(block,w)*mass*multiplier//denominator)) > 0}
