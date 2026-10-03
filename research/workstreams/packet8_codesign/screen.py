"""Bounded floating original-order screens for eight-bit routed packets.

Exact finite formulas define the local moments and hypergeometric placement.
All numerical evaluations are floating proposals, never certified endpoints.
Eight finite birth families retain their full 65535-state shapes before the
next fresh transitive update. State persists through every step and region.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, isfinite, log
from pathlib import Path
import sys
from time import monotonic

import numpy as np
from scipy.special import gammaln, logsumexp

import maps

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'k16_design'))
import packet_regional_log as logarithmic
import rs_outer
import rs_uniform_envelope


def walsh(values):
    result = np.asarray(values, dtype=np.longdouble).copy()
    length = 1
    while length < len(result):
        shape = (-1, 2*length, *result.shape[1:])
        blocks = result.reshape(shape)
        left, right = blocks[:, :length].copy(), blocks[:, length:].copy()
        blocks[:, :length], blocks[:, length:] = left+right, left-right
        length *= 2
    return result


def profile_polynomials(profiles, inactive, active):
    W = int(profiles[0].sum())
    if np.any(profiles.sum(axis=1) != W):
        raise ValueError('one complete packet profile per state required')
    result = []
    for profile in profiles:
        polynomial = np.array([1.], dtype=np.longdouble)
        for r, count in enumerate(profile):
            for _ in range(int(count)):
                polynomial = np.convolve(polynomial, [inactive[r], active[r]])
        result.append(polynomial / np.array([comb(W, j) for j in range(W+1)], dtype=np.longdouble))
    return np.asarray(result)


def moments(data, z):
    """Return W_j(s)=E[z^wt X 1{CX=s}] and M_j(a)=E[z^wt(X+Aa)]."""
    z = np.longdouble(z)
    if not 0 < z <= 1:
        raise ValueError('weight variable in (0,1] required')
    S, W, b = 1 << data['bits'], data['windows'], data['packet_bits']
    labels = (1 << b)-1
    factors = [((1+z)**(b-r)*(1-z)**r-1)/labels for r in range(b+1)]
    character = profile_polynomials(data['character_profiles'], [1.]*(b+1), factors)
    weighted = walsh(character[data['character_indices']])/S
    expected = np.array([(((1+z)**b-1)/labels)**j for j in range(W+1)])
    tolerance = 5e-12*expected
    if not np.isfinite(weighted).all() or np.any(weighted.min(axis=0) < -tolerance):
        raise ArithmeticError('unstable weighted syndrome inversion')
    negative_mass = np.maximum(-weighted, 0).sum(axis=0)
    weighted = np.maximum(weighted, 0)
    # Structural zero-occupancy law is exact. Single-packet laws are small
    # enough to recompute positively, avoiding cancellation outside support.
    weighted[:, 0] = 0
    weighted[0, 0] = 1
    direct = np.zeros(S, dtype=np.longdouble)
    for h in range(W):
        columns = data['columns'][b*h:b*(h+1)]
        for label in range(1, labels+1):
            direct[maps.apply(columns, label)] += z**label.bit_count()/(W*labels)
    direct_discrepancy = float(np.max(np.abs(direct-weighted[:, 1])))
    if direct_discrepancy > float(tolerance[1]):
        raise ArithmeticError('Walsh and positive single-packet census disagree')
    weighted[:, 1] = direct
    for j in data['zero_feedback_forbidden']:
        weighted[0, j] = 0
    if np.any(np.abs(weighted.sum(axis=0)-expected) > tolerance):
        raise ArithmeticError('weighted syndrome mass does not match input moment')
    expansion = profile_polynomials(data['expansion_profiles'],
        [z**r for r in range(b+1)], [((1+z)**b-z**r)/labels for r in range(b+1)])
    emission = expansion[data['expansion_indices']]
    if not np.isfinite(emission).all() or emission.min() < 0 or emission.max() > 1+1e-12:
        raise ArithmeticError('invalid finite emission polynomial')
    diagnostic = dict(walsh_clipped_negative_mass=list(map(float, negative_mass)),
        single_packet_walsh_discrepancy=direct_discrepancy,
        birth_mass_relative_error=list(map(float, np.abs(weighted.sum(axis=0)-expected)/expected)))
    return np.asarray(weighted, dtype=float), np.asarray(emission, dtype=float), diagnostic


def operators(weighted, emission):
    """Positive closure on delta0, uniform nonzero U, and normalized births.

    Each nonzero-source weighted output is dominated by T*U+(T/L)*delta0.
    The delta0 term is omitted at zero occupancy. Fresh M is independent of
    both the entering state and the current input. This is not an exact
    Markov kernel: excluded destinations are deliberately filled in.
    """
    if weighted.shape != emission.shape or weighted.ndim != 2:
        raise ValueError('matching state-by-occupancy moments required')
    S, count = weighted.shape
    W, L = count-1, S-1
    born = weighted.copy()
    born[0, :] = 0
    mass = born.sum(axis=0)
    if mass[0] != 0 or np.any(mass[1:] <= 0):
        raise ArithmeticError('positive nonzero birth families required')
    basis = born[:, 1:]/mass[None, 1:]
    means = basis.T @ emission
    uniform = emission[1:].mean(axis=0)
    family = np.zeros((W+1, W+2, W+2))
    for j in range(W+1):
        family[j, 0, 0] = weighted[0, j]
        if j:
            family[j, 0, j+1] = mass[j]
        family[j, 1, 1] = uniform[j]
        family[j, 2:, 1] = means[:, j]
        if j:
            family[j, 1:, 0] = family[j, 1:, 1]/L
    if not np.isfinite(family).all() or np.any(family < 0):
        raise ArithmeticError('invalid finite operator family')
    return family


def scaled_placement(local, degree, *, epochs, windows):
    """Ordered conditional placement, scaled separately at every occupancy.

    R_(e,j)=sum_k C(W,k)C((e-1)W,j-k)/C(eW,j) R_(e-1,j-k) T_k.
    This is exact without-replacement geometry, evaluated in floating point.
    Numerical underflow raises; the caller can retry entirely in log space.
    """
    local = np.asarray(local, dtype=float)
    if (local.ndim != 3 or local.shape[1] != local.shape[2]
            or len(local) != windows+1 or not 0 <= degree <= epochs*windows
            or not np.isfinite(local).all() or np.any(local < 0)):
        raise ValueError('complete positive family and feasible placement required')
    size = local.shape[1]
    scales = local.max(axis=(1, 2))
    if np.any(scales <= 0):
        raise ValueError('every occupancy requires a positive operator')
    with np.errstate(over='raise', invalid='raise', divide='raise', under='raise'):
        normalized = local/scales[:, None, None]
        logs = np.log(scales)
        current, current_logs = np.eye(size)[None, :, :], np.zeros(1)
        for epoch in range(1, epochs+1):
            limit, total, previous = min(degree, epoch*windows), epoch*windows, (epoch-1)*windows
            terms, base = [], np.full(limit+1, -np.inf)
            for k in range(min(windows, limit)+1):
                js = np.arange(k, min(limit, k+len(current)-1)+1)
                old = js-k
                weight = (log(comb(windows, k))+gammaln(previous+1)-gammaln(old+1)
                    -gammaln(previous-old+1)-gammaln(total+1)+gammaln(js+1)+gammaln(total-js+1))
                exponent = weight+current_logs[old]+logs[k]
                base[js] = np.maximum(base[js], exponent)
                terms.append((k, js, old, exponent))
            nxt = np.zeros((limit+1, size, size))
            for k, js, old, exponent in terms:
                nxt[js] += np.exp(exponent-base[js])[:, None, None]*(current[old] @ normalized[k])
            scale = nxt.max(axis=(1, 2))
            if not np.isfinite(scale).all() or np.any(scale <= 0):
                raise FloatingPointError('nonpositive scaled placement')
            current, current_logs = nxt/scale[:, None, None], base+np.log(scale)
    result = np.full(current.shape, -np.inf)
    positive = current > 0
    result[positive] = np.log(current[positive])
    return result+current_logs[:, None, None]


def placement(local, degree, *, epochs, windows, force_log=False):
    reason = None
    if not force_log and np.all(np.max(local, axis=(1, 2)) > 0):
        try:
            return scaled_placement(local, degree, epochs=epochs, windows=windows), 'scaled floating'
        except FloatingPointError as error:
            reason = str(error)
    value = logarithmic.log_placement(logarithmic.array_logs(local), degree, epochs=epochs, windows=windows)
    return value, 'logarithmic floating'+(f'; scaled fallback: {reason}' if reason else '')


def uniform_mixture(regional, q, *, packet_bits=8):
    js = np.arange(q+1)
    labels, alphabet = (1 << packet_bits)-1, 1 << packet_bits
    weights = gammaln(q+1)-gammaln(js+1)-gammaln(q-js+1)+js*log(labels)-q*log(alphabet)
    if abs(float(logsumexp(weights))) > 1e-8:
        raise ArithmeticError('binomial weights failed normalization')
    return logsumexp(regional[:q+1]+weights[:, None, None], axis=0)


def q1_support_logs(regional, *, regions, tilt, cutoff):
    """Average ordered supports using both regional boundary-state operators."""
    size = regional.shape[1]
    coefficients = np.full((regions+1, size), -np.inf)
    coefficients[0, 0] = 0
    for region in range(regions):
        following = np.full_like(coefficients, -np.inf)
        following[:region+1] = logsumexp(coefficients[:region+1, :, None]+regional[0], axis=1)
        hit = logsumexp(coefficients[:region+1, :, None]+regional[1], axis=1)
        following[1:region+2] = np.logaddexp(following[1:region+2], hit)
        coefficients = following
    value = logsumexp(coefficients, axis=1)-np.array([log(comb(regions, v)) for v in range(regions+1)])
    return np.minimum(0., value+float(tilt)*cutoff)


def geometry(data):
    if (data['bits'] != 16 or data['packet_bits'] != 8
            or data['width'] != 8*data['windows'] or 512 % data['windows']):
        raise ValueError('actual16-state maps and complete byte-region steps required')
    return dict(K=65536, N=131072, outer_groups=512, group_dimension=128,
        group_output_bits=256, regions=32, packet_bits=8, slots_per_region=512,
        physical_t=data['width'], physical_packet_slots=data['windows'],
        physical_steps_per_region=512//data['windows'],
        state_bits=16, zero_initial_state=True, final_flush=False,
        continuous_state_across_all_steps_and_regions=True, cutoff=13107,
        order='y=X+A*a; a_next=M*a+C*X')


def sources(name):
    result = maps.sources(name)
    for module in (sys.modules[__name__], logarithmic, rs_outer, rs_uniform_envelope):
        path = Path(module.__file__).resolve()
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    tilts = tuple(Fraction(value) for value in args.tilts)
    if (args.output.exists() or not tilts or min(tilts) <= 0 or len(set(tilts)) != len(tilts)
            or not 2 <= args.min_q <= args.max_q <= 512):
        raise ValueError('fresh output, distinct positive tilts and occupancy interval2..512 required')
    start = monotonic()
    data, record = maps.prepare(args.map)
    pins = sources(args.map)
    shape = geometry(data)
    envelope = rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta = envelope.beta
    beta_log = log(beta.numerator)-log(beta.denominator)
    counts = rs_outer.expected_group_support_counts(16, 8, 8, 2)
    count_logs = np.array([log(x.numerator)-log(x.denominator) if x else -np.inf for x in counts])
    best_q1 = np.zeros(33)
    result = dict(schema='packet8-finite-birth-family-proposal-1', geometry=shape,
        map_record=record, proposal_only=True, whole_code_certificate=False,
        has_outward_endpoints=False, local_arithmetic='floating exact finite formulas',
        global_arithmetic='scaled/logarithmic ordered hypergeometric placement',
        occupancy_range=[args.min_q, args.max_q], q1_checked=True,
        q2_exact_shells_checked=False, envelope=envelope.metadata(),
        source_sha256=pins, source_pins_verified_at_finish=False,
        input_tilts=list(map(str, tilts)), trials=[], best={}, tilt_choices={})
    print(f'{args.map}: fresh map census in {monotonic()-start:.2f}s; '
          f'profiles A={len(data["expansion_profiles"])} C={len(data["character_profiles"])}; '
          f'ranks={record["packet_restriction_rank_counts"]}', flush=True)
    for tilt in tilts:
        weighted, emission, diagnostics = moments(data, np.exp(-float(tilt)))
        local = operators(weighted, emission)
        regional, backend = placement(local, args.max_q,
            epochs=shape['physical_steps_per_region'], windows=data['windows'], force_log=args.force_log)
        # This restricted path is a NONNEGATIVE CONTRIBUTION to the same
        # comparison expression, not a stronger upper bound or failure lower bound.
        zero_regional, zero_backend = placement(local[:, :1, :1], args.max_q,
            epochs=shape['physical_steps_per_region'], windows=data['windows'], force_log=args.force_log)
        witnesses = {}
        for q in range(args.min_q, args.max_q+1):
            moment = logarithmic.log_power_matrix(uniform_mixture(regional, q), 32)
            margin = -(log(comb(512, q))+q*beta_log+float(tilt)*13107+moment)/log(2)
            zero_moment = 32*float(uniform_mixture(zero_regional, q)[0, 0])
            zero_margin = -(log(comb(512, q))+q*beta_log+float(tilt)*13107+zero_moment)/log(2)
            if not isfinite(margin):
                raise FloatingPointError('invalid global moment')
            witnesses[str(q)] = dict(estimated_margin_bits=margin, log_moment=moment,
                zero_state_only_margin_bits=zero_margin,
                full_over_zero_moment_log2=(moment-zero_moment)/log(2))
            if str(q) not in result['best'] or margin > result['best'][str(q)]:
                result['best'][str(q)] = margin
                result['tilt_choices'][str(q)] = str(tilt)
        best_q1 = np.minimum(best_q1, q1_support_logs(regional, regions=32, tilt=tilt, cutoff=13107))
        q1_terms = log(512)+count_logs+best_q1
        result['q1_margin_bits'] = -float(logsumexp(q1_terms))/log(2)
        result['q1_support_log_uppers_floating'] = best_q1.tolist()
        result['q1_dominant_supports'] = np.argsort(q1_terms)[::-1][:5].tolist()
        result['trials'].append(dict(tilt=str(tilt), backend=backend, zero_path_backend=zero_backend,
            local_diagnostics=diagnostics,
            q1_accumulated_margin_bits=result['q1_margin_bits'], witnesses=witnesses))
        if pins != sources(args.map):
            raise ArithmeticError('source changed during proposal')
        result['elapsed_seconds'] = monotonic()-start
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        worst = min(result['best'], key=result['best'].get)
        print(f'tilt={tilt}: weakest q={worst}, {result["best"][worst]:.6f} bits; '
              f'q1={result["q1_margin_bits"]:.6f}; {backend}; elapsed={monotonic()-start:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--map', choices=maps.NAMES, default='byte_native')
    parser.add_argument('--tilts', nargs='+', default=['.00256', '.00512', '.01', '.0256', '.0512', '.10', '.20', '.40'])
    parser.add_argument('--min-q', type=int, default=10)
    parser.add_argument('--max-q', type=int, default=256)
    parser.add_argument('--force-log', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
