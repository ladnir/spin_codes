"""Floating exact-formula gate for the actual24-bit byte-native candidate.

Large weighted-syndrome vectors are inverted one occupancy at a time and
immediately compressed by the expansion profile. No16-bit moment is reused.
"""
from __future__ import annotations

import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '1'

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from math import comb, log
from pathlib import Path
import sys
from time import monotonic

import numpy as np
from scipy.special import logsumexp

import maps24

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import screen as prior


def profile_polynomials(data, z, *, character=False):
    b, W = data['packet_bits'], data['windows']
    labels = (1 << b)-1
    z = np.longdouble(z)
    if character:
        inactive = [1.]*(b+1)
        active = [((1+z)**(b-r)*(1-z)**r-1)/labels for r in range(b+1)]
    else:
        inactive = [z**r for r in range(b+1)]
        active = [((1+z)**b-z**r)/labels for r in range(b+1)]
    return np.asarray(prior.profile_polynomials(data['profiles'], inactive, active), dtype=float)


def walsh_in_place(values, *, block_bits=18):
    """Unnormalized float64 FWHT, with at most2^block_bits scratch entries."""
    if values.dtype != np.float64 or values.ndim != 1 or len(values) & (len(values)-1):
        raise ValueError('one power-of-two float64 vector required')
    S, cap, length = len(values), 1 << block_bits, 1
    while length < S:
        if length < cap:
            for start in range(0, S, 2*cap):
                pairs = values[start:min(S, start+2*cap)].reshape(-1, 2*length)
                left = pairs[:, :length].copy()
                pairs[:, :length] += pairs[:, length:]
                pairs[:, length:] *= -1
                pairs[:, length:] += left
        else:
            for start in range(0, S, 2*length):
                for offset in range(0, length, cap):
                    left_view = values[start+offset:start+min(length, offset+cap)]
                    right_view = values[start+length+offset:start+length+min(length, offset+cap)]
                    left = left_view.copy()
                    left_view += right_view
                    right_view *= -1
                    right_view += left
        length *= 2
    return values


def direct_birth(data, z, j=1):
    if j != 1:
        raise ValueError('only single-packet census is bounded here')
    S, b, W = 1 << data['bits'], data['packet_bits'], data['windows']
    mass = np.zeros(len(data['profiles']))
    for h in range(W):
        columns = data['columns'][b*h:b*(h+1)]
        states = np.array([maps24.apply(columns, label) for label in range(1, 1 << b)], dtype=np.uint32)
        indices = maps24.profile_indices(data, states)
        weights = np.array([float(z)**label.bit_count()/(W*((1 << b)-1)) for label in range(1, 1 << b)])
        mass += np.bincount(indices, weights=weights, minlength=len(mass))
    return mass


def local_operators(data, counts, z, *, only_q1=False, progress=False):
    """Return the exact finite birth closure evaluated in float64.

    The only analytic relaxation is the same dominated-uniform outgoing
    measure as the small-state screen. Floating clipping is explicitly logged.
    """
    S, W, b = 1 << data['bits'], data['windows'], data['packet_bits']
    L, labels = S-1, (1 << b)-1
    expansion = profile_polynomials(data, z)
    uniform = (counts['histogram'] @ expansion-expansion[maps24.profile_indices(data, [0])[0]])/L
    J = 1 if only_q1 else W
    zero = np.zeros(J+1)
    zero[0] = 1.
    born = np.zeros((J, len(data['profiles'])))
    born[0] = direct_birth(data, z)
    diagnostic = dict(arithmetic='float64 exact finite formulas; not outward', clipping=[],
        state_count=S, profile_count=len(data['profiles']), inverse_transform_count=0)
    if not only_q1:
        if counts['character_indices'] is None:
            raise ValueError('character census needed for the middle gate')
        character = profile_polynomials(data, z, character=True)
        packet = ((1+float(z))**b-1)/labels
        for j in range(2, W+1):
            started = monotonic()
            values = character[:, j][counts['character_indices']].copy()
            walsh_in_place(values)
            values /= S
            expected = packet**j
            minimum = float(values.min())
            negative = -float(values[values < 0].sum())
            if minimum < -2e-10*expected:
                raise ArithmeticError('unstable weighted syndrome inversion')
            np.maximum(values, 0., out=values)
            zero[j] = float(values[0]) if j > 3 else 0.
            values[0] = 0.
            # np.bincount converts all16M indices to int64 internally;
            # chunk it to keep the peak memory bounded.
            for start in range(0, S, 1 << 18):
                end = min(S, start+(1 << 18))
                born[j-1] += np.bincount(counts['expansion_indices'][start:end],
                    weights=values[start:end], minlength=len(data['profiles']))
            total = float(born[j-1].sum()+zero[j])
            if abs(total-expected) > 2e-9*expected:
                raise ArithmeticError('weighted syndrome total failed')
            diagnostic['clipping'].append(dict(j=j, minimum=minimum, negative_mass=negative,
                relative_mass_error=abs(total-expected)/expected, seconds=monotonic()-started))
            diagnostic['inverse_transform_count'] += 1
            if progress:
                print(f'  local j={j}: {monotonic()-started:.2f}s, clippedmass={negative:.3g}', flush=True)
            del values
    masses = born.sum(axis=1)
    if np.any(masses <= 0):
        raise ArithmeticError('all retained births must have positive mass')
    means = (born @ expansion)/masses[:, None]
    local = np.zeros((J+1, J+2, J+2))
    for j in range(J+1):
        local[j, 0, 0] = zero[j]
        if j:
            local[j, 0, j+1] = masses[j-1]
        local[j, 1, 1] = uniform[j]
        local[j, 2:, 1] = means[:, j]
        if j:
            local[j, 1:, 0] = local[j, 1:, 1]/L
    if np.any(local < 0) or not np.isfinite(local).all():
        raise ArithmeticError('invalid local comparison matrices')
    return local, diagnostic


def sources():
    result = maps24.sources()
    for module in (sys.modules[__name__], prior, prior.maps, prior.logarithmic,
                   prior.rs_outer, prior.rs_uniform_envelope):
        path = Path(module.__file__).resolve()
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    tilts = tuple(Fraction(x) for x in args.tilts)
    qs = tuple(args.q)
    if args.output.exists() or not tilts or min(tilts) <= 0 or any(q < 2 or q > 512 for q in qs):
        raise ValueError('fresh output, positive tilts and occupancies2..512 required')
    pins, started = sources(), monotonic()
    data = maps24.make_maps()
    counts = maps24.census(data, include_character=not args.q1_only)
    census_seconds = monotonic()-started
    print(f'actual24-bit census: {census_seconds:.2f}s; '
          f'Aprofiles={np.count_nonzero(counts["histogram"])}, '
          f'minweight={np.flatnonzero(counts["spectrum"])[1]}', flush=True)
    envelope = prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta = envelope.beta
    logbeta = log(beta.numerator)-log(beta.denominator)
    outer = prior.rs_outer.expected_group_support_counts(16, 8, 8, 2)
    logs = np.array([log(x.numerator)-log(x.denominator) if x else -np.inf for x in outer])
    best_q1 = np.zeros(33)
    result = dict(schema='actual24-bit-width8-floating-1', map_record=data['record'],
        geometry=dict(K=65536, N=131072, groups=512, regions=32, slots_per_region=512,
            physical_bits=64, physical_slots=8, steps_per_region=64, state_bits=24,
            zero_initial_state=True, final_flush=False, continuous_state=True, cutoff=13107),
        proposal_only=True, whole_code_certificate=False, has_outward_endpoints=False,
        actual_larger_state_maps=True, q1_only=args.q1_only, occupancy_values=list(qs),
        source_sha256=pins, source_pins_verified_at_finish=False, census_seconds=census_seconds,
        expansion_spectrum={str(i): int(n) for i, n in enumerate(counts['spectrum']) if n},
        expansion_profile_count=int(np.count_nonzero(counts['histogram'])),
        envelope=envelope.metadata(), trials=[], best={}, tilt_choices={})
    for tilt in tilts:
        local_started = monotonic()
        local, diagnostic = local_operators(data, counts, np.exp(-float(tilt)),
            only_q1=args.q1_only, progress=True)
        local_seconds = monotonic()-local_started
        # At q1, no physical step receives two active packets. Pad exact zeros
        # to satisfy the generic eight-slot placement interface, degree1 only.
        if args.q1_only:
            padded = np.zeros((9, 3, 3))
            padded[:2] = local
            regional, backend = prior.placement(padded, 1, epochs=64, windows=8)
        else:
            regional, backend = prior.placement(local, max(qs), epochs=64, windows=8)
        best_q1 = np.minimum(best_q1, prior.q1_support_logs(regional, regions=32, tilt=tilt, cutoff=13107))
        result['q1_margin_bits'] = -float(logsumexp(log(512)+logs+best_q1))/log(2)
        result['q1_support_log_uppers_floating'] = best_q1.tolist()
        witnesses = {}
        if not args.q1_only:
            for q in qs:
                moment = prior.logarithmic.log_power_matrix(prior.uniform_mixture(regional, q), 32)
                margin = -(log(comb(512, q))+q*logbeta+float(tilt)*13107+moment)/log(2)
                witnesses[str(q)] = dict(estimated_margin_bits=margin, log_moment=moment)
                if str(q) not in result['best'] or margin > result['best'][str(q)]:
                    result['best'][str(q)] = margin
                    result['tilt_choices'][str(q)] = str(tilt)
        result['trials'].append(dict(tilt=str(tilt), local_seconds=local_seconds,
            local_diagnostics=diagnostic, backend=backend, witnesses=witnesses,
            q1_accumulated_margin_bits=result['q1_margin_bits']))
        if sources() != pins:
            raise ArithmeticError('source changed during screen')
        result['elapsed_seconds'] = monotonic()-started
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(f'tilt={tilt}: q1={result["q1_margin_bits"]:.6f}; '
              f'middle={result["best"]}; local={local_seconds:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts', nargs='+', default=['.00256', '.00384', '.00512', '.00768'])
    parser.add_argument('--q', nargs='+', type=int, default=[64, 119, 128, 256])
    parser.add_argument('--q1-only', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
