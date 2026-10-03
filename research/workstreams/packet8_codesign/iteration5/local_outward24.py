"""Outward local comparison matrices from exact profile polynomials.

The large Walsh transform uses an explicit IEEE-binary64 error budget, not
directed rounding. Endpoint conversions and the final profile arithmetic
are checked with exact Fractions. See LOCAL_ARITHMETIC.md for the contract.
"""
from __future__ import annotations

import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '1'

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, isfinite, nextafter, inf
from pathlib import Path
import sys
from time import monotonic

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'larger_state'))
import maps24
import screen24

ROUNDING_UNIT = Fraction(1, 1 << 53)
CONTRACT = (
    'IEEE binary64, round-to-nearest; NumPy add/subtract/negate/ldexp are correctly rounded; '
    'no reassociation or fast-math; each weighted bincount performs at most one binary64 '
    'addition per input atom; no exceptional or subnormal intermediate in the audited circuit. '
    'Exact Fraction comparisons certify conversions, not the platform arithmetic contract.'
)


def fraction(value):
    return value if isinstance(value, Fraction) else Fraction.from_float(float(value))


def enclose(value):
    """Check dyadic endpoints by exact comparisons; conversion is not trusted."""
    value = Fraction(value)
    rounded = float(value)
    if not isfinite(rounded):
        raise OverflowError('rational endpoint outside binary64 range')
    lower = upper = rounded
    while fraction(lower) > value:
        lower = nextafter(lower, -inf)
    while fraction(upper) < value:
        upper = nextafter(upper, inf)
    if not (isfinite(lower) and isfinite(upper)):
        raise OverflowError('finite outward endpoints required')
    return lower, upper


def gamma(count):
    if type(count) is not int or not 0 <= count < 1 << 53:
        raise ValueError('nonnegative integer rounding count below2^53 required')
    return Fraction(count, (1 << 53)-count)


def validate(data, z):
    z = Fraction(z)
    if not 0 < z < 1 or z.denominator & (z.denominator-1):
        raise ValueError('z must be a dyadic rational strictly between0 and1')
    b, W = data['packet_bits'], data['windows']
    restrictions = data['record']['packet_restriction_rank_counts']
    for j in range(1, min(3, W)+1):
        if restrictions[j] != {j*b: comb(W, j)}:
            raise ValueError('exact full-rank restrictions through three packets required')
    return z


def coefficient_product(inactive, active):
    values = [1]
    for a, b in zip(inactive, active):
        following = [0]*(len(values)+1)
        for j, value in enumerate(values):
            following[j] += value*a
            following[j+1] += value*b
        values = following
    return values


def exact_profiles(data, z, selected, *, character=False):
    """Integer numerators by profile and common denominator by occupancy."""
    z = validate(data, z)
    b, W, u, v = data['packet_bits'], data['windows'], z.numerator, z.denominator
    labels = (1 << b)-1
    if character:
        inactive = [1]*(b+1)
        active = [(v+u)**(b-r)*(v-u)**r-v**b for r in range(b+1)]
        denominators = [comb(W, j)*(labels*v**b)**j for j in range(W+1)]
    else:
        inactive = [u**r*v**(b-r) for r in range(b+1)]
        active = [(v+u)**b-value for value in inactive]
        denominators = [comb(W, j)*labels**j*v**(b*W) for j in range(W+1)]
    numerators = {}
    for index in map(int, selected):
        weights = [r for r, count in enumerate(data['profiles'][index]) for _ in range(int(count))]
        if len(weights) != W:
            raise ArithmeticError('profile has wrong packet count')
        numerators[index] = coefficient_product([inactive[r] for r in weights], [active[r] for r in weights])
    return numerators, denominators


def index_histogram(indices, size, *, chunk=1 << 18):
    result = np.zeros(size, dtype=np.int64)
    for start in range(0, len(indices), chunk):
        result += np.bincount(indices[start:start+chunk], minlength=size)
    return result


def walsh_enclosure(profile_numerators, denominator, character_indices, bits, *, block_bits=18):
    """One bounded float FWHT, returning its center and a uniform exact error."""
    S = 1 << bits
    if len(character_indices) != S:
        raise ValueError('one character index per state required')
    lookup = np.zeros(max(profile_numerators)+1, dtype=np.float64)
    delta, largest, quantum_bits = Fraction(0), Fraction(0), 0
    for index, numerator in profile_numerators.items():
        exact = Fraction(numerator, denominator)
        rounded = float(exact)
        if not isfinite(rounded):
            raise FloatingPointError('nonfinite character input')
        stored = fraction(rounded)
        delta = max(delta, abs(stored-exact))
        largest = max(largest, abs(stored))
        quantum_bits = max(quantum_bits, stored.denominator.bit_length()-1)
        lookup[index] = rounded
    # Every input is an integer multiple of2^-quantum_bits. Addition and
    # rounding preserve that lattice. Scaling by2^-bits is also exact.
    if quantum_bits+bits > 1022:
        raise FloatingPointError('Walsh lattice does not exclude subnormal intermediates')
    if largest*S/(1-gamma(bits)) > Fraction(1 << 1022):
        raise FloatingPointError('Walsh magnitude does not exclude overflow')
    values = lookup[character_indices].copy()
    with np.errstate(over='raise', under='raise', invalid='raise', divide='raise'):
        screen24.walsh_in_place(values, block_bits=block_bits)
        np.ldexp(values, -bits, out=values)
    if not np.isfinite(values).all():
        raise FloatingPointError('nonfinite Walsh output')
    error = delta+gamma(bits)*largest
    return values, error, dict(
        input_rounding_error=str(delta), transform_rounding_error=str(gamma(bits)*largest),
        uniform_absolute_error=str(error), rounded_input_max_abs=str(largest),
        common_input_quantum_bits=quantum_bits, scaled_quantum_bits=quantum_bits+bits,
        rounding_stages=bits, no_subnormal_or_overflow_by_lattice_and_magnitude=True)


def compress_enclosure(values, error, expansion_indices, histogram, *, chunk=1 << 18):
    """Bounds sums over nonzero states, including clipping and bincount errors."""
    S, count = len(values), len(histogram)
    if len(expansion_indices) != S or int(histogram.sum()) != S:
        raise ValueError('inconsistent exact state census')
    zero_value = fraction(values[0])
    zero_interval = (max(Fraction(0), zero_value-error), max(Fraction(0), zero_value+error))
    # Projection onto the nonnegative half-line does not increase error from
    # the genuine, nonnegative weighted syndrome mass.
    np.maximum(values, 0., out=values)
    values[0] = 0.
    accumulated = np.zeros(count, dtype=np.float64)
    chunks = (S+chunk-1)//chunk
    with np.errstate(over='raise', under='raise', invalid='raise'):
        for start in range(0, S, chunk):
            accumulated += np.bincount(expansion_indices[start:start+chunk],
                                      weights=values[start:start+chunk], minlength=count)
    if not np.isfinite(accumulated).all():
        raise FloatingPointError('nonfinite profile sum')
    rounding = gamma(S+chunks)
    populations = histogram.copy()
    populations[int(expansion_indices[0])] -= 1
    result = {}
    for index in np.flatnonzero(populations):
        center = fraction(accumulated[index])
        uncertainty = int(populations[index])*error
        result[int(index)] = (max(Fraction(0), center/(1+rounding)-uncertainty),
                             center/(1-rounding)+uncertainty)
    return result, zero_interval, dict(positive_sum_rounding_count=S+chunks,
        positive_sum_gamma=str(rounding), projection_preserves_absolute_error=True,
        nonzero_state_count=int(populations.sum()))


def direct_birth_one(data, z):
    """Exact2040-atom birth census (or its smaller-field analogue)."""
    b, W, u, v = data['packet_bits'], data['windows'], z.numerator, z.denominator
    labels = (1 << b)-1
    counts = {}
    for h in range(W):
        states = [maps24.apply(data['columns'][b*h:b*(h+1)], value) for value in range(1, labels+1)]
        indices = maps24.profile_indices(data, states)
        for value, state, index in zip(range(1, labels+1), states, indices):
            if state == 0:
                raise ArithmeticError('single-packet restriction unexpectedly has a kernel')
            weight = value.bit_count()
            counts[int(index)] = counts.get(int(index), 0)+u**weight*v**(b-weight)
    denominator = W*labels*v**b
    return {index: (Fraction(value, denominator),)*2 for index, value in counts.items()}


def local_operators(data, census, z, *, progress=False):
    """Enclose the filled comparison family in its genuine physical birth basis.

    The lower array is a lower endpoint for that comparison family, not a
    lower bound on the physical transition kernel.
    """
    z = validate(data, z)
    bits, b, W = data['bits'], data['packet_bits'], data['windows']
    S, L, labels = 1 << bits, (1 << bits)-1, (1 << b)-1
    histogram = census['histogram']
    if int(histogram.sum()) != S or census['character_indices'] is None:
        raise ValueError('complete literal expansion and character census required')
    expansion_selected = np.flatnonzero(histogram)
    character_hist = index_histogram(census['character_indices'], len(data['profiles']))
    characters, character_den = exact_profiles(data, z, np.flatnonzero(character_hist), character=True)
    emissions, emission_den = exact_profiles(data, z, expansion_selected)
    zero_profile = int(census['expansion_indices'][0])
    if int(histogram[zero_profile]) != 1 or data['profile_weights'][zero_profile] != 0:
        raise ArithmeticError('full-rank expansion must have exactly one zero image')
    nonzero_hist = histogram.copy()
    nonzero_hist[zero_profile] -= 1
    uniform = [Fraction(sum(int(nonzero_hist[p])*emissions[int(p)][j] for p in expansion_selected),
                        L*emission_den[j]) for j in range(W+1)]
    births = [None]+[direct_birth_one(data, z)]
    zeros = [(Fraction(1), Fraction(1)), (Fraction(0), Fraction(0))]
    packet = ((1+z)**b-1)/labels
    diagnostics = []
    for j in range(2, W+1):
        started = monotonic()
        values, error, transform = walsh_enclosure(
            {p: row[j] for p, row in characters.items()}, character_den[j], census['character_indices'], bits)
        birth, zero, compression = compress_enclosure(values, error, census['expansion_indices'], histogram)
        if j <= 3:
            zero = (Fraction(0), Fraction(0))
        if zero[0] > packet**j:
            raise ArithmeticError('zero-state lower mass exceeds total input moment')
        zero = (zero[0], min(zero[1], packet**j))
        births.append(birth)
        zeros.append(zero)
        diagnostics.append(dict(j=j, transform=transform, compression=compression, seconds=monotonic()-started))
        if progress:
            print(f'outward local j={j}: epsilon<={float(error):.3g}, elapsed={monotonic()-started:.2f}s', flush=True)
    birth_masses = [(Fraction(0), Fraction(0))]
    for i in range(1, W+1):
        total = packet**i
        lower, upper = total-zeros[i][1], total-zeros[i][0]
        if lower <= 0:
            raise ArithmeticError('birth normalization lacks a positive lower endpoint')
        birth_masses.append((lower, upper))
    low, high = np.zeros((W+1, W+2, W+2)), np.zeros((W+1, W+2, W+2))
    def assign(j, source, dest, lower, upper=None):
        upper = lower if upper is None else upper
        if not 0 <= lower <= upper:
            raise ArithmeticError('invalid positive rational interval')
        low[j, source, dest] = enclose(lower)[0]
        high[j, source, dest] = enclose(upper)[1]
    for j in range(W+1):
        assign(j, 0, 0, *zeros[j])
        if j:
            assign(j, 0, j+1, *birth_masses[j])
        assign(j, 1, 1, uniform[j])
        if j:
            assign(j, 1, 0, uniform[j]/L)
        for i in range(1, W+1):
            numerator_low = sum(interval[0]*emissions[p][j] for p, interval in births[i].items())
            numerator_high = sum(interval[1]*emissions[p][j] for p, interval in births[i].items())
            mean_low = numerator_low/(emission_den[j]*birth_masses[i][1])
            mean_high = numerator_high/(emission_den[j]*birth_masses[i][0])
            # Every exact emitted Hamming weight lies in0..bW.
            mean_low = max(z**(b*W), mean_low)
            mean_high = min(Fraction(1), mean_high)
            assign(j, i+1, 1, mean_low, mean_high)
            if j:
                assign(j, i+1, 0, mean_low/L, mean_high/L)
    positive = high > 0
    widths = high-low
    return low, high, dict(arithmetic_contract=CONTRACT, z=str(z), state_count=S,
        character_profile_count=len(characters), expansion_profile_count=len(emissions),
        exact_uniform_means=True, exact_single_packet_birth=True, exact_zero_birth_through_occupancy=3,
        normalization='genuine physical birth family; exact total moment minus enclosed zero mass',
        lower_endpoints_are_for_comparison_not_physical_kernel=True,
        max_absolute_width=float(widths.max()),
        max_relative_width=float(np.max(widths[positive]/high[positive])), transforms=diagnostics)


def sources():
    pins = screen24.sources()
    for path in (Path(__file__).resolve(), HERE/'test_local_outward24.py'):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


def run(args):
    if args.output.exists():
        raise ValueError('fresh output path required')
    pins, started = sources(), monotonic()
    data = maps24.make_maps()
    census = maps24.census(data)
    lower, upper, diagnostic = local_operators(data, census, Fraction(args.z), progress=True)
    if sources() != pins:
        raise ArithmeticError('source changed during local certificate run')
    exact_z = Fraction(args.z)
    receipt = dict(schema='packet8-actual24-local-outward-1', z=str(exact_z),
        z_numerator=exact_z.numerator, z_denominator=exact_z.denominator,
        z_binary64_hex=float(exact_z).hex() if fraction(float(exact_z)) == exact_z else None,
        map_record=data['record'],
        local_operator_lower=lower.tolist(), local_operator_upper=upper.tolist(),
        local_operator_lower_hex=[[[float(x).hex() for x in row] for row in matrix] for matrix in lower],
        local_operator_upper_hex=[[[float(x).hex() for x in row] for row in matrix] for matrix in upper],
        source_sha256=pins, source_pins_verified_at_finish=True, diagnostics=diagnostic,
        outward_under_stated_arithmetic_contract=True, whole_code_certificate=False,
        elapsed_seconds=monotonic()-started)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(f'local outward receipt saved: max relative width={diagnostic["max_relative_width"]:.3g}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--z', required=True, help='exact dyadic rational, not a floating weight tilt')
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args())
