"""Outward G4 macros from authenticated local upper comparison operators.

The selected dyadic upper family is thinned, rebased with exact rational
normalization, and divided by the known occupancy factor a^j. Conditional
tuple weights stay outside the entrywise fractional power.
"""
from __future__ import annotations

import os
for _name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_name] = '1'

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
from time import monotonic

import numpy as np

import local_outward24 as local
import outward_positive as op

HERE = Path(__file__).resolve().parent


def checked_pins(pins):
    return all(Path(path).is_file() and hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
               for path, digest in pins.items())


def load_local(path):
    path = Path(path).resolve()
    receipt = json.loads(path.read_text(encoding='utf-8'))
    if (receipt.get('schema') != 'packet8-actual24-local-outward-1'
            or not receipt.get('source_pins_verified_at_finish')
            or not receipt.get('outward_under_stated_arithmetic_contract')
            or not checked_pins(receipt['source_sha256'])):
        raise ValueError('authenticated outward local receipt required')
    if receipt['diagnostics']['arithmetic_contract'] != local.CONTRACT:
        raise ValueError('unrecognized local arithmetic contract')
    z = Fraction(receipt['z'])
    if z != Fraction(receipt['z_numerator'], receipt['z_denominator']):
        raise ValueError('inconsistent rational weight parameter')
    # JSON canonicalizes integer keys in the restriction-rank tables to text.
    expected = json.loads(json.dumps(local.maps24.make_maps()['record']))
    if receipt['map_record'] != expected:
        raise ValueError('literal baseline24 map identity required')
    lower = np.asarray(receipt['local_operator_lower'], dtype=float)
    upper = np.asarray(receipt['local_operator_upper'], dtype=float)
    if (lower.shape != (9, 10, 10) or upper.shape != lower.shape
            or not np.isfinite(lower).all() or not np.isfinite(upper).all()
            or np.any(lower < 0) or np.any(upper < lower)):
        raise ValueError('invalid local endpoint arrays')
    for key, array in (('local_operator_lower_hex', lower), ('local_operator_upper_hex', upper)):
        stored = np.asarray(receipt[key], dtype=object)
        if stored.shape != array.shape or any(float(value).hex() != text
                for value, text in zip(array.flat, stored.flat)):
            raise ValueError('decimal and hexadecimal endpoints disagree')
    if receipt['z_binary64_hex'] != float(z).hex() or Fraction.from_float(float(z)) != z:
        raise ValueError('local receipt must identify an exactly representable dyadic z')
    validate_family(upper)
    return upper, z, receipt


def exact_array(value):
    array = np.asarray(value)
    result = np.empty(array.shape, dtype=object)
    for index in np.ndindex(array.shape):
        item = array[index]
        if isinstance(item, Fraction):
            result[index] = item
        elif isinstance(item, (int, np.integer)):
            result[index] = Fraction(int(item))
        else:
            result[index] = Fraction.from_float(float(item))
    return result


def zero_array(shape):
    result = np.empty(shape, dtype=object)
    result.fill(Fraction(0))
    return result


def upper_array(exact):
    result = np.empty(exact.shape, dtype=float)
    for index in np.ndindex(exact.shape):
        result[index] = op.upper_float(exact[index])
    return result


def validate_family(value):
    array = np.asarray(value)
    if (array.ndim != 3 or array.shape[1] != array.shape[2]
            or array.shape[1] != len(array)+1 or not 1 <= len(array)-1 <= 8):
        raise ValueError('delta0/U/one-birth-per-occupancy family required')
    if array.dtype != object and not np.isfinite(array).all():
        raise ValueError('finite entries required')
    if (np.any(array < 0) or np.any(array[:, 1:, 2:] != 0)
            or np.any(array[:, 0, 1] != 0) or np.any(array[0, 0, 1:] != 0)
            or np.any(array[0, 1:, 0] != 0) or array[0, 0, 0] != 1):
        raise ValueError('comparison family violates required structural zeros')
    return array


def thin_exact(physical, *, packet_bits=8):
    """Exact rational thinning of selected upper dyadics, before conversion."""
    validate_family(physical)
    if type(packet_bits) is not int or not 1 <= packet_bits <= 8:
        raise ValueError('packet width1..8 required')
    old = exact_array(physical)
    W, size, labels = len(old)-1, old.shape[1], (1 << packet_bits)-1
    result = zero_array(old.shape)
    for j in range(W+1):
        for k in range(j+1):
            weight = Fraction(comb(j, k)*labels**k, 1 << (packet_bits*j))
            for source in range(size):
                for target in range(size):
                    result[j, source, target] += weight*old[k, source, target]
    return result


def rebase_exact(selected_potential):
    """Exactly rebase a chosen positive upper family; no inverse is needed."""
    validate_family(selected_potential)
    old = exact_array(selected_potential)
    W, size = len(old)-1, old.shape[1]
    change = zero_array((size, size))
    change[0, 0] = change[1, 1] = Fraction(1)
    mass = [sum(old[j, 0, 2:]) for j in range(1, W+1)]
    if any(value <= 0 for value in mass):
        raise ValueError('every nonempty potential occupancy needs positive birth mass')
    for i in range(W):
        for k in range(W):
            change[i+2, k+2] = old[i+1, 0, k+2]/mass[i]
    new = zero_array(old.shape)
    for j in range(W+1):
        new[j, 0, 0] = old[j, 0, 0]
        if j:
            new[j, 0, j+1] = mass[j-1]
        for source in range(1, size):
            for target in (0, 1):
                new[j, source, target] = sum(change[source, k]*old[j, k, target] for k in range(size))
    if any(sum(row) != 1 for row in change):
        raise ArithmeticError('exact change of basis is not stochastic')
    for j in range(W+1):
        if np.any(new[j]@change != change@old[j]):
            raise ArithmeticError('exact comparison intertwining failed')
    return new, change


def prepare_physical_upper(physical, z, *, packet_bits=8):
    """Select upper dyadic thinning, rebase exactly, and extract a^j."""
    z = Fraction(z)
    if not 0 < z < 1:
        raise ValueError('weight parameter must lie in(0,1)')
    thinned = upper_array(thin_exact(physical, packet_bits=packet_bits))
    exact, change = rebase_exact(thinned)
    a = ((1+z)/2)**packet_bits
    normalized = exact.copy()
    for j in range(len(exact)):
        normalized[j] /= a**j
    upper = upper_array(normalized)
    validate_family(upper)
    metadata = dict(packet_bits=packet_bits, a=str(a), a_numerator=a.numerator, a_denominator=a.denominator,
        exact_stochastic_change_of_basis=True, exact_intertwining_before_final_rounding=True,
        birth_basis='artificial normalized mixtures defined from selected upper thinned dyadics',
        normalization='operator j divided upward by exact a^j',
        change_of_basis_rationals=[[str(x) for x in row] for row in change],
        selected_thinned_operators=thinned.tolist())
    return upper, a, metadata


def macro_operators(potential, alpha, *, group_steps=4, tangent_bins=256):
    """Outward chronological products, powers, and conditional tuple average."""
    op.check_runtime()
    potential = np.asarray(validate_family(potential), dtype=float)
    alpha = Fraction(alpha)
    if not 0 < alpha <= 1 or group_steps not in (1, 2, 4):
        raise ValueError('alpha in(0,1] and group size1/2/4 required')
    W, size = len(potential)-1, potential.shape[1]
    choices = np.array([comb(W, j) for j in range(W+1)], dtype=np.int64)
    products = potential.copy()
    totals = np.arange(W+1, dtype=np.int64)
    multiplicities = choices.copy()
    for _ in range(1, group_steps):
        products = op.matmul_upper(products[:, None], potential[None]).reshape(-1, size, size)
        totals = (totals[:, None]+np.arange(W+1)[None]).reshape(-1)
        multiplicities = (multiplicities[:, None]*choices[None]).reshape(-1)
    weights = np.array([op.upper_float(Fraction(int(n), comb(W*group_steps, int(j))))
                        for n, j in zip(multiplicities, totals)])
    # The integer tuple multiplicity and conditional denominator remain
    # outside alpha; only the chronological matrix product is powered.
    powered = op.fractional_power_upper(products, alpha, bins=tangent_bins)
    weighted = op.multiply_upper(powered, weights[:, None, None])
    macro = np.zeros((W*group_steps+1, size, size))
    for total in range(W*group_steps+1):
        macro[total] = op.sum_upper(weighted[totals == total], axis=0)
    return macro, dict(group_steps=group_steps, physical_windows=W, macro_windows=W*group_steps,
        tuple_count=(W+1)**group_steps, alpha=str(alpha), tangent_bins=tangent_bins,
        exact_hypergeometric_weights_rounded_up=True, tuple_weights_outside_fractional_power=True,
        structural_zeros_preserved=True, positive_arithmetic_runtime_check_passed=True,
        arithmetic_contract=op.__doc__.strip())


def source_pins():
    pins = local.sources()
    for path in (Path(__file__).resolve(), HERE/'test_macro_outward.py', Path(op.__file__).resolve(),
                 HERE/'test_outward_positive.py'):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


def run(args):
    if args.output.exists():
        raise ValueError('fresh macro receipt path required')
    physical, z, saved = load_local(args.local)
    pins = source_pins()
    pins.update(saved['source_sha256'])
    pins[str(args.local.resolve())] = hashlib.sha256(args.local.read_bytes()).hexdigest()
    if not checked_pins(pins):
        raise ArithmeticError('initial source authentication failed')
    started = monotonic()
    potential, a, algebra = prepare_physical_upper(physical, z)
    macro, diagnostic = macro_operators(potential, Fraction(args.alpha), tangent_bins=args.bins)
    if not checked_pins(pins):
        raise ArithmeticError('source changed during outward macro construction')
    result = dict(schema='packet8-actual24-G4-outward-1', z=str(z), alpha=str(Fraction(args.alpha)),
        z_numerator=z.numerator, z_denominator=z.denominator, z_binary64_hex=float(z).hex(),
        a=str(a), a_numerator=a.numerator, a_denominator=a.denominator,
        map_record=saved['map_record'], local_receipt=str(args.local.resolve()),
        normalized_potential_operators=potential.tolist(), macro_operator_upper=macro.tolist(),
        macro_operator_upper_hex=[[[float(x).hex() for x in row] for row in matrix] for matrix in macro],
        source_sha256=pins, source_pins_verified_at_finish=True,
        algebra=algebra, diagnostics=diagnostic, elapsed_seconds=monotonic()-started,
        local_arithmetic_contract=saved['diagnostics']['arithmetic_contract'],
        outward_under_stated_arithmetic_contract=True, whole_code_certificate=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(f'outward G4 macro saved: alpha={args.alpha}, {diagnostic["tuple_count"]} tuples, '
          f'elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local', type=Path, default=HERE/'local_theta003_v1.json')
    parser.add_argument('--alpha', default='1/2')
    parser.add_argument('--bins', type=int, default=256)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
