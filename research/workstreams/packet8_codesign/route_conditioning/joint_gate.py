"""Bounded joint H1/H2 route-conditioned gate for unchanged16-bit maps.

For each119-group subset, require H1>=1483 and H2>=2559, where H1 counts
occupied potential steps and H2 sums min(potential occupancy,2). The bad
route probability is bounded by the SUM of two exact rational bounds;
independence is neither asserted nor used. Good-message evaluation is
floating only and does not certify the complete code.
"""
import argparse
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import comb, isfinite, log
from pathlib import Path
from time import monotonic

import numpy as np

import gate
import cap_gate
import cap_counts
import route_counts


def marked(local, nu1, nu2, *, packet_bits=8):
    if not isfinite(nu1) or not isfinite(nu2) or min(nu1, nu2) < 0:
        raise ValueError('finite nonnegative route tilts required')
    result = gate.potential_operators(local, 0., packet_bits=packet_bits)
    occupancies = np.arange(len(local))
    exponent = float(nu1)*(occupancies > 0)+float(nu2)*np.minimum(occupancies, 2)
    result *= np.exp(exponent)[:, None, None]
    if not np.isfinite(result).all():
        raise FloatingPointError('marked operator overflow')
    return result


def exact_bad_union(q=119, h1=1483, h2=2559, x1=Fraction(3, 20), x2=Fraction(9, 35)):
    one = route_counts.rational_bad_route_bound(q, h1, x1)
    two = cap_counts.rational_bad_route_bound(q, h2, x2, cap=2)
    if (one.numerator << 60) > one.denominator or (two.numerator << 60) > two.denominator:
        raise ArithmeticError('both individual bounds must certify60bits')
    bound = one+two
    if bound <= 0 or bound > 1:
        raise ArithmeticError('positive probability upper bound required')
    floor = bound.denominator.bit_length()-bound.numerator.bit_length()
    if bound.numerator << floor > bound.denominator:
        floor -= 1
    identity = hex(bound.numerator)+'/'+hex(bound.denominator)
    return bound, dict(q=q, h1=h1, h2=h2, x1=str(x1), x2=str(x2),
        combination='exact rational sum; no independence assumption',
        individual_60bit_checks=[True, True], certified_margin_bits_floor=floor,
        exact_59bit_sum_check=(bound.numerator << 59) <= bound.denominator,
        rational_sum_sha256=hashlib.sha256(identity.encode()).hexdigest(),
        numerator_bit_length=bound.numerator.bit_length(), denominator_bit_length=bound.denominator.bit_length(),
        estimated_sum_margin_bits=(log(bound.denominator)-log(bound.numerator))/log(2),
        whole_code_certificate=False)


def sources():
    pins = cap_gate.sources()
    for path in (Path(__file__).resolve(), gate.HERE/'test_joint_gate.py'):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    pins, started = sources(), monotonic()
    q, h1, h2, theta = 119, 1483, 2559, Fraction(2, 5)
    bad, witness = exact_bad_union(q, h1, h2)
    logbad = log(bad.numerator)-log(bad.denominator)
    data, record = gate.maps.prepare('byte_native')
    weighted, emission, diagnostics = gate.prior.moments(data, np.exp(-float(theta)))
    local = gate.prior.operators(weighted, emission)
    beta = gate.prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2).beta
    fixed = log(comb(512, q))+q*(log(beta.numerator)-log(beta.denominator))+float(theta)*13107
    reference_local = cap_gate.marked(local, 3.5, cap=2)
    control_local = marked(local, 0., 3.5)
    np.testing.assert_allclose(control_local, reference_local, rtol=2e-15, atol=0.)
    reference, _ = gate.prior.placement(reference_local, q, epochs=64, windows=8)
    reference_moment = gate.prior.logarithmic.log_power_matrix(reference[q], 32)
    result = dict(schema='packet8-joint-H1-H2-floating-1', q=q, tilt=str(theta),
        map_record=record, geometry=gate.prior.geometry(data), unchanged_baseline_construction=True,
        proposal_only=True, whole_code_certificate=False, has_outward_endpoints=False,
        bad_route_bound_exact=True, good_message_bound_outward=False,
        bad_route_union=witness, local_diagnostics=diagnostics, source_sha256=pins,
        source_pins_verified_at_finish=False, trials=[], best_margin_bits=-float('inf'),
        q1_checked=False, all_occupancies_checked=False)
    choices = [(Fraction(0), Fraction(7, 2))]+list(product(
        (Fraction(1, 2), Fraction(1), Fraction(3, 2)),
        (Fraction(5, 2), Fraction(3), Fraction(7, 2))))
    for nu1, nu2 in choices:
        potential = marked(local, float(nu1), float(nu2))
        regional, backend = gate.prior.placement(potential, q, epochs=64, windows=8)
        moment = gate.prior.logarithmic.log_power_matrix(regional[q], 32)
        exponent = fixed-float(nu1)*h1-float(nu2)*h2+moment
        good = -exponent/log(2)
        combined = -float(np.logaddexp(exponent, logbad))/log(2)
        trial = dict(nu1=str(nu1), nu2=str(nu2), backend=backend,
            marked_log_moment=moment, good_message_margin_bits=good, combined_margin_bits=combined)
        if nu1 == 0:
            error = abs(moment-reference_moment)
            if error > 2e-8:
                raise ArithmeticError('H2-only regression failed')
            trial['H2_control_log_moment_discrepancy'] = error
            trial['H2_control_good_margin_bits'] = -(fixed-3.5*h2+reference_moment)/log(2)
        result['trials'].append(trial)
        if combined > result['best_margin_bits']:
            result['best_margin_bits'], result['best_choice'] = combined, trial
        if sources() != pins:
            raise ArithmeticError('source changed during gate')
        result['elapsed_seconds'] = monotonic()-started
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(f'joint nu1={nu1} nu2={nu2}: good={good:.6f}, combined={combined:.6f}, '
              f'elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
