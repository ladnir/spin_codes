"""Bounded actual24 local-slack diagnostics; no new state census.

The shadow and no-return operators are deliberately NOT probability upper
bounds. The shadow embeds below the true weighted physical kernel; deleting
all returns is a more optimistic, generally nonphysical ablation.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, exp, log, log1p
from pathlib import Path
import sys
from time import monotonic

import numpy as np

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent/'iteration3'
sys.path.insert(0, str(PREVIOUS))
import capped_gate
import potential_birth_gate as birth

DEFAULT_CACHE = PREVIOUS/'cache24_v2/theta_3_50.json'


def analytic_bounds(theta=Fraction(3, 50), alpha=Fraction(2, 5), q=16, regions=32):
    """Exact rational checks plus floating evaluations of analytic bounds."""
    theta, alpha = Fraction(theta), Fraction(alpha)
    if not 0 < theta < 1 or not 0 < alpha <= 1 or type(q) is not int or q < 1:
        raise ValueError('small positive tilt, fractional power and positive occupancy required')
    z = exp(-float(theta))
    delta = (1+z)**-8
    eta = 1/((1+z)**8-1)
    slots = regions*q
    # e^-theta >= 1-theta: both upper atom bounds are exact rationals.
    normalizer_lower = (2-theta)**8
    delta_upper = 1/normalizer_lower
    eta_upper = 1/(normalizer_lower-1)
    epsilon_upper = delta_upper**3
    returns_permitted_steps = slots//4
    # -ln(1-epsilon)<=epsilon/(1-epsilon), 1/ln2<3/2.
    return_fractional_upper = alpha*returns_permitted_steps*Fraction(3, 2)*epsilon_upper/(1-epsilon_upper)
    check_15 = None
    if theta == Fraction(3, 50) and alpha == Fraction(2, 5) and slots == 512:
        check_15 = (1-delta_upper)**2048 > Fraction(1, 2**15)
        if not check_15:
            raise ArithmeticError('claimed strict1.5-bit conservative bound failed')
    return dict(theta=str(theta), alpha=str(alpha), q=q, regions=regions,
        potential_packets=slots, delta=delta, eta=eta,
        unpowered_all_local_slack_bits=-slots*log1p(-delta)/log(2),
        same_basis_fractional_slack_bits=-float(alpha)*slots*log1p(-delta)/log(2),
        false_return_only_fractional_slack_bits=-float(alpha)*returns_permitted_steps*log1p(-delta**3)/log(2),
        rational_delta_upper=str(delta_upper), rational_eta_upper=str(eta_upper),
        strict_1_5_bit_bound_exact_integer_check=check_15,
        conservative_return_fractional_upper=str(return_fractional_upper),
        return_gain_below_1e_minus5_exact_check=return_fractional_upper < Fraction(1, 100000),
        interpretation='analytic local-kernel sandwich; not a bound on fractional representation slack')


def shadow_potential(potential, z, *, return_only=False):
    """Lower comparison from excluded-atom bounds; not a certificate kernel."""
    source = np.asarray(potential, dtype=float)
    if (source.ndim != 3 or source.shape[1] != source.shape[2] or len(source) != 9
            or not np.isfinite(source).all() or np.any(source < 0) or not 0 < z <= 1
            or np.any(source[0, 1:, 0] != 0)):
        raise ValueError('actual24 eight-slot family with zero-input/no-return preservation required')
    delta = (1+z)**-8
    result = source.copy()
    for j in range(1, 9):
        if not return_only:
            result[j, 1:, 1] *= 1-delta**min(j, 3)
        if j >= 4:
            result[j, 1:, 0] *= 1-delta**3
    return result


def run(args):
    if args.output.exists():
        raise ValueError('fresh diagnostic receipt required')
    start = monotonic()
    local, cache_record = capped_gate.load_local(args.cache)
    if cache_record['tilt'] != '3/50':
        raise ValueError('this bounded diagnostic is for cached theta=.06 only')
    pins = birth.source_pins()
    pins.update(cache_record['source_sha256'])
    for path in (Path(__file__).resolve(), HERE/'test_local_slack.py',
                 Path(capped_gate.__file__).resolve(), args.cache.resolve()):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    potential, _, algebra = birth.rebase_potential(birth.marked.potential_operators(local, 0.))
    z, theta, alpha, q = exp(-.06), Fraction(3, 50), Fraction(2, 5), 16
    no_return = potential.copy()
    no_return[:, 1:, 0] = 0.
    variants = [('unchanged', potential),
        ('false_return_shadow', shadow_potential(potential, z, return_only=True)),
        ('all_local_shadow', shadow_potential(potential, z)),
        ('delete_all_returns_NONPHYSICAL', no_return)]
    envelope = birth.prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    logbeta = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    result = dict(schema='actual24-local-slack-diagnostic-1', proposal_only=True,
        whole_code_certificate=False, has_outward_endpoints=False,
        construction_changed=False, diagnostic_altered_operators_are_not_upper_bounds=True,
        q=q, tilt=str(theta), alpha=str(alpha), group_steps=4,
        analytic=analytic_bounds(theta, alpha, q), algebra=algebra,
        source_sha256=pins, source_pins_verified_at_finish=False, variants=[])
    for name, operators in variants:
        grouped = birth.fine.fine_operators(birth.fine.tuple_products(operators, 4), alpha)
        margin, moment, backend = birth.evaluate(grouped, q, theta, alpha, logbeta, 13107, force_log=True)
        baseline = result['variants'][0]['margin_bits'] if result['variants'] else margin
        result['variants'].append(dict(name=name, margin_bits=margin,
            gain_over_unchanged_bits=margin-baseline, log_moment=moment, backend=backend))
        print(f'{name}: {margin:.12f} bits; gain={margin-baseline:.12f}', flush=True)
    if any(hashlib.sha256(Path(p).read_bytes()).hexdigest() != h for p, h in pins.items()):
        raise ArithmeticError('diagnostic source or cache changed')
    result['source_pins_verified_at_finish'] = True
    result['elapsed_seconds'] = monotonic()-start
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=DEFAULT_CACHE)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
