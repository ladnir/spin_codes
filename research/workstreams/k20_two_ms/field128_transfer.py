"""Bounded exact transfer to a systematic field128 outer at K=2^20.

This changes only the outer, retaining the authenticated small-outer s22
inner and routing. It authenticates the complete old certificate, then
uses cached sparse operators and the old, no-worse uniform-envelope tail.
No local census, numerical witness search, or encoder benchmark is run.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import comb, log2
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REFERENCE = ROOT / 'tmp/k20-small-outer-whole-v1.json'
REFERENCE_HASH = '90492fc15834d2e28a9ddc8289a68dbf7bd17984146524016f9361a9ec54fe65'
Q1 = ROOT / 'tmp/k20-small-outer-q1-p192-v1.json'
Q1_HASH = '27f1b3586b7c40b7f9429466f3baf6e3aa5bf751ea86cb5aed1d73f7d400ed29'
GROUPS = 8192
FIELD = 1 << 128


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def dyadic(pair):
    n, e = pair
    if type(n) is not int or n < 0 or type(e) is not int:
        raise ValueError('nonnegative exact dyadic required')
    return F(n << e) if e >= 0 else F(n, 1 << -e)


def ceil_dyadic(value, precision=256):
    """An outward endpoint using integer arithmetic only."""
    value = F(value)
    if not value:
        return [0, 0]
    if value < 0:
        raise ValueError('nonnegative rational required')
    exponent = value.numerator.bit_length() - value.denominator.bit_length() - precision
    if exponent >= 0:
        numerator, denominator = value.numerator, value.denominator << exponent
    else:
        numerator, denominator = value.numerator << -exponent, value.denominator
    mantissa = (numerator + denominator - 1) // denominator
    result = [mantissa, exponent]
    if dyadic(result) < value:
        raise ArithmeticError('endpoint rounded inward')
    return result


def margin(value):
    """Display only; acceptance uses exact rational comparison."""
    value = F(value)
    return log2(value.denominator) - log2(value.numerator)


def field_shells(half_packets=32, packet_bits=4):
    labels = (1 << packet_bits) - 1
    half = [0] + [comb(half_packets, v) * labels**v for v in range(1, half_packets + 1)]
    denominator = (1 << (packet_bits * half_packets)) - 1
    counts = [F(0)] * (2 * half_packets + 1)
    for a, left in enumerate(half):
        for b, right in enumerate(half):
            counts[a + b] += F(left * right, denominator)
    if counts[0] or counts[1] or sum(counts) != denominator:
        raise ArithmeticError('incorrect field outer measure')
    return counts


def transfer(reference, sparse, rs_counts):
    counts = field_shells()
    if len(rs_counts) != 65 or sum(rs_counts) != FIELD - 1:
        raise ArithmeticError('invalid reference shell measure')
    missing = [v for v in range(65) if counts[v] and not rs_counts[v]]
    if missing != list(range(2, 9)):
        raise ArithmeticError('unexpected missing reference supports')
    # F is the global event that any group has a word with packet support < 9.
    # This bound counts every local message, without assuming independence
    # between messages or conditioning the remaining group setups on F^c.
    rare = GROUPS * sum((counts[v] for v in missing), F())
    ratios = {v: counts[v] / rs_counts[v] for v in range(9, 65)}
    witness = max(ratios, key=ratios.__getitem__)
    ratio = ratios[witness]
    if any(counts[v] > ratio * rs_counts[v] for v in range(9, 65)):
        raise ArithmeticError('regular shell domination failed')
    probabilities = list(map(dyadic, sparse['support_probability_uppers']))
    if len(probabilities) != 65 or any(not 0 <= p <= 1 for p in probabilities):
        raise ArithmeticError('invalid sparse support endpoints')
    q1 = GROUPS * sum((a * p for a, p in zip(counts, probabilities)), F())
    q2 = ratio**2 * dyadic(reference['occupancy_uppers']['2'])
    # Retaining old q>=3 endpoints is conservative. The new uniform-input
    # comparison has a smaller beta; no fractional exponent needs replay.
    beta_new = F(1 << 256, FIELD - 1)
    beta_old = F(1 << 256, 65535**8)
    if not beta_new < beta_old:
        raise ArithmeticError('new uniform envelope must be no worse')
    tail = sum((dyadic(value) for q, value in reference['occupancy_uppers'].items()
                if int(q) >= 3), F())
    values = dict(rare=rare, q1=q1, q2_regular=q2, tail=tail)
    endpoints = {name: ceil_dyadic(value) for name, value in values.items()}
    total = sum(map(dyadic, endpoints.values()), F())
    upper = ceil_dyadic(total)
    return dict(exact_expected_shells=list(map(str, counts)),
                excluded_packet_supports=missing,
                regular_shell_ratio=str(ratio), regular_ratio_witness=witness,
                beta_new=str(beta_new), beta_old=str(beta_old),
                component_uppers=endpoints,
                component_margin_bits_display={name: margin(value) for name, value in values.items()},
                union_upper=upper, margin_bits_display=margin(dyadic(upper)),
                target_40_met=dyadic(upper) < F(1, 1 << 40))


def run():
    if digest(REFERENCE) != REFERENCE_HASH or digest(Q1) != Q1_HASH:
        raise ValueError('retained certificate identity changed')
    audit = load_module('k20_two_ms_old_audit', HERE.parent / 'k20_small_outer_screen/verify_whole.py')
    authentication = audit.verify(REFERENCE)
    reference = json.loads(REFERENCE.read_text(encoding='utf-8'))
    sparse = json.loads(Q1.read_text(encoding='utf-8'))
    if not any(Path(component['path']).resolve() == Q1.resolve() for component in reference['components']):
        raise ValueError('sparse receipt is not authenticated by the complete receipt')
    outer_path = HERE.parent / 'k16_design/rs_outer.py'
    inherited_pins = {str(Path(p).resolve()): h for p, h in reference['source_sha256'].items()}
    if inherited_pins.get(str(outer_path.resolve())) != digest(outer_path):
        raise ValueError('reference RS enumerator must already be authenticated')
    outer = load_module('k20_two_ms_rs_outer', outer_path)
    result = transfer(reference, sparse, outer.expected_group_support_counts(
        n=16, k=8, packet_bits=4, packets_per_symbol=4))
    result.update(schema='k20-systematic-field128-s22-transfer-1', K=1048576, N=2097152,
                  distance='1/10', minimum_distance=209716, bad_cutoff=209715,
                  outer_forward='x -> (x,L_c^T x), independent uniform nonzero c in GF(2^128) per group',
                  outer_transpose='(u,v) -> u+L_c v', groups=GROUPS,
                  inner_and_route='unchanged exact certified22 maps and independent uniform four-bit-packet routing',
                  map_sha256=reference['map_sha256'], continuous_state=True,
                  zero_initial_state=True, final_flush=False,
                  reference=dict(path=str(REFERENCE), sha256=REFERENCE_HASH),
                  sparse_reference=dict(path=str(Q1), sha256=Q1_HASH),
                  reference_audit=authentication,
                  proof_method='global rare-outer event, exact q1 shell fold, restricted expected-measure q2 domination, and unchanged no-worse uniform-envelope q>=3 endpoints',
                  all_occupancies_covered=True, local_operators_recomputed=False,
                  arithmetic='exact rationals; outward integer dyadic conversion; exact acceptance comparison',
                  source_sha256={str(Path(__file__).resolve()): digest(__file__),
                                 str(Path(audit.__file__).resolve()): digest(audit.__file__),
                                 str(outer_path.resolve()): digest(outer_path)},
                  scope='Proof-only transfer under ideal independent setup. No implementation or performance claim. Requires independent mathematical review before promotion.')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('fresh output required')
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps({name: result[name] for name in
                     ('margin_bits_display', 'component_margin_bits_display',
                      'regular_ratio_witness', 'target_40_met')}, indent=2))
