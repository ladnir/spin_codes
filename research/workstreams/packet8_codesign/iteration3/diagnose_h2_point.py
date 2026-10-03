"""Replay H2 trajectory marks at a chosen cached actual24 point.

Counts refer to a tilted positive comparison expression, not to failing codes.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log
from pathlib import Path
from time import monotonic

import conditioned_diagnostics as cd
import capped_gate as cg


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    local, saved = cg.load_local(args.receipt)
    theta = float(Fraction(saved['tilt']))
    nu = float(Fraction(args.nu))
    if not 2 <= args.q <= 512 or nu < 0:
        raise ValueError('q in2..512 and nonnegative route tilt required')
    pins = dict(saved['source_sha256'])
    for path in (Path(__file__).resolve(), Path(cd.__file__).resolve(),
                 Path(cg.__file__).resolve(), args.receipt.resolve(),
                 Path(cg.cap_counts.__file__).resolve()):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    epsilon = 1e-4
    started = monotonic()
    families, labels = cd.make_perturbations(local, nu, epsilon=epsilon)
    moments, terminals = cd.batch_moments(families, args.q)
    counts = {label: float((moments[2+2*i]-moments[1+2*i])/(2*epsilon))
              for i, label in enumerate(labels)}
    checks = [0, 1+2*labels.index('active_0'), 2+2*labels.index('active_0')]
    errors = {str(i): float(moments[i]-cd.logarithmic_moment(families[i], args.q))
              for i in checks}
    if max(map(abs, errors.values())) > 1e-7:
        raise ArithmeticError('scaled diagnostic disagrees with all-log replay')
    witness = cg.search_threshold(args.q, 2)
    beta = cg.gate.prior.rs_uniform_envelope.UniformInputEnvelope(16,8,8,2).beta
    exponent = (log(comb(512,args.q))+args.q*(log(beta.numerator)-log(beta.denominator))
                +theta*13107-nu*witness['h']+float(moments[0]))
    if not cd.checked_pins(pins):
        raise ArithmeticError('source pins changed')
    result = dict(schema='packet8-H2-point-diagnostic-1', q=args.q, theta=theta,
        route_tilt=nu, route_witness=witness, map_record=saved['map_record'],
        counts=counts, expected_active_packets=sum(k*counts[f'active_{k}'] for k in range(9)),
        potential_packet_count=32*args.q, minimum_outer_nonzero_symbols=9*args.q,
        good_message_margin_bits=-exponent/log(2),
        source_sha256=pins, source_pins_verified_at_finish=True,
        proposal_only=True, whole_code_certificate=False, actual_failure_distribution=False,
        log_replay_errors=errors, elapsed_seconds=monotonic()-started)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k: result[k] for k in ('q','theta','good_message_margin_bits',
        'expected_active_packets','potential_packet_count','minimum_outer_nonzero_symbols','counts')},indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q',type=int,required=True)
    parser.add_argument('--nu',default='3')
    parser.add_argument('--receipt',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
