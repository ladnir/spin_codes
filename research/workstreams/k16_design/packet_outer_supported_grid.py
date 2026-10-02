"""Proposal driver for the support-checked regional recurrence."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
from time import monotonic

from packet_outer_cost_search import Model, prepare
from packet_regional_supported import estimate
from rs_uniform_envelope import UniformInputEnvelope


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, required=True)
    parser.add_argument('--K', type=int, default=1 << 20)
    parser.add_argument('--outer', choices=('256', '512'), default='512')
    parser.add_argument('--occupancies', help='Explicit comma-separated q list, else use q-min..q-max')
    parser.add_argument('--q-min', type=int, default=1)
    parser.add_argument('--q-max', type=int, default=224)
    parser.add_argument('--tilts', required=True)
    parser.add_argument('--banded', action='store_true',
        help='At each tilt evaluate only q between0.35*lambda*L and0.9*lambda*L')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    if output.exists():
        raise ValueError('output must be fresh')
    qs = list(map(int, args.occupancies.split(','))) if args.occupancies else list(range(args.q_min, args.q_max+1))
    tilts = list(map(Q, args.tilts.split(',')))
    env = UniformInputEnvelope(16, 8, 4, 8) if args.outer == '256' else UniformInputEnvelope(32, 16, 4, 8)
    start = monotonic()
    print(f'fresh census s={args.bits}', flush=True)
    model = Model(*prepare(args.bits))
    best, trials = {}, []
    for tilt in tilts:
        scale = tilt*(args.K//env.message_bits)
        trial_qs = [q for q in qs if Q(7, 20)*scale <= q <= Q(9, 10)*scale] if args.banded else qs
        if not trial_qs:
            continue
        result = estimate(model.local(tilt), K=args.K, envelope=env, occupancies=trial_qs, tilt=tilt)
        model.authenticate()
        trials.append(result)
        for q, value in result['witnesses'].items():
            if q not in best or value['estimated_margin_bits'] > best[q]['margin']:
                best[q] = dict(margin=value['estimated_margin_bits'], tilt=str(tilt))
        weakest = sorted(best, key=lambda q: best[q]['margin'])[:8]
        print(json.dumps(dict(elapsed=monotonic()-start, tilt=str(tilt), backend=result['backend'],
            worst={q: best[q] for q in weakest}, passing=sum(v['margin'] > 40 for v in best.values()))), flush=True)
    receipt = dict(proposal_only=True, whole_code_certificate=False,
        K=args.K, state_bits=args.bits, envelope=env.metadata(), map_record=model.map_record,
        source_sha256=model.sources, best=best, trials=trials, elapsed=monotonic()-start)
    if set(best) != {str(q) for q in qs}:
        raise ValueError('the tilt bands did not cover every requested occupancy')
    with output.open('x') as handle:
        json.dump(receipt, handle, indent=2)
        handle.write('\n')


if __name__ == '__main__':
    main()
