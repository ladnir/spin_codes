"""Run the explicit larger-outer conditioned-iid component, not a benchmark."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

from flint import arb
from packet_outer_cost_search import prepare
from packet_outer_cost_dense import run
from rs_uniform_envelope import UniformInputEnvelope


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, required=True)
    parser.add_argument('--K', type=int, default=1 << 20)
    parser.add_argument('--outer', choices=('256', '512'), default='512')
    parser.add_argument('--q-min', type=int, default=1)
    parser.add_argument('--q-max', type=int)
    parser.add_argument('--tilts', default='.0008,.0032,.0128,.0512,.1024,.2048,.4096,.8192,1.6384,2.1972246')
    parser.add_argument('--markers', default='1/1024,1/256,1/64,1/16,1/4,1/2,3/4,1')
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    if output.exists():
        raise ValueError('output must be fresh')
    envelope = UniformInputEnvelope(16, 8, 4, 8) if args.outer == '256' else UniformInputEnvelope(32, 16, 4, 8)
    last = args.q_max if args.q_max is not None else args.K//envelope.message_bits

    def progress(tilt, best, choices, elapsed):
        passing = [q for q, value in best.items() if value is not None and value < arb(2)**-40]
        passed = set(passing)
        holes = [q for q in best if q not in passed]
        print(json.dumps(dict(tilt=str(tilt), elapsed=elapsed, passing=len(passing),
            unresolved_range=[min(holes), max(holes)] if holes else None,
            selected_margins={str(q): str(-best[q].log()/arb(2).log()) for q in
                (1, 2, 32, 64, 96, 128, 160, 192, 256, 384, 512, last) if q in best})), flush=True)

    print(f'fresh census s={args.bits}', flush=True)
    data, record = prepare(args.bits)
    result = run(data, record, K=args.K, envelope=envelope,
        occupancies=range(args.q_min, last+1), tilts=map(Q, args.tilts.split(',')),
        marker_probabilities=map(Q, args.markers.split(',')),
        precision=args.precision, progress=progress)
    with output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')


if __name__ == '__main__':
    main()
