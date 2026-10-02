"""Small command-line proof-cost screen; receipts are explicitly noncertificates."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
from time import monotonic

from packet_outer_cost_search import Model, prepare
from rs_uniform_envelope import UniformInputEnvelope


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, required=True)
    parser.add_argument('--K', type=int, default=1 << 20)
    parser.add_argument('--outer', choices=['256', '512', 'both'], default='both')
    parser.add_argument('--occupancies', default='1,2,32,64,96,128,160,192,256')
    parser.add_argument('--tilts', default='1/1024,1/256,1/128,1/64,1/32,3/64,1/16,3/32,1/8,3/16,1/4')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    output = Path(args.output)
    if output.exists():
        raise ValueError('output must be fresh')
    qs = list(map(int, args.occupancies.split(',')))
    tilts = list(map(Q, args.tilts.split(',')))
    started = monotonic()
    print(f'fresh census s={args.bits}', flush=True)
    model = Model(*prepare(args.bits))
    print(f'census prepared in {monotonic()-started:.2f}s', flush=True)
    outers = {256: UniformInputEnvelope(16, 8, 4, 8), 512: UniformInputEnvelope(32, 16, 4, 8)}
    if args.outer != 'both':
        outers = {int(args.outer): outers[int(args.outer)]}
    records = []
    best = {name: {} for name in outers}
    for tilt in tilts:
        for name, envelope in outers.items():
            record = model.estimate(K=args.K, envelope=envelope, occupancies=qs, tilt=tilt, backend='scaled-or-log')
            records.append(record)
            for q, witness in record['witnesses'].items():
                value = witness['estimated_margin_bits']
                if q not in best[name] or value > best[name][q]['margin']:
                    best[name][q] = dict(margin=value, tilt=str(tilt))
            print(json.dumps(dict(elapsed=monotonic()-started, bits=args.bits,
                outer=name, tilt=str(tilt), backend=record['backend'], margins={q: round(w['estimated_margin_bits'], 4)
                for q, w in record['witnesses'].items()})), flush=True)
    receipt = dict(proposal_only=True, whole_code_certificate=False,
        K=args.K, bits=args.bits, best=best, trials=records, elapsed=monotonic()-started)
    with output.open('x') as handle:
        json.dump(receipt, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(best=best, elapsed=receipt['elapsed'])), flush=True)


if __name__ == '__main__':
    main()
