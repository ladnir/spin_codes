"""Fresh q1 bounds and floating tail proposals for explicit monomial maps."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
from time import monotonic

import monomial_maps as maps
import screen_t128 as screen
from flint import arb, ctx
from packet_outer_geometry_proposal import estimate, upper_arrays
from rs_uniform_envelope import UniformInputEnvelope


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--map', choices=('minrank', 'bipartite', 'prism', 'lex', 'disjoint'), required=True)
    parser.add_argument('--q1-tilts', nargs='+', default=['.00256', '.00512', '.01024', '.02048'])
    parser.add_argument('--tail-tilts', nargs='*', default=[])
    parser.add_argument('--max-q', type=int, default=512)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or not 3 <= args.max_q <= 512:
        raise ValueError('fresh output and bounded q range required')
    start = monotonic()
    ctx.prec = 192
    if args.map == 'disjoint':
        import disjoint_pair_maps
        data, record = disjoint_pair_maps.prepare()
    else:
        data, record = maps.prepare(args.map)
    _, counts = screen.tail.uniform_envelope(16,8,4)
    best_one = [arb(1)] * 65
    local_cache = {}
    for tilt in args.q1_tilts:
        local = screen.q1.kernel_t64.local_operators(data, Q(tilt), Q(1,2))
        local_cache[Q(tilt)] = local
        regional = screen.q1.placement(local, epochs=16, windows=32,
            rounding=screen.q1.rounded, maximum_groups=1)
        factor = (screen.aq(Q(tilt)) * 13107).exp()
        moments = screen.q1.support_moments(regional[0], regional[1], 64)
        for v,m in enumerate(moments):
            best_one[v] = min(best_one[v], screen.up(factor*m))
        upper = screen.up(512*sum((screen.aq(c)*m for c,m in zip(counts,best_one)), arb(0)))
        print(f'{args.map} tilt={tilt} q1={-upper.log()/arb(2).log()}', flush=True)
    best, choices, trials = {}, {}, []
    result = dict(schema='monomial-t64-screen-1', map_record=record,
        q1_upper=screen.q1.endpoint(upper), q1_margin_bits=str(-upper.log()/arb(2).log()),
        q1_fresh_outward=True, q1_precision=192, q1_tilts=args.q1_tilts,
        tail_proposal_only=True, whole_code_certificate=False,
        K=65536, N=131072, threshold=13107, distance='1/10',
        zero_initial_state=True, final_flush=False, continuous_state=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    for tilt in args.tail_tilts:
        local = local_cache.get(Q(tilt))
        if local is None:
            local = screen.q1.kernel_t64.local_operators(data, Q(tilt), Q(1,2))
        proposal = estimate(upper_arrays(local), K=65536,
            envelope=UniformInputEnvelope(16,8,4,4),
            occupancies=tuple(range(3,args.max_q+1)), tilt=tilt)
        for q,w in proposal['witnesses'].items():
            margin=w['estimated_margin_bits']
            if q not in best or margin>best[q]:
                best[q],choices[q]=margin,tilt
        trials.append(proposal)
        worst=sorted(best,key=best.get)[:6]
        print(f'{args.map} tilt={tilt} tailworst={[(q,round(best[q],4)) for q in worst]}',flush=True)
        result.update(tail_estimated_margin_bits=best, tail_tilt_choices=choices,
            tail_trials=trials, elapsed_seconds=monotonic()-start)
        args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
