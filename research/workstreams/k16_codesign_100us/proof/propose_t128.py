"""Floating tilt proposals for the explicit t128/s16 screen; no certificate."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as Q

import screen_t128 as screen
from flint import ctx
from packet_outer_geometry_proposal import estimate, upper_arrays
from rs_uniform_envelope import UniformInputEnvelope


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', choices=('prefix16', 'repeat64'), required=True)
    parser.add_argument('--tilts', nargs='+', required=True)
    parser.add_argument('--max-q', type=int, default=64)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or not 3 <= args.max_q <= 512:
        raise ValueError('fresh output and q=3..512 range required')
    ctx.prec = 192
    data, record = screen.prepare(args.candidate)
    best, choices, trials = {}, {}, []
    for tilt in args.tilts:
        operators = screen.local_operators(data, Q(tilt))
        proposal = estimate(upper_arrays(operators), K=65536,
            envelope=UniformInputEnvelope(16, 8, 4, 4),
            occupancies=tuple(range(3, args.max_q + 1)), tilt=tilt)
        for q, witness in proposal['witnesses'].items():
            margin = witness['estimated_margin_bits']
            if q not in best or margin > best[q]:
                best[q], choices[q] = margin, tilt
        trials.append(proposal)
        worst = sorted(best, key=best.get)[:8]
        print(f'tilt={tilt} worst={[(q, round(best[q], 4)) for q in worst]}', flush=True)
        args.output.write_text(json.dumps(dict(proposal_only=True,
            whole_code_certificate=False, candidate=args.candidate, map_record=record,
            best_margin_bits=best, tilt_choices=choices, trials=trials), indent=2) + '\n')


if __name__ == '__main__':
    main()
