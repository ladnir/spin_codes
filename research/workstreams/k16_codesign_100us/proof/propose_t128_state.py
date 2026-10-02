"""Bounded floating state-size proposals; all output is noncertificate."""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import screen_t128 as screen
from flint import ctx
from packet_outer_geometry_proposal import estimate, upper_arrays
from rs_uniform_envelope import UniformInputEnvelope


def prepare(bits):
    if not 16 <= bits <= 20:
        raise ValueError('a prefix of sixteen through twenty rows required')
    source = screen.t128.BASE_MAP
    if hashlib.sha256(source.read_bytes()).hexdigest() != screen.t128.BASE_SHA256:
        raise ArithmeticError('pinned t128/s20 map changed')
    rows = tuple(int(x, 16) for x in json.loads(source.read_bytes())['generator_rows_hex'][:bits])
    columns = screen.t128.columns_from_rows(rows)
    maps = screen.q1.kernel_t64.s16_maps
    if (maps.binary_rank(rows) != bits or maps.binary_rank(columns) != bits
            or any((a & b).bit_count() & 1 for a in rows for b in rows)
            or any(maps.binary_rank(columns[i:i+4]) != 4 for i in range(0, 128, 4))):
        raise ArithmeticError('rank, CA=0, or full-rank-packet check failed')
    images = tuple(maps.images_from_rows(rows))
    spectrum = dict(sorted(Counter(x.bit_count() for x in images).items()))
    data = screen.q1.kernel_t64.kernel_maps.prepare_maps(images, columns, bits=bits,
        distribution='uniform_gl', birth_density='capped')
    return data, dict(physical_t=128, state_bits=bits,
        expansion_rows_hex=list(map(hex, rows)), feedback_columns=columns,
        feedback_definition='C=A^T', full_state_census=True,
        expansion_spectrum={str(w): n for w, n in spectrum.items()},
        minimum_expansion_weight=min(w for w in spectrum if w),
        map_sha256=data['map_sha256'])


def local(data, tilt):
    z = (-screen.aq(Q(tilt))).exp()
    family = screen.q1.kernel_t64.sparse_kernel.outward_at_z(data, z)
    return screen.q1.kernel_t64.kernel_birth_density.refine_local(data, family, z, Q(1, 2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, required=True)
    parser.add_argument('--tilts', nargs='+', required=True)
    parser.add_argument('--max-q', type=int, default=128)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or not 3 <= args.max_q <= 512:
        raise ValueError('fresh output and q=3..512 range required')
    ctx.prec = 192
    data, record = prepare(args.bits)
    print(f'Prepared t128/s{args.bits}, minimum expansion weight={record["minimum_expansion_weight"]}', flush=True)
    best, choices, trials = {}, {}, []
    for tilt in args.tilts:
        proposal = estimate(upper_arrays(local(data, tilt)), K=65536,
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
            whole_code_certificate=False, map_record=record,
            best_margin_bits=best, tilt_choices=choices, trials=trials), indent=2) + '\n')


if __name__ == '__main__':
    main()
