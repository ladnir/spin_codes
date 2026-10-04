"""Outward routing-conditioned bounds for every q in a requested interval.

Uses the fixed certified s22 physical operators at K20. This component is
not a complete certificate. Fractional exponents and tilts are proof witnesses,
not code parameters. Each recorded endpoint already includes the group union.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from pathlib import Path
import json
from time import monotonic

from flint import ctx
import outward_screen as screen


def run(local, output, minimum, maximum, alpha):
    output = Path(output)
    alpha = F(alpha)
    if output.exists() or not 3 <= minimum <= maximum <= 8192 or not 0 < alpha <= 1:
        raise ValueError('fresh output, q interval in3..8192, and alpha in(0,1] required')
    ctx.prec = 256
    physical, tilt, source = screen.replay.load_local(local)
    pins = dict(source['source_sha256'])
    for path in (Path(local), Path(__file__), Path(screen.__file__),
                 Path(screen.replay.__file__), Path(screen.replay.placement.__file__)):
        pins[str(path.resolve())] = screen.replay.digest(path)
    start = monotonic()
    blocks = screen.replay.fractional_blocks(physical, alpha)
    print(f'K20 range q{minimum}..{maximum}: tilt={tilt} alpha={alpha}; local ready', flush=True)
    regional = screen.replay.placement.placement_power(blocks, epochs=256,
        windows=32, maximum_groups=maximum)
    values = {str(q): screen.replay.endpoint(screen.scalar(regional[q], q, tilt, alpha, 1 << 20))
              for q in range(minimum, maximum+1)}
    margins = {q: screen.replay.bits(screen.replay.fraction(v)) for q, v in values.items()}
    good = [int(q) for q, value in margins.items() if value >= 65]
    screen.replay.authenticate(pins)
    result = dict(schema='k20-small-outer-fractional-range-1', K=1048576,
        N=2097152, threshold=209715, distance='1/10', state_bits=22,
        group_dimension=128, group_output_bits=256, groups=8192,
        regions=64, slots_per_region=8192, physical_t=64, macro_slots=32,
        macros_per_region=256, physical_steps_per_macro=2,
        zero_initial_state=True, final_flush=False, continuous_state=True,
        map_record=source['map_record'], precision=256,
        local_receipt=str(Path(local).resolve()), local_sha256=screen.replay.digest(local),
        occupancy_covered=list(range(minimum, maximum+1)), occupancy_uppers=values,
        tilt=str(tilt), alpha=str(alpha), whole_code_certificate=False,
        source_sha256=pins, elapsed_seconds=monotonic()-start,
        scope='Outward setup-failure contribution for every listed active-group occupancy; '
              'ideal independent setup ensemble, not a whole-code certificate.')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(output=str(output), good65_range=[min(good),max(good)] if good else [],
        good65_count=len(good), minimum_margin=min(margins.values()),
        elapsed_seconds=result['elapsed_seconds'])), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--local', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--min-q', type=int, default=3)
    parser.add_argument('--max-q', type=int, required=True)
    parser.add_argument('--alpha', type=F, default=F(4,5))
    args = parser.parse_args()
    run(args.local, args.output, args.min_q, args.max_q, args.alpha)
