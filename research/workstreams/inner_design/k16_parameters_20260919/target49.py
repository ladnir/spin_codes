"""Bounded Q1 feasibility screen for a 49-bit K16 target, not a certificate."""
import argparse
import json
from pathlib import Path
import time
import model


def run(args):
    if args.output.exists():
        raise FileExistsError(args.output)
    saved = json.loads(args.maps.read_text())
    subspaces = {r['s']: r['inner'] for r in saved['rows'] if r['found']}
    cases = [('selected12819', model.study.record_for(256,128))]
    cases += [(f'sub64s{s}', subspaces[s]) for s in (12,13)]
    cases += [(f't{t}s{s}', model.study.core.inner(t,s))
              for t,s in [(16,10),(32,11),(32,12),(32,13),(32,14),(32,15),(64,16),(64,19)]]
    rows = []
    for name, record in cases:
        start = time.monotonic()
        diagnostic = model.study.core.evaluate(record,256,16,sharp=True)
        row = dict(name=name, inner=record, diagnostic=diagnostic,
                   elapsed_seconds=time.monotonic()-start)
        rows.append(row)
        print(name, round(diagnostic['q1_margin_bits'],6),
              'tilt', diagnostic['dominant_tilt'], flush=True)
    result = dict(status='K16_TARGET49_BINARY64_Q1_SCREEN', target_margin_bits=49,
                  full_distance_proved=False, rows=rows, source_sha256=model.base.sources())
    result['source_sha256'][Path(__file__).resolve().relative_to(model.base.ROOT).as_posix()] = model.base.base.sha(Path(__file__))
    result['source_sha256'][args.maps.resolve().relative_to(model.base.ROOT).as_posix()] = model.base.base.sha(args.maps)
    model.base.base.write_new(args.output,result)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maps',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args())
