"""Separate smaller state, shorter step, and map choices in the length study."""
import argparse
from pathlib import Path

import adaptive_step as core
import adaptive_length_study as study


def run(output):
    if output.exists(): raise FileExistsError('Use a fresh output')
    cells = []
    for b,m in ((128,12),(128,16),(256,16),(256,18)):
        for t,s in ((32,15),(64,15),(64,19 if b==256 else 20)):
            row = core.evaluate(core.inner(t,s),b,m,sharp=True)
            cells.append(row)
            print(row,flush=True)
    core.grid.ladder.model.base.write_new(output,dict(status='BINARY64_IMT_ADAPTIVE_STATE_CONTROLS',
        cells=cells,full_distance_proved=False,source_sha256={**study.sources(),
            Path(__file__).resolve().relative_to(core.grid.ROOT).as_posix():core.grid.maps.sha(Path(__file__))}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
