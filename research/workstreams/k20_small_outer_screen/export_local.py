"""Export fresh certified-s22 physical bounds at two tilts from one census.

Physical operators are independent of global K. The frozen physical-export
schema is retained so the authenticated outward replay can consume them.
No regional sweep or whole-code certificate is produced here.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from pathlib import Path
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'k18_margin60'))
import tune

from flint import ctx


def authenticate(pins):
    for path, expected in pins.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
            raise ArithmeticError(f'source changed during physical export: {path}')


def run(output_prefix):
    output_prefix = Path(output_prefix)
    paths = [output_prefix.with_name(output_prefix.name + suffix)
             for suffix in ('-t002.json', '-t004.json')]
    if any(path.exists() for path in paths):
        raise ValueError('both physical output paths must be fresh')
    ctx.prec = 192
    started = monotonic()
    print('Preparing one fresh certified s22 census for physical exports', flush=True)
    data, map_record = tune.certified.prepare(22, birth_density='capped')
    if map_record['map_sha256'] != tune.MAP_SHA256:
        raise ArithmeticError('exact certified s22 map required')
    tune.certified.authenticate(data, map_record)
    print(f'Fresh census prepared in {monotonic()-started:.3f}s', flush=True)
    model = tune.ActivityModel(data, map_record, activity=F(1, 2), precision=192)
    pins = dict(model.sources)
    for path in (Path(__file__), Path(tune.__file__), Path(tune.dense.__file__)):
        pins[str(path.resolve())] = hashlib.sha256(path.read_bytes()).hexdigest()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for path, tilt in zip(paths, (F(1, 50), F(1, 25))):
        model.physical(tilt)
        authenticate(pins)
        tune.export_physical(path, model, tilt, pins)
        authenticate(pins)
        print(f'Tilt {tilt}: elapsed {monotonic()-started:.3f}s; '
              f'SHA256 {hashlib.sha256(path.read_bytes()).hexdigest()}', flush=True)
    return paths


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-prefix', type=Path, required=True)
    args = parser.parse_args()
    run(args.output_prefix)
