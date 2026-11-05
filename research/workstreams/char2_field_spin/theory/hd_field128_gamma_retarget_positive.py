"""Copy pinned positive proof parameters to a new original weight; no bounds."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import shutil

import hd_field128_gamma_midpoint_adaptive as adaptive


HERE = Path(__file__).resolve().parent


def retarget(source_path, expected_sha256, source_weight, weight):
    if (type(weight) is not int or not 2 <= weight < 4096
            or type(source_weight) is not int or not 2 <= source_weight < 4096):
        raise ValueError('Expected explicit original weights2..4095')
    path = Path(source_path).resolve()
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Frozen positive-literal source changed')
    source = json.loads(raw)
    construction = adaptive.provider.gamma.construction(30)
    if (source.get('active_groups') != source_weight
            or source.get('state_dimension') != 30
            or source.get('construction') != construction
            or F(source['target_distance']) != F(2, 5)):
        raise ValueError('Unexpected source weight or gamma30/d40 construction')
    old = source['witness']
    markers = [[str(F(v)) for v in pair] for pair in old['markers']]
    z, alpha = F(old['z']), F(old['alpha'])
    if (not 1 <= len(markers) <= 128
            or any(len(row) != 2 or any(F(v) <= 0 for v in row) for row in markers)
            or not 0 < z < 1 or not 0 < alpha <= 1
            or old.get('grid', 128) not in (64, 128)
            or old.get('state_dimension', 30) != 30
            or old.get('inner_width', 64) != 64
            or int(old.get('field_order', source.get('field_order', 1 << 128))) != 1 << 128):
        raise ValueError('Invalid exact positive witness or field dimensions')
    pins = adaptive.source_pins()
    for name, value in source['source_sha256'].items():
        if name in pins and pins[name] != value:
            raise ValueError('Source closure conflict: '+name)
        pins[name] = value
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    for name, value in pins.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != value:
            raise RuntimeError('Pinned helper changed: '+name)
    if path.read_bytes() != raw:
        raise RuntimeError('Source proposal changed')
    return dict(schema='gamma-original-weight-positive-literal-proposal-v2',
        construction=construction, active_groups=weight, state_dimension=30,
        field_order=1 << 128, target_distance='2/5',
        source=dict(path=str(path), sha256=expected_sha256), source_weight=source_weight,
        source_sha256=dict(sorted(pins.items())),
        witness=dict(z=str(z), alpha=str(alpha), markers=markers,
                     grid=old.get('grid', 128), state_dimension=30),
        marker_count=len(markers),
        source_role='Exact positive literals only; every numerical cap must be rebuilt at the new weight',
        old_numerical_bounds_reused=False, no_weight_interpolation=True,
        directed=False, selected_weight_certificate=False, whole_code_certificate=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--source-sha256', required=True)
    parser.add_argument('--source-weight', type=int, required=True)
    parser.add_argument('--weight', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = retarget(args.source, args.source_sha256, args.source_weight, args.weight)
    if shutil.disk_usage(args.output.resolve().parent).free < 750 << 20:
        raise RuntimeError('Less than750MiB free before proposal publication')
    pin = adaptive.provider.search.write_new(args.output, result)
    print(json.dumps(pin), flush=True)


if __name__ == '__main__':
    main()
