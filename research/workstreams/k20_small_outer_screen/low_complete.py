"""Outward q2 and consecutive low-occupancy components for K20/s22.

One fresh census serves every component. Frozen exact-shell and regional
placement helpers supply the bounds; this driver adds immutable receipts
and its source pin. These are components, never a whole-code certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
from time import monotonic

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'k16_design'))
sys.path.insert(0, str(HERE.parent / 'k18_margin60'))
import packet_inner_quadratic_extension as maps
import packet_rs_state_sparse as sparse
import packet_rs_state_ladder as ladder
import tune
import replay_fractional as replay

from flint import ctx

K = 1 << 20
MAP = 'e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581'
GEOMETRY = sparse.q1.Geometry(8192, 64, 128)
Q2_TILTS = ('.0002', '.0004', '.0006', '.0008')
TAIL_TILTS = ('.0008', '.0016', '.0032', '.006', '.01')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fraction(pair):
    if (not isinstance(pair, (list, tuple)) or len(pair) != 2 or
            any(type(value) is not int for value in pair) or pair[0] <= 0):
        raise ValueError('positive dyadic endpoint required')
    return F(pair[0]) * F(2) ** pair[1]


def tail_endpoints(record, maximum):
    expected = list(range(3, maximum + 1))
    if (record.get('occupancy_covered') != expected or
            set(record.get('occupancy_uppers', {})) != set(map(str, expected)) or
            record.get('whole_code_certificate') is not False or
            (record.get('K'), record.get('N'), record.get('threshold'), record.get('groups')) !=
            (K, 2*K, (2*K)//10, 8192) or record.get('state_bits') != 22 or
            record.get('map_record', {}).get('map_sha256') != MAP):
        raise ValueError('explicit consecutive K20/s22 component required')
    return {q: fraction(record['occupancy_uppers'][str(q)]) for q in expected}


def write_component(path, record, driver_pins):
    pins = dict(record['source_sha256'])
    for name, value in driver_pins.items():
        if name in pins and pins[name] != value:
            raise ArithmeticError('conflicting source pin')
        pins[name] = value
    replay.authenticate(pins)
    record['source_sha256'] = pins
    record['component_driver'] = str(Path(__file__).resolve())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')
    return digest(path)


def run(prefix, *, maximum=63, q2_tilts=Q2_TILTS, tail_tilts=TAIL_TILTS,
        precision=256, target_bits=65, export_if_needed=True):
    prefix = Path(prefix)
    if type(maximum) is not int or not 3 <= maximum <= 128:
        raise ValueError('bounded consecutive low range 3..maximum<=128 required')
    if type(target_bits) is not int or target_bits < 1:
        raise ValueError('positive target margin required')
    q2_tilts, _ = sparse._options(2, GEOMETRY, q2_tilts, precision, None)
    tail_tilts, _ = sparse._options(1, GEOMETRY, tail_tilts, precision, None)
    q2_path = prefix.with_name(prefix.name + '-q2.json')
    paths = [prefix.with_name(prefix.name + '-tail-t' + str(F(t)).replace('/', '_') + '.json')
             for t in tail_tilts]
    exports = [prefix.with_name(prefix.name + '-physical-t' + str(F(t)).replace('/', '_') + '.json')
               for t in ('.0032', '.006', '.01')]
    if any(path.exists() for path in (q2_path, *paths, *exports)):
        raise ValueError('all component and export paths must be fresh')
    driver_pins = {str(path.resolve()): digest(path)
                   for path in (Path(__file__), Path(tune.__file__), Path(replay.__file__))}
    started = monotonic()
    ctx.prec = precision
    print('Preparing one fresh K20/s22 census for q2 and consecutive low q', flush=True)
    data, map_record = maps.prepare(22, birth_density='capped')
    if map_record['map_sha256'] != MAP:
        raise ArithmeticError('exact certified s22 map required')
    print(f'Census prepared in {monotonic()-started:.3f}s', flush=True)
    q2_record = sparse.run_low(2, data=data, map_record=map_record, geometry=GEOMETRY,
                               tilts=q2_tilts, precision=precision)
    write_component(q2_path, q2_record, driver_pins)
    print(f'Q2 ready: {q2_path}; margin {q2_record["margin_bits"]}', flush=True)
    best, witnesses = {}, {}
    for path, tilt in zip(paths, tail_tilts):
        record = ladder.run_tail(K=K, bits=22, occupancies=range(3, maximum+1),
            tilts=[tilt], precision=precision, method='exact', placement_backend='binary',
            data=data, map_record=map_record)
        values = tail_endpoints(record, maximum)
        write_component(path, record, driver_pins)
        for q, value in values.items():
            if q not in best or value < best[q]:
                best[q], witnesses[q] = value, str(tilt)
        gaps = [q for q, value in best.items() if value >= F(1, 1 << target_bits)]
        worst = max(best, key=best.get)
        print(json.dumps(dict(tilt=tilt, receipt=str(path),
            worst_q=worst, worst_bits=replay.bits(best[worst]), gaps=gaps,
            margins={q: replay.bits(best[q]) for q in (3,4,8,16,32,48,63) if q in best},
            elapsed_seconds=monotonic()-started)), flush=True)
        if not gaps:
            break
    if gaps and export_if_needed:
        # The export helper explicitly records 192-bit arithmetic.
        ctx.prec = 192
        model = tune.ActivityModel(data, map_record, activity=F(1, 2), precision=192)
        pins = dict(model.sources)
        pins.update(driver_pins)
        pins[str(Path(tune.dense.__file__).resolve())] = digest(tune.dense.__file__)
        for path, tilt in zip(exports, (F('.0032'), F('.006'), F('.01'))):
            model.physical(tilt)
            replay.authenticate(pins)
            tune.export_physical(path, model, tilt, pins)
            replay.authenticate(pins)
    return dict(q2=str(q2_path), covered=list(range(2, maximum+1)),
                target_bits=target_bits, tail_gaps=gaps, whole_code_certificate=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-prefix', type=Path, required=True)
    parser.add_argument('--maximum', type=int, default=63)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--target-bits', type=int, default=65)
    parser.add_argument('--q2-tilts', nargs='+', default=Q2_TILTS)
    parser.add_argument('--tail-tilts', nargs='+', default=TAIL_TILTS)
    parser.add_argument('--no-export', action='store_true')
    args = parser.parse_args()
    print(json.dumps(run(args.output_prefix, maximum=args.maximum, precision=args.precision,
        target_bits=args.target_bits, q2_tilts=args.q2_tilts, tail_tilts=args.tail_tilts,
        export_if_needed=not args.no_export)), flush=True)
