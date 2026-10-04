"""Fresh outward q1 screen for the K20 small-outer, certified-s22 ensemble.

This bounds occupancy one only, including the union over all 8192 groups.
It reuses the exact expected RS16 shell law, not a K18 numerical endpoint.
The actual 2^22-state census is regenerated once per invocation.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

K16 = Path(__file__).resolve().parents[1] / 'k16_design'
sys.path.insert(0, str(K16))

from flint import ctx
import packet_inner_quadratic_extension as maps
import packet_rs_state_sparse as sparse


MAP_SHA256 = 'e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581'
GEOMETRY = sparse.q1.Geometry(8192, 64, 128)
DEFAULT_TILTS = ('.0002', '.0004', '.0006', '.0008')


def run(output, *, tilts=DEFAULT_TILTS, precision=192):
    """Certify only the fresh q1 first-moment bound at the selected tilts."""
    tilts, output = sparse._options(1, GEOMETRY, tilts, precision, output)
    ctx.prec = precision
    _, _, identity = maps.construction(22)
    if identity['map_sha256'] != MAP_SHA256:
        raise ArithmeticError('the selected certified K18 s22 map changed')
    print('Preparing fresh certified s22 census for K20 q1', flush=True)
    data, record = maps.prepare(22, birth_density='capped')
    maps.authenticate(data, record)
    return sparse.run_low(1, data=data, map_record=record, geometry=GEOMETRY,
                          tilts=tilts, precision=precision, output=output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--tilts', nargs='+', default=DEFAULT_TILTS)
    parser.add_argument('--precision', type=int, default=192)
    args = parser.parse_args()
    run(args.output, tilts=args.tilts, precision=args.precision)


if __name__ == '__main__':
    main()
