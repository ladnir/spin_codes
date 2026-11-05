"""Explicit metadata projection of the actual r3072 bank's exact proposal."""

from hashlib import sha256
import json
from pathlib import Path

import hd_field128_gamma_midpoint_neighbors as neighbors


HERE = Path(__file__).resolve().parent
SOURCE = HERE/'hd_field128_gamma30_d40_r3072_midpoint_directed256_20261010/proposal.json'
OUTPUT = HERE/'hd_field128_gamma30_d40_r3072_exact_neighbor_seed_20261010.json'


def main():
    raw = SOURCE.read_bytes()
    actual = json.loads(raw)
    parent_path = Path(actual['source']['path'])
    parent_raw = parent_path.read_bytes()
    if sha256(parent_raw).hexdigest() != actual['source']['sha256']:
        raise RuntimeError('Original positive proposal changed')
    parent = json.loads(parent_raw)
    expected = dict(active_groups=3072, state_dimension=30, target_distance='2/5')
    if any(actual[key] != value or parent[key] != value for key, value in expected.items()):
        raise ValueError('Calibration metadata mismatch')
    a, p = actual['witness'], parent['witness']
    if a['markers'] != p['exact_positive_marker_proposals'] or any(a[key] != p[key] for key in ('z', 'alpha')):
        raise ValueError('Actual-bank literals differ from the original positive proposal')
    construction = neighbors.provider.gamma.construction(30)
    if p['construction'] != construction or p['grid'] != 128:
        raise ValueError('Original bank construction or grid mismatch')
    pins = neighbors.enrich.source_pins()
    for source in (actual['source_sha256'], parent['source_sha256']):
        for name, digest in source.items():
            if name in pins and pins[name] != digest:
                raise ValueError('Source closure conflict: '+name)
            pins[name] = digest
    for path in (Path(neighbors.__file__), Path(__file__)):
        pins[path.name] = sha256(path.read_bytes()).hexdigest()
    for name, digest in pins.items():
        if sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Changed source: '+name)
    if SOURCE.read_bytes() != raw or parent_path.read_bytes() != parent_raw:
        raise RuntimeError('Proposal source changed')
    result = dict(schema='gamma-r3072-actual-exact-proposal-metadata-projection-v1', **expected,
                  construction=construction, witness=dict(a, grid=p['grid']),
                  source=dict(path=str(SOURCE), sha256=sha256(raw).hexdigest()),
                  metadata_source=dict(path=str(parent_path), sha256=sha256(parent_raw).hexdigest()),
                  source_sha256=dict(sorted(pins.items())), exact_marker_count=len(a['markers']),
                  projection_rule='Exact actual-bank proposal literals; construction/grid copied from pinned original positive proposal',
                  old_numerical_bounds_reused=False, directed=False, whole_code_certificate=False, proposal_only=True)
    print(json.dumps(neighbors.provider.search.write_new(OUTPUT, result)), flush=True)


if __name__ == '__main__':
    main()
