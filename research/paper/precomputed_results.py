"""Accepted optimized precomputed campaign; never runs benchmarks.

The convenience-library rerun is a separate measurement, not a replacement
for this retained campaign. See PRECOMPUTED_PERFORMANCE.md for provenance.
"""
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import imt_results

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / 'workstreams/spin_optimized/measurements/forward_schedule'
RAW_PIN = 'fe4afcc33fec10d1a802500394991ad6091b03fef646bbd0f4f62f68c540b087'
R2_PIN = ('k16_parameters_20260919/measurements/target49_sub12_r2_full.json',
          'dd171397b8ba21860aea9566f2c155b6410947d195d76eea45747d1cfb355c15')
INPLACE18_PINS = (
    'd4594ea42c45912e95e54e71b4c4004f356edd3e0e6c6228f56efe0d949e4982',
    'c0c585f4147f9ecb3b75151c0e702003682214a5ea5a4970d03c6f83995cae5a',
    '74db2bc62c71298feeb702da0f0d1271997a7921a99db1f9ee6e4fe592898536',
)


def inplace18():
    """Authenticate the separately reported in-place timing and derived rate."""
    require = imt_results.require
    medians = []
    for run, pin in enumerate(INPLACE18_PINS, 1):
        raw = (DIRECTORY / f'confirm-m18-inplace-r{run}.json').read_bytes()
        require(hashlib.sha256(raw).hexdigest() == pin, 'In-place timing receipt changed')
        row = json.loads(raw)
        require(row['m'] == 18 and row['direction'] == 'inplace'
                and row['configuration'] == 't128_s19_weight5_seed0_r1', 'Wrong in-place map')
        samples = row['samples_ms']
        require(len(samples) == 101 and all(math.isfinite(x) and x > 0 for x in samples),
                'Invalid in-place samples')
        require(statistics.median(samples) == row['median_ms'], 'Incorrect in-place median')
        medians.append(row['median_ms'])
    latency = statistics.median(medians)
    return latency, 2**18 / (latency * 1000)  # million output blocks / second


def load():
    require = imt_results.require
    paths = sorted(p for p in DIRECTORY.glob('confirm-m*-*.json') if '-inplace-' not in p.name)
    require(len(paths) == 18, 'Expected eighteen retained precomputed runs')
    digest = hashlib.sha256()
    groups = defaultdict(list)
    for path in paths:
        raw = path.read_bytes()
        digest.update(path.name.encode() + b'\n' + raw)
        row = json.loads(raw)
        key = row['m'], row['direction']
        require(key[0] in (16, 18, 20) and key[1] in ('forward', 'transpose'), 'Wrong timing cell')
        samples = row['samples_ms']
        require(len(samples) == 101 and all(math.isfinite(x) and x > 0 for x in samples), 'Invalid samples')
        require(statistics.median(samples) == row['median_ms'], 'Incorrect process median')
        expected = ('imt_t64_s12_subspace_v1_r2' if key[0] == 16 else
                    'imt_t128_s19_weight5' if key[1] == 'forward' else 't128_s19_weight5_seed0_r1')
        require(row['configuration'] == expected, 'Wrong precomputed configuration')
        groups[key].append(row)
    require(digest.hexdigest() == RAW_PIN, 'Precomputed timing receipt changed')
    require(len(groups) == 6 and all(len(v) == 3 for v in groups.values()), 'Incomplete timing grid')
    result = {}
    for key, rows in groups.items():
        require(len({(r['tile_rows'], r['setup_bytes'], r['workspace_bytes']) for r in rows}) == 1,
                'Inconsistent measured setup')
        result[key] = dict(median_ms=statistics.median(r['median_ms'] for r in rows),
                           process_medians=[r['median_ms'] for r in rows],
                           setup_bytes=rows[0]['setup_bytes'], workspace_bytes=rows[0]['workspace_bytes'])
    # Authenticate the additional K16 point separately from the paper's original
    # one-round ladder; do not transfer that ladder's certificate to a new map.
    r2 = imt_results.authenticate(R2_PIN, {})
    require(r2['status'] == 'VERIFIED_K16_TWO_ROUND_FULL_DISTANCE_BOUND'
            and r2['full_distance_proved'] and r2['covered_occupancies'] == [1, 512], 'Incomplete R2 certificate')
    instance = r2['instance']
    require((instance['message_bits'], instance['output_bits'], instance['cutoff']) == (65536, 131072, 13107),
            'Wrong R2 geometry')
    inner = instance['inner']
    require((inner['t'], inner['s'], inner['transvection_rounds']) == (64, 12, 2), 'Wrong R2 inner')
    bound = imt_results.exact_union(r2['component_upper'], r2['union_upper'], 49)
    require(abs(bound - r2['margin_bits']) < 1e-8, 'Wrong R2 margin')
    source = (ROOT.parent / 'spin/src/kernels/Map64S12.h').read_text()
    for field, name in [('expansion_columns', 'columns'), ('feedback_columns', 'feedbackColumns')]:
        match = re.search(r'\b' + name + r'\{([^}]+)\}', source)
        require(match is not None, 'Missing selected map')
        require([int(x, 0) for x in match[1].split(',')] == inner[field], 'R2 implementation map changed')
    return result


if __name__ == '__main__':
    for key, cell in sorted(load().items()):
        print(key, cell['median_ms'])
    latency, rate = inplace18()
    paper = Path(__file__).resolve().parent
    for name in ('applications.tex', 'implementation_appendix.tex'):
        prose = (paper / name).read_text()
        imt_results.require(f'{latency:.3f}' in prose and f'{rate:.0f}' in prose,
                            'Stale in-place latency or throughput')
    print('K18 in-place:', latency, 'ms;', rate, 'million output blocks/s (encoding only)')
