"""Bind the fresh serial half-rate timings to sources, maps, and full proofs."""
import argparse
import json
from pathlib import Path
import re
import statistics

import ladder
import verify_ladder
import verify_full

model = ladder.model
HERE = Path(__file__).resolve().parent
IMPLEMENTATION = HERE.parent / 'asymmetric/bch256/weight5/implementation'


def hashes(path):
    result = {}
    for line in path.read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        assert name not in result
        result[name] = digest
    return result


def verify(directory):
    proof18 = verify_ladder.verify(18, HERE / 'Q1_LOWER.json', HERE / 'SPARSE_M18_v3.json',
                                   HERE / 'DENSE_M18_short_v3.json')
    proof20 = verify_full.verify()
    manifest = model.base.read(IMPLEMENTATION / 'IMPLEMENTATION.json')
    model.authenticate(manifest)
    assert manifest['instance'] == proof20['instance']
    for name, digest in manifest['generated_sha256'].items():
        assert model.base.sha(IMPLEMENTATION / name) == digest, name
    maps = re.findall(r'columns\{([^}]+)\};', (IMPLEMENTATION / 'AsymmetricMap.h').read_text())
    parsed = [[int(c, 0) for c in row.split(',')] for row in maps]
    assert parsed == [proof20['instance']['inner']['expansion_columns'], proof20['instance']['inner']['feedback_columns']]
    assert proof18['instance']['inner'] == proof20['instance']['inner']
    assert proof18['instance']['outer_manifest_sha256'] == proof20['instance']['outer_manifest_sha256']
    measured_sources = hashes(directory / 'sources.sha256')
    for name, digest in measured_sources.items():
        assert model.base.sha(model.ROOT / name) == digest, name
    binaries = hashes(directory / 'binaries.sha256')
    old_binaries = hashes(IMPLEMENTATION / 'measurements/binaries.sha256')
    variants = ('baseline', 'baseline_tuned', 'sparse_pages')
    for name in variants:
        assert binaries[f'build/{name}_bench'] == old_binaries[f'build/{name}_bench']
        assert binaries[f'build/{name}_test'] == old_binaries[f'build/{name}_test']
        flags = (directory / f'flags-{name}.txt').read_text()
        assert '-O3' in flags and '-march=znver4' in flags
    assert (directory / 'measure_half.sh').read_bytes() == (HERE / 'measure_half.sh').read_bytes()
    log = (directory / 'correctness.log').read_text()
    assert '100% tests passed, 0 tests failed out of 12' in log
    assert all(f'{name}_test' in log for name in variants)
    cells = []
    for exponent in (16, 18, 20):
        summaries = {}
        for name in variants:
            rows = [model.base.read(directory / f'{name}-m{exponent}-r{i}.jsonl') for i in (1, 2, 3)]
            expected = 't128_s19_weight5_seed0_r1' if name == 'sparse_pages' else 't128_s19'
            for row in rows:
                assert row['configuration'] == expected
                assert (row['m'], row['outer_length'], row['outer_dimension'], row['trials'],
                        row['layout'], row['inplace']) == (exponent, 256, 128, 101, 'packed24', True)
                assert row['tile_rows'] == (2048 if exponent == 20 else 256)
                assert 0 < row['p10_ms'] <= row['median_ms'] <= row['p90_ms']
            for key in ('output_hash', 'retained_setup_bytes', 'workspace_bytes'):
                assert len({r[key] for r in rows}) == 1
            medians = [r['median_ms'] for r in rows]
            summaries[name] = dict(median_ms=statistics.median(medians), process_medians_ms=medians,
                **{key: rows[0][key] for key in ('output_hash', 'retained_setup_bytes', 'workspace_bytes')})
        assert summaries['baseline']['output_hash'] == summaries['baseline_tuned']['output_hash']
        assert len({x['workspace_bytes'] for x in summaries.values()}) == 1
        cells.append(dict(message_exponent=exponent, summaries=summaries,
            distance_certificate_available=exponent in (18, 20),
            certified_margin_bits={18: proof18['margin_bits'], 20: proof20['margin_bits']}.get(exponent),
            reduction_vs_baseline_percent=100*(1-summaries['sparse_pages']['median_ms']/summaries['baseline']['median_ms'])))
    paths = [Path(__file__).resolve(), HERE / 'measure_half.sh', IMPLEMENTATION / 'IMPLEMENTATION.json',
             IMPLEMENTATION / 'measurements/binaries.sha256'] + list(directory.iterdir())
    sources = {**proof18['source_sha256'], **proof20['source_sha256'], **measured_sources}
    sources.update({p.relative_to(model.ROOT).as_posix(): model.base.sha(p) for p in paths if p.is_file()})
    return dict(status='VERIFIED_IMT_HALF_RATE_SERIAL_TIMINGS', selected='sparse_pages',
        cells=cells, new_serial_measurements=True, exact_map_columns_checked=True,
        matched_existing_binaries=True, correctness_tests=12,
        policy='Peach CPU 15; 3 process medians per cell, 101 in-place calls after 3 warmups; no reset/copy; setup and workspace preparation excluded.',
        full_certificate_exponents=[18, 20], source_sha256=sources)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, default=HERE / 'measurements/half_20260916')
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    result = verify(a.directory.resolve())
    if a.output:
        model.base.write_new(a.output.resolve(), result)
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}, indent=2))
