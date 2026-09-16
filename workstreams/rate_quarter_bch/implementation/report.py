"""Validate the retained timing receipts and rebuild PERFORMANCE.{json,md}."""
import hashlib
import json
from pathlib import Path
from statistics import median
import sys

HERE = Path(__file__).resolve().parent


def main():
    receipts = Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'receipts'
    rows = [json.loads(line) for line in (receipts/'final.jsonl').read_text().splitlines() if line.strip()]
    assert len(rows) == 12
    for row in rows:
        assert row['trials']==101 and row['inplace'] and row['layout']=='packed24'
        assert row['configuration']=='t128_s19'
    summaries = []
    for length, dimension, exponent, tile in ((128,32,16,256),(128,32,18,2048),(128,32,20,4096),(256,128,20,2048)):
        group = [r for r in rows if r['outer_length']==length and r['m']==exponent]
        assert len(group)==3 and all(r['outer_dimension']==dimension and r['tile_rows']==tile for r in group)
        assert len({r['output_hash'] for r in group})==1
        item = dict(outer_length=length,outer_dimension=dimension,m=exponent,tile_rows=tile,
                    median_ms=median(r['median_ms'] for r in group),
                    run_medians_ms=[r['median_ms'] for r in group],
                    setup_median_ms=median(r['setup_ms'] for r in group),
                    retained_setup_bytes=group[0]['retained_setup_bytes'],workspace_bytes=group[0]['workspace_bytes'],
                    inplace_buffer_bytes=(1<<exponent)*(length//dimension)*16)
        summaries.append(item)
    manifest = HERE/'generated/MANIFEST.json'
    manifest_hash = hashlib.sha256(manifest.read_bytes()).hexdigest()
    assert manifest_hash in (receipts/'binary_hashes.txt').read_text()
    identity = json.loads(manifest.read_text())
    root = HERE.parents[2]
    for name, expected in identity['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected, name
    for name, expected in identity['generated_sha256'].items():
        assert hashlib.sha256((HERE/'generated'/name).read_bytes()).hexdigest()==expected, name
    for name in ('correctness.log','sanitizer.log'):
        assert '100% tests passed, 0 tests failed out of 2' in (receipts/name).read_text()
    payload = dict(date='2026-09-11',cpu='AMD Ryzen 9 7950X',cpu_affinity=15,
                   compiler='GCC 15.2.0',flags='-O3 -DNDEBUG -std=c++20 -march=znver4 -mavx2 -mpclmul -mvpclmulqdq',
                   method='Median of three run medians; 101 online encodes per run, three warmups; all timing cells sequential',
                   scope='128-way bitsliced complete transposed encoder; in-place; setup and allocation excluded; no input copy',
                   implementation_manifest_sha256=manifest_hash,results=summaries,
                   receipt_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(receipts.iterdir()) if p.is_file()})
    (HERE/'PERFORMANCE.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')
    q=summaries[2]; half=summaries[3]
    report = ['# Quarter-rate transposed encoder performance', '',
              f'The fixed BCH-derived [128,32,32] outer with RM2Sub `(t,s)=(128,19)` takes **{q["median_ms"]:.3f} ms** at `K=2^20`.',
              'This is the complete in-place transposed encoder, not an isolated inner kernel.',
              'Each 128-bit block represents 128 parallel binary instances. The operation reads `N=4K` blocks and replaces the first `K` blocks.', '',
              '| Message size K | Encoded size N | In-place time | Three run medians (ms) | Tile rows |',
              '|---:|---:|---:|---|---:|']
    for r in summaries[:3]:
        runs=', '.join(f'{x:.3f}' for x in r['run_medians_ms'])
        report.append(f'| 2^{r["m"]} | 2^{r["m"]+2} | {r["median_ms"]:.3f} ms | {runs} | {r["tile_rows"]} |')
    report += ['',f'The same build measures the existing rate-half [256,128] / t128_s19 encoder at **{half["median_ms"]:.3f} ms** for `K=2^20`.',
               f'The quarter-rate encoder takes **{q["median_ms"]/half["median_ms"]:.2f}x** as long at the same message size, while processing twice as many encoded coordinates.',
               'The rate-half measurement is a throughput reference, not an equal-distance comparison.', '',
               '## Distance certificates', '',
               'The implemented outer spans exactly the certified fixed [128,32,32] code. Its systematic message basis does not change its spectrum.',
               'The inner reuses the exact selected columns, field multiplication, and optimized kernels from the existing t128_s19 implementation.',
               'For the setup distribution in [the outward certificate](../SMALLER_OUTWARD_CERTIFICATE.md), the two operating points at `K=2^20` are:', '',
               '| Relative distance target | Certified setup-failure probability | Margin diagnostic |',
               '|---:|---:|---:|', '| 16.5% | < 2^-40 | 41.083485 bits |', '| 19% | < 2^-30 | 30.052513 bits |', '',
               'Both rows use the same encoder and have the same runtime. The performance rows at smaller K do not add new distance certificates.', '',
               '## Optimizations and memory', '',
               'The generated outer uses 664 XORs per row, versus 1,568 for dense evaluation. AVX2 evaluates two independent rows together.',
               'The implementation retains the shared unrolled inner, optimized nibble grouping, pruned zeta transform, and two-stage routing.',
               'There are no online allocations or runtime callback wrappers. Outer selection dispatches once before the specialized hot loop.', '',
               'A bounded sweep compared packed24 and indices32 routing, with 128 through 32,768 rows per tile at K=2^20.',
               'A follow-up in-place sweep selected packed24 with 4,096 rows. Representative screening medians were:', '',
               '| K=2^20, in-place | Time |', '|---|---:|',
               '| Packed24, 2,048 rows | 19.343 ms |', '| Packed24, 4,096 rows | 18.221 ms |',
               '| Packed24, 8,192 rows | 19.100 ms |', '| Indices32, 4,096 rows | 19.587 ms |', '',
               'These screening runs used 31 trials. The headline result uses the separate three-run measurement above.',
               'Smaller tile sweeps selected 256 rows at K=2^16 and 2,048 rows at K=2^18. These choices are now the defaults.', '',
               f'At K=2^20, the in-place buffer occupies {q["inplace_buffer_bytes"]/2**20:.0f} MiB; reusable workspace occupies {q["workspace_bytes"]/2**20:.0f} MiB.',
               f'Retained setup occupies {q["retained_setup_bytes"]/2**20:.3f} MiB after discarding diagnostics and unused routing.',
               f'Setup itself took a median {q["setup_median_ms"]:.1f} ms, outside the online interval.', '',
               '## Measurement and verification', '',
               'Measurements used Peach, Ryzen 9 7950X, one thread pinned to CPU 15, GCC 15.2.0, and the release flags in PERFORMANCE.json.',
               'The recorded machine state had frequency boost disabled. Runs share the existing benchmark lock and reject other visible benchmark executables.',
               'Each timing run contains three warmups and 101 measured calls. The summary is the median of the three run medians.',
               'Setup, allocation, and input initialization are outside the timer. In-place calls reuse the buffer without resetting or copying it between calls.',
               'The timed entry point is `encodeUnchecked` with identical input and output pointers; the checked `encodeInplace` API adds argument validation.',
               'Thus later calls see the preceding output prefix; the encoder uses the same fixed operation schedule for every input value.',
               'One separate out-of-place check measured 17.701 ms at K=2^20. It is not the in-place headline measurement.', '',
               'Release and AddressSanitizer/UndefinedBehaviorSanitizer tests pass for both outers. The quarter-rate tests cover:', '',
               '- Every outer coordinate impulse in both AVX2 lanes.',
               '- Full independent dense oracles at K=2^16, 2^18, and 2^20.',
               '- Both routing layouts, alternate tile sizes, a second setup, boundary impulses, zero input, and linearity.',
               '- In-place equality, unchanged buffer suffix, compaction, and rejected invalid API arguments.', '',
               'Three Python identity tests verify source hashes, equality of the outer row spaces, and the exact selected inner.',
               'The existing 28 quarter-rate proof tests also pass. The code generator verifies its XOR circuit symbolically.', '',
               'See [README.md](README.md) for build commands and the API. [receipts](receipts/) retains compact timing rows, test logs, environment data, and binary hashes.',
               'Rebuild this summary with `python workstreams/rate_quarter_bch/implementation/report.py`.', '']
    (HERE/'PERFORMANCE.md').write_text('\n'.join(report),encoding='utf-8',newline='\n')
    print(json.dumps(summaries,indent=2))


if __name__ == '__main__':
    main()
