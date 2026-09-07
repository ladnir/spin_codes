"""Validate serial benchmark receipts and produce a compact performance report."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re
import statistics

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
NAMES = ('t64_s20', 't128_s19', 't64_s16', 't256_s14')


def build(directory):
    checks = {}
    log = (directory/'correctness-final.txt').read_text()
    if not log.rstrip().endswith('correctness=PASS'):
        raise ValueError('Missing complete correctness result')
    for name, m, digest in re.findall(r'(t\d+_s\d+) m=(\d+) hash=([a-f0-9]+) route=[a-f0-9]+ PASS', log):
        key = (name, int(m))
        if key in checks:
            raise ValueError('Duplicate correctness cell')
        checks[key] = digest
    expected = {(n,m) for n in NAMES for m in (16,18,20)}
    if set(checks) != expected:
        raise ValueError('Incomplete correctness grid')
    source_hashes = {}
    for line in (directory/'sources-final.sha256').read_text().splitlines():
        digest, relative = line.split(maxsplit=1)
        path = ROOT/relative
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'Final source hash mismatch: {relative}')
        source_hashes[relative] = digest
    samples = defaultdict(list)
    raw = []
    for run in (1,2,3):
        rows = [json.loads(s) for s in (directory/f'final-r{run}.jsonl').read_text().splitlines()]
        if len(rows) != 12 or {(r['configuration'],r['m']) for r in rows} != expected:
            raise ValueError('Incomplete or duplicate performance grid')
        for row in rows:
            key = (row['configuration'], row['m'])
            if row['output_hash'] != checks[key] or row['trials'] != 101 or row['layout'] != 'packed24':
                raise ValueError('Performance receipt mismatch')
            if row['tile_rows'] != (256 if row['m'] <= 18 else 2048):
                raise ValueError('Wrong tile policy')
            if not 0 < row['p10_ms'] <= row['median_ms'] <= row['p90_ms']:
                raise ValueError('Invalid timing quantiles')
            samples[key].append(row)
            raw.append(dict(run=run, **row))
    summary = []
    for name in NAMES:
        for m in (16,18,20):
            rows = samples[name,m]
            medians = [r['median_ms'] for r in rows]
            evidence = (f'full certificate, { {16:53.944367,18:52.346388,20:50.448203}[m]:.6f} bits' if name=='t128_s19' else
                        'full reference, 50.487298 bits' if name=='t64_s20' and m==20 else
                        'shortlisted, not fully certified' if name in ('t128_s19','t64_s16') else
                        'exploratory watch' if name=='t256_s14' else 'reference parameters; no certificate replay here')
            summary.append(dict(configuration=name,m=m,median_ms=statistics.median(medians),
                run_min_median_ms=min(medians),run_max_median_ms=max(medians),
                p10_ms=min(r['p10_ms'] for r in rows),p90_ms=max(r['p90_ms'] for r in rows),
                setup_median_ms=statistics.median(r['setup_ms'] for r in rows),
                retained_setup_bytes=rows[0]['retained_setup_bytes'],workspace_bytes=rows[0]['workspace_bytes'],
                evidence=evidence))
    result = dict(status='MEASURED_SERIAL_PERFORMANCE_NOT_NEW_DISTANCE_CERTIFICATES',
        machine='Peach, AMD Ryzen 9 7950X, Linux, CPU 15, GCC 15.2.0',
        compiler_flags='Release (-O3 -DNDEBUG), -march=znver4 -mavx2 -mpclmul -mvpclmulqdq',
        method='Three serial runs per cell; three warmups then 101 calls per run; median of run medians.',
        scope='Transposed block encoder; each block carries 128 parallel binary instances; no fanout.',
        source_sha256=source_hashes,
        binary_sha256=(directory/'binaries-final.sha256').read_text(),
        summary=summary,raw_final_samples=raw)
    (HERE/'PERFORMANCE.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    lookup = {(r['configuration'],r['m']):r for r in summary}
    lines = ['# Bare BCH-256/RM2Sub performance', '',
        'Measured on Peach on 2026-09-07. All 12 cells pass the independent dense oracle;',
        'the GCC and MSVC output hashes agree. No fanout is included.', '',
        '## Online encoding time', '',
        'Times are milliseconds for the transposed block encoder. Each 128-bit block',
        'carries 128 parallel binary instances. Values are medians of three run medians.', '',
        '| Inner (t,s) | K=2^16 | K=2^18 | K=2^20 | K=2^20 interpretation |',
        '|---|---:|---:|---:|---|']
    for n in NAMES:
        cells = ' | '.join(f"{lookup[n,m]['median_ms']:.3f}" for m in (16,18,20))
        lines.append(f"| {n} | {cells} | {lookup[n,20]['evidence']} |")
    baseline = lookup['t64_s20',20]['median_ms']
    s19 = lookup['t128_s19',20]['median_ms']
    lines += ['', f'S19 takes {100*(1-s19/baseline):.1f}% less online time than the S20 reference at K=2^20.',
        'The wider-step candidate is faster, but its timing does not resolve its weaker proof evidence.',
        'Use the selected (128,19) map as the main certified performance point; retain (64,20) as a reference.', '',
        '## Variation and setup costs at K=2^20', '',
        '| Inner | Range of run medians (ms) | Setup (ms) | Retained setup (MiB) | Workspace (MiB) |',
        '|---|---:|---:|---:|---:|']
    for n in NAMES:
        r=lookup[n,20]
        lines.append(f"| {n} | {r['run_min_median_ms']:.3f}–{r['run_max_median_ms']:.3f} | {r['setup_median_ms']:.2f} | {r['retained_setup_bytes']/2**20:.2f} | {r['workspace_bytes']/2**20:.2f} |")
    lines += ['', 'Setup is measured separately and includes creating oracle and alternative routing data.',
        'The benchmark discards that data before online timing. The setup column excludes this final',
        'deallocation and workspace construction. Input/output buffers add 48 MiB at K=2^20.',
        'Those buffers are separate from the listed workspace. Setup storage before compaction is',
        'approximately 36–39 MiB; PERFORMANCE.json records exact byte counts and all final observations.', '',
        '## Measurement method', '',
        '- Ryzen 9 7950X, Linux, one thread pinned to CPU 15; GCC 15.2.0.',
        '- CMake Release: `-O3 -DNDEBUG -march=znver4 -mavx2 -mpclmul -mvpclmulqdq`.',
        '- Three independent serial processes, three warmups, 101 timed calls per cell per process.',
        '- Identical inputs, routing seeds, and coefficient seeds across repetitions.',
        '- Setup, workspace allocation, correctness checks, and output hashing are outside online timing.',
        '- Every cell checks for other benchmark executables; a process lock prevents concurrent copies.',
        '- The same buffers and workspace are reused. This is steady-state latency, not cold-cache latency.',
        '- The tuning sweep precedes the final runs; final numbers are not minima selected from that sweep.', '',
        '## Optimizations retained', '',
        'All configurations share routing, setup, the paired BCH transpose, and the RM2Sub template.',
        'Each selected map supplies compile-time constants and an exactly checked quadratic XOR circuit.',
        'The kernels use pruned SIMD zeta stages, unrolled emission, fixed nibble tables, and',
        'a templated emitter. The hot call performs no allocation or type-erased callback dispatch.',
        'Configuration and layout dispatch occur once per encode.', '',
        'The bounded tuning effort selected these choices:', '',
        '- Packed 24-bit routing rather than 32-bit indices at K=2^20.',
        '- 256-row tiles at K=2^16 and K=2^18; 2048-row tiles at K=2^20.',
        '- Optimized nibble grouping for S19 only. Three A/B repetitions showed about 3% lower latency.',
        '  Other maps keep sequential grouping because fewer nominal lookups did not give a reliable win.',
        '- Setup compaction removes the dense route, scalar coefficient list, and unused routing format.', '',
        'The configuration constructor exposes tile size, and CMake exposes `SPIN_GROUPED_A`.',
        'These are measured defaults on this host, not a global optimum for every processor.', '',
        '## Construction and correctness', '',
        'The BCH matrix is not the imported Eq3 matrix. Let p and q denote the length-255 BCH',
        'generator polynomials for designed distances 37 and 39. Let P and Q denote their',
        'parity-extended codes. The implemented code spans Q and the parity extensions of',
        'five representatives `p*x^j`, for `j=0..4`, and uses a systematic basis.',
        'Generation checks rank 128, Q containment, P containment, and inclusion of the all-one word.',
        'This is the intermediate-code family used by the existing bounds.', '',
        'The generated paired transpose uses 3029 XORs in its scalar circuit and processes two rows',
        'with AVX2. Symbolic checks verify its matrix and transpose. Separate symbolic checks verify',
        'the pruned zeta schedule and the selected inner maps. Runtime tests compare complete outputs',
        'against dense BCH multiplication and a separately implemented scalar RM2Sub recurrence.',
        'Tests also cover both route layouts, tile changes, compaction, multiple setups, boundary',
        'impulses, zero input, and invalid or overlapping buffers.', '',
        'The known 50.487298-bit label applies to the selected S20 reference at K=2^20.',
        'Subsequent proof work certified the selected S19 map at K=2^16, 2^18, 2^20, 2^22, and 2^24.',
        'See [the paper handoff](../bch_rm2sub_bridge/PAPER_HANDOFF.md) for the exact scope and evidence.',
        'Timings and implementation support remain limited to the sizes described above.',
        'Historical Eq3/fanout timings concern a different map and are not used as the comparison baseline.', '',
        '## Reproduction and next step', '',
        'See README.md for build commands and the API. PERFORMANCE.json includes source hashes,',
        'binary hashes, exact memory sizes, quantiles, and the 36 final cell observations.',
        'Local tuning and correctness logs are in `out/bare_spin_receipts/`.', '',
        'Next, integrate S20 and S19 into the calling application and measure end-to-end throughput.',
        'The S19 distance/setup proof is complete at the listed sizes; retain S16 and T256 as comparison points.', '']
    (HERE/'PERFORMANCE.md').write_text('\n'.join(lines),encoding='utf-8')
    print('\n'.join(lines[:22]))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    build(p.parse_args().directory)
