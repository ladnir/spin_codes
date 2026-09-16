"""Quarter-rate diagnostic bounds with the existing optimized t128_s19 map.

No inner search, sampling, benchmarks, or outward certification. Counts and
geometry are exact; transfer calculations use nearest binary64 arithmetic.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
import sys

import numpy as np
import reconstruct as outer

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LANDSCAPE = ROOT / 'workstreams/finite_asymptotic_theory/landscape_db'
BRIDGE = ROOT / 'workstreams/bch_rm2sub_bridge'
sys.path[:0] = [str(LANDSCAPE), str(BRIDGE)]
import activation_q1 as q1
import activation_occupation as general
import composition_occupation as composition
import audit_imported_rm2sub as maps


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_inner():
    directory = BRIDGE / 'generated/larger_state_inputs_v1'
    manifest_path = directory / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    paths = [directory / f't128_s19_{suffix}.json' for suffix in
             ('selection', 'a_spectrum', 'b_kernel_spectrum')]
    for path in paths:
        outer.require(sha(path) == manifest['snapshot_sha256'][path.name], 'inner snapshot changed')
    selection, a, kernel = [json.loads(p.read_text()) for p in paths]
    columns = [int(v, 16) for v in selection['selected']['B_columns_hex']]
    rows = maps.generators(columns, 19)
    outer.require(rows == [int(v, 16) for v in selection['selected']['A_generator_words_hex']], 'A != transpose(B)')
    outer.require(len(set(columns)) == 128 and 0 not in columns, 'invalid columns')
    outer.require(not any((x & y).bit_count() % 2 for x in rows for y in rows), 'BA != 0')
    exact = maps.spectrum(rows)
    outer.require(exact == {r['weight']: r['count'] for r in a['spectrum'] if r['count']}, 'A spectrum mismatch')
    dual = maps.dual_spectrum(exact, 128, 19)
    claimed = {r['total_weight']: r['kernel_words'] for r in kernel['by_total_weight'] if r['kernel_words']}
    outer.require(dual == claimed, 'kernel spectrum mismatch')
    implementation = ROOT / 'workstreams/bare_bch_rm2sub/generated/manifest.json'
    record = json.loads(implementation.read_text())['t128_s19']
    outer.require(record['selection_sha256'] == sha(paths[0]) and record['columns'] == columns,
                  'map differs from optimized implementation')
    return {w: n for w, n in exact.items() if w}, [dual.get(w, 0) for w in range(129)], paths + [manifest_path, implementation]


def geometry(exponent, delta):
    bits = 1 << exponent
    outer.require(bits % 64 == 0 and (bits // 64) % 128 == 0, 'incomplete rows or epochs')
    outer.require(0 < delta < Fraction(1, 2), 'invalid distance target')
    length = bits // 64
    size = 256 * length
    return dict(message_exponent=exponent, message_bits=bits, output_bits=size,
                outer_rows=length, epochs_per_region=length // 128,
                distance_target=str(delta), bad_weight=size * delta.numerator // delta.denominator)


def aggregate_q1(moments, tilts, counts, geo):
    lam = np.exp(tilts)
    values = np.minimum(0., moments + geo['bad_weight'] * lam[:, None])
    witnesses = np.argmin(values, axis=0)
    best = values[witnesses, np.arange(257)]
    weights = sorted(counts)
    terms = [math.log(geo['outer_rows']) + math.log(counts[w]) + best[w] for w in weights]
    dominant = weights[int(np.argmax(terms))]
    return dict(occupation=1, margin_bits=-float(np.logaddexp.reduce(terms))/math.log(2),
                dominant_weight=dominant, log_surprisal=float(tilts[witnesses[dominant]]),
                witness_at_grid_edge=bool(witnesses[dominant] in (0, len(tilts)-1)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exponents', type=int, nargs='+', default=[16, 18, 20])
    parser.add_argument('--distances', type=Fraction, nargs='+', default=[Fraction(1, 10), Fraction(19, 100)])
    parser.add_argument('--maximum-occupation', type=int, default=64)
    parser.add_argument('--output', type=Path, default=HERE / 'FIXED_INNER_SCREEN.json')
    args = parser.parse_args()
    outer.require(args.maximum_occupation >= 1, 'invalid occupation limit')
    a, kernel, sources = load_inner()
    spectrum_path = HERE / 'BCH256_64.wd'
    spectrum = outer.read_spectrum(spectrum_path)
    outer.audit_spectrum(spectrum, 64, 62)
    small = outer.read_spectrum(HERE / 'sources/EBCH256_63.wd')
    parent = outer.read_spectrum(HERE / 'sources/EBCH256_71.wd')
    outer.require(spectrum == outer.reconstruct(small, parent), 'outer spectrum mismatch')
    counts = {w: n for w, n in enumerate(spectrum) if w and n}
    tilts = np.arange(-180, 1, dtype=float)/10
    sparse_tilts = np.arange(-14., .01, .5)
    shifts = [-.5, 0., .5]
    results = []
    for exponent in args.exponents:
        geometries = [geometry(exponent, delta) for delta in args.distances]
        length = geometries[0]['outer_rows']
        maximum = min(args.maximum_occupation, length)
        moments = q1.coefficient_logs(*q1.region_logs(*q1.epoch_logs(128, 19, a, np.exp(tilts)), length//128), 256)
        rows = [dict(g, q1=aggregate_q1(moments, tilts, counts, g)) for g in geometries]
        print(json.dumps([dict(exponent=exponent, distance=r['distance_target'], **r['q1']) for r in rows]), flush=True)
        best = np.full((len(rows), maximum), np.inf)
        witnesses = np.zeros((len(rows), maximum, 2))
        compositions = {q: composition.SparseComposition(counts, 256, q) for q in range(2, min(4, maximum)+1)}
        component_best = {q: np.full((len(rows), len(model.indices)), np.inf) for q, model in compositions.items()}
        for log_lam in sparse_tilts if maximum > 1 else []:
            lam = math.exp(log_lam)
            epoch = general.epoch_logs(128, 19, a, kernel, lam, min(128, maximum))
            regions = general.region_logs(epoch, 128, length, maximum)
            for i, row in enumerate(rows):
                for shift in shifts:
                    values = general.occupation_bounds(regions, counts, 256, length, row['bad_weight'], lam, shift, 8)
                    improved = values < best[i]
                    best[i, improved] = values[improved]
                    witnesses[i, improved] = [log_lam, shift]
                    for q, model in compositions.items():
                        values = model.components(regions, row['bad_weight'], lam, shift)
                        np.minimum(component_best[q][i], values, out=component_best[q][i])
        for i, row in enumerate(rows):
            row['composition_refinements'] = []
            adaptive = best[i].copy()
            for q, model in compositions.items():
                value = model.aggregate(component_best[q][i], length)
                row['composition_refinements'].append(dict(occupation=q, margin_bits=-value/math.log(2)))
                best[i, q-1] = min(best[i, q-1], value)
            row['higher_occupations'] = [dict(occupation=q, margin_bits=-float(best[i, q-1])/math.log(2),
                method='composition' if best[i, q-1] < adaptive[q-1] else 'adaptive',
                adaptive_margin_bits=-float(adaptive[q-1])/math.log(2),
                adaptive_log_surprisal=float(witnesses[i, q-1, 0]),
                adaptive_shift=float(witnesses[i, q-1, 1])) for q in range(2, maximum+1)]
            logs = [-row['q1']['margin_bits']*math.log(2)] + list(best[i, 1:])
            row['covered_occupations'] = [1, maximum]
            row['partial_union_margin_bits'] = -float(np.logaddexp.reduce(logs))/math.log(2)
            row['full_occupation_coverage'] = maximum == length
        results.extend(rows)
        print(f'e{exponent}: Q1..Q{maximum} evaluated', flush=True)
    sources += [Path(__file__), Path(outer.__file__), Path(q1.__file__), Path(general.__file__),
                Path(maps.__file__), Path(composition.__file__), spectrum_path, HERE/'sources/EBCH256_63.wd', HERE/'sources/EBCH256_71.wd']
    payload = dict(status='BINARY64_PARTIAL_OCCUPATION_DIAGNOSTIC', inner='existing optimized t128_s19',
        inner_reoptimized=False, outer=[256, 64, 62], q1_log_surprisals=tilts.tolist(),
        sparse_log_surprisals=sparse_tilts.tolist(), sparse_shifts=shifts, sparse_band_count=8,
        source_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in sources}, results=results)
    args.output.write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8', newline='\n')


if __name__ == '__main__':
    main()
