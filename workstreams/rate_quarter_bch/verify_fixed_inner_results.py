"""Replay retained dense witnesses and check their exact integer coverage.

Numerical replay remains binary64. Exact partition checks do not convert
floating-point transfer evaluations into an outward-rounded certificate.
"""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import evaluate_fixed_inner as fixed
import typed_dense_boxes as typed


def check_coverage(boxes, length, minimum):
    for box in boxes:
        lo, hi = box['lower'], box['upper']
        fixed.outer.require(len(lo) == len(hi) == 3 and all(type(x) is int for x in lo+hi), 'invalid type box')
        fixed.outer.require(all(0 <= a <= b <= length for a,b in zip(lo,hi)), 'invalid box endpoints')
        fixed.outer.require(hi[0] <= length-minimum, 'box includes omitted sparse occupations')
    # For fixed total active Q, the all-one count parametrizes every type.
    for q in range(minimum, length+1):
        intervals = []
        for box in boxes:
            lo, hi = box['lower'], box['upper']
            if lo[0] <= length-q <= hi[0]:
                a, b = max(lo[2], q-hi[1], 0), min(hi[2], q-lo[1], q)
                if a <= b:
                    intervals.append((a,b))
        following = 0
        for a,b in sorted(intervals):
            fixed.outer.require(a == following, f'gap or overlap at occupation {q}')
            following = b+1
        fixed.outer.require(following == q+1, f'incomplete occupation {q}')


def verify(path, a, kernel):
    payload = json.loads(path.read_text())
    for name, digest in payload['source_sha256'].items():
        fixed.outer.require(fixed.sha(fixed.ROOT/name) == digest, 'changed dependency: '+name)
    geo, result = payload['geometry'], payload['result']
    fixed.outer.require(geo == fixed.geometry(geo['message_exponent'],Fraction(geo['distance_target'])), 'geometry or cutoff changed')
    fixed.outer.require(payload['inner'] == 'existing optimized t128_s19' and not payload['inner_reoptimized'], 'inner changed')
    length = geo['outer_rows']
    fixed.outer.require(geo['output_bits'] == 4*geo['message_bits'] == 256*length, 'wrong rate')
    fixed.outer.require(result['occupation_max'] == length, 'incomplete dense endpoint')
    boxes = result['selected_boxes']
    check_coverage(boxes, length, result['occupation_min'])
    spectrum = fixed.outer.read_spectrum(fixed.HERE/'BCH256_64.wd')
    counts = {w:n for w,n in enumerate(spectrum) if w and n}
    fixed.outer.require(result['bands'] == [sorted(w for w in counts if w != 256), [256]], 'wrong type bands')
    cache, logs = {}, []
    for box in boxes:
        witness = box['witness']
        z = witness['log_surprisal']
        if z not in cache:
            cache[z] = typed.general.epoch_logs(128,19,a,kernel,math.exp(z),128)
        bank = result['probability_banks'][witness['probability_bank']]
        probabilities, costs = np.array(bank['probabilities']), np.array(bank['log_density_costs'])
        fixed.outer.require(probabilities[0] == 0 and probabilities[2] == 1 and costs[0] == costs[2] == 0, 'wrong endpoint type')
        for g, band in enumerate(result['bands'],1):
            if band == [256]:
                continue
            p = probabilities[g]
            expected = max(math.log(counts[w])-math.log(math.comb(256,w))-w*math.log(p)-(256-w)*math.log1p(-p) for w in band)
            fixed.outer.require(abs(expected-costs[g]) < 1e-10, 'density cost changed')
        proposal = np.array(witness['proposal'])
        theta = float(proposal@probabilities)
        moment, _ = typed.moment_and_density(cache[z],128,256*length,theta)
        corners = typed.vertices(box['lower'],box['upper'],length)
        bound = max(typed.point_logs(corners,length,256,costs,proposal,moment,geo['bad_weight'],math.exp(z)))
        bound += typed.lattice_log_count(box['lower'],box['upper'])
        fixed.outer.require(abs(bound-box['own_log_bound']) < 2e-6, 'dense witness replay differs')
        logs.append(float(bound))
    total = float(np.logaddexp.reduce(logs))
    fixed.outer.require(abs(total-result['log_union_upper']) < 2e-6, 'dense union differs')
    return payload, total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipts', type=Path, nargs='+')
    parser.add_argument('--output', type=Path, help='optional combined diagnostic summary')
    args = parser.parse_args()
    a, kernel, _ = fixed.load_inner()
    sparse = json.loads((fixed.HERE/'FIXED_INNER_SCREEN.json').read_text())
    for name,digest in sparse['source_sha256'].items():
        fixed.outer.require(fixed.sha(fixed.ROOT/name) == digest, 'changed sparse dependency: '+name)
    summaries = []
    for path in args.receipts:
        payload, dense_log = verify(path,a,kernel)
        geo = payload['geometry']
        row = next(r for r in sparse['results'] if r['message_exponent'] == geo['message_exponent'] and r['distance_target'] == geo['distance_target'])
        fixed.outer.require(row['covered_occupations'][1]+1 == payload['result']['occupation_min'], 'gap between sparse and dense')
        total = float(np.logaddexp(-row['partial_union_margin_bits']*math.log(2),dense_log))
        summary = dict(file=path.name, message_exponent=geo['message_exponent'], distance_target=geo['distance_target'],
            integer_coverage='all nonzero messages', sparse_margin_bits=row['partial_union_margin_bits'],
            dense_margin_bits=-dense_log/math.log(2),
            dense_witness_replay=True, combined_margin_bits=-total/math.log(2),
            status='BINARY64_DIAGNOSTIC_NOT_OUTWARD_CERTIFICATE')
        summaries.append(summary)
        print(json.dumps(summary),flush=True)
    if args.output:
        dependencies = [Path(__file__),fixed.HERE/'FIXED_INNER_SCREEN.json']+args.receipts
        args.output.write_text(json.dumps(dict(results=summaries,
            source_sha256={p.resolve().relative_to(fixed.ROOT).as_posix():fixed.sha(p) for p in dependencies}),indent=2)+'\n',
            encoding='utf-8',newline='\n')


if __name__ == '__main__':
    main()
