"""Fast all-grid zero-state diagnostics using a two-weight coefficient term.

For each regional weight q, retain only epochs of the two adjacent even
weights bracketing q/E. This gives an explicit positive coefficient lower
bound and avoids a large polynomial table at the longest message lengths.
"""
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import gammaln

import bch_zero_state_lower_v1 as base
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


def two_weight_logs(kernel, epochs, maximum):
    length = (len(kernel)-1)*epochs
    q = np.arange(0, maximum+1, 2, dtype=np.int64)
    a = 2*(q//(2*epochs)); h = (q-a*epochs)//2
    log_kernel = np.array([math.log(n) if n else -np.inf for n in kernel]+[-np.inf]*2)
    first = np.zeros(len(q)); second = np.zeros(len(q))
    np.multiply(epochs-h, log_kernel[a], out=first, where=(epochs-h) != 0)
    np.multiply(h, log_kernel[a+2], out=second, where=h != 0)
    return (gammaln(epochs+1)-gammaln(h+1)-gammaln(epochs-h+1)+first+second
            -gammaln(length+1)+gammaln(q+1)+gammaln(length-q+1))


def evaluate(row, counts, config):
    block, t, s, exponent = (int(row[k]) for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
    length = (1 << exponent)//(block//2); cutoff = block*length//10
    queries = []
    for weight, count in counts.items():
        if weight == block:
            continue
        maximum_q = min(length, cutoff//weight)
        for fraction in (.5, .75, 1.):
            q = 2*int(maximum_q*fraction/2)
            if q < 2:
                continue
            lo, hi, good = base.concentrated_range(q, weight/block, block)
            first, last = lo+lo % 2, hi-hi % 2
            if first >= 4 and first <= last:
                queries.append((q, weight, count, first, last, good))
    if not queries:
        return dict(log_first_moment_lower=None, evidence='NO_USABLE_CONCENTRATION_INTERVAL')
    logs = two_weight_logs(config['kernel_counts'], length//t, max(r[4] for r in queries))
    best = None
    for q, weight, count, lo, hi, good in queries:
        value, segment = base.convex_lower_value(np.arange(lo, hi+1, 2), logs[lo//2:hi//2+1], q*weight/block)
        if not math.isfinite(value):
            continue
        choose = float(gammaln(length+1)-gammaln(q+1)-gammaln(length-q+1))
        bound = choose+q*math.log(count)+good+block*value
        if best is None or bound > best['log_first_moment_lower']:
            best = dict(log_first_moment_lower=bound, first_moment_lower_bits=bound/LN2,
                        occupation=q, outer_weight=weight, output_weight=q*weight,
                        regional_weight_range=[lo, hi], lower_convex_segment=segment,
                        good_parity_log_lower=good, evidence='TWO_WEIGHT_KERNEL_COEFFICIENT_LOWER')
    return best or dict(log_first_moment_lower=None, evidence='TWO_WEIGHT_SUPPORT_INSUFFICIENT')


def main():
    rows, spectra, maps, dependencies = study.load_inputs()
    for path in (Path(__file__), Path(base.__file__)):
        dependencies[path.relative_to(study.ROOT).as_posix()] = study.sha(path)
    results = []
    for index, row in enumerate(rows, 1):
        b, t, s, e = (int(row[k]) for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
        value = evaluate(row, spectra[b], maps[t, s])
        result = dict(block_bits=b, step_bits=t, state_bits=s, message_exponent=e, **value)
        results.append(result)
        if index % 10 == 0:
            print(f'zero-state grid {index}/{len(rows)}', flush=True)
    payload = dict(status='BINARY64_FIRST_MOMENT_LOWER_GRID', source_sha256=dependencies,
                   rows=results, obstructed=sum(r.get('first_moment_lower_bits', -math.inf) > 0 for r in results),
                   limitations=['A nonpositive or unavailable lower bound does not establish a useful upper bound.',
                                'A positive lower bound obstructs this first-moment certificate, not every possible proof.',
                                'Nearest binary64; selected high-precision replay is required before interpretation.'])
    (HERE/'bch_zero_state_lower_grid_v2.json').write_text(json.dumps(payload, indent=2)+'\n')
    fields = ['block_bits', 'step_bits', 'state_bits', 'message_exponent', 'first_moment_lower_bits', 'occupation', 'outer_weight', 'evidence']
    with (HERE/'bch_zero_state_lower_grid_v2.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore'); writer.writeheader(); writer.writerows(results)
    print('obstructed', payload['obstructed'], 'of', len(rows), flush=True)


if __name__ == '__main__':
    main()
