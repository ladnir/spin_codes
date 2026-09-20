"""Exact kernel coefficient diagnostics on all BCH t/s choices at K=2^20."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import bch_zero_state_lower_v1 as base
import kernel_integer_coefficients_v1 as exact
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


def evaluate(block, t, s, counts, config, exponent=20):
    length = (1 << exponent)//(block//2); cutoff = block*length//10
    queries = []
    for weight, count in counts.items():
        maximum_q = min(length, cutoff//weight)
        for fraction in (.35, .5, .65, .75, .8, .9, 1.):
            q = 2*int(maximum_q*fraction/2)
            if q < 2 or weight == block:
                continue
            lo, hi, good = base.concentrated_range(q, weight/block, block)
            first, last = lo+lo % 2, hi-hi % 2
            if first >= 4 and first <= last:
                queries.append((q, weight, count, first, last, good))
    maximum = max(r[4] for r in queries)
    logs = exact.regional_logs(config['kernel_counts'], length//t, maximum)
    best = None
    for q, weight, count, lo, hi, good in queries:
        value, segment = base.convex_lower_value(np.arange(lo, hi+1, 2), logs[lo//2:hi//2+1], q*weight/block)
        if not math.isfinite(value):
            continue
        bound = math.log(math.comb(length, q))+q*math.log(count)+good+block*value
        if best is None or bound > best['log_first_moment_lower']:
            best = dict(log_first_moment_lower=bound, first_moment_lower_bits=bound/LN2,
                        occupation=q, outer_weight=weight, output_weight=q*weight,
                        regional_weight_range=[lo, hi], lower_convex_segment=segment,
                        good_parity_log_lower=good)
    return dict(block_bits=block, step_bits=t, state_bits=s, message_exponent=exponent,
                kernel_coefficient_maximum=maximum,
                **(best or dict(log_first_moment_lower=None, first_moment_lower_bits=None,
                                 evidence='EVEN_PROFILE_SUPPORT_INSUFFICIENT')))


def main():
    _, spectra, maps, dependencies = study.load_inputs()
    for path in (Path(__file__), Path(base.__file__), Path(exact.__file__)):
        dependencies[path.relative_to(study.ROOT).as_posix()] = study.sha(path)
    fingerprint = hashlib.sha256(json.dumps(dependencies, sort_keys=True).encode()).hexdigest()
    directory = HERE/'bch_zero_state_exact_v3'; directory.mkdir(exist_ok=True)
    rows = []
    for block in (64, 128):
        for t, s in sorted(maps):
            path = directory/f'b{block}_t{t}_s{s}_e20.json'
            if path.exists():
                data = json.loads(path.read_text())
                if data['input_fingerprint'] != fingerprint:
                    raise ValueError('stale exact zero-state checkpoint')
                result = data['row']
            else:
                result = evaluate(block, t, s, spectra[block], maps[t, s])
                path.write_text(json.dumps(dict(input_fingerprint=fingerprint, row=result), indent=2)+'\n')
            rows.append(result)
            value = result['first_moment_lower_bits']
            print(f'B{block} t{t} s{s}: lower exponent {value}', flush=True)
    payload = dict(status='BINARY64_EXACT_KERNEL_LOWER_GRID', rows=rows, source_sha256=dependencies,
                   input_fingerprint=fingerprint, obstructed=sum(r['first_moment_lower_bits'] is not None
                                                                 and r['first_moment_lower_bits'] > 0 for r in rows),
                   limitations=['Kernel coefficients are exact integers; logarithms use nearest binary64 arithmetic.',
                                'Expected bad-word counts are not failure probabilities.',
                                'Nonpositive lower bounds do not certify closure.'])
    (HERE/'bch_zero_state_exact_grid_v3.json').write_text(json.dumps(payload, indent=2)+'\n')


if __name__ == '__main__':
    main()
