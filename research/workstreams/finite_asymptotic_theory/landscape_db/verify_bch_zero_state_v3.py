"""90-digit replay of exact-coefficient zero-state witnesses."""
import json
import math
from pathlib import Path

import mpmath as mp

import bch_zero_state_lower_v1 as concentration
import kernel_integer_coefficients_v1 as exact
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent


def main():
    _, spectra, maps, _ = study.load_inputs()
    path = HERE/'bch_zero_state_exact_grid_v3.json'
    payload = json.loads(path.read_text())
    for name, digest in payload['source_sha256'].items():
        if study.sha(study.ROOT/name) != digest:
            raise ValueError(f'changed lower-bound dependency: {name}')
    selected = [r for r in payload['rows'] if r.get('first_moment_lower_bits') is not None and r['first_moment_lower_bits'] > 0]
    checks = []; cache = {}
    with mp.workdps(90):
        def log_choose(n, k):
            key = (n, k)
            if key not in cache:
                cache[key] = mp.loggamma(n+1)-mp.loggamma(k+1)-mp.loggamma(n-k+1)
            return cache[key]
        for row in selected:
            block, t, s, exponent = (row[k] for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
            length = (1 << exponent)//(block//2); epochs = length//t
            q, w = row['occupation'], row['outer_weight']
            if q % 2 or row['output_weight'] != q*w or q*w > block*length//10:
                raise ValueError('invalid restricted bad-word class')
            p = mp.mpf(w)/block
            lo, hi, _ = concentration.concentrated_range(q, w/block, block)
            if row['regional_weight_range'] != [lo+lo % 2, hi-hi % 2]:
                raise ValueError('regional count interval changed')
            def divergence(a):
                if a == 0:
                    return -mp.log1p(-p)
                if a == 1:
                    return -mp.log(p)
                return a*mp.log(a/p)+(1-a)*mp.log((1-a)/(1-p))
            outside = block*((mp.exp(-q*divergence(mp.mpf(lo-1)/q)) if lo else 0)
                             +(mp.exp(-q*divergence(mp.mpf(hi+1)/q)) if hi < q else 0))
            parity = mp.power(2, 1-block)
            good = parity*mp.mpf(31)/32
            if parity-outside < good:
                raise ArithmeticError('concentration charge is too small')
            kernel = maps[t, s]['kernel_counts']
            integer_coefficients = exact.coefficients(kernel, epochs, row['regional_weight_range'][1])
            regional_cache = {}
            def regional(weight):
                if weight not in regional_cache:
                    coefficient = integer_coefficients[weight//2]
                    regional_cache[weight] = (mp.log(coefficient)-log_choose(length, weight)
                                              if coefficient else mp.ninf)
                return regional_cache[weight]
            segment = row['lower_convex_segment']
            a, b = segment[0][0], segment[-1][0]
            va, vb = regional(a), regional(b)
            slope = (vb-va)/(b-a) if a != b else mp.mpf(0)
            first, last = row['regional_weight_range']
            # Check the proposed supporting line against every regional
            # weight, correcting even a tiny numerical hull discrepancy.
            correction = min(mp.mpf(0), min(regional(j)-(va+(j-a)*slope) for j in range(first, last+1, 2)))
            moment_lower = va+(mp.mpf(q*w)/block-a)*slope+correction
            bound = log_choose(length, q)+q*mp.log(spectra[block][w])+mp.log(good)+block*moment_lower
            error = abs(float(bound)-row['log_first_moment_lower'])
            if error > 2e-5 or bound <= 0:
                raise ArithmeticError(f'positive lower bound failed replay: {row} error={error}')
            check = dict(block_bits=block, step_bits=t, state_bits=s, message_exponent=exponent,
                         absolute_log_error=error, lower_bits=float(bound/mp.log(2)),
                         supporting_line_correction=str(correction), regional_weights_checked=(last-first)//2+1)
            checks.append(check)
            print(f'B{block} t{t} s{s}: {check["lower_bits"]:.3f} lower bits, error {error:.3g}', flush=True)
    result = dict(status='ALL_POSITIVE_GRID_WITNESSES_REPLAYED', decimal_digits=90, checks=checks,
                  maximum_log_error=max(c['absolute_log_error'] for c in checks),
                  source_sha256={str(Path(__file__)): study.sha(Path(__file__)), str(path): study.sha(path), str(Path(exact.__file__)): study.sha(Path(exact.__file__))},
                  limitations=['High-precision diagnostics, not outward interval certificates.',
                               'A first-moment obstruction is not a lower bound on failure probability.'])
    (HERE/'bch_zero_state_verification_v3.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()

