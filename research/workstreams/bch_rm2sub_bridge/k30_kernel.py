"""Kernel-aware four-state transfer; certificates cover only listed occupancies.

The kernel spectrum gives the exact probability of zero syndrome. Separate
positive bounds on its emission moment and the nonzero-syndrome moment avoid
assuming independence between syndrome and emitted weight.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import numpy as np
from flint import arb, arb_poly, ctx
import bridge as base
import larger_state_maps as maps
import certify_refresh_q1 as matrix
import k30_sparse as sparse
from certify_k30_q1 import ROWS, CUTOFF
from audit_bch_q1_full_arb import rational
from migrate_legacy_workspace import restore_or_verify


def epochs(t, s, spectrum, kernel, z, maximum, number):
    assert 0 <= maximum <= t
    assert 0 < z < 1
    assert sum(spectrum.values()) == (1 << s) - 1
    result = [matrix.epoch(t, s, spectrum, z, number)[0]]
    m = (1 << s) - 1
    kappa = number(m) / (m - 1)
    upper = (lambda x: x.upper()) if number is arb else (lambda x: x)
    for j in range(1, maximum + 1):
        beta = number(kernel.get(j, 0)) / math.comb(t, j)
        uniform = sum((number(n) * sum((
            number(math.comb(w, v) * math.comb(t-w, j-v)) * z**(w+j-2*v)
            for v in range(max(0, j-t+w), min(w, j)+1)), number(0))
            / math.comb(t, j) for w, n in spectrum.items()), number(0)) / m
        arbitrary = z**min(abs(w-j) for w in spectrum)
        row = [number(0)] * 16
        row[0] = beta * z**j
        row[1] = (1-beta) * z**j
        for state, moment in ((1, arbitrary), (2, uniform), (3, kappa*uniform)):
            zero = min(upper(beta*arbitrary), upper(moment))
            nonzero = min(upper((1-beta)*arbitrary), upper(moment))
            row[4*state+2] = zero
            row[4*state] = nonzero / m
            row[4*state+3] = nonzero * (m-1) / m
        result.append(tuple(row))
    return result


def regions(t, s, spectrum, kernel, z, maximum, length=ROWS, reverse=False):
    assert length % t == 0 and 0 <= maximum <= length
    epoch = epochs(t, s, spectrum, kernel, z, min(maximum, t), arb)
    power = tuple(arb_poly([row[k]*math.comb(t, j) for j, row in enumerate(epoch)])
                  for k in range(16))
    current = tuple(arb_poly([int(i == j)]) for i in range(4) for j in range(4))
    def mul(a, b):
        return tuple(sum((a[4*i+k]*b[4*k+j] for k in range(4)), arb_poly()).truncate(maximum+1)
                     for i in range(4) for j in range(4))
    count = length // t
    if reverse:
        for bit in bin(count)[2:]:
            current = mul(current, current)
            if bit == '1':
                current = mul(current, power)
    else:
        while count:
            if count & 1:
                current = mul(current, power)
            count >>= 1
            if count:
                power = mul(power, power)
    return [tuple(max(arb(0), (p[j]/math.comb(length, j)).upper()) for p in current)
            for j in range(maximum+1)]


def run(output, first, last, verify=False):
    output = output.resolve()
    saved = base.read(output) if verify else None
    if not verify:
        assert not output.exists()
    restore_or_verify(base.ROOT)
    if verify:
        for name, digest in saved['source_sha256'].items():
            assert base.sha(base.ROOT/name) == digest
        first, last = saved['occupancy_range']
        assert [r['occupation'] for r in saved['rows']] == list(range(first, last+1))
    assert 8 <= first <= last <= 128
    ctx.prec = 512 if verify else 256
    t, s, spectrum, kernel = maps.load('t64_s20')
    counts = sparse.caps()
    results = []
    for q in range(first, last+1):
        if verify:
            old = saved['rows'][q-first]
            candidates = [(base.decode(old['scaled_tilt']), [base.decode(p) for p in old['p']])]
        else:
            candidates = [(F(q*n, 2), None) for n in (5, 8, 12, 18)]
        best = None
        for a, ps in candidates:
            lam = (arb(a.numerator)/a.denominator)/ROWS
            region = regions(t, s, spectrum, kernel, (-lam).exp(), q, reverse=verify)
            if ps is None:
                approx = np.array([[float(v) for v in m] for m in region]).reshape(q+1, 4, 4)
                ps = sparse.choose_ps(approx, q, counts)
            assert len(ps) == len(sparse.BANDS) and all(0 < p < 1 for p in ps[:-1]) and ps[-1] == 1
            m = sparse.adaptive(region, q, ps, counts, arb)
            for _ in range(8):
                m = matrix.product(m, m)
            upper = math.comb(ROWS, q)*len(sparse.BANDS)**q*rational(
                (sum(m[:4], arb(0))*(CUTOFF*lam).exp()).upper())
            if best is None or upper < best[0]:
                best = upper, a, ps
        bound, a, ps = best
        assert 0 < bound < F(1, 1 << 64), f'Q{q} misses budget'
        if verify:
            assert bound <= base.decode(old['upper'])
        results.append(dict(occupation=q, upper=base.encode(bound), scaled_tilt=base.encode(a),
            p=[base.encode(p) for p in ps], margin_bits=math.log2(bound.denominator)-math.log2(bound.numerator)))
        print('replay' if verify else 'producer', 'Q', q, 'margin', results[-1]['margin_bits'], flush=True)
    if verify:
        assert sum((base.decode(r['upper']) for r in saved['rows']), F(0)) == base.decode(saved['range_upper'])
        base.write_new(output.with_name(output.stem+'_replay.json'), dict(
            status='K30_KERNEL_RANGE_512_BIT_REPLAY_PASSED', producer_sha256=base.sha(output),
            occupancy_range=[first, last], full_distance_proved=False))
    else:
        paths = [Path(__file__), Path(sparse.__file__), Path(matrix.__file__),
            base.HERE/'certify_k30_q1.py', base.HERE/'general_batch_certificate.py',
            base.HERE/'occupation_three.py', base.HERE/'MIGRATION_MANIFEST.json']
        base.write_new(output, dict(status='K30_KERNEL_OUTWARD_CERTIFICATE',
            parameters=dict(message_bits=1 << 30, outer_rows=ROWS, cutoff=CUTOFF, step_bits=t, state_bits=s),
            occupancy_range=[first, last], rows=results,
            range_upper=base.encode(sum((base.decode(r['upper']) for r in results), F(0))),
            full_distance_proved=False, source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in paths}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--first', type=int, default=8)
    parser.add_argument('--last', type=int, default=32)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    run(args.output, args.first, args.last, args.verify)
