"""Replay a scaled partition, replacing every inner moment and retained bound.

The old tree supplies geometry and candidate witnesses only. Its numerical
bounds and old-inner claims are never imported.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
from flint import arb, ctx
import dense_search as search
import search_ladder_dense_tilt as geometry
model = search.model

SEED = model.base.HERE/'generated/t128_s19_m20_ladder_fixed_reference_v2/cover_0000.json'


def mapped_leaves(seed, exponent, minimum):
    geometry.check_partition(seed['leaves'], seed['splits'], 512, 8192)
    scale = 1 << (20-exponent)
    assert minimum >= (512+scale-1)//scale
    output = []
    for key, node in seed['leaves'].items():
        lo, hi, a, b = geometry.geometry(node)
        lo, hi = max(minimum, (lo+scale-1)//scale), hi//scale
        if lo <= hi:
            output.append((key, lo, hi, a, b, node['witness']))
    # Every new integer Q maps to the old integer scale*Q. Closed density
    # boundaries may overlap, but no occupancy or density interval is omitted.
    for q in range(minimum, (1 << (exponent-7))+1):
        intervals = sorted((a, b) for _, lo, hi, a, b, _ in output if lo <= q <= hi)
        end = F(0)
        assert intervals and intervals[0][0] == 0
        for a, b in intervals:
            assert a <= end and a < b
            end = max(end, b)
        assert end == 1
    return output


def power_bound(value):
    bits = model.exact((value/arb(2).log()).upper())
    return max(-200, -(-bits.numerator//bits.denominator))


def run(output, exponent, minimum, retunes, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        exponent, minimum = saved['message_exponent'], saved['minimum']
        assert saved['seed_sha256'] == model.base.sha(SEED)
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    seed = model.base.read(SEED)
    checker = search.Checker(exponent, [model.base.decode(p) for p in seed['probabilities']])
    mapped = mapped_leaves(seed, exponent, minimum)
    results = []
    used = 0
    if saved:
        assert len(saved['leaves']) == len(mapped) and saved['instance'] == checker.engine.identity()
    for i, (key, lo, hi, a, b, witness) in enumerate(mapped):
        old = saved['leaves'][i] if saved else None
        if old:
            assert old['key'] == key and old['geometry'] == [lo, hi, str(a), str(b)]
            witness = old['witness']
        value = checker.bound(lo, hi, a, b, witness)
        power = power_bound(value)
        if not saved and power > -80 and used < retunes:
            alternative = checker.witness(lo, hi, a, b)
            proposed = power_bound(checker.bound(lo, hi, a, b, alternative))
            used += 1
            if proposed < power:
                power, witness = proposed, alternative
        if old:
            assert power <= old['power'], ('replay', key, power, old['power'])
            power = old['power']
        results.append(dict(key=key, geometry=[lo, hi, str(a), str(b)], witness=witness, power=power))
        if i%100 == 0:
            print(ctx.prec, 'dense', exponent, i+1, '/', len(mapped),
                  'weak leaves', sum(r['power'] > -80 for r in results), flush=True)
    # Very weak leaves can have million-bit positive exponents. Keep their
    # powers for diagnosis without allocating a huge irrelevant rational sum.
    unresolved = sum(r['power'] > -80 for r in results)
    total = sum((F(2)**r['power'] for r in results), F(0)) if not unresolved else None
    if saved:
        assert unresolved == saved['unresolved']
        if total is not None:
            assert total == model.base.decode(saved['upper'])
        model.base.write_new(output.with_name(output.stem+'_replay.json'),
                             dict(status='BCH256_ASYMMETRIC_DENSE_512_BIT_REPLAY_PASSED',
                                  producer_sha256=model.base.sha(output), unresolved=unresolved,
                                  complete_dense_certificate=not unresolved, full_distance_proved=False))
    else:
        sources = model.sources()
        sources[SEED.relative_to(model.ROOT).as_posix()] = model.base.sha(SEED)
        result = dict(status='OUTWARD_BCH256_ASYMMETRIC_DENSE_COVER' if not unresolved else 'OUTWARD_DENSE_COVER_UNRESOLVED',
                      message_exponent=exponent, minimum=minimum, instance=checker.engine.identity(),
                      seed_sha256=model.base.sha(SEED), leaves=results, unresolved=unresolved,
                      retunes=used, full_distance_proved=False, source_sha256=sources)
        if total is not None:
            result.update(upper=model.base.encode(total), margin_bits=math.log2(total.denominator)-math.log2(total.numerator))
        model.base.write_new(output, result)
    print('dense completed', exponent, 'leaves', len(results), 'unresolved', unresolved, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=(16, 18, 20), default=20)
    p.add_argument('--minimum', type=int, default=512)
    p.add_argument('--retunes', type=int, default=16)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    assert 1 <= a.minimum <= 1 << (a.m-7) and a.retunes >= 0
    run(a.output.resolve(), a.m, a.minimum, a.retunes, a.verify)
