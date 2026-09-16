"""Authenticate replays and aggregate partial coverage, never imply full closure."""
from fractions import Fraction as F
import math
from pathlib import Path
import bch_model as model
import certify_dense as dense
import christoffel_caps


def receipt(name, status):
    path = model.HERE/name
    record = model.base.read(path)
    replay = model.base.read(path.with_name(path.stem+'_replay.json'))
    model.authenticate(record)
    assert replay['status'] == status and replay['producer_sha256'] == model.base.sha(path)
    return record


def outer_caps():
    counts = {w:c for w,c in christoffel_caps.deterministic_caps().items() if w}
    first = model.base.read(model.base.HERE/'generated/christoffel_oa29_caps.json')
    assert counts == {r['weight']:r['cap'] for r in first['rows']}
    dependencies = {}
    for path in [model.base.HERE/'generated/christoffel_oa29_caps.json']+sorted((model.base.HERE/'generated').glob('joint_shell_*/cap.json')):
        record = model.base.read(path)
        for root, key in ((model.base.HERE, 'local_sha256'), (model.base.BCH, 'outer_sha256')):
            for name, digest in record[key].items():
                dependency = (root/name).resolve()
                assert model.base.sha(dependency) == digest, dependency
                dependencies[dependency.relative_to(model.ROOT).as_posix()] = digest
        if 'weight' in record:
            assert record['rational_primal_dual_checks_passed']
            w = record['weight']
            assert record['cap'] == math.floor(model.base.decode(record['physical_upper']))
            counts[w] = counts[256-w] = min(counts[w], record['cap'])
    assert counts == model.outer.inputs.caps_module.caps()
    return dependencies


def main():
    dependencies = outer_caps()
    q1 = receipt('Q1_CERTIFICATE.json', 'BCH256_ASYMMETRIC_Q1_512_BIT_LINEAR_REPLAY_PASSED')
    sparse = receipt('SPARSE_Q64.json', 'BCH256_ASYMMETRIC_SPARSE_512_BIT_REPLAY_PASSED')
    assert sparse['occupations'] == list(range(2, 65))
    results = []
    for exponent, q1row in zip((16,18,20), q1['results']):
        engine = model.Engine(exponent)
        assert q1row['instance'] == engine.identity()
        coefficients = {int(w):model.exact(model.unpack(v)) for w, v in q1row['coefficient_upper'].items()}
        assert set(coefficients) == set(model.base.WEIGHTS)
        upper, factor, rest = model.base.bch_bound(coefficients)
        assert (upper, factor, rest) == tuple(model.base.decode(q1row[key]) for key in ('upper','factor','rest'))
        best = {1:upper}
        selected = [r for r in sparse['results'] if r['instance']['message_bits'] == 1 << exponent]
        assert [r['occupation'] for r in selected] == list(range(2, 65))
        for row in selected:
            assert row['instance'] == engine.identity()
            best[row['occupation']] = model.exact(model.unpack(row['upper']))
        bank = receipt(f'SPARSE_M{exponent}.json', 'ASYMMETRIC_SPARSE_RANGE_512_BIT_REPLAY_PASSED')
        assert bank['instance'] == engine.identity()
        reconstructed = {}
        for witness in bank['witnesses']:
            lo, hi = witness['interval']
            assert 2 <= lo <= hi <= bank['last']
            assert [r['occupation'] for r in witness['rows']] == list(range(lo, hi+1))
            for row in witness['rows']:
                q, power = row['occupation'], row['power']
                assert type(power) is int
                if power <= -60:
                    reconstructed[q] = min(reconstructed.get(q, power), power)
        assert reconstructed == {int(q):v for q,v in bank['best_powers'].items()}
        assert sum((F(2)**v for v in reconstructed.values()), F(0)) == model.base.decode(bank['upper'])
        for q, power in reconstructed.items():
            best[q] = min(best.get(q, F(2)**power), F(2)**power)
        last = max(best)
        assert set(best) == set(range(1, last+1))
        total = sum(best.values(), F(0))
        assert total < F(1, 1 << 40)
        results.append(dict(message_exponent=exponent, covered_occupancies=[1,last], outer_rows=engine.length,
                            covered_union_upper=model.base.encode(total),
                            covered_union_margin_bits=math.log2(total.denominator)-math.log2(total.numerator),
                            q1_margin_bits=q1row['margin_bits'], full_distance_proved=False))
    checkpoint = receipt('DENSE_ADAPTIVE_M20.json', 'ASYMMETRIC_ADAPTIVE_DENSE_512_BIT_REPLAY_PASSED')
    dense.geometry.check_partition(checkpoint['leaves'], checkpoint['splits'], 512, 8192)
    assert checkpoint['unresolved'] == sum(v['power'] > -80 for v in checkpoint['leaves'].values())
    result = dict(status='VERIFIED_PARTIAL_BCH256_ASYMMETRIC_PROGRESS', results=results,
                  dense_m20=dict(leaves=len(checkpoint['leaves']), unresolved=checkpoint['unresolved']),
                  full_distance_proved=False, source_sha256={**model.sources(), **dependencies})
    for name in ('Q1_CERTIFICATE','SPARSE_Q64','SPARSE_M16','SPARSE_M18','SPARSE_M20','DENSE_ADAPTIVE_M20'):
        for suffix in ('.json','_replay.json'):
            path = model.HERE/(name+suffix)
            result['source_sha256'][path.relative_to(model.ROOT).as_posix()] = model.base.sha(path)
    path = model.HERE/'PROGRESS_VERIFIED.json'
    # This summary is regenerated from immutable numerical receipts.
    import json
    path.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k != 'source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
