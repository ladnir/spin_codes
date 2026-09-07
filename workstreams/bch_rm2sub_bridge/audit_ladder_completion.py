"""All-rung completion audit with fresh 768-bit dense evaluation.

Reconstructs sparse coverage from authenticated numerical receipts and
recomputes every dense leaf. This is not a new witness search.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

from flint import ctx

import ladder_instance as identity
import search_ladder_dense_fixed as dense
import poisson_density_factor as comparison

core,base = identity.core,identity.core.base
PRIOR = {
    't128_s19_m16_full_split_coverage_v1.json':'3684d355322ee59d9cffdeb3ac5edd29ff3bdc22ca37b73fb723a4cbe95b0b18',
    't128_s19_m18_full_split_coverage_v1.json':'acac24114eb25d658ffcd86ee3a7fe2c012f0c156e0fa73fb50cbfa4f5819550',
    't128_s19_m20_ladder_full_v1.json':'657198cb67f805ec5fda4c9caa48e95063b67090badbb1eea5a38d6343ff03cd',
}


def full_fields(saved,m,sparse_values,dense_value):
    spec = identity.instance(m)
    core.require(saved['instance'] == spec,'Changed instance or selected map')
    core.require(saved['status'] == 'LADDER_FULL_DISTANCE_SETUP_CERTIFIED_40_BITS' and
                 saved['complete_coverage'] and saved['covered_intervals'] == [[1,spec['rows']]],'Not a full certificate')
    first,last = saved['dense_interval']
    core.require(last == spec['rows'] and saved['sparse_interval'] == [1,first-1] and
                 set(sparse_values) == set(range(1,first)),'Occupancy coverage gap or overlap')
    core.require(all(v > 0 for v in sparse_values.values()) and dense_value > 0,'Nonpositive bound')
    total = sum(sparse_values.values(),F(0))+dense_value
    core.require(total == base.decode(saved['union_upper']) and total <= F(2)**-40,'Wrong full union or missed target')
    core.require(sparse_values[1] == base.decode(saved['q1_upper']) and
                 sum((v for q,v in sparse_values.items() if q != 1),F(0)) == base.decode(saved['sparse_rest_upper']) and
                 dense_value == base.decode(saved['dense_upper']),'Wrong component union')
    return total


def run(output):
    core.require(not output.exists(),'Use a fresh audit output')
    hashes = {}
    def digest(path):
        path = path.resolve()
        core.require(path.is_relative_to(base.ROOT),'Dependency outside repository')
        name = path.relative_to(base.ROOT).as_posix()
        if name not in hashes:
            hashes[name] = base.sha(path)
        return hashes[name]
    def authenticated(path):
        value = base.read(path)
        digest(path)
        for name,expected in value['source_sha256'].items():
            core.require(digest(base.ROOT/name) == expected,'Changed pinned dependency: '+name)
        return value
    for name,expected in PRIOR.items():
        path = base.HERE/'generated'/name
        core.require(digest(path) == expected,'Earlier certificate changed: '+name)
        authenticated(path)
    rows = []
    common = None
    for m in (20,22,24):
        ledger_path = base.HERE/'generated'/f't128_s19_m{m}_ladder_full_v1.json'
        ledger = authenticated(ledger_path)
        spec = identity.instance(m)
        core.require(ledger['instance'] == spec,'Cross-size ledger')
        shared = {k:v for k,v in spec.items() if k not in ('message_exponent','rows','cutoff','fingerprint')}
        if common is None:
            common = shared
        core.require(common == shared,'Construction or setup changed across ladder')
        reports = [base.ROOT/name for name in ledger['source_sha256']
                   if Path(name).name.startswith('report_') and Path(name).suffix == '.json']
        core.require(len(reports) == 1,'Ambiguous sparse report')
        report_path = reports[0]
        report = base.read(report_path)
        core.require(report['instance'] == spec and not report['reused_legacy_certificate'],'Wrong sparse instance or reused old-size bound')
        sparse_values = {}
        max_degree = 0
        for item in report['accepted_certificates']:
            cp,rp = report_path.parent/item['certificate'],report_path.parent/item['replay']
            core.require(cp.resolve().is_relative_to(report_path.parent) and rp.resolve().is_relative_to(report_path.parent),'Escaped sparse receipt path')
            cert,replay = authenticated(cp),authenticated(rp)
            task,result = cert['task'],cert['result']
            core.require(task['instance'] == spec,'Cross-size sparse producer')
            core.require(replay['status'] == 'NUMERICAL_512_BIT_REPLAY_PASSED' and
                         replay['certificate_sha256'] == digest(cp) and replay['instance_fingerprint'] == spec['fingerprint'],
                         'Missing matching sparse replay')
            expected = [1] if task['kind'] == 'q1' else list(range(task['interval'][0],task['interval'][1]+1))
            core.require([r['occupation'] for r in result['bounds']] == expected,'Sparse producer gap or duplicate')
            selected = item['occupancies']
            core.require(len(selected) == len(set(selected)) and set(selected) <= set(expected),'Invalid sparse selection')
            if task['kind'] == 'q1':
                coefficients = {int(w):base.decode(v) for w,v in result['coefficient_upper'].items()}
                core.require(set(coefficients) == set(base.WEIGHTS),'Incomplete Q1 spectrum coefficients')
                core.require(base.bch_bound(coefficients)[0] == base.decode(result['bounds'][0]['upper']),'Wrong exact Q1 BCH aggregation')
            else:
                core.require(task['kind'] == 'range','Unknown sparse worker')
                max_degree = max(max_degree,task['interval'][1])
            core.merge_rows(sparse_values,[r for r in result['bounds'] if r['occupation'] in selected],spec['rows'])
        core.require(core.summarize(sparse_values,spec['rows'],report['target_bits']) == report['coverage'],'Sparse report union mismatch')
        refined_path = base.HERE/'generated'/f't128_s19_m{m}_ladder_dense_refined_v1.json'
        refined = authenticated(refined_path)
        exact_replay = authenticated(refined_path.with_name(refined_path.stem+'_replay.json'))
        core.require(refined['instance'] == spec and refined['status'] == 'REFINED_DENSITY_DENSE_RANGE_CERTIFIED','Wrong dense instance or status')
        core.require(exact_replay['status'] == 'REFINED_DENSITY_EXACT_REPLAY_PASSED' and exact_replay['producer_sha256'] == digest(refined_path),
                     'Missing exact refinement replay')
        source_path = base.ROOT/refined['source_cover']
        source = authenticated(source_path)
        source_replay_path = source_path.with_name('replay_'+source_path.stem+'.json')
        source_replay = authenticated(source_replay_path)
        core.require(source['instance'] == spec and source['policy'] == dense.POLICY and
                     digest(source_path) == refined['source_cover_sha256'],'Changed dense cover')
        core.require(source_replay['status'] == 'TWO_TILT_512_BIT_REPLAY_PASSED' and source_replay['producer_sha256'] == digest(source_path) and
                     exact_replay['source_cover_replay_sha256'] == digest(source_replay_path),'Missing dense numerical replay')
        coverage = dense.old.check_partition(source['leaves'],source['splits'],source['minimum'],spec['rows'])
        core.require(coverage == refined['coverage'] and refined['interval'] == [source['minimum'],spec['rows']],
                     'Wrong dense domain or partition')
        ctx.prec = 768
        checker = dense.dense.Checker(spec,[base.decode(p) for p in source['probabilities']])
        powers = []
        for index,node in enumerate(source['leaves'].values()):
            power = node['power']
            core.require(type(power) is int and -200 <= power <= 100000,'Invalid dyadic bound')
            core.require(dense.old.retained_power(checker,node) <= power,'Fresh 768-bit dense check failed')
            powers.append(power)
            if index % 250 == 0:
                print('m',m,'dense 768-bit leaf',index+1,'/',len(source['leaves']),flush=True)
        old_sum = sum((F(2)**p for p in powers),F(0))
        D = comparison.density_factor(spec['rows'])
        ratio = F(D,spec['rows']+1)**256
        dense_value = old_sum*ratio
        core.require(refined['density_factor'] == D and refined['reference_density_factor'] == spec['rows']+1 and
                     old_sum == base.decode(refined['source_union_upper']) and ratio == base.decode(refined['exact_ratio']) and
                     dense_value == base.decode(refined['union_upper']),'Wrong density-factor conversion')
        total = full_fields(ledger,m,sparse_values,dense_value)
        rows.append(dict(message_exponent=m,rows=spec['rows'],cutoff=spec['cutoff'],instance_fingerprint=spec['fingerprint'],
                         ledger_sha256=digest(ledger_path),union_upper=base.encode(total),
                         margin_bits=math.log2(total.denominator)-math.log2(total.numerator),
                         sparse_occupancies=len(sparse_values),sparse_maximum_polynomial_degree=max_degree,
                         dense_rectangles=len(powers),dense_precision_bits=768,all_occupancies_covered=True))
        print('AUDITED full m',m,'margin',rows[-1]['margin_bits'],flush=True)
    # Pin the small completion records, not the mutable global status document.
    for name in ('T128_S19_LADDER_TO_M24.md','T128_S19_M20_CLOSURE.md','T128_S19_M22_M24_CLOSURE.md'):
        digest(base.HERE/name)
    hashes.update(dense.old.provenance.sources())
    result = dict(status='FULL_M20_M22_M24_LADDER_COMPLETION_AUDIT_PASSED',target_bits=40,
                  instances=rows,common_construction=common,preserved_prior_sha256=PRIOR,
                  fresh_dense_precision_bits=768,sparse_replays='authenticated existing 512-bit numerical replays',
                  source_sha256=hashes)
    base.write_new(output,result)
    print(result['status'],flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    run(args.output)
