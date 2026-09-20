"""Exact full occupancy union from sparse receipts and a refined dense cover."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import refine_ladder_density_cover as refined
import verify_certificate_search as sparse

base,core = refined.base,refined.core


def full_union(rows,values,interval,dense_upper):
    first,last = interval
    core.require(2 <= first <= last == rows and set(values) == set(range(1,first)),
                 'Sparse/dense occupancy gap, overlap, or wrong endpoint')
    core.require(dense_upper > 0 and all(v > 0 for v in values.values()),'Invalid component upper bound')
    return sum(values.values(),F(0))+dense_upper


def run(report_path,dense_path,output):
    core.require(not output.exists(),'Use fresh ledger path')
    sparse.verify_report(report_path)
    report = base.read(report_path)
    spec = report['instance']
    density = base.read(dense_path)
    replay_path = dense_path.with_name(dense_path.stem+'_replay.json')
    replay = base.read(replay_path)
    core.authenticate(density,base.ROOT)
    core.authenticate(replay,base.ROOT)
    core.require(density['status'] == 'REFINED_DENSITY_DENSE_RANGE_CERTIFIED' and density['instance'] == spec,
                 'Wrong dense instance or status')
    core.require(replay['status'] == 'REFINED_DENSITY_EXACT_REPLAY_PASSED' and
                 replay['producer_sha256'] == base.sha(dense_path),'Missing exact dense replay')
    source_path = base.ROOT/density['source_cover']
    source = base.read(source_path)
    source_replay_path = source_path.with_name('replay_'+source_path.stem+'.json')
    source_replay = base.read(source_replay_path)
    core.authenticate(source,base.ROOT)
    core.authenticate(source_replay,base.ROOT)
    core.require(base.sha(source_path) == density['source_cover_sha256'] and source['instance'] == spec and
                 source['policy'] == refined.search.POLICY,'Wrong dense source geometry or policy')
    core.require(source_replay['status'] == 'TWO_TILT_512_BIT_REPLAY_PASSED' and
                 source_replay['producer_sha256'] == base.sha(source_path) and
                 replay['source_cover_replay_sha256'] == base.sha(source_replay_path),'Missing dense numerical replay')
    first = density['interval'][0]
    core.require(first == source['minimum'],'Changed dense lower endpoint')
    coverage = refined.search.old.check_partition(source['leaves'],source['splits'],first,spec['rows'])
    core.require(coverage == density['coverage'],'Wrong dense coverage record')
    old_sum = sum((F(2)**node['power'] for node in source['leaves'].values()),F(0))
    ratio = F(refined.factor.density_factor(spec['rows']),spec['rows']+1)**256
    dense_upper = old_sum*ratio
    core.require(dense_upper == base.decode(density['union_upper']),'Wrong dense exact union')
    values,pins = {},refined.search.old.provenance.sources()
    def pin(path):
        pins[path.resolve().relative_to(base.ROOT).as_posix()] = base.sha(path)
    for item in report['accepted_certificates']:
        path = report_path.parent/item['certificate']
        receipt = report_path.parent/item['replay']
        entries = sparse.load_verified(path,receipt,spec)
        for entry in entries:
            q = entry['occupation']
            if q in item['occupancies'] and q < first:
                value = base.decode(entry['upper'])
                values[q] = min(values.get(q,value),value)
        pin(path)
        pin(receipt)
    total = full_union(spec['rows'],values,density['interval'],dense_upper)
    core.require(total <= F(2)**-40,'Full union does not reach 40 bits')
    for path in (report_path,dense_path,replay_path,source_path,source_replay_path):
        pin(path)
    result = dict(status='LADDER_FULL_DISTANCE_SETUP_CERTIFIED_40_BITS',instance=spec,
                  covered_intervals=[[1,spec['rows']]],complete_coverage=True,union_upper=base.encode(total),
                  margin_bits=math.log2(total.denominator)-math.log2(total.numerator),
                  q1_upper=base.encode(values[1]),sparse_rest_upper=base.encode(sum((v for q,v in values.items() if q != 1),F(0))),
                  sparse_interval=[1,first-1],dense_interval=density['interval'],dense_upper=base.encode(dense_upper),
                  source_sha256=pins,numerical_replays_rerun_by_this_ledger=False)
    base.write_new(output,result)
    print(result['status'],'m',spec['message_exponent'],'margin',result['margin_bits'],flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sparse-report',type=Path,required=True)
    p.add_argument('--dense',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    run(a.sparse_report,a.dense,a.output)
