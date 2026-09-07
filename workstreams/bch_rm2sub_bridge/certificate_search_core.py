"""Instance identity, exact union budgets, and authenticated legacy coverage."""
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import inner_candidate_screen as inputs

base=inputs.base
LEGACY_PINS={
    'generated/larger_coverage_full.json':'27b579f9f786ad377acebe84db1e7b2efabde88f0444b334ad79144147ccf28e',
    'generated/refresh_q1_outward_v1.json':'889417722fe2d6ef33696da9297e4f61caa63d14db2d0f3b1ad60613ebbe939c',
}


def require(condition,message):
    if not condition: raise ValueError(message)


def outer_dependencies():
    """Authenticate the actual frozen outer inputs, not only their manifest."""
    result={}
    for row in base.read(base.HERE/'MIGRATION_MANIFEST.json')['files']:
        name=row['path']
        if not name.startswith('bch_spectrum_work/'):continue
        path=(base.ROOT/name).resolve()
        require(path.is_relative_to(base.ROOT),'Outer dependency outside repository')
        with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
        require(digest==row['sha256'],'Changed frozen outer dependency: '+name)
        result[name]=digest
    require(bool(result),'No frozen outer dependencies')
    return result


def instance(name,m):
    require(type(m) is int and 13<=m<=20,'Initial driver supports message exponents 13..20')
    outer_dependencies()
    t,s,spectrum,kernel=inputs.load(name)
    length=1<<(m-7)
    require(length>=t and length%t==0,'Complete epochs per region required')
    directory=inputs.maps.DIRECTORY if name in inputs.maps.NAMES else base.HERE/'inputs'
    value=dict(configuration=name,message_exponent=m,rows=length,cutoff=256*length//10,
        step_bits=t,state_bits=s,outer='fixed-p37-syndrome-subspace-0..31-GF256-0x14D',
        outer_bits=256,outer_dimension=128,minimum_outer_weight=38,
        setup='independent-row-and-region-permutations;fresh-nonzero-scalars;zero-start;output-before-update;no-flush',
        selection_sha256=base.sha(directory/f'{name}_selection.json'),
        outer_inputs_sha256=base.sha(base.HERE/'MIGRATION_MANIFEST.json'))
    value['fingerprint']=hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()
    return value


def authenticate(data,root,key='source_sha256'):
    for name,digest in data[key].items():
        require(base.sha(root/Path(name))==digest,'Changed dependency: '+name)


def intervals(points):
    result=[]
    for q in sorted(set(points)):
        if result and result[-1][1]+1==q: result[-1][1]=q
        else: result.append([q,q])
    return result


def summarize(best,rows,bits):
    require(type(bits) is int and 1<=bits<=256,'Invalid target margin')
    require(all(type(q) is int and 1<=q<=rows and v>0 for q,v in best.items()),'Invalid occupancy bounds')
    total=sum(best.values(),F(0));target=F(1,1<<bits)
    complete=len(best)==rows
    return dict(status='CERTIFIED' if complete and total<=target else 'UNRESOLVED',
        complete_coverage=complete,covered_occupancies=len(best),covered_intervals=intervals(best),
        uncovered_intervals=intervals(q for q in range(1,rows+1) if q not in best),
        covered_union_upper=base.encode(total),target_upper=base.encode(target),
        covered_margin_bits=math.log2(total.denominator)-math.log2(total.numerator) if total else None,
        covered_margin_is_full_margin=complete,remaining_budget=base.encode(target-total))


def allocation(best,rows,bits):
    """Equal remaining per-occupancy budget rounded DOWN to a dyadic value."""
    remaining=F(1,1<<bits)-sum(best.values(),F(0))
    count=rows-len(best)
    if remaining<=0 or count<=0: return F(0)
    value=remaining/count
    exponent=value.numerator.bit_length()-value.denominator.bit_length()
    candidate=F(2)**exponent
    return candidate if candidate<=value else candidate/2


def merge_rows(best,values,rows):
    seen=set()
    for row in values:
        q=row['occupation'];upper=base.decode(row['upper'])
        require(type(q) is int and 1<=q<=rows and q not in seen and upper>0,'Invalid or duplicate bound')
        seen.add(q);best[q]=min(best.get(q,upper),upper)


def legacy_baseline(spec):
    """Authenticate and reconstruct prior coverage, without rerunning producers."""
    if spec['configuration']!='t64_s20' or spec['message_exponent']!=20: return None
    require(spec==instance('t64_s20',20),'Wrong baseline instance')
    for name,digest in LEGACY_PINS.items(): require(base.sha(base.HERE/name)==digest,'Legacy pin mismatch')
    ledger=base.read(base.HERE/'generated/larger_coverage_full.json')
    authenticate(ledger,base.HERE,'local_sha256')
    old_q1=base.read(base.HERE/'generated/larger_t64_s20_q1_outward.json')
    authenticate(old_q1,base.HERE,'local_sha256')
    best={1:base.decode(old_q1['Q1_upper'])}
    for name in ledger['local_sha256']:
        path=base.HERE/Path(name)
        if not (path.name.startswith('larger_range_') and path.name.endswith('_outward.json')): continue
        saved=base.read(path);authenticate(saved,base.HERE,'local_sha256')
        replay=base.read(path.with_name(path.name.replace('_outward.json','_replay.json')))
        require(saved['configuration']=='t64_s20','Wrong legacy map')
        require(replay['status']=='LARGER_STATE_512_BIT_RANGE_REPLAY_PASSED' and
                replay['producer_sha256']==base.sha(path),'Missing legacy replay')
        lo,hi=saved['occupancy_range']
        require(len(saved['rows'])==hi-lo+1==replay['rows_checked'],'Legacy coverage length mismatch')
        for q,row in zip(range(lo,hi+1),saved['rows']):
            require(q==row['occupation'],'Legacy occupancy mismatch')
            upper=base.decode(row['upper'])
            require(upper>0,'Nonpositive legacy bound')
            if upper<F(1,1<<60): best[q]=min(best.get(q,upper),upper)
    require(set(best)==set(range(1,8193)),'Legacy coverage gap')
    require(sum(best.values(),F(0))==base.decode(ledger['covered_union_upper']),'Legacy union mismatch')
    improved=base.read(base.HERE/'generated/refresh_q1_outward_v1.json')
    authenticate(improved,base.ROOT)
    replay=base.read(base.HERE/'generated/refresh_q1_outward_v1_replay.json')
    require(replay['producer_sha256']==LEGACY_PINS['generated/refresh_q1_outward_v1.json'] and
            replay['status']=='REFRESH_Q1_512_BIT_LINEAR_EPOCH_REPLAY_PASSED','Missing refresh replay')
    co={int(w):base.decode(v) for w,v in improved['coefficient_upper'].items()}
    require(base.bch_bound(co)==tuple(base.decode(improved[k]) for k in ('Q1_upper','factor','rest')),
            'Legacy BCH aggregation mismatch')
    require(improved['message_bits']==1<<20 and improved['cutoff']>=spec['cutoff'],'Wrong legacy cutoff')
    best[1]=base.decode(improved['Q1_upper'])
    require(sum(best.values(),F(0))==base.decode(improved['combined_full_upper']),'Refresh union mismatch')
    return best
