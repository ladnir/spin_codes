"""Exact partial coverage ledger; not a replacement for numerical replays."""
from fractions import Fraction as F
import bridge as base
from verify_coverage_1024 import verify as previous


def verify():
    previous()
    floor=base.read(base.HERE/'generated/region_floor_d1200_h513_outward.json')
    endpoint=base.read(base.HERE/'generated/constant_split_ordinary8192_outward.json')
    assert floor['status']=='OUTWARD_DIRECT_REGION_FLOOR_CLASS'
    assert endpoint['status']=='OUTWARD_LISTED_CONSTANT_SPLIT_CELLS'
    for receipt in (floor,endpoint):
        assert receipt['configuration']=='t128_s15'
        assert receipt['parameters']==dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,
            outer_length=256,step_bits=128,state_bits=15,distance_cutoff=209716)
        assert receipt['all_occupations_certified'] is False
        for name,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    assert floor['scope']==dict(ordinary_rows_maximum=1200,all_one_rows_minimum=513)
    assert len(endpoint['rows'])==1
    row=endpoint['rows'][0]
    assert (row['ordinary_rows'],row['all_one_rows'],row['occupation'])==(8192,0,8192)
    a=base.decode(floor['failure_contribution_upper']);b=base.decode(row['upper'])
    assert 0<a<F(1,1<<5390) and 0<b<F(1,1<<9382)
    assert b==base.decode(endpoint['listed_cells_upper'])
    shared=[base.read(base.HERE/'generated'/n) for n in
        ('oa29_shared_q129_q1024_outward.json','oa29_gap339_outward.json')]
    bounds={}
    for receipt in shared:
        for row in receipt['rows']:
            q=row['occupation'];v=base.decode(row['upper']);bounds[q]=min(bounds.get(q,v),v)
    low=(base.decode(base.read(base.HERE/'generated/adaptive_q17_q64_outward.json')['Q1_through_Q64_upper'])
        +base.decode(base.read(base.HERE/'generated/tightened_q65_q128_extension_outward.json')['range_upper'])
        +sum(bounds.values(),F(0)))
    extension=base.read(base.HERE/'generated/exactcap_shared_q1025_q1280_outward.json')
    assert extension['status']=='OUTWARD_EXACTCAP_SHARED_RANGE' and extension['configuration']=='t128_s15'
    assert extension['parameters']==floor['parameters'] and extension['occupancy_range']==[1025,1280]
    for name,digest in extension['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    assert [r['occupation'] for r in extension['rows']]==list(range(1025,1281))
    assert sum((base.decode(r['upper']) for r in extension['rows']),F(0))==base.decode(extension['range_upper'])
    added=sum((base.decode(r['upper']) for r in extension['rows'] if r['occupation']<=1251),F(0))
    assert 0<added<F(1,1<<64)
    extension2=base.read(base.HERE/'generated/full_exactcaps_q1252_q1792_outward.json')
    assert extension2['status']=='OUTWARD_FULL_EXACTCAPS_RANGE' and extension2['configuration']=='t128_s15'
    assert extension2['parameters']==floor['parameters'] and extension2['occupancy_range']==[1252,1792]
    for name,digest in extension2['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    assert [r['occupation'] for r in extension2['rows']]==list(range(1252,1793))
    assert sum((base.decode(r['upper']) for r in extension2['rows']),F(0))==base.decode(extension2['range_upper'])
    added2=sum((base.decode(r['upper']) for r in extension2['rows'] if r['occupation']<=1655),F(0))
    assert 0<added2<F(1,1<<79)
    assert low+added+added2+a+b<F(1,1<<49)
    # Count compositions, not messages; counts are diagnostic only.
    total=8193*8194//2-1
    covered=0
    for d in range(8193):
        low_hi=min(8192-d,1655-d)
        low_count=max(0,low_hi+1)-(1 if d==0 else 0)
        tail_count=max(0,8192-d-max(513,low_hi+1)+1) if d<=1200 else 0
        covered+=low_count+tail_count+(1 if d==8192 else 0)
    assert covered<total
    print('Verified Q1..1655 plus d<=1200,h>=513 and (8192,0): partial union <2^-49.')
    print('Q1025..1251 subtotal <2^-64; Q1252..1655 subtotal <2^-79.')
    print('Nonpassing rows Q1656..1792 are retained but excluded.')
    print('Composition cells covered:',covered,'of',total,'; open:',total-covered)
    print('Full occupancy range NOT closed. Numerical replay commands remain separate.')


if __name__=='__main__':verify()
