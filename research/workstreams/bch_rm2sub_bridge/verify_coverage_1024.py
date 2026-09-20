"""Exact aggregation of the overlapping outward certificates through Q1024.

This ledger verifies sums and source identities. Separate --verify commands
recompute the numerical certificates; the ledger does not replace them.
"""
import math
from fractions import Fraction as F
import bridge as base
from verify_coverage_128 import verify as verify_previous


def verify():
    verify_previous()
    receipts=[base.read(base.HERE/'generated'/name) for name in
              ('oa29_shared_q129_q1024_outward.json','oa29_gap339_outward.json')]
    parameters=base.read(base.HERE/'generated/tightened_q65_q128_extension_outward.json')['parameters']
    bounds={}
    for receipt in receipts:
        assert receipt['parameters']==parameters and receipt['configuration']=='t128_s15'
        for name,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
        for name,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
        lo,hi=receipt['occupancy_range']
        assert [r['occupation'] for r in receipt['rows']]==list(range(lo,hi+1))
        subtotal=F(0)
        for row in receipt['rows']:
            q=row['occupation'];value=base.decode(row['upper']);assert value>0
            subtotal+=value
            bounds[q]=min(bounds.get(q,value),value)
        assert subtotal==base.decode(receipt['range_upper'])
    assert set(bounds)==set(range(129,1025))
    total=sum(bounds.values(),F(0));assert total<F(1,1<<77)
    partial=(base.decode(base.read(base.HERE/'generated/adaptive_q17_q64_outward.json')['Q1_through_Q64_upper'])
             +base.decode(base.read(base.HERE/'generated/tightened_q65_q128_extension_outward.json')['range_upper'])+total)
    assert partial<F(1,1<<49)
    print('Exact coverage Q1..1024 <2^-49; added Q129..1024 <2^-77; diagnostic margin',
          math.log2(total.denominator)-math.log2(total.numerator),flush=True)
    print('Q1025..8192 remains open. This is not the full distance theorem.',flush=True)


if __name__=='__main__':verify()
