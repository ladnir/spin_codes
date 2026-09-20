"""Check the exact coverage ledger and source hashes through occupancy 128."""
from fractions import Fraction as F
import bridge as base


def verify():
    previous=base.read(base.HERE/'generated/adaptive_q17_q64_outward.json')
    extension=base.read(base.HERE/'generated/tightened_q65_q128_extension_outward.json')
    for receipt in (previous,extension):
        for name,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
        for name,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
    assert extension['parameters']==previous['parameters']
    assert [row['occupation'] for row in extension['rows']]==list(range(65,129))
    total=sum((base.decode(row['upper']) for row in extension['rows']),F(0))
    assert total==base.decode(extension['range_upper']) and 0<total<F(1,1<<148)
    partial=base.decode(previous['Q1_through_Q64_upper'])+total
    assert 0<partial<F(1,1<<49)
    print('Exact coverage ledger Q1..128 <2^-49; Q129..8192 remains open.',flush=True)


if __name__=='__main__':verify()
