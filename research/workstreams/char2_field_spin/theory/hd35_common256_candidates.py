"""Directed audit of transferred search hints for RS256/common-h2.

The outer has twice as many regions and half as many rows per region.
Halving the active-row count preserves potential input density. This only
proposes rational witnesses: every bound is recomputed for the new maps.
Rejected witnesses and uncovered classes remain explicit in the result.
"""
from hd35_common128_witnesses import WITNESSES as SOURCE
from hd35_common_certificate import replay


def candidates():
    result=[]
    for kind,values in SOURCE:
        if kind=='paired':
            first,last,*parameters=values
            first,last=(first+1)//2,last//2
            if first<=last:result.append([kind,[first,last]+parameters])
        else:
            center,*parameters=values
            # Keep paired composition for the cheaper feedback's harder
            # middle classes. Dense classes can use the single-packet bound.
            choice='global' if center<6500 else 'global-single'
            result.append([choice,[(center+1)//2]+parameters])
    return result


if __name__=='__main__':
    import argparse
    import json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--precision',type=int,default=256)
    p.add_argument('--workers',type=int,default=4)
    args=p.parse_args()
    witnesses=candidates()
    def progress(item):
        item=dict(item)
        receipt=item.pop('receipt')
        item.update(intervals=receipt.get('certified_intervals'),
                    bound=receipt.get('union_log2_upper',receipt.get('log2_upper')))
        print(json.dumps(item),flush=True)
    result=replay(witnesses,args.precision,progress,outer_length=256,
                  local_rank=2,workers=args.workers,allow_rejected=True)
    intervals=[]
    for r in result.pop('uncovered_rows'):
        if intervals and r==intervals[-1][1]+1:intervals[-1][1]=r
        else:intervals.append([r,r])
    result['uncovered_intervals']=intervals
    print(json.dumps(result),flush=True)
