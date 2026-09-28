"""Refine the largest support-box contribution against a global budget.

The heap always partitions the full support domain. After the floating
search passes with slack, every retained leaf is replayed outward.
"""
import argparse
from heapq import heappush,heappop
from math import comb,log,prod
from scipy.special import logsumexp
from flint import arb,ctx

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from group_rank_one_verify import up
from occupancy_model import local_data
from occupancy_memory import prepare,build
from occupancy_memory_verify import TILTS,replay
from occupancy_adaptive import volume,split,geometry_test,witness_selector


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=8,choices=range(3,33))
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--target-bits',type=int,default=55)
    parser.add_argument('--max-splits',type=int,default=20000)
    parser.add_argument('--screen-only',action='store_true')
    parser.add_argument('--pair-conditioned',action='store_true')
    parser.add_argument('--balanced-witness',action='store_true')
    parser.add_argument('--tilts',nargs='+',default=list(TILTS)+['.0032','.004','.005','.0064','.008','.01'])
    args=parser.parse_args()
    assert args.precision>=128
    geometry_test()
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    counts=[sum(row[u] for row in caps) for u in range(257)]
    prepared=prepare(local_data(min(args.groups,4)))
    ctx.prec=args.precision
    operators={}
    for tilt in args.tilts:
        operators[tilt]=build(prepared,tilt,args.groups,args.pair_conditioned)
        print('Built best-first operators',tilt,flush=True)
    witness=witness_selector(operators,args.groups,args.balanced_witness)
    location_count=comb(2048,args.groups)
    heap=[]
    serial=0
    def push(box,mult):
        nonlocal serial
        value,tilt,ps=witness(box)
        count=location_count*mult*prod(counts[hi] for lo,hi in box)
        score=min(0.,value)+log(count)
        serial+=1
        heappush(heap,(-score,serial,box,mult,tilt,ps,count))
    push(((38,256),)*args.groups,1)
    threshold=-(args.target_bits+2)*log(2)
    steps=0
    total_score=-heap[0][0]
    while steps<args.max_splits and total_score>threshold:
        item=heappop(heap)
        _,_,box,mult,_,_,_=item
        if all(lo==hi for lo,hi in box):
            heappush(heap,item)
            print('Cannot refine dominating singleton',box,flush=True)
            break
        for child,child_mult in split(box,mult):
            push(child,child_mult)
        steps+=1
        if steps%50==0:
            # Recompute the full positive sum; no subtractive update or
            # truncation of small contributions is used to certify stopping.
            total_score=logsumexp([-item[0] for item in heap])
            if steps%500==0:
                print('Best-first splits',steps,'leaves',len(heap),'union log2',total_score/log(2),flush=True)
    total_score=logsumexp([-item[0] for item in heap])
    assert sum(volume(item[2],item[3]) for item in heap)==219**args.groups
    print('BINARY64 full-domain cover:',steps,'splits;',len(heap),'leaves; log2 union',total_score/log(2),flush=True)
    if total_score>threshold:
        print('TARGET NOT CLOSED by this cover; no outward certificate.',flush=True)
        return
    if args.screen_only:
        print('Screen only; every leaf still requires outward replay.',flush=True)
        return
    total=arb(0)
    for index,(_,_,box,_,tilt,ps,count) in enumerate(sorted(heap),1):
        rectangles=[(lo,hi,counts[hi]) for lo,hi in box]
        total=up(total+replay(operators[tilt][0],tilt,ps,rectangles)*count)
        if index%2000==0:
            print('Outward leaves',index,'/',len(heap),flush=True)
    assert 0<total<arb(2)**-args.target_bits
    print('VERIFIED best-first occupancy',args.groups,'precision',args.precision,'upper',total,
          'margin',-total.log()/arb(2).log(),flush=True)
    print('Full-code certificate: NO. Other occupancies remain.')


if __name__=='__main__':
    main()
