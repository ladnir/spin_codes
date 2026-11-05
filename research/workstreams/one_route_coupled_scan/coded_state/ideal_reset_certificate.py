"""Ideal MDS-reset bridge, reusing the frozen persistent-state bound.

H_t[r,j]=xi_(bt+j)^r on globally distinct GF(2^32) points, and uniform
full-field feedback. This is NOT the running-sum update or a native kernel.
Source-only runs: one precision per invocation, bounded cover and replay.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'theory'))
import persistent_mds_certificate as retained


def partition(boxes,n,cmax,deadline):
    canonical=[]
    for vl,vh,cl,ch,kind,marker,_ in boxes:
        tau=Fraction(*marker)
        if not (1<=vl<=vh<=n and 1<=cl<=ch<=cmax and
                kind in ('forward','inverse') and 1<tau<=2**128):
            raise ArithmeticError('Invalid box')
        canonical.append((vl,vh,cl,ch,kind,tau.numerator,tau.denominator))
    total=0
    for c in range(1,cmax+1):
        if time.monotonic()>deadline: raise RuntimeError('Partition time cap')
        lower=2*c-1
        spans=sorted((max(vl,lower),vh) for vl,vh,cl,ch,*_ in canonical
                     if cl<=c<=ch and vh>=lower)
        cursor=lower
        for lo,hi in spans:
            if lo!=cursor: raise ArithmeticError('Partition gap/overlap')
            total+=hi-lo+1
            cursor=hi+1
        if cursor!=n+1: raise ArithmeticError('Missing partition suffix')
    if total!=cmax*(n+1-cmax): raise ArithmeticError('Wrong pair count')
    digest=hashlib.sha256(json.dumps(sorted(canonical),separators=(',',':')).encode()).hexdigest()
    return total,digest


def controls(boxes,n,cmax,deadline):
    changed=list(boxes)
    item=list(changed[0]); item[0]+=1; changed[0]=tuple(item)
    for bad in (boxes[1:],boxes+boxes[:1],changed):
        try: partition(bad,n,cmax,deadline)
        except ArithmeticError: continue
        raise ArithmeticError('Corrupt-cover control accepted')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--batch',type=int,choices=(256,512,1024),default=512)
    ap.add_argument('--outer',type=int,choices=(256,512),default=256)
    ap.add_argument('--precision',type=int,choices=(80,120))
    ap.add_argument('--distance',type=Fraction,default=Fraction(491,1000))
    ap.add_argument('--seconds',type=float,default=50,help='Cap EACH of proposal and replay phases')
    ap.add_argument('--fine-markers',action='store_true')
    ap.add_argument('--check-controls',action='store_true')
    ap.add_argument('--replay-part',type=int,nargs=2,metavar=('INDEX','PARTS'))
    args=ap.parse_args()
    if not 0<args.seconds<=120 or not 0<args.distance<Fraction(1,2):
        ap.error('Invalid distance/time cap')
    if args.replay_part and (args.precision is None or
            not 0<=args.replay_part[0]<args.replay_part[1]<=8):
        ap.error('Require a precision and 0<=INDEX<PARTS<=8')
    N,B,b,bits=2**21,args.outer,args.batch,128
    e=int(args.distance*N); n=N//b
    cmax=retained.c_limit(n,b,e)
    term_bits=40+(n*cmax-1).bit_length()
    if args.fine_markers:
        retained.MARKERS=tuple((32+i,32) for i in range(1,33))+((4,1),(8,1))
    started=time.monotonic()
    boxes=retained.cover(N,B,b,bits,args.distance,term_bits,args.seconds)
    deadline=time.monotonic()+args.seconds
    pairs,digest=partition(boxes,n,cmax,deadline)
    if args.check_controls: controls(boxes,n,cmax,deadline)
    print(dict(kind='exact_partition_only',pairs=pairs,sha256=digest,
               corrupted_controls=args.check_controls),flush=True)
    if args.precision is None: return
    index,parts=args.replay_part or (0,1)
    first,last=len(boxes)*index//parts,len(boxes)*(index+1)//parts
    bound=retained.PersistentUpward(args.precision)
    worst=Decimal('-Infinity')
    with localcontext(bound.up):
        for box in boxes[first:last]:
            if time.monotonic()>deadline: raise RuntimeError('Replay time cap; NO certificate')
            value=bound.rectangle(N,B,b,bits,e,box)
            if not value<=-term_bits: raise ArithmeticError(f'Failed directed box {box}: {value}')
            worst=max(worst,value)
    if time.monotonic()>deadline: raise RuntimeError('Final time cap; NO certificate')
    print(dict(kind=('directed_ideal_mds_reset_certificate' if parts==1 else
                     'directed_partial_replay_NOT_whole_certificate'),state=b,batch=b,
               outer=f'RS[{B},{B//2}]',field_bits=bits,N=N,threshold=e,
               precision=args.precision,boxes=len(boxes),pairs=pairs,
               partition_sha256=digest,worst_log2_upper=str(worst),
               term_bits=term_bits,margin_bits=(40 if parts==1 else None),
               verified_box_range=[first,last],part=index,parts=parts,
               seconds=round(time.monotonic()-started,3)),flush=True)


if __name__=='__main__': main()
