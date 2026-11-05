"""FLOAT proposal using grouped MDS coset weight enumerators. No cert claim.

For 0<=z<=1, a linear code's zero coset maximizes its weight enumerator.
Shorten each closed group to two batches; sum remaining coordinates freely.
This sums zero/nonzero group choices BEFORE the fractional-moment step.
"""
import argparse
import math
import time
from group32_joint import old, logsum
from striped_reset_probe import box as coarse_box

LNQ = 128*math.log(2)


def volume(b,d,v,c,W,parameters=False):
    r, free = c-1, b*(v-2*(c-1))
    terms = [(w, math.lgamma(2*d+1)-math.lgamma(w+1)-math.lgamma(2*d-w+1)
              +(w-d)*LNQ) for w in range(d+1,2*d+1)]
    def evaluate(lz):
        lp = logsum([0.0]+[coef+w*lz for w,coef in terms])
        ls = logsum([0.0, math.log1p(-math.exp(-LNQ))+LNQ+lz])
        return -W*lz+(b//d)*r*lp+free*ls
    lo,hi = -2*LNQ, 0.0
    for _ in range(65):
        left,right = (2*lo+hi)/3,(lo+2*hi)/3
        if evaluate(left) < evaluate(right):
            hi=right
        else:
            lo=left
    lz=(lo+hi)/2
    value=evaluate(lz)/LNQ
    return (value,lz) if parameters else min(value,float(b*(v-r)))


def rectangle(B,b,d,item):
    N,W,n = 2**21,491*2**21//1000,2**21//b
    vl,vh,cl,ch = item
    coarse = coarse_box(N,B,b,d,W,vl,vh,cl,ch)
    if coarse[0] == -math.inf:
        return coarse
    n1,n2=vh-cl,n-vl+ch
    count=(old.entropy_logchoose(n1,min(ch-1,n1//2))
           +old.entropy_logchoose(n2,min(ch,n2//2)))
    # The count is increasing in v and decreasing in c: at every z,
    # P_2d(z) <= S(z)^(2d) for the ACTUAL code enumerator. Our coefficient
    # envelope may exceed S^(2d), so explicitly take their minimum below.
    _,lz=volume(b,d,vh,cl,W,True)
    terms=[0.0]+[math.lgamma(2*d+1)-math.lgamma(w+1)-math.lgamma(2*d-w+1)
                   +(w-d)*LNQ+w*lz for w in range(d+1,2*d+1)]
    ls=logsum([0.0,math.log1p(-math.exp(-LNQ))+LNQ+lz])
    lp=min(logsum(terms),2*d*ls)
    out=(-W*lz+(b//d)*(cl-1)*lp+b*(vh-2*(cl-1))*ls)/LNQ
    out=min(out,float(b*(vh-cl+1)))
    al,ah=max(0,b*(vl-ch)),b*(vh-cl)
    selected=min(ah,N//2) if al<=N//2 else N-al
    answers=[coarse]
    for marker in old.MARKERS:
        theta=math.log(marker[0]/marker[1])
        excess=old.excess_log(B,B//2,ch,selected*B/N,theta)
        routed=(N//B)*(max(0,excess)+math.log1p(math.exp(-abs(excess))))
        value=(count+theta*(out-min(al,N//2))+routed)/old.LN2
        answers.append((value,'coset_inverse',marker))
    return min(answers)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outer',type=int,choices=(256,512),default=512)
    ap.add_argument('--group',type=int,choices=(1,8,32,64,128),default=32)
    ap.add_argument('--cover',action='store_true')
    args=ap.parse_args()
    B,d=args.outer,args.group
    b=512 if B==256 else 256
    n,W=2**21//b,491*2**21//1000
    start=time.monotonic()
    if not args.cover:
        for v,c in ((n//2,32),(n//2+128,128),(n//2+256,256),(n//2+512,512)):
            print(dict(v=v,c=c,volume=volume(b,d,v,c,W),
                       bound=rectangle(B,b,d,(v,v,c,c))),flush=True)
        return
    cmax=min((n+1)//2,(W-1)//(d+1)+1)
    bits=40+(n*cmax-1).bit_length()
    queue=[(1,n,1,cmax)]
    checked=accepted=0
    worst=-math.inf
    while queue:
        if time.monotonic()-start>50 or checked>=10000:
            print(dict(status='CAP_NOT_CERTIFICATE',checked=checked,accepted=accepted,
                       pending=len(queue),next_box=queue[-1]),flush=True)
            return
        vl,vh,cl,ch=queue.pop()
        vl,ch=max(vl,2*cl-1),min(ch,(vh+1)//2)
        if vl>vh or cl>ch:
            continue
        value,kind,marker=rectangle(B,b,d,(vl,vh,cl,ch))
        checked+=1
        if value<=-bits:
            accepted+=1
            worst=max(worst,value)
        elif vl==vh and cl==ch:
            print(dict(status='FAILED_BOUND_NOT_COUNTEREXAMPLE',v=vl,c=cl,
                       log2_upper=value,kind=kind,marker=marker,checked=checked),flush=True)
            return
        elif vh-vl>=2*(ch-cl) and vl<vh:
            mid=(vl+vh)//2
            queue.extend(((vl,mid,cl,ch),(mid+1,vh,cl,ch)))
        else:
            mid=(cl+ch)//2
            queue.extend(((vl,vh,cl,mid),(vl,vh,mid+1,ch)))
    print(dict(status='FLOAT_COMPLETE_NOT_CERTIFICATE',checked=checked,accepted=accepted,
               worst=worst,term_bits=bits,seconds=round(time.monotonic()-start,3)),flush=True)


if __name__=='__main__':
    main()
