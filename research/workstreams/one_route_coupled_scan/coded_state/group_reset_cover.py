"""Uniform FLOAT rectangle proposal for grouped resets; NOT a certificate.

Three coordinates (v,c,A), retaining the joint length/mask generating
function. Any completed cover needs independent directed replay. Default45s,
maximum180s,
10000 rectangle checks, no files, O(number of frontier rectangles) memory.
"""
import argparse
import math
import time
from functools import lru_cache
from group32_joint import old, weighted_lengths, logsum
from striped_reset_probe import box as coarse_box


@lru_cache(maxsize=4096)
def cached_parameters(G,d,r,active,extra,theta):
    return weighted_lengths(G,d,r,active,extra,theta,parameters=True)[1:]


def proposal_parameters(G,d,r,active,extra,theta):
    # Reuse nearby optimization proposals. Their original objective value
    # is NEVER used as a bound: the rectangle evaluates its actual theta,
    # actual endpoints, and the returned positive x,z afresh.
    rr = max(1,8*round(r/8)) if r else 0
    aa = min(G*rr,max(rr,round(G*rr*round(32*active/(G*r))/32))) if r else 0
    ee = max(0,128*round(extra/128))
    tt = max(1e-6,round(theta,4))
    return cached_parameters(G,d,rr,aa,ee,tt)


def linear_max(coeffs, ranges, G):
    """Maximize over the feasible integer polytope, not impossible corners.

    r<=A<=Gr and v>=2r+1. For fixed r, choose v and A endpoints.
    The resulting function is piecewise affine in r; test every breakpoint.
    """
    (vl,vh),(rl,rh),(al,ah) = ranges
    rl, rh = max(rl,(al+G-1)//G), min(rh,ah,(vh-1)//2)
    if rl > rh:
        return -math.inf
    candidates = {rl,rh}
    for value in ((vl-1)/2, al, ah/G):
        candidates.update((max(rl,min(rh,math.floor(value))),
                           max(rl,min(rh,math.ceil(value)))))
    cv,cr,ca = coeffs
    return max(cv*(vh if cv>=0 else max(vl,2*r+1))+cr*r
               +ca*(min(ah,G*r) if ca>=0 else max(al,r)) for r in candidates)


def route(B, c, lo, hi, theta, projective):
    N = 2**21
    closest = min(hi, N//2) if lo <= N//2 else N-lo
    excess = old.excess_log(B, B//2, c, closest*B/N, theta)
    answer = (N//B)*(max(0, excess)+math.log1p(math.exp(-abs(excess))))
    if projective and hi <= N//2:
        corrected = -math.inf if excess == -math.inf else (
            math.log(N//B)+excess+(N//B)*math.exp(excess)-math.log(math.expm1(theta)))
        answer = min(answer, corrected)
    return answer


def rectangle(B, b, d, item):
    N, W, n, G = 2**21, 491*2**21//1000, 2**21//b, b//d
    vl, vh, cl, ch, al, ah = item
    rl, rh = cl-1, ch-1
    ranges = ((vl, vh), (rl, rh), (al, ah))
    # Deterministic support exclusion inherited from the retained proof.
    if b*vh*B+2*ch*N <= (B//2)*N:
        return -math.inf, 'outer_excludes'
    gap_n = n-vl+ch
    gap = old.entropy_logchoose(gap_n, min(ch, gap_n//2))
    rc = (rl+rh)//2
    ac = min(G*rc, max(rc, (al+ah)//2))
    vc = max(2*rc+1, (vl+vh)//2)
    Slo, Shi = b*(vl-rh)+d*al, b*(vh-rl)+d*ah
    # Impossible corners can lie outside [0,N]; clip only to widen the
    # feasible cardinality interval, never evaluate a negative complement.
    Slo, Shi = max(0, Slo), min(N, Shi)
    alo, ahi = max(0, b*(vl-ch)), min(N, b*(vh-cl))
    Ehi = min(N, b*(vh-2*rl)+2*d*ah)
    tangents = [(0.0, 1.0), (float(W), 1/128)]
    if Ehi > W:
        for fraction in (1, .25):
            anchor = W+(Ehi-W)*fraction
            slope = math.log(anchor/(anchor-W))/(128*old.LN2)
            constant = W+old.entropy_logchoose(anchor, W)/(128*old.LN2)-slope*anchor
            tangents.append((constant, slope))
    best = (math.inf, None)
    for marker in ((17,16), (9,8), (5,4), (3,2), (2,1)):
        theta = math.log(marker[0]/marker[1])
        for forward, constant, slope in [(True, 0, 1)]+[(False, *t) for t in tangents]:
            lx, llz = proposal_parameters(G,d,rc,ac,vc-2*rc-1,theta*slope)
            lz = -math.exp(llz)
            logF = logsum([math.log(math.comb(G,a))+a*lx
                           -math.log(-math.expm1(lz-theta*slope*d*(G-a)))
                           for a in range(1, G+1)])
            base_constant = lz-math.log(-math.expm1(lz))
            if forward:
                small = linear_max((-lz, logF+2*lz-theta*b, -lx), ranges,G)
                large = linear_max((-lz+theta*b, logF+2*lz-2*theta*b,
                                    -lx+theta*d), ranges,G)-theta*N/2
                options = ([small] if Shi <= N//2 else
                           [large] if Slo >= N//2 else [small, large])
                value = gap+base_constant+max(options)+route(B,ch,Slo,Shi,theta,True)
            else:
                small = linear_max((-lz+theta*b*(slope-1),
                                    logF+2*lz+theta*b*(1-2*slope),
                                    -lx+theta*d*(2*slope-1)), ranges,G)+theta*b
                large = linear_max((-lz+theta*b*slope,
                                    logF+2*lz-2*theta*b*slope,
                                    -lx+theta*d*(2*slope-1)), ranges,G)-theta*N/2
                options = ([small] if ahi <= N//2 else
                           [large] if alo >= N//2 else [small, large])
                value = (gap+base_constant+theta*constant+max(options)
                         +route(B,ch,alo,ahi,theta,False))
            best = min(best, (value/old.LN2, ('forward' if forward else 'inverse', marker)))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outer', type=int, choices=(256,512), required=True)
    ap.add_argument('--group', type=int, choices=(32,64,128), default=64)
    ap.add_argument('--seconds', type=int, default=45, choices=range(1,181))
    args = ap.parse_args()
    B, d = args.outer, args.group
    b = 512 if B == 256 else 256
    n, G, W = 2**21//b, b//d, 491*2**21//1000
    cmax = min((n+1)//2, (W-1)//(d+1)+1)
    amax = min(G*(cmax-1), (W-1)//(d+1))
    bits = 40+(n*cmax*(amax+1)-1).bit_length()
    queue = [(1,n,1,cmax,0,amax)]
    checked, accepted = 0, 0
    start = time.monotonic()
    worst = -math.inf
    while queue:
        if time.monotonic()-start > args.seconds or checked >= 10000:
            print(dict(status='CAP_NOT_CERTIFICATE', checked=checked, accepted=accepted,
                       pending=len(queue), term_bits=bits, next_box=queue[-1]), flush=True)
            return
        vl,vh,cl,ch,al,ah = queue.pop()
        vl, ch = max(vl,2*cl-1), min(ch,(vh+1)//2)
        al, ah = max(al,cl-1), min(ah,G*(ch-1))
        cl, ch = max(cl,(al+G-1)//G+1), min(ch,ah+1)
        if vl>vh or cl>ch or al>ah:
            continue
        item = vl,vh,cl,ch,al,ah
        coarse = coarse_box(2**21,B,b,d,W,vl,vh,cl,ch)
        if coarse[0] <= -bits:
            value, method = coarse[0], ('coarse', coarse[1:])
        else:
            value, method = rectangle(B,b,d,item)
        checked += 1
        if value <= -bits:
            accepted += 1
            worst = max(worst,value)
        elif vl==vh and cl==ch and al==ah:
            print(dict(status='FAILED_BOUND_NOT_COUNTEREXAMPLE', item=item,
                       log2_upper=value, method=method, checked=checked,
                       accepted=accepted, term_bits=bits), flush=True)
            return
        else:
            widths = (b*(vh-vl), b*(ch-cl), d*(ah-al))
            axis = max(range(3), key=lambda i: widths[i])
            bounds = [[vl,vh],[cl,ch],[al,ah]]
            lo,hi = bounds[axis]
            mid = (lo+hi)//2
            for l,h in ((lo,mid),(mid+1,hi)):
                new = [x[:] for x in bounds]
                new[axis] = [l,h]
                queue.append(tuple(x for pair in new for x in pair))
    print(dict(status='FLOAT_COMPLETE_NOT_CERTIFICATE', checked=checked,
               accepted=accepted, worst=worst, term_bits=bits,
               seconds=round(time.monotonic()-start,3)), flush=True)


if __name__ == '__main__':
    main()
