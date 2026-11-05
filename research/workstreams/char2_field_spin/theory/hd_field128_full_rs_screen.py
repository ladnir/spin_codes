"""Isolated FLOAT prototype for one RS[64,16] emission/feedback block.

This is a different construction from the frozen four-block literal lift.
All local caps depend only on total occupancy. Positive expansions retain
birth/live terms even far below binary64's value range. No result here is
directed arithmetic, a whole-code certificate, or an implementation timing.
"""

from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb,isfinite,log2
from pathlib import Path

import numpy as np

import hd_field128_lift_screen as old


CONSTRUCTION = 'toeplitz-fullfield-full-rs64-16-alternating-swapped-pairs-k20-rate-half-v1'
FIELD_ORDER = 1<<128
NEG = float('-inf')


def parameters(q,n,h):
    bits = old.checked_q(q)
    if (type(n) is not int or type(h) is not int or not 1<=h<=16 or not 2*h<=n<=64
            or n&(n-1) or q<n):
        raise ValueError('Use a power-of-two width<=64, rank<=16, and q>=width')
    return bits


@lru_cache(256)
def exact_moments(j,z,q=FIELD_ORDER,n=64,h=16,joint=True):
    """Tiny rational oracle: K, T-K, return cap, uniform-nonzero live V."""
    parameters(q,n,h)
    if isinstance(z,float):
        raise ValueError('The reference mark must be exact')
    z = F(z)
    if type(j) is not int or not 0<=j<=n or not 0<z<1:
        raise ValueError('Invalid occupancy or exact mark')
    w,t,rho = 1-z,z+(1-z)/q,min(h,j)
    k = sum((F(comb(j,c),q**(c+min(h,j-c)))*z**(j-c)*w**c for c in range(j+1)),F(0))
    birth = sum((F(comb(j,c),q**c)*z**(j-c)*w**c*(1-F(1,q**min(h,j-c)))
                 for c in range(j)),F(0))
    inactive = n-j
    excess = sum((F(comb(inactive,e),q**e)*z**(inactive-e)*w**e*(1-F(1,q**(h-e)))
                  for e in range(min(h,inactive+1))),F(0))
    live = t**j*excess/(1-F(1,q**h))
    ret = min(birth,(1-F(1,q**rho))*z**max(0,n-h+1-j),(q**h-1)*live)
    if joint:
        joint_bound = F(0)
        for c in range(j+1):
            for e in range(min(h,inactive+1)):
                probability = min(1-F(1,q**rho),
                                  F(q**min(rho,h-e)-1,q**(j-max(0,j-c-h))))
                joint_bound += comb(j,c)*comb(inactive,e)*z**(n-c-e)*w**(c+e)*probability
        ret = min(ret,joint_bound)
    return k,birth,ret,live


def exact_packet(j,z,q=FIELD_ORDER,n=64,h=16,joint=True):
    k,birth,ret,live = exact_moments(j,z,q,n,h,joint)
    d = q**h-1
    density = F(d,(q-1)**h)
    return k,density*birth,ret/d,min(density*live,F(z)**(n-h+1)) if not j else density*live


@lru_cache(32)
def local_logs(z,q=FIELD_ORDER,n=64,h=16,joint=True):
    """One matrix for every occupancy0..n; toy n,h do not claim a full encoder."""
    bits = parameters(q,n,h)
    if not isfinite(z) or not 0<z<1:
        raise ValueError('Expected finite0<z<1')
    lz,lw = log2(z),log2(1-z)
    tilt = float(np.logaddexp2(lz,lw-bits))
    nonzero = old.log1mexp2(-h*bits)
    log_d = h*bits+nonzero
    density = nonzero-h*old.log1mexp2(-bits)
    quiet,birth,excess = [],[],[]
    for j in range(n+1):
        quiet.append(float(np.logaddexp2.reduce([
            log2(comb(j,c))+(j-c)*lz+c*lw-(c+min(h,j-c))*bits for c in range(j+1)])))
        birth.append(float(np.logaddexp2.reduce([
            log2(comb(j,c))+(j-c)*lz+c*lw-c*bits+old.log1mexp2(-min(h,j-c)*bits)
            for c in range(j)])) if j else NEG)
        excess.append(float(np.logaddexp2.reduce([
            log2(comb(j,e))+(j-e)*lz+e*lw-e*bits+old.log1mexp2(-(h-e)*bits)
            for e in range(min(h,j+1))])))
    result = []
    for j in range(n+1):
        rho,inactive = min(h,j),n-j
        live = j*tilt+excess[inactive]-nonzero
        ret = min(birth[j],old.log1mexp2(-rho*bits)+max(0,n-h+1-j)*lz,log_d+live)
        if joint and j:
            terms = []
            for c in range(j+1):
                for e in range(min(h,inactive+1)):
                    dimension = min(rho,h-e)
                    count = (dimension+max(0,j-c-h)-j)*bits+old.log1mexp2(-dimension*bits)
                    probability = min(old.log1mexp2(-rho*bits),count)
                    terms.append(log2(comb(j,c)*comb(inactive,e))+(n-c-e)*lz+(c+e)*lw+probability)
            ret = min(ret,float(np.logaddexp2.reduce(terms)))
        live += density
        if not j:
            live = min(live,(n-h+1)*lz)
        result.append((quiet[j],density+birth[j],ret-log_d,live))
    entries = np.array(result).reshape((n+1,2,2))
    if np.any(np.isnan(entries)) or np.any(np.isposinf(entries)):
        raise ArithmeticError('Nonfinite full-RS local log operator')
    entries.setflags(write=False)
    return entries


def geometry(maximum,*,good=False):
    """Fair orientation swaps cancel2**s; only total occupancy remains."""
    return tuple((s,f,comb(32,f)*comb(32-f,s))
                 for f in range(33) for s in range(min(1 if good else 32,32-f)+1)
                 if s+2*f<=maximum)


def single_coefficients(z,alpha,maximum,q,budget,joint=True):
    logs = local_logs(z,q,joint=joint).reshape((65,4))
    result = old._empty(maximum)
    for s,f,count in geometry(maximum):
        budget.check()
        result[:,s,f] = log2(count)+alpha*logs[s+2*f]
    return result


def good_pair_coefficients(z,alpha,maximum,q,budget,joint=True):
    logs = local_logs(z,q,joint=joint)
    good = geometry(maximum,good=True)
    result = old._empty(maximum)
    for s,f,count in good:
        budget.check()
        left = logs[s+2*f]
        for t,g,other in good:
            if s+t+2*(f+g)>maximum:
                continue
            right = logs[t+2*g]
            composed = np.logaddexp2(left[:,0,None]+right[None,0,:],left[:,1,None]+right[None,1,:])
            result[:,s+t,f+g] = np.logaddexp2(result[:,s+t,f+g],log2(count*other)+alpha*composed.ravel())
    return result


def regional_screen(z,alpha,maximum,*,q=FIELD_ORDER,paired=True,joint=True,max_seconds=60.,_budget=None):
    old.checked_q(q,construction=True)
    if type(maximum) is not int or not 1<=maximum<=64 or not 0<z<1 or not 0<alpha<=1:
        raise ValueError('Expected maximum1..64 and valid marks')
    budget = _budget or old.Budget(max_seconds)
    single = single_coefficients(z,alpha,maximum,q,budget,joint)
    if paired and alpha!=1:
        good,dirty = single.copy(),single.copy()
        good[:,2:,:],dirty[:,:2,:] = NEG,NEG
        operator = good_pair_coefficients(z,alpha,maximum,q,budget,joint)
        for left,right in ((good,dirty),(dirty,good),(dirty,dirty)):
            operator = np.logaddexp2(operator,old.polynomial_product(left,right,maximum,budget))
        run = old.polynomial_power(operator,64,maximum,budget)
    else:
        run = old.polynomial_power(single,128,maximum,budget)
    result = {}
    for r in range(maximum+1):
        for f in range(r//2+1):
            s = r-2*f
            row = run[:,s,f]-old.log_multinomial((old.PAIRS-s-f,s,f))
            result[r,f] = tuple(map(tuple,row.reshape((2,2))))
    budget.check()
    return result


def weight_screen(r,z,alpha,*,q=FIELD_ORDER,distance=F(7,20),paired=True,joint=True,
                  outer_measure='invertible-toeplitz',max_seconds=60.):
    bits = old.checked_q(q,construction=True)
    distance = F(distance)
    if not 0<distance<F(1,2):
        raise ValueError('Expected target distance between0and1/2')
    budget = old.Budget(max_seconds)
    regions = regional_screen(z,alpha,r,q=q,paired=paired,joint=joint,_budget=budget)
    forbidden = int(distance*old.OUTPUT_LENGTH)
    prefix = alpha*(r*old.outer_log_density(q,outer_measure)-(bits+old.log1mexp2(-bits))-forbidden*log2(z))
    profiles = []
    for a in range(r//2+1):
        budget.check()
        for b in range(r//2+1):
            count = old.profile_count_upper(old.PAIRS,r,a,b)
            if count:
                run = old.power(old.multiply(regions[r,a],regions[r,b]),old.REGIONS//2)
                profiles.append((float(log2(count)+prefix+old.total(run[0])),a,b))
    value = float(old.total(row[0] for row in profiles))
    if not isfinite(value):
        raise ArithmeticError('Nonfinite selected-weight prototype result')
    budget.check()
    return dict(construction=CONSTRUCTION,active_groups=r,z=z,alpha=alpha,field_order=q,field_bits=bits,
                target_distance=str(distance),forbidden_weight=forbidden,outer_measure=outer_measure,
                inner_emission='single-RS[64,16]',feedback='A transpose',CA_zero=True,
                fixed_nodes='GF256 additive subspace0..63; not GF64',state_nodes='128..143',
                pair_swaps=True,composition_steps=2 if paired else 1,joint_rank_return_cap=joint,
                log2_upper_screen=value,worst_profile=max(profiles),profile_count=len(profiles),
                covers_all_supports_at_selected_weight=True,covers_all_even_odd_splits=True,
                directed=False,whole_code_certificate=False,arithmetic='FLOAT',
                scope='new structural candidate; selected weight only, no implementation or global-distance claim')


def source_pins():
    pins = old.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return pins


def main():
    import argparse,json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r',type=int,required=True)
    parser.add_argument('--z',type=float,required=True)
    parser.add_argument('--alpha',type=float,required=True)
    parser.add_argument('--distance',default='7/20')
    parser.add_argument('--field-bits',type=int,default=128)
    parser.add_argument('--simple-return',action='store_true')
    parser.add_argument('--unpaired',action='store_true')
    parser.add_argument('--max-seconds',type=float,default=60.)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output already exists')
    if not 8<=args.field_bits<=128 or args.field_bits%8:
        parser.error('Field bits must be a multiple of8in8..128')
    result = weight_screen(args.r,args.z,args.alpha,q=1<<args.field_bits,distance=F(args.distance),
                            paired=not args.unpaired,joint=not args.simple_return,max_seconds=args.max_seconds)
    result['source_sha256'] = source_pins()
    if args.output:
        with args.output.open('x',encoding='utf-8') as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write('\n')
    print(json.dumps(result,allow_nan=False))


if __name__=='__main__':
    main()
