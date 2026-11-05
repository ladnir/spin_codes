"""Directed selected-weight bounds for the literal extension-field lift.

This standalone verifier keeps rational packet formulas, positive polynomial
coefficients, and Arb enclosures throughout. It certifies only explicitly
requested source weights. Neither a passing row nor their partial union is
a whole-code distance certificate. Existing GF256 artifacts are unchanged.
"""

import ast
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from importlib.metadata import version
from math import comb, isfinite
from pathlib import Path
import sys
from time import monotonic

from flint import arb, arb_poly, ctx

from hd_common_refresh_principal import nullities
from hd_grouped_quiet import quiet_cap
from hd_shifted_pairs_count import profile_count_upper


FIELD_ORDER = 1 << 128
PAIRS, GROUPS, REGIONS = 4096, 8192, 256
DIMENSION, OUTPUT_LENGTH, MAXIMUM = 1 << 20, 1 << 21, 64
CONSTRUCTION = 'toeplitz-fullfield-h4-alternating-swapped-pairs-k20-rate-half-v1'
OUTER_MEASURES = ('invertible-toeplitz', 'inherited-diagonal')


@dataclass
class Budget:
    max_seconds: float = 60.

    def __post_init__(self):
        if not isfinite(self.max_seconds) or not 0 < self.max_seconds <= 3600:
            raise ValueError('max_seconds must lie in (0,3600]')
        self.started = monotonic()

    def check(self):
        if monotonic()-self.started >= self.max_seconds:
            raise TimeoutError('Selected-weight directed replay exhausted max_seconds')


def exact(value):
    if isinstance(value, float):
        raise ValueError('Witnesses must be exact rationals, not binary floats')
    return F(value)


def checked_q(q):
    if (not isinstance(q, int) or q < 256 or q & (q-1)
            or (q.bit_length()-1) % 8 or q > FIELD_ORDER):
        raise ValueError('The fixed maps require field bits divisible by eight, in 8..128')
    return q.bit_length()-1


def rational(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def upper(value):
    endpoint = value.upper()
    if endpoint < 0 or not endpoint.is_finite():
        raise ArithmeticError('Invalid directed positive quantity')
    return endpoint


def round_up(value, bits=80):
    """Exact upward dyadic rounding with a relative mantissa budget."""
    value = F(value)
    if value < 0:
        raise ValueError('Negative packet entry')
    if not value:
        return value
    exponent = value.numerator.bit_length()-value.denominator.bit_length()
    quantum = F(2)**(exponent-bits)
    scaled = value/quantum
    return ((scaled.numerator+scaled.denominator-1)//scaled.denominator)*quantum


def rational_ceiling(value):
    if not value.is_finite():
        raise ArithmeticError('Nonfinite logarithmic bound')
    return F(int((value.upper()*10**12).upper().ceil().fmpz()), 10**12)


@lru_cache(1)
def defects():
    values = nullities(4)
    result = [0]*17
    for mask, value in enumerate(values):
        size = mask.bit_count()
        result[size] = max(result[size], int(value))
    return tuple(result)


@lru_cache(8)
def local_moments(z, q=FIELD_ORDER):
    """Exact local h4 moments, using positive differences throughout."""
    z = exact(z)
    checked_q(q)
    if not 0 < z < 1:
        raise ValueError('Expected 0<z<1')
    w, nu = 1-z, defects()
    quiet, excess, returns, births = [], [], [], []
    for length in range(17):
        quiet.append(sum(F(comb(length,c), q**(c+min(4,length-c)))
                         *z**(length-c)*w**c for c in range(length+1)))
        births.append(sum(F(comb(length,c),q**c)*z**(length-c)*w**c
                          *(1-F(1,q**min(4,length-c))) for c in range(length)))
        excess.append(sum(F(comb(length,c))*z**(length-c)*w**c
                          *(F(1,q**c)-F(1,q**4)) for c in range(min(4,length+1))))
        rho, ret = min(4,length), F(0)
        for c in range(length+1):
            for e in range(min(4,17-length)):
                rank = max(0,c-min(nu[c],max(0,4-length+c)))
                count = F(q**min(rho,4-e)-1,q**(length-max(0,length-c-4)))
                probability = min(1-F(1,q**rho),F(1,q**rank),count)
                ret += comb(length,c)*comb(16-length,e)*z**(16-c-e)*w**(c+e)*probability
        returns.append(min(births[-1],(1-F(1,q**rho))*z**max(0,13-length),ret))
    return tuple(quiet),tuple(excess),tuple(returns),tuple(births),z+w/q


def packet(occupancies, z, q=FIELD_ORDER, *, tighten_quiet=True):
    """Exact rational four-entry cap, before upward mantissa rounding."""
    z = exact(z)
    if len(occupancies)!=4 or any(not isinstance(n,int) or not 0<=n<=16 for n in occupancies):
        raise ValueError('Four integer occupancies in 0..16 are required')
    k,e,j,b,t = local_moments(z,q)
    zero, birth, running = F(1),F(0),F(1)
    for n in occupancies:
        birth = birth*k[n]+running*b[n]
        running *= t**n
        zero *= k[n]
    nonzero, running, base = F(0),F(1),F(1,q**4)
    for n in occupancies:
        nonzero = nonzero*base+running*e[16-n]
        running *= base+e[16-n]
    live = t**sum(occupancies)*nonzero/(1-F(1,q**16))
    ret, running = F(0),F(1)
    for n in occupancies:
        ret = ret*k[n]+running*j[n]
        running *= k[n]+j[n]
    d = q**16-1
    ret = min(ret,d*live)
    density = F(d,(q-1)**16)
    live *= density
    if tighten_quiet and not sum(occupancies):
        live = min(live,quiet_cap(z,q=q,state=16,width=64,h=4))
    return zero,density*birth,ret/d,live


@lru_cache(8)
def sparse_geometry(maximum):
    """Labeled pair arrangements; averaging swaps later divides by 2**s."""
    halves = []
    for a in range(17):
        for b in range(17-a):
            for f in range(17-a-b):
                if a+b+2*f<=maximum:
                    halves.append((a,b,f,comb(16,a)*comb(16-a,b)*comb(16-a-b,f)))
    counts = {}
    for i,(a,b,f,n) in enumerate(halves):
        for index in range(i,len(halves)):
            u,v,g,m = halves[index]
            s,full = a+b+u+v,f+g
            if s+2*full<=maximum:
                key = s,full,tuple(sorted((a+f,b+f,u+g,v+g)))
                counts[key] = counts.get(key,0)+n*m*(1 if i==index else 2)
    return tuple((s,f,n,occ) for (s,f,occ),n in counts.items())


def good_geometry(maximum):
    """Swap-averaged geometries with zero or one singleton per step."""
    result = []
    for a in range(17):
        for b in range(a,17):
            if 2*(a+b)<=maximum:
                result.append((0,a+b,comb(16,a)*comb(16,b)*(1 if a==b else 2),(a,a,b,b)))
    for a in range(16):
        for b in range(17):
            if 1+2*(a+b)<=maximum:
                result.append((1,a+b,32*comb(15,a)*comb(16,b),tuple(sorted((a,a+1,b,b)))))
    return tuple(result)


def valid_degrees(maximum):
    return tuple((s,f) for f in range(maximum//2+1) for s in range(maximum-2*f+1))


def matrix_product(left,right):
    return tuple(upper(left[2*r]*right[c]+left[2*r+1]*right[2+c])
                 for r in range(2) for c in range(2))


def matrix_power(base,exponent):
    result = tuple(arb(int(i in (0,3))) for i in range(4))
    while exponent:
        if exponent&1:
            result = matrix_product(result,base)
        exponent >>= 1
        if exponent:
            base = matrix_product(base,base)
    return result


def directed_parts(z,alpha,maximum,q,bits,budget,paired):
    geometry = sparse_geometry(maximum)
    budget.check()
    table = {}
    for _,_,_,occ in geometry:
        if occ not in table:
            budget.check()
            table[occ] = tuple(rational(round_up(x,bits)) for x in packet(occ,z,q))
    exponent = rational(alpha)
    powered = {occ:tuple(upper(x**exponent) if not x.is_zero() else arb(0) for x in values)
               for occ,values in table.items()}
    single = {}
    for s,f,n,occ in geometry:
        values = single.setdefault((s,f),[arb(0) for _ in range(4)])
        coefficient = rational(F(n,2**s))
        for i in range(4):
            values[i] += coefficient*powered[occ][i]
    if not paired or alpha==1:
        return single,None
    composed = {}
    good = good_geometry(maximum)
    for s,f,n,a in good:
        budget.check()
        for t,g,m,b in good:
            if s+t+2*(f+g)>maximum:
                continue
            values = composed.setdefault((s+t,f+g),[arb(0) for _ in range(4)])
            for i,value in enumerate(matrix_product(table[a],table[b])):
                if not value.is_zero():
                    values[i] += n*m*value**exponent
    return single,composed


def encode(coefficients,maximum,marker):
    stride = 2*maximum+1
    degree = stride*(maximum//2)+(maximum%2)
    result = [[arb(0) for _ in range(degree+1)] for _ in range(4)]
    marker = rational(marker)
    for (s,f),values in coefficients.items():
        if s+2*f<=maximum:
            for i in range(4):
                result[i][s+stride*f] = upper(values[i]*marker**(s+f))
    return tuple(arb_poly(row) for row in result)


def polynomial_product(left,right,maximum,budget):
    """Carry-free radix and downward-closed positive truncation."""
    budget.check()
    stride = 2*maximum+1
    degree = stride*(maximum//2)+(maximum%2)
    indices = tuple(s+stride*f for s,f in valid_degrees(maximum))
    result = []
    for r in range(2):
        for c in range(2):
            value = (left[2*r]*right[c]+left[2*r+1]*right[2+c]).truncate(degree+1)
            kept = [arb(0) for _ in range(degree+1)]
            for index in indices:
                kept[index] = upper(value[index])
            result.append(arb_poly(kept))
    return tuple(result)


def polynomial_power(base,exponent,maximum,budget):
    result = tuple(arb_poly([int(i in (0,3))]) for i in range(4))
    while exponent:
        if exponent&1:
            result = polynomial_product(result,base,maximum,budget)
        exponent >>= 1
        if exponent:
            base = polynomial_product(base,base,maximum,budget)
    return result


def regional_directed(z,alpha,maximum,*,q=FIELD_ORDER,bits=80,paired=True,budget=None):
    z,alpha = exact(z),exact(alpha)
    checked_q(q)
    if (not isinstance(maximum,int) or not 1<=maximum<=MAXIMUM or not 0<z<1
            or not 0<alpha<=1 or not isinstance(bits,int) or not 32<=bits<=512):
        raise ValueError('Invalid bounded regional scope or exact witness')
    budget = budget or Budget()
    marker = F(maximum,PAIRS)
    single,composed = directed_parts(z,alpha,maximum,q,bits,budget,paired)
    if composed is None:
        run = polynomial_power(encode(single,maximum,marker),128,maximum,budget)
    else:
        good = encode({k:v for k,v in single.items() if k[0]<=1},maximum,marker)
        dirty = encode({k:v for k,v in single.items() if k[0]>=2},maximum,marker)
        parts = (encode(composed,maximum,marker),polynomial_product(good,dirty,maximum,budget),
                 polynomial_product(dirty,good,maximum,budget),polynomial_product(dirty,dirty,maximum,budget))
        base = tuple(sum((part[i] for part in parts),arb_poly()) for i in range(4))
        run = polynomial_power(base,64,maximum,budget)
    stride,regions = 2*maximum+1,{}
    for s,f in valid_degrees(maximum):
        normalizer = comb(PAIRS,s)*comb(PAIRS-s,f)*rational(marker)**(s+f)
        regions[s+2*f,f] = tuple(upper(run[i][s+stride*f]/normalizer) for i in range(4))
    budget.check()
    return regions


def outer_density(q=FIELD_ORDER,mode='invertible-toeplitz'):
    checked_q(q)
    if mode not in OUTER_MEASURES:
        raise ValueError('Unknown proved outer measure')
    return F(q**160,(q-1)**32) if mode=='invertible-toeplitz' else F(q**256,(q-1)**128)


def source_pins():
    """Hash the complete recursive static local-import closure, including lazy imports."""
    directory = Path(__file__).resolve().parent
    pending,seen = [Path(__file__).name],{}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        path = directory/name
        raw = path.read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw,filename=str(path))):
            modules = ([alias.name for alias in node.names] if isinstance(node,ast.Import)
                       else [node.module] if isinstance(node,ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def weight_receipt(r,z,alpha,*,q=FIELD_ORDER,distance=F(7,20),precision=256,bits=80,
                   paired=True,outer_measure='invertible-toeplitz',max_seconds=60.,_pins=None):
    z,alpha,distance = exact(z),exact(alpha),exact(distance)
    field_bits = checked_q(q)
    if (not isinstance(r,int) or not 1<=r<=MAXIMUM or not 0<distance<F(1,2)
            or not isinstance(precision,int) or not 128<=precision<=4096):
        raise ValueError('Invalid selected-weight scope or precision')
    pins = source_pins() if _pins is None else _pins
    budget = Budget(max_seconds)
    forbidden = int(distance*OUTPUT_LENGTH)
    with ctx.workprec(precision):
        regions = regional_directed(z,alpha,r,q=q,bits=bits,paired=paired,budget=budget)
        union,profiles = arb(0),0
        for a in range(r//2+1):
            budget.check()
            for b in range(r//2+1):
                count = profile_count_upper(PAIRS,r,a,b)
                if count:
                    run = matrix_power(matrix_product(regions[r,a],regions[r,b]),REGIONS//2)
                    union += count*(run[0]+run[1])
                    profiles += 1
        log2arb = lambda value:value.log()/arb(2).log()
        prefix = rational(alpha)*(r*log2arb(rational(outer_density(q,outer_measure)))
                                  -log2arb(arb(q-1))-forbidden*log2arb(rational(z)))
        value = prefix+log2arb(union)
        ceiling = rational_ceiling(value)
        enclosure = str(value)
    budget.check()
    if source_pins()!=pins:
        raise RuntimeError('A pinned transitive helper changed during the replay')
    return dict(construction=CONSTRUCTION,field_order=q,field_bits=field_bits,
                active_groups=r,profiles=profiles,covers_all_supports_at_weight=True,
                covers_all_even_odd_splits=True,outer_groups=GROUPS,source_dimension=DIMENSION,
                output_length=OUTPUT_LENGTH,forbidden_weight=forbidden,target_distance=str(distance),
                pair_swaps=True,state_persists_across_regions=True,full_field_random_coefficients=True,
                outer_measure=outer_measure,projective_divisor=str(q-1),
                log2_upper_rational=str(ceiling),log2_upper_enclosure=enclosure,
                meets_40_bit_selected_weight_target=ceiling<=-40,
                witness=dict(z=str(z),alpha=str(alpha),precision=precision,packet_mantissa_bits=bits,
                             paired=paired,coefficient_maximum=r,coefficient_marker=str(F(r,PAIRS))),
                coefficient_method='positive-carry-free-Arb-truncated-bivariate',
                source_sha256=pins,pin_scope='complete recursive static local-import closure, including lazy imports',
                runtime=dict(python=sys.version,python_flint=version('python-flint'),numpy=version('numpy')),
                elapsed_seconds=monotonic()-budget.started,max_seconds=max_seconds,
                directed=True,whole_code_certificate=False,
                scope='every nonempty indexed support with exactly the selected active-group weight')


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r',type=int,required=True)
    parser.add_argument('--z',required=True,help='Exact rational or decimal witness')
    parser.add_argument('--alpha',required=True,help='Exact rational or decimal witness')
    parser.add_argument('--field-bits',type=int,default=128)
    parser.add_argument('--distance',default='7/20')
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--bits',type=int,default=80)
    parser.add_argument('--max-seconds',type=float,default=60.)
    parser.add_argument('--unpaired',action='store_true')
    parser.add_argument('--outer-measure',choices=OUTER_MEASURES,default='invertible-toeplitz')
    parser.add_argument('--output',help='Unused output path; existing files are never replaced')
    args = parser.parse_args()
    if not 8<=args.field_bits<=128 or args.field_bits%8:
        parser.error('--field-bits must be a multiple of eight in 8..128')
    if args.output and Path(args.output).exists():
        parser.error('--output already exists')
    row = weight_receipt(args.r,F(args.z),F(args.alpha),q=1<<args.field_bits,distance=F(args.distance),
                         precision=args.precision,bits=args.bits,paired=not args.unpaired,
                         outer_measure=args.outer_measure,max_seconds=args.max_seconds)
    encoded = json.dumps(row,indent=2,allow_nan=False)
    if args.output:
        with open(args.output,'x',encoding='utf-8') as stream:
            stream.write(encoded+'\n')
    print(encoded)


if __name__=='__main__':
    main()
