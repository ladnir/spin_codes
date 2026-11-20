"""Outward discrete-Fourier L1 bound for one actual shared region pair type.

Roots are enclosed by Arb. A point-dependent absolute error enclosure covers each
binary64 evaluation of the degree-128 pair enumerator; see the proof note.
The discrete Fourier sum includes nonnegative aliases, so no truncation of
the coefficient polynomial and no tail approximation is used.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from flint import arb,arb_poly,ctx
import bridge as base
import tightened_occupancy as tight
from screen_pair_type_bound import load_pairs,SIGNS
import pair_fourier_error as tracked


def up(x):return np.nextafter(x,np.inf)


def roots(length):
    result=[]
    for j in range(length):
        angle=2*arb.pi()*j/length;c=angle.cos();s=angle.sin()
        re=float(c);im=float(s)
        assert abs(arb(re)-c).upper()<=arb(2)**(-52)
        assert abs(arb(im)-s).upper()<=arb(2)**(-52)
        result.append(complex(re,im))
    return np.array(result,dtype=np.complex128)


def evaluate(values,powers,counts):
    linear=SIGNS@values;squared=linear*linear;fourth=squared*squared
    table=np.empty((4,33,values.shape[1]),dtype=np.complex128);table[:,0]=1
    for n in range(1,33):table[:,n]=table[:,n-1]*fourth
    result=np.zeros(values.shape[1],dtype=np.complex128)
    for exponents,count in zip(powers,counts):
        result+=count*table[0,exponents[0]]*table[1,exponents[1]]*table[2,exponents[2]]*table[3,exponents[3]]
    return result


def discrete_upper(p,powers,counts,grid):
    assert np.finfo(np.float64).eps==2**-52 and np.all(4*powers.sum(axis=1)==128)
    assert len(counts)==155 and float(counts.sum())==1.
    assert all(F.from_float(float(c)).denominator<=1<<30 for c in counts)
    tables=[roots(n) for n in grid];M=math.prod(grid);total=F(0)
    for start in range(0,M,2048):
        indices=np.arange(start,min(start+2048,M));i=indices//(grid[1]*grid[2]);j=indices//grid[2]%grid[1];k=indices%grid[2]
        values=np.vstack([np.full(len(indices),p[0],dtype=np.complex128),
            p[1]*tables[0][i],p[2]*tables[1][j],p[3]*tables[2][k]])
        value,error=tracked.evaluate(values,powers,counts)
        magnitude=up(np.sqrt(up(up(value.real*value.real)+up(value.imag*value.imag))))
        # Scale by 2^30 before the 64th power, avoiding tiny absolute values.
        scaled=up(up(magnitude*2**30)+up(error*2**30))
        for _ in range(6):scaled=up(scaled*scaled)
        # At most 2048 nonnegative terms. The factor exceeds the worst
        # gamma_2047 summation correction, regardless of reduction order.
        chunk=float(up(float(np.sum(scaled))*(1+2**-38)))
        assert math.isfinite(chunk) and chunk>0;total+=F.from_float(chunk)
    return total/M/F(1<<1920)


def run(tag,verify=False):
    assert tag.isidentifier();path=base.HERE/'generated'/f'pair_type_{tag}_outward.json'
    old=base.read(path) if verify else None
    if not verify:assert not path.exists()
    if old:
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    screen_path=base.HERE/'generated/pair_type_central01_screen.json';screen=base.read(screen_path)
    for name,digest in screen['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    row=screen['rows'][0];m=row['type'];grid=screen['grid']
    assert m==[6638,736,736,82] and grid==[128,128,64]
    denominator=1<<40;integers=[round(p*denominator) for p in row['tilt_probabilities'][1:]]
    integers=[denominator-sum(integers)]+integers
    assert min(integers)>0;ps=[F(x,denominator) for x in integers]
    p=np.array([float(x) for x in ps]);assert sum(ps)==1 and all(F.from_float(float(x))==y for x,y in zip(p,ps))
    ctx.prec=512 if verify else 256;powers,counts=load_pairs()
    averaged=discrete_upper(p,powers,counts,grid)
    multinomial=math.factorial(8192)//math.prod(math.factorial(v) for v in m)
    tilted_type_mass=arb(multinomial)
    for value,power in zip(ps,m):tilted_type_mass*=(arb(value.numerator)/value.denominator)**power
    from audit_bch_q1_full_arb import rational
    bound=rational(((arb(averaged.numerator)/averaged.denominator)/tilted_type_mass).upper())
    t,s,spectrum=base.load_map(tight.NAME);kernel=tight.kernel_spectrum()
    poly=arb_poly([arb(kernel.get(j,0)) for j in range(129)])**64
    single=poly[818]/math.comb(8192,818);product_lower=rational((single*single).lower())
    ratio=bound/product_lower;assert 0<bound and ratio<F(11,10)
    if old:
        assert bound<=base.decode(old['pair_region_upper']) and ratio<=base.decode(old['ratio_to_single_product_upper'])
        print('512-bit root/final replay passed the central shared-type Fourier bound',flush=True);return
    local=[Path(__file__),Path(tracked.__file__),screen_path,base.HERE/'screen_pair_type_bound.py',base.HERE/'generated/t128_s15_pair_spectrum.json',
        Path(tight.__file__),base.HERE/'general_occupancy.py',base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']
    base.write_new(path,dict(status='OUTWARD_ONE_SHARED_REGION_PAIR_TYPE',configuration=tight.NAME,
        region_type=m,grid=grid,tilt_probabilities=[base.encode(v) for v in ps],evaluation_error_rule='2^-40 A + 2^-49 D + 2^-1000, outward',
        pair_region_upper=base.encode(bound),single_region_product_lower=base.encode(product_lower),
        ratio_to_single_product_upper=base.encode(ratio),ratio_below_eleven_tenths=True,
        full_second_moment_certified=False,local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local}))
    print('Outward central pair ratio <=',float(ratio),'log2',math.log2(ratio.numerator)-math.log2(ratio.denominator),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tag',required=True);p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.tag,a.verify)


