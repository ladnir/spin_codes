"""Independent checks of the count-box majorants; not sampled certificates."""
from itertools import product
import random
import time
import math
import numpy as np
from cover_layered_location_regions import (joint_dimension,region_env,location_box,location_coeff,
    local_logs,plain_box,cover_box,POWER,BANK,M,b,N,inverse_binomial,tube)
from layered_general_components import mixed_dimension
from layered_outer_locations import hardest,rows,start,objective
from layered_two_layer_components import exact_replay
from check_layered_component_count import log_upper_units
from fractions import Fraction as F


def run():
    began=time.monotonic();rng=random.Random(904);center=hardest()[1:6];tested=0;strict=0
    # Random affine coefficients exercise both signs. Corners AND interior
    # points test the formula; the proof of interior coverage is algebraic.
    for R in (190,191,192,193):
        P=(R*M-b)//2;env=region_env(R)
        for trial in range(16):
            lo=tuple(v-rng.randrange(1,257) for v in center)
            hi=tuple(v+rng.randrange(1,257) for v in center)
            a=tuple(rng.randrange(-40*POWER,40*POWER) for _ in range(5));rate=POWER*rng.randrange(1,33)
            majorant=joint_dimension(R,lo,hi,a,rate)
            dim=mixed_dimension(R,P,lo[0],hi[1],hi[2],lo[3],env)
            old=rate*dim+sum(c*(h if c>=0 else l) for c,l,h in zip(a,lo,hi))
            strict+=majorant<old
            points=list(product(*[(l,h) for l,h in zip(lo,hi)]))
            points += [tuple(rng.randrange(l,h+1) for l,h in zip(lo,hi)) for _ in range(16)]
            for x in points:
                d=mixed_dimension(R,P,*x[:4],env)
                assert sum(c*v for c,v in zip(a,x))+rate*d<=majorant
                tested+=1
    # An exhaustive small box, rather than just its corners.
    lo=tuple(v-1 for v in center);hi=tuple(v+1 for v in center)
    for i in (1,3,5,7):
        _,a,rate,_=location_coeff(192,i);bound=joint_dimension(192,lo,hi,a,rate)
        for x in product(*(range(l,h+1) for l,h in zip(lo,hi))):
            d=mixed_dimension(192,(192*M-b)//2,*x[:4],region_env(192))
            assert sum(c*v for c,v in zip(a,x))+rate*d<=bound;tested+=1
    # Independent floating evaluation of the original full expression, using
    # actual point dimensions, must lie below the integer box upper bound.
    # Floats here are only regression checks, never acceptance criteria.
    direct_checks=0
    for i in (1,3,5,7):
        _,m,two,p,markers=BANK[i];x=start(markers)
        x=np.append(x,math.log2(p/(1-p)));rr,ee=rows(x,M,m,two)
        x[14:18]=rr.max(axis=1);x[18:20]=ee.max(axis=1)
        lo=tuple(v-32 for v in center);hi=tuple(v+32 for v in center)
        bound=location_box(192,lo,hi,i,7900,8100)
        for point in product(*[(l,h) for l,h in zip(lo,hi)]):
            d=mixed_dimension(192,(192*M-b)//2,*point[:4],region_env(192))
            profile=(192,)+point+(d,0)
            for s in (7900,7950,8000,8050,8100):
                value=objective(x,profile,s,m,two)*(192*M)
                assert value<=bound+1e-5,(i,point,s,value,bound);direct_checks+=1
    for s in (5901,6903,7636,7990,8061,8125,8192):
        tight=inverse_binomial(s)
        old=(POWER//2048)*log_upper_units(F(1,math.comb(M,s)),2048)
        assert old>=tight and old-tight<2*(POWER//2048)
    # The generic smaller-active-row box agrees conservatively with the
    # original complete per-profile expression and its exact interpolation.
    reference=exact_replay(hardest(),power=2048)-4*64*(M-5900)+M
    assert 0<=reference-plain_box(192,center,center)<20000
    # A failed cell must survive a subdivision budget, never become coverage.
    failed=tube(193,49168,128,budget=1)
    assert not failed['complete'] and failed['log2_union_upper'] is None
    assert failed['failed_cells']==1 and failed['failed'][0][:2]==(0,49168)
    passed=tube(192,48000,128)
    assert passed['complete'] and passed['log2_union_upper']<=-417
    print(dict(joint_affine_dimension_checks=tested,strict_improvements=strict,
               direct_expression_checks=direct_checks,binomial_rounding_checks=7,
               honest_failed_region=True,gapless_region_replay=True,
               elapsed_seconds=time.monotonic()-began,whole_code_certificate=False))
    assert time.monotonic()-began<60


if __name__=='__main__':run()
