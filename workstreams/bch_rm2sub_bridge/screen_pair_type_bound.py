"""Discovery for actual shared-permutation region pair types.

Positive coefficient extraction first gives a Chernoff bound. A discrete
Fourier L1 bound additionally bounds the tilted type atom, including aliases.
Floating computations here are diagnostics, not outward certificates.
"""
import argparse
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
from flint import arb,arb_poly,ctx
import bridge as base
import tightened_occupancy as tight


def load_pairs():
    path=base.HERE/'generated/t128_s15_pair_spectrum.json';receipt=base.read(path)
    for name,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    powers=[];counts=[]
    for r in receipt['rows']:
        a,b,c=r['first_weight'],r['second_weight'],r['sum_weight']
        powers.append([128-(a+b+c)//2,(b+c-a)//2,(a+c-b)//2,(a+b-c)//2])
        counts.append(r['ordered_pairs']/2**30)
    powers=np.array(powers,dtype=np.int64);assert np.all(powers%4==0)
    return powers//4,np.array(counts)


SIGNS=np.array([[1,1,1,1],[1,-1,1,-1],[1,1,-1,-1],[1,-1,-1,1]])


def evaluate(values,powers,counts):
    linear=SIGNS@values
    # All exponents are multiples of four. Store the 33 powers once per
    # linear form and batch; no per-monomial exponentiation or heap callbacks.
    fourth=linear**4;table=np.empty((4,33,values.shape[1]),dtype=values.dtype);table[:,0]=1
    for n in range(1,33):table[:,n]=table[:,n-1]*fourth
    out=np.zeros(values.shape[1],dtype=values.dtype)
    for exponents,count in zip(powers,counts):
        out+=count*table[0,exponents[0]]*table[1,exponents[1]]*table[2,exponents[2]]*table[3,exponents[3]]
    return out


def probabilities(theta):
    values=np.exp(np.r_[0.,theta]-max(0,float(max(theta))))
    return values/values.sum()


def fourier_l1(p,powers,counts,grid):
    total=0.;peak=float(evaluate(p[:,None],powers,counts)[0]);M=math.prod(grid)
    for start in range(0,M,2048):
        indices=np.arange(start,min(start+2048,M));i=indices//(grid[1]*grid[2]);j=indices//grid[2]%grid[1];k=indices%grid[2]
        values=np.vstack([np.full(len(indices),p[0]),p[1]*np.exp(2j*np.pi*i/grid[0]),
            p[2]*np.exp(2j*np.pi*j/grid[1]),p[3]*np.exp(2j*np.pi*k/grid[2])])
        characteristic=evaluate(values,powers,counts)/peak
        total+=float(np.sum(np.abs(characteristic)**64))
    return total/M


def run(types,grid,tag):
    powers,counts=load_pairs();ctx.prec=192;kernel=tight.kernel_spectrum()
    poly=arb_poly([arb(kernel.get(j,0)) for j in range(129)])**64
    rows=[]
    for entry in types:
        m=np.array(entry,dtype=np.int64);assert len(m)==4 and sum(m)==8192 and min(m)>0
        logmulti=float(gammaln(8193)-np.sum(gammaln(m+1)))
        def objective(theta):
            p=probabilities(theta);value=float(evaluate(p[:,None],powers,counts)[0]);assert value>0
            return 64*math.log(value)-float(m@np.log(p))-logmulti
        initial=np.log(m[1:]/m[0]);opt=minimize(objective,initial,method='BFGS',options={'gtol':1e-6})
        p=probabilities(opt.x);value=objective(opt.x)
        x=int(m[2]+m[3]);y=int(m[1]+m[3]);assert x%2==y%2==0
        single=float((poly[x]/math.comb(8192,x)).log()+(poly[y]/math.comb(8192,y)).log())
        l1=fourier_l1(p,powers,counts,grid)
        row=dict(type=m.tolist(),tilt_probabilities=p.tolist(),chernoff_upper_log2=value/math.log(2),
            chernoff_ratio_to_single_product_bits=(value-single)/math.log(2),
            fourier_atom_upper_screen=l1,fourier_ratio_to_single_product_bits=(value+math.log(l1)-single)/math.log(2))
        rows.append(row);print('Pair type bound',row,flush=True)
    base.write_new(base.HERE/'generated'/f'pair_type_{tag}_screen.json',dict(status='PAIR_TYPE_FOURIER_SCREEN_ONLY',grid=grid,rows=rows,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),
            base.HERE/'generated/t128_s15_pair_spectrum.json',Path(tight.__file__)]}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--types',nargs='+',required=True,help='00:01:10:11')
    p.add_argument('--grid',nargs=3,type=int,default=[128,128,64]);p.add_argument('--tag',required=True);a=p.parse_args()
    run([list(map(int,t.split(':'))) for t in a.types],a.grid,a.tag)
