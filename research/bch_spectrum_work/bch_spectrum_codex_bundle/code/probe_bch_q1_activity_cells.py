"""Positive binary64 diagnostic for a lower bound on true one-row tails.

Charge a complete cell whenever state enters active or the cell contains the
unique input one. Charged time dominates actual active time. This is not an
outward arithmetic certificate; a separate checker is required for that claim.
"""
from __future__ import annotations
import argparse
import json
import math
import sys
import time
from pathlib import Path
sys.dont_write_bytecode = True
import numpy as np
from scipy.special import bdtr
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'


def cell_matrices(length, memory):
    epsilon=2.**-memory
    decay=math.log1p(-epsilon)
    survive=math.exp(length*decay)
    after_one=(1-epsilon)*(-math.expm1(length*decay))/(length*epsilon)
    zero=np.zeros((2,2,2))
    one=np.zeros_like(zero)
    zero[0,0,0]=1
    zero[1,0,1]=1-survive
    zero[1,1,1]=survive
    one[:,0,1]=1-after_one
    one[:,1,1]=after_one
    return zero,one


def region_matrices(length, cells, memory):
    assert length%cells==0
    z,a=cell_matrices(length//cells,memory)
    zero=np.zeros((2,2,cells+1))
    one=np.zeros_like(zero)
    zero[0,0,0]=zero[1,1,0]=1
    for count in range(cells):
        nz=np.zeros_like(zero)
        na=np.zeros_like(one)
        for start in range(2):
            for end in range(2):
                for mid in range(2):
                    for charge in (0,1):
                        width=count+1
                        nz[start,end,charge:charge+width]+=zero[start,mid,:width]*z[mid,end,charge]
                        na[start,end,charge:charge+width]+=(one[start,mid,:width]*z[mid,end,charge]
                                                           +zero[start,mid,:width]*a[mid,end,charge])
        zero,one=nz,na
    one/=cells
    assert np.max(np.abs(zero.sum(axis=(1,2))-1))<1e-12
    assert np.max(np.abs(one.sum(axis=(1,2))-1))<1e-12
    return zero,one


def charge_distribution(zero, one, regions, max_weight, max_charge):
    cells=zero.shape[-1]-1
    current=np.zeros((max_weight+1,2,max_charge+1))
    current[0,0,0]=1
    for completed in range(regions):
        n=completed+1
        maximum=min(max_weight,n)
        oldmax=min(max_weight,completed)
        updated=np.zeros_like(current)
        skip=((n-np.arange(oldmax+1))/n)[:,None]
        take=(np.arange(1,maximum+1)/n)[:,None]
        prior_width=min(completed*cells,max_charge)+1
        for start in range(2):
            for end in range(2):
                for charge in range(cells+1):
                    width=min(prior_width,max_charge+1-charge)
                    if width<=0:
                        continue
                    z=zero[start,end,charge]
                    a=one[start,end,charge]
                    if z:
                        updated[:oldmax+1,end,charge:charge+width]+=(current[:oldmax+1,start,:width]*skip)*z
                    if a:
                        updated[1:maximum+1,end,charge:charge+width]+=(current[:maximum,start,:width]*take)*a
        current=updated
    return current.sum(axis=1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cells',type=int,required=True)
    args=parser.parse_args()
    output=GEN/f'bch256_q1_activity_cells_{args.cells}_diagnostic.json'
    assert not output.exists()
    start=time.perf_counter()
    cell_length=8192//args.cells
    maximum=(2*209716)//cell_length+2
    z,a=region_matrices(8192,args.cells,22)
    distribution=charge_distribution(z,a,256,42,maximum)
    n=np.arange(maximum+1)*cell_length
    cdf=np.array([1. if h<=209716 else bdtr(209716,int(h),.5) for h in n])
    # Hoeffding gives a simpler analytic lower bound without a binomial CDF.
    hoeffding=np.array([1. if h<=209716 else -math.expm1(-2*(209717-h/2)**2/h)
                        if h<2*209717 else 0. for h in n])
    values=distribution@cdf
    analytic=distribution@hoeffding
    old=json.loads((GEN/'bch256_q1_full_arb_transfer.json').read_text())
    rows=[]
    for w in (38,40,42):
        log=math.log2(8192*values[w])
        rows.append(dict(weight=w,point_probability_lower_diagnostic=float(values[w]),
                         row_coefficient_lower_log2_diagnostic=log,
                         hoeffding_row_coefficient_lower_log2_diagnostic=math.log2(8192*analytic[w]),
                         old_chernoff_coefficient_log2=old['coefficient_rows'][str(w)]['coefficient_upper']['log2_diagnostic'],
                         maximum_possible_chernoff_improvement_bits=old['coefficient_rows'][str(w)]['coefficient_upper']['log2_diagnostic']-log))
    result=dict(classification='Binary64 diagnostic of mathematically valid lower-tail bound; not an arithmetic certificate',
                cells_per_region=args.cells,cell_length=cell_length,maximum_charge=maximum,
                memory_bits=22,outer_regions=256,outer_rows=8192,cutoff=209716,
                rows=rows,elapsed_seconds=time.perf_counter()-start,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                               (Path(__file__),GEN/'bch256_q1_full_arb_transfer.json')})
    write_new(output,result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
