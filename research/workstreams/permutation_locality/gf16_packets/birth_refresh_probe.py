"""Floating diagnostic for the density of newly activated states.

This module only proposes bounds; birth_refresh.py implements the separate
outward evaluation used in certificate searches.
Low packet occupancies remain arbitrary; a Fourier bound controls the
weighted feedback density of the complementary nonnegative input measure.
"""
import argparse
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
import numpy as np
import rank_return
import scalar_cover as sc


def prepare(data,degree=8):
    W=data['windows'];degree=min(degree,W)
    return dict(data,birth_records_float=np.asarray(data['records'],dtype=float),birth_degree=degree)


def elementary(records,values,degree):
    """Truncated product coefficients via Newton identities, proposals only."""
    sums=records@(values[:,None]**np.arange(1,degree+1)[None,:])
    coefficients=np.zeros((len(records),degree+1));coefficients[:,0]=1
    for j in range(1,degree+1):
        for i in range(1,j+1):
            coefficients[:,j]+=(-1)**(i-1)*coefficients[:,j-i]*sums[:,i-1]
        coefficients[:,j]/=j
    return coefficients


def candidates(data,p,z,original):
    W=data['windows'];S=1<<data['bits'];degree=data['birth_degree']
    packet=((1+z)**4-1)/15
    values=np.array([((1+z)**(4-r)*(1-z)**r-1)/15 for r in range(5)])
    factors=1-p+p*values
    fourier=np.prod(factors[None,:]**data['records'],axis=1)
    coefficients=elementary(data['birth_records_float'],values,degree)
    zero=float(data['multiplicities']@fourier/S)
    result={'original':original};tail=fourier.copy();arbitrary=0.
    for cutoff in range(1,degree+2):
        j=cutoff-1
        tail-=((1-p)**(W-j)*p**j)*coefficients[:,j]
        if j:arbitrary+=comb(W,j)*(1-p)**(W-j)*p**j*packet**j
        order=np.argsort(tail);mass=np.cumsum(data['multiplicities'][order])
        center=tail[order[np.searchsorted(mass,S/2)]]
        tail_mass=sum(comb(W,k)*(1-p)**(W-k)*p**k*packet**k for k in range(cutoff,W+1))
        cap=min(tail_mass,float(data['multiplicities']@np.abs(tail-center)/S))
        candidate=original.copy()
        candidate[0,0]=min(candidate[0,0],max(0.,zero))
        candidate[0,1]=arbitrary
        candidate[0,2]=(S-1)*cap
        result[f'birth_from_j{cutoff}']=candidate
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    source=json.loads(args.input.read_text());data=prepare(rank_return.actual());rows=[]
    for row in source['rows']:
        p=float(Q(row['activity']));lam=float(Q(row['witness']['parameters'][0]));z=np.exp(-lam)
        original=rank_return.floating(data,sc.probabilities(Q(row['activity'])),lam)
        base=sc.log_power(original)
        scores={name:float(row['proposal']+(sc.log_power(matrix)-base)/np.log(2))
                for name,matrix in candidates(data,p,z,original).items()}
        result=dict(x=row['x'],scores=scores);rows.append(result);print(json.dumps(result),flush=True)
    if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,rows=rows),indent=2)+'\n')


if __name__=='__main__':main()
