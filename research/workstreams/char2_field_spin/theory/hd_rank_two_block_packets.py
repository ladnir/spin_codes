"""Exact rational packet operators for block emission and rank-two refresh.

Four RS[16,4] emission blocks; feedback uses four RS[16,4] or two RS[32,8]
blocks. The larger-feedback bound takes the worst labelling of an emission
profile. This supplies directed-verifier inputs, not a whole-code claim.
"""
from functools import lru_cache
from fractions import Fraction as F
from math import comb,prod
from hd_bucket_feedback import profiles


@lru_cache(None)
def local_data(z):
    z=F(z);q=65536
    quiet=lambda n,h:sum(F(comb(n,c),q**(c+min(h,n-c)))*z**(n-c)*(1-z)**c for c in range(n+1))
    k4=tuple(quiet(n,4) for n in range(17));k8=tuple(quiet(n,8) for n in range(33))
    flat=tuple(sum(F(comb(n,c),q**min(4,c))*z**(n-c)*(1-z)**c for c in range(n+1)) for n in range(17))
    powers=tuple((z+(1-z)/q)**j for j in range(65))
    return k4,k8,flat,powers


def matrix(occupancy,z,feedback_dimension=8):
    z=F(z);q=65536;d=q**16-1;epsilon=F(1,q**2)
    if len(occupancy)!=4 or any(not 0<=n<=16 for n in occupancy):raise ValueError('Invalid block occupancy')
    if feedback_dimension not in (4,8):raise ValueError('Feedback must use four or eight checks')
    k4,k8,flat,powers=local_data(z);j=sum(occupancy);t=powers[j]
    if feedback_dimension==4:
        k=prod(k4[n] for n in occupancy);quiet_lower=k;rho=sum(min(4,n) for n in occupancy)
    else:
        pairings=tuple((occupancy[a]+occupancy[b],occupancy[c]+occupancy[e])
                       for a,b,c,e in ((0,1,2,3),(0,2,1,3),(0,3,1,2)))
        quiet_values=tuple(k8[l]*k8[r] for l,r in pairings)
        k=max(quiet_values);quiet_lower=min(quiet_values)
        rho=min(min(8,l)+min(8,r) for l,r in pairings)
    v=t*(q**16*prod(flat[16-n] for n in occupancy)-1)/d
    gamma=min(k,z**max(0,13-max(occupancy))/q**rho,d*v)
    return k,t-quiet_lower,v/d if j else F(0),(v if j else (1-epsilon)*v)+epsilon*gamma


def grouped_packets(z,feedback_dimension=8,bits=12,maximum=64,cutoff=64):
    from hd_sumrank_paired_certificate import round_up
    if not 0<=maximum<=64:raise ValueError('Invalid degree cutoff')
    groups={};tails={}
    for occupancy,count in profiles(16,64,4):
        j=sum(occupancy)
        if j>maximum:continue
        entries=matrix(occupancy,z,feedback_dimension)
        if j>cutoff:
            if j in tails:
                n,old=tails[j];tails[j]=(n+count,tuple(max(a,b) for a,b in zip(old,entries)))
            else:tails[j]=(count,entries)
        else:
            key=(j,tuple(round_up(v,bits) for v in entries));groups[key]=groups.get(key,0)+count
    for j,(count,entries) in tails.items():
        key=(j,tuple(round_up(v,bits) for v in entries));groups[key]=groups.get(key,0)+count
    result=tuple((j,count,entries) for (j,entries),count in groups.items())
    for j in range(maximum+1):
        if sum(n for k,n,_ in result if k==j)!=comb(64,j):raise ArithmeticError('Incorrect profile multiplicity')
    return result
