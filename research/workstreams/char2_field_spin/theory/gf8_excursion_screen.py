"""Excursion diagnostics for the existing positive sparse transfer.

Marks count paths in the upper-bound transfer, not actual sampled codes.
The regional transfer already retains all zero/live transitions; adding
these counters alone is not a new bound. Source coupling is still needed.
"""
from __future__ import annotations
import argparse
import json
import math
import numpy as np
from scipy.special import logsumexp

from cauchy_firstmoment_screen import geometry, log_choose
from gf8_feedthrough_mds_screen import weighted_matrices, search
from reduced_distance_certificate import log_polynomial_product
from sparse_active_interval_screen import regional_from_coefficients


MARKS=('births','returns','live_steps','quiet_live_steps','occupied_live_steps','input_weight')


def marked_local(y,tau,mark=None,amount=0.,state=8):
    result=np.asarray(weighted_matrices(y,tau,False,True,state,True),dtype=float)
    if mark is None or not amount:
        return result
    value=math.exp(amount)
    if mark=='births':
        result[:,0,1]*=value
    elif mark=='returns':
        result[:,1,0]*=value
    elif mark=='live_steps':
        result[:,1,:]*=value
    elif mark=='quiet_live_steps':
        result[0,1,:]*=value
    elif mark=='occupied_live_steps':
        result[1:,1,:]*=value
    elif mark=='input_weight':
        result*=np.exp(np.arange(9)*amount)[:,None,None]
    elif mark=='balance':
        result[:,0,1]*=value
        result[:,1,0]/=value
    else:
        raise ValueError('unknown excursion mark')
    return result


def log_terminal_vector(matrix,steps):
    vector=np.array([0.,-math.inf])
    for _ in range(steps):
        vector=logsumexp(vector[:,None]+matrix,axis=0)
    return vector


def path_mass(active,point,mark=None,amount=0.,state=8,dimension=96,k_log=20,terminal=False):
    g=geometry(k_log,128,dimension)
    logtilt,logtau=map(float,point)
    logy=-math.exp(logtilt)*active/g['n']
    local=marked_local(math.exp(logy),min(255.,math.exp(logtau)),mark,amount,state)
    with np.errstate(divide='ignore'):
        factor=np.log(local[:min(8,active)+1])
    coefficients=np.array([[[0.,-math.inf],[-math.inf,0.]]])
    steps=g['rows']//8
    while steps:
        if steps&1:
            coefficients=log_polynomial_product(coefficients,factor,active)
        steps >>= 1
        if steps:
            factor=log_polynomial_product(factor,factor,active)
    regional=regional_from_coefficients(coefficients,active,g['rows'],logtau)
    vector=log_terminal_vector(regional,128)
    return vector if terminal else float(logsumexp(vector))


def diagnose(active,point,state=8,dimension=96,k_log=20,epsilon=1e-4):
    g=geometry(k_log,128,dimension)
    vector=path_mass(active,point,state=state,dimension=dimension,k_log=k_log,terminal=True)
    total=float(logsumexp(vector))
    means={}
    for mark in MARKS:
        plus=path_mass(active,point,mark,epsilon,state,dimension,k_log)
        minus=path_mass(active,point,mark,-epsilon,state,dimension,k_log)
        means[mark]=(plus-minus)/(2*epsilon)
    logtilt,logtau=map(float,point)
    logy=-math.exp(logtilt)*active/g['n']
    bits=(log_choose(g['rows'],active)-((128-dimension)*active+1)*logtau
          -g['threshold']*logy+total)/math.log(2)
    final_live=math.exp(float(vector[1])-total)
    return dict(outer=[128,dimension],width=8,state=state,k_log=k_log,active=active,
                log_parameters=list(point),margin_bits=-bits,tilted_cover_means=means,
                final_live_mass_fraction=final_live,
                balance_residual=means['births']-means['returns']-final_live,
                time_residual=means['live_steps']-means['quiet_live_steps']-means['occupied_live_steps'],
                mean_inputs_per_completed_excursion=means['input_weight']/means['returns'],
                mean_live_packets_per_excursion=means['live_steps']/means['births'],
                scope='positive-upper-cover path diagnostic, not a code distribution',
                rigorous=False,whole_code=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active',type=int,nargs='+',default=[128,256])
    parser.add_argument('--iterations',type=int,default=40)
    args=parser.parse_args()
    previous=None
    for active in args.active:
        candidate=search(active,True,iterations=args.iterations,start=previous,state=8,rank_caps=False)
        previous=candidate['log_parameters']
        print(json.dumps(diagnose(active,previous)),flush=True)
