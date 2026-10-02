"""Diagnostic only: couple shuffle variance to the outer counting bound.

At a fixed base-tilted mean, partition compositions by their average
Bernoulli variance. Each interval gets its own outer Chernoff witness
and variance-sensitive shuffle factor. Floating scores are not proofs.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
import scalar_cover as sc
import shuffle_variance as sv
from variance_partition import outer_witness


def split_point(model,row,bins):
    if type(bins) is not int or bins<1:raise ValueError('positive bin count required')
    x=Q(row['x']);witness=row['witness'];tilt=Q(witness['tilt'])
    if not 0<x<1 or tilt!=model.tilt:raise ValueError('interior point and base-tilt witness required')
    _,_,logs=model.family(tilt);features=np.array(list(map(float,model.features)))
    active=np.array(model.active);eta,mu=map(lambda v:float(Q(v)),witness['parameters'][1:])
    outside=sc.G*(logsumexp(logs+eta*features+mu*active)-eta*float(x)-mu*model.q_min/sc.G)/np.log(2)
    loss=model.comparison_loss((x,x),tilt,witness.get('variance_dual'))
    offset=row['proposal']-outside-sc.REGIONS*sc.logq(loss)/np.log(2)
    upper=min(x,1-x);parts=[]
    for j in range(bins):
        interval=(upper*j/bins,upper*(j+1)/bins)
        count,dual=outer_witness(logs,features,active,float(x),model.q_min/sc.G,interval)
        factor=min(float(loss),float(sv.density_upper(sc.G,sc.G*x,sc.G*interval[0]))) if interval[0]<=x*(1-x) else float(loss)
        # Intervals above the possible variance remain in the sum. They
        # retain the universal bound; the counting dual can suppress them.
        score=float(offset+sc.G*count/np.log(2)+sc.REGIONS*np.log2(factor))
        parts.append(dict(interval=list(map(str,interval)),score=score,dual=list(map(str,dual)),factor=factor))
    score=float(logsumexp([p['score']*np.log(2) for p in parts])/np.log(2))
    return dict(x=str(x),original=row['proposal'],split_score=score,parts=parts,
                note='Fixed inner/output witness; diagnostic only, no outward replay')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--bins',type=int,default=16)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();source=json.loads(args.input.read_text())
    if not source.get('diagnostic_only'):parser.error('diagnostic frontier input required')
    # Older birth-class diagnostics omitted this field but always used
    # parity. This probe only handles that same fixed row envelope.
    if source.get('row_parity',True) is not True:parser.error('row-parity frontier input required')
    model=sc.Model(sc.actual_components(129,Q(21,50),Q(1001,1000),True),{},
        int(Q(source['distance'])*sc.N),source['minimum_groups'],tilt=Q(1,8),
        variance_shuffle=source.get('variance_shuffle',False))
    rows=[]
    for row in source['rows']:
        if Q(row['witness']['tilt'])!=model.tilt:
            print('SKIP alternate tilt at',row['x'],flush=True);continue
        result=split_point(model,row,args.bins);rows.append(result)
        print(json.dumps({k:v for k,v in result.items() if k!='parts'}),flush=True)
    if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,source=str(args.input),
        distance=source['distance'],minimum_groups=source['minimum_groups'],bins=args.bins,rows=rows),indent=2)+'\n')


if __name__=='__main__':main()
