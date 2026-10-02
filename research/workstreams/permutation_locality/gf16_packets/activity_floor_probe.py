"""Diagnostic: impose the actual outer's minimum total packet activity.

This is a floating proposal only. No certificate evaluator uses its output.
An active four-row BCH group has union support at least 38, hence J>=38q.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
import birth_classes as inner
import frontier_probe
import scalar_cover as sc


def marked_inner(model,p,lam,gamma):
    """Log moment after marking every active packet by exp(gamma)."""
    if not 0<=p<=1 or lam<=0 or gamma<0:raise ValueError('valid activity and nonnegative activity mark required')
    normalizer=1+float(p)*np.expm1(gamma)
    marked=float(p)*np.exp(gamma)/normalizer
    matrix=model.inner.floating(model.data,sc.probabilities(Q(marked)),lam)
    return sc.log_power(matrix)+sc.PACKETS*np.log(normalizer)+lam*model.threshold


def probe(model,p,gammas):
    x=model.tilt*p/(1-p+model.tilt*p);cell=(x,x)
    baseline,witness=model.proposal(cell)
    tilt=Q(witness['tilt']);weights,_=model.weights(cell,tilt,witness.get('weights_dual'))
    reference=weights[1]/sum(weights)
    _,outside=frontier_probe.outer_logs(model,cell,witness)
    offset=outside+sc.PACKETS*sc.logq(sum(weights));rows=[]
    for gamma in gammas:
        grid=[.000001,.0001,.001,.004,.016,.064,.256,1.,4.]
        objective=lambda lam:marked_inner(model,reference,lam,float(gamma))
        scores=[objective(lam) for lam in grid];i=min(range(len(grid)),key=scores.__getitem__)
        fit=minimize_scalar(objective,bounds=(grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]),method='bounded')
        best=min((float(fit.fun),float(fit.x)),(scores[i],grid[i]))
        _,marked_outside=frontier_probe.outer_logs(model,cell,witness,activity_penalty=38*gamma)
        rows.append(dict(gamma=str(gamma),output_tilt=best[1],
            minimum_q_proposal=float((offset+best[0]-38*model.q_min*float(gamma))/np.log(2)),
            per_group_proposal=float((marked_outside+sc.PACKETS*sc.logq(sum(weights))+best[0])/np.log(2))))
    return dict(activity=str(p),x=str(x),baseline=baseline,outer_witness=witness,rows=rows)


def marked_counts(model,cell,parts,gamma,floor=38):
    """Reoptimize each outer dual after charging every active group."""
    import variance_partition as variance
    if gamma<0 or floor<0:raise ValueError('nonnegative activity mark and packet floor required')
    _,_,logs=model.family(model.tilt)
    features=np.array(list(map(float,model.features)));active=np.array(model.active)
    marked=logs-floor*float(gamma)*active
    return [sc.G*variance.outer_witness(marked,features,active,cell,model.q_min/sc.G,interval)[0]
            for interval,_ in parts]


def regional_probe(model,p,gammas,tilt_scales):
    """Diagnostic only: regional counts plus a reoptimized packet-floor mark."""
    import occupancy_birth_classes
    import regional_count as regional
    import regional_count_probe as floating
    import variance_partition as variance
    x=model.tilt*p/(1-p+model.tilt*p);cell=(x,x)
    _,base=variance.propose(model,cell,model.propose_with(cell,model.tilt))
    base.update(regional_tilted_atom=True,regional_fine_tilts=True,regional_tilted_variance=True)
    parts,witness=regional.prepare_witness(model,cell,base)
    ratios=[]
    for (interval,_),part in zip(parts,witness['regional_count_parts']):
        cap=variance.factor(model,cell,interval[0],witness['variance_dual'])
        values=regional.count_ratios(model.features,model.active,cell,interval,model.q_min,sc.G,
            regional.aq(cap),part['mgf_witnesses'])
        ratios.append(np.array([float(value.log()) for value in values]))
    weights,_=model.weights(cell,model.tilt);counts=np.arange(sc.G+1)
    families=[]
    for gamma in gammas:
        values=np.array(list(map(float,weights)));values[1]*=np.exp(float(gamma))
        scale=values.sum();reference=values[1]/scale
        log_binomial=(floating.log_choose(sc.G,counts)+counts*np.log(reference)
            +(sc.G-counts)*np.log1p(-reference))
        outer=marked_counts(model,cell,parts,gamma)
        families.append((gamma,scale,log_binomial,outer))
        print('REGIONAL FLOOR outer duals prepared for mark',str(gamma),flush=True)
    rows=[]
    for tilt_scale in tilt_scales:
        if tilt_scale<=0:raise ValueError('positive output tilt scale required')
        lam=Q(witness['parameters'][0])*tilt_scale
        local=occupancy_birth_classes.outward(model.data,lam);size=local[0].nrows()
        arrays=np.array([[[float(m[i,j]) for j in range(size)] for i in range(size)] for m in local])
        region=floating.placement(arrays,sc.G//model.data['windows'])
        for gamma,scale,log_binomial,outer in families:
            terms=[]
            for ratio,outside in zip(ratios,outer):
                matrix,shift=floating.weighted_region(region,log_binomial+ratio)
                terms.append(outside+sc.log_power(matrix,sc.REGIONS)+sc.REGIONS*shift)
            value=float((logsumexp(terms)+sc.PACKETS*np.log(scale)+float(lam)*model.threshold)/np.log(2))
            if not np.isfinite(value):raise ArithmeticError('nonfinite marked regional proposal')
            row=dict(gamma=str(gamma),tilt_scale=str(tilt_scale),output_tilt=str(lam),proposal=value)
            rows.append(row);print('REGIONAL FLOOR',json.dumps(row),flush=True)
    return dict(activity=str(p),x=str(x),regional_count=True,rows=rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distance',default='19/200')
    parser.add_argument('--minimum-groups',type=int,default=97)
    parser.add_argument('--base-tilt',default='3/16')
    parser.add_argument('--activities',nargs='+',default=['.15','.2'])
    parser.add_argument('--gammas',nargs='+',default=['0','1/32','1/16','1/8','1/4','1/2','1'])
    parser.add_argument('--regional',action='store_true',help='Use regional count bounds and reoptimize every marked outer dual')
    parser.add_argument('--variance-bins',type=int,default=8)
    parser.add_argument('--tilt-scales',nargs='+',default=['1','.975'])
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    model=sc.Model(sc.actual_components(129,Q(21,50),Q(1001,1000),True),inner.actual(),
        int(Q(args.distance)*sc.N),args.minimum_groups,tilt=Q(args.base_tilt),inner=inner,
        variance_shuffle=True,variance_bins=args.variance_bins)
    rows=[]
    for p in map(Q,args.activities):
        row=(regional_probe(model,p,list(map(Q,args.gammas)),list(map(Q,args.tilt_scales))) if args.regional
             else probe(model,p,list(map(Q,args.gammas))))
        rows.append(row)
        print(json.dumps({k:v for k,v in row.items() if k!='outer_witness'}),flush=True)
        if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,
            distance=args.distance,minimum_groups=args.minimum_groups,base_tilt=args.base_tilt,
            packet_floor_per_group=38,regional_count=args.regional,variance_bins=args.variance_bins,rows=rows),indent=2)+'\n')


if __name__=='__main__':main()
