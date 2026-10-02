"""Diagnostic at individual compositions of the positive outer envelope.

These are comparison-measure components, not actual BCH codewords.
Floating direct Poisson-binomial counts isolate a source of proof slack;
neither a positive nor a negative score is a complete certificate.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import numpy as np
from scipy.special import gammaln,logsumexp
from flint import ctx
import mass_density_potential as potential
import regional_potential


def round_counts(probabilities,total):
    p=np.asarray(probabilities,dtype=float)
    if (p.ndim!=1 or not len(p) or not np.isfinite(p).all() or (p<0).any()
            or abs(p.sum()-1)>1e-10 or type(total) is not int or total<1):
        raise ValueError('probability vector and positive integer total required')
    expected=total*p/p.sum();counts=np.floor(expected).astype(np.int64)
    deficit=total-int(counts.sum())
    counts[np.argsort(-(expected-counts),kind='stable')[:deficit]]+=1
    if int(counts.sum())!=total:raise ArithmeticError('rounding changed the group count')
    return counts.tolist()


def poisson_binomial_logs(features,counts):
    if (len(features)!=len(counts) or any(not 0<=p<=1 for p in features)
            or any(type(n) is not int or n<0 for n in counts)):
        raise ValueError('probabilities and nonnegative integer multiplicities required')
    result=np.array([0.]);shift=0;total=sum(counts)
    for probability,count in zip(features,counts):
        p=Q(probability)
        if not count or not p:continue
        if p==1:shift+=count;continue
        positive=float(np.log(float(p)));negative=float(np.log(float(1-p)))
        for _ in range(count):
            row=np.full(len(result)+1,-np.inf);row[:-1]=result+negative
            row[1:]=np.logaddexp(row[1:],result+positive);result=row
    output=np.full(total+1,-np.inf);output[shift:shift+len(result)]=result
    if abs(logsumexp(output))>1e-9:raise ArithmeticError('direct count recurrence lost normalization')
    return output


def proposed_composition(model,part,groups):
    _,_,logs=model.family(model.tilt);fs=np.array(list(map(float,model.features)));active=np.array(model.active)
    eta,mu,gamma=map(lambda x:float(Q(x)),part['dual'])
    scores=logs+eta*fs+mu*active+gamma*fs*(1-fs)
    counts=round_counts(np.exp(scores-logsumexp(scores)),groups)
    mean=sum(n*f for n,f in zip(counts,model.features))/groups
    variance=sum(n*f*(1-f) for n,f in zip(counts,model.features))/groups
    log_count=float(gammaln(groups+1)-sum(gammaln(n+1) for n in counts)+np.array(counts)@logs)
    return dict(counts=counts,mean=str(mean),variance=str(variance),
        active_groups=sum(n*a for n,a in zip(counts,model.active)),log_outer_count=log_count)


def zero_path_logs(zero_entries,epochs):
    """Conditional bounds for paths staying at zero after every step.

    Only the zero-to-zero entries are used. They do not depend on the
    transvection count, because every transvection fixes zero.
    The supplied zero entries may themselves be loose upper bounds.
    This is a floating diagnostic of the positive comparison measure.
    """
    import regional_count_probe as placement
    entries=np.asarray(zero_entries,dtype=float)
    if entries.ndim!=1 or not np.isfinite(entries).all() or (entries<0).any():
        raise ValueError('finite nonnegative zero-to-zero moments required')
    values,logs=placement.scaled_placement(entries[:,None,None],epochs)
    with np.errstate(divide='ignore'):
        return logs+np.log(values[:,0,0])


def main():
    import birth_classes
    import fiber_density
    import scalar_cover as sc
    import regional_count as regional
    import regional_count_probe as placement
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('--activities',nargs='+',default=['.4'])
    parser.add_argument('--parts',nargs='+',type=int,default=[4,5])
    parser.add_argument('--iterations',type=int,default=4);parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--tilts',nargs='+',help='Absolute output tilts; compositions stay fixed across this scan')
    parser.add_argument('--updates',type=int,help='Diagnostic inner-update override; does not change production parameters')
    parser.add_argument('--linear-only',action='store_true',help='Fast fixed-composition screen without nonlinear potential iteration')
    parser.add_argument('--zero-path-only',action='store_true',help='Only the comparison paths staying at zero; not a certificate or an actual-code lower bound')
    parser.add_argument('--exact-zero',action='store_true',help='Replace the zero-to-zero entry by its weighted Fourier count')
    parser.add_argument('--output',type=Path);args=parser.parse_args();source=json.loads(args.input.read_text())
    if source.get('schema')!='gf16-regional-count-screen-1' or not source.get('diagnostic_only'):
        parser.error('regional floating screen required')
    if args.precision<128:parser.error('precision at least 128 required')
    updates=source['updates'] if args.updates is None else args.updates
    if not 1<=updates<=32 or args.tilts and any(Q(t)<=0 for t in args.tilts):parser.error('1..32 updates and positive tilts required')
    data=fiber_density.attach(birth_classes.actual(updates));ctx.prec=args.precision
    model=sc.Model(sc.actual_components(129,Q(21,50),Q(1001,1000),True),data,source['threshold'],
        source['minimum_groups'],tilt=Q(3,16),inner=birth_classes,variance_shuffle=True,variance_bins=16,regional_count=True)
    rows=[]
    for p in map(Q,args.activities):
        matching=[r for r in source['rows'] if Q(r['activity'])==p]
        if not matching:parser.error('requested activity missing')
        row=min(matching,key=lambda r:r['proposal'])
        for lam in map(Q,args.tilts or [row['witness']['parameters'][0]]):
            witness=dict(row['witness'],parameters=[str(lam),*row['witness']['parameters'][1:]])
            if args.exact_zero:witness['regional_exact_zero']=True
            baseline=regional.local_operators(model,witness)
            raw={k:v for k,v in witness.items() if k not in (
                'regional_lazy_density_through','regional_feedback_classes_from','regional_feedback_classes_through',
                'regional_feedback_uniform_classes','regional_feedback_uniform_replace')}
            original=regional.local_operators(model,raw)
            bound=potential.Bound(data,original,(-potential.aq(lam)).exp(),[baseline])
            zero_logs=zero_path_logs(bound.float_linear[1,:,0,0],64)
            if not args.zero_path_only:region,logs=placement.scaled_placement(bound.float_linear[1],64)
            for part_index in args.parts:
                if not 0<=part_index<len(witness['regional_count_parts']):parser.error('variance part out of range')
                part=witness['regional_count_parts'][part_index];composition=proposed_composition(model,part,sc.G)
                mean=Q(composition['mean']);reference=mean/(model.tilt+(1-model.tilt)*mean)
                count_logs=poisson_binomial_logs(model.features,composition['counts'])
                weights=count_logs-np.arange(sc.G+1)*sc.logq(model.tilt)
                offset=float(lam)*model.threshold
                zero_score=float((composition['log_outer_count']+offset+sc.REGIONS*logsumexp(zero_logs+weights))/np.log(2))
                linear_score=None
                if not args.zero_path_only:
                    matrix,shift=placement.weighted_region(region,logs+weights)
                    linear_score=float((composition['log_outer_count']+offset+sc.log_power(matrix,sc.REGIONS)+sc.REGIONS*shift)/np.log(2))
                score=None;proposals=[None]
                if not args.linear_only and not args.zero_path_only:
                    initial=potential.iid_proposal(bound,float(reference))['potential']
                    score,proposals=regional_potential.propose(bound,[(weights,composition['log_outer_count'])],offset,
                        initial,iterations=args.iterations)
                result=dict(activity=str(p),variance_part=part_index,tilt=str(lam),composition=composition,
                    reference_activity=str(reference),cell=row['cell'],variance_interval=part['interval'],
                    mean_in_source_cell=Q(row['cell'][0])<=mean<=Q(row['cell'][1]),
                    variance_in_source_part=Q(part['interval'][0])<=Q(composition['variance'])<=Q(part['interval'][1]),
                    linear_score=linear_score,shared_score=score,zero_path_score=zero_score,potential=proposals[0])
                rows.append(result);print('FIXED COMPOSITION',json.dumps(result),flush=True)
                if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,
                    scope='individual positive-envelope compositions, not actual codewords or a complete bound',
                    threshold=source['threshold'],updates=updates,source_updates=source['updates'],
                    precision=args.precision,linear_only=args.linear_only,zero_path_only=args.zero_path_only,
                    exact_zero=args.exact_zero,rows=rows),indent=2)+'\n')
    print('Comparison-measure diagnostic only; no full-domain or actual-codeword claim.',flush=True)


if __name__=='__main__':main()
