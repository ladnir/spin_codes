"""Floating point diagnostics for the dense frontier, not certificates."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
import scalar_cover as sc
import shape_return as shape


def counterfactual_matrices(matrix,data):
    """Deliberately optimistic modifications; never certificate operators."""
    L=(1<<data['bits'])-1
    result={name:matrix.copy() for name in ('no_returns','no_lazy_returns',
        'no_arbitrary_lazy_mass','no_uniform_lazy_mass','no_lazy_mass_to_arbitrary')}
    # Every nonzero source, including the added birth classes, can return.
    result['no_returns'][1:,0]=0
    result['no_lazy_returns'][1:,0]=np.minimum(matrix[1:,0],matrix[1:,2]/L)
    result['no_arbitrary_lazy_mass'][1,1]=0
    result['no_uniform_lazy_mass'][2,1]=0
    result['no_lazy_mass_to_arbitrary'][1:,1]=0
    return result


def outer_logs(model,cell,witness,activity_penalty=0):
    """Count and count-times-shuffle logs under the same saved duals."""
    tilt=Q(witness['tilt']);_,_,logs=model.family(tilt)
    features=np.array(list(map(float,model.features)));active=np.array(model.active)
    logs=logs-float(activity_penalty)*active
    if 'variance_partition' in witness:
        import variance_partition as variance
        counts=[];weighted=[]
        for interval,dual in variance.validate(cell,witness['variance_partition']):
            eta,mu,gamma=map(float,dual);lo,hi=map(float,interval)
            count=(sc.G*logsumexp(logs+eta*features+mu*active+gamma*features*(1-features))
                -sc.G*(min(eta*float(x) for x in cell)+min(gamma*lo,gamma*hi))-mu*model.q_min)
            loss=variance.factor(model,cell,interval[0],witness['variance_dual'])
            counts.append(count);weighted.append(count+sc.REGIONS*sc.logq(loss))
        return float(logsumexp(counts)),float(logsumexp(weighted))
    eta,mu=map(lambda v:float(Q(v)),witness['parameters'][1:])
    count=(sc.G*logsumexp(logs+eta*features+mu*active)
        -sc.G*min(eta*float(x) for x in cell)-mu*model.q_min)
    loss=model.comparison_loss(cell,tilt,witness.get('variance_dual'))
    return float(count),float(count+sc.REGIONS*sc.logq(loss))


def probe(model,x):
    cell=(x,x);score,witness=model.proposal(cell)
    tilt=Q(witness['tilt']);lam=Q(witness['parameters'][0])
    weights,_=model.weights(cell,tilt,witness.get('weights_dual'))
    p=weights[1]/sum(weights);probs=sc.probabilities(p)
    matrix=model.inner.floating(model.data,probs,float(lam))
    base=sc.log_power(matrix);tests={}
    for name,other in counterfactual_matrices(matrix,model.data).items():
        tests[name]=float(score+(sc.log_power(other)-base)/np.log(2))
    for updates in (3,4,8):
        data=dict(model.data,updates=updates)
        other=model.inner.floating(data,probs,float(lam))
        tests['updates_'+str(updates)]=float(score+(sc.log_power(other)-base)/np.log(2))
    count,weighted=outer_logs(model,cell,witness)
    reconstructed=(base+float(lam)*model.threshold+sc.PACKETS*sc.logq(sum(weights))+weighted)/np.log(2)
    if abs(reconstructed-score)>1e-5:raise ArithmeticError('diagnostic decomposition disagrees with proposed score')
    shuffle_bits=(weighted-count)/np.log(2)
    tests['no_shuffle_comparison_loss']=float(score-shuffle_bits)
    tests['no_returns_or_shuffle_comparison_loss']=float(tests['no_returns']-shuffle_bits)
    tests['no_lazy_mass_or_shuffle_comparison_loss']=float(tests['no_lazy_mass_to_arbitrary']-shuffle_bits)
    data=model.data;W=data['windows'];S=1<<data['bits'];pf=float(p);z=np.exp(-float(lam))
    c=1-16*pf/15;b=pf/15;shapes=[]
    if pf<=15/16:
        factors=np.array([[c*z**w+b*(1+z)**(4-r)*(1-z)**r for r in range(5)] for w in range(5)])
        for row in data['shapes']:
            value=float(data['multiplicities']@shape.float_products(factors[row['weights']].ravel(),row['records'],W))/S-c**W*z**row['weight']
            histogram=[0]*5
            for w,n in zip(row['weights'],row['records'][0].reshape(-1,5).sum(axis=1)):histogram[w]=int(n)
            shapes.append(dict(histogram=histogram,multiplicity=row['multiplicity'],bound=value))
    return dict(x=str(x),proposal=score,activity=str(p),witness=witness,matrix=matrix.tolist(),
                shuffle_comparison_bits=float(shuffle_bits),
                counterfactual_scores=tests,shapes=sorted(shapes,key=lambda row:-row['bound']))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points',nargs='+',default=['1/20','3/50','7/100'])
    parser.add_argument('--distance',default='2/25')
    parser.add_argument('--minimum-groups',type=int,default=97)
    parser.add_argument('--kernel',choices=('shape-return','conditioned-return','trimmed-return','rank-return','birth-refresh','birth-classes'),default='shape-return')
    parser.add_argument('--row-parity',action='store_true')
    parser.add_argument('--variance-shuffle',action='store_true')
    parser.add_argument('--variance-bins',type=int,default=0)
    parser.add_argument('--base-tilt',default='1/8')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.kernel=='birth-classes':
        import birth_classes as inner
    elif args.kernel=='birth-refresh':
        import birth_refresh as inner
    elif args.kernel=='rank-return':
        import rank_return as inner
    elif args.kernel=='trimmed-return':
        import trimmed_return as inner
    elif args.kernel=='conditioned-return':
        import conditioned_return as inner
    else:inner=shape
    model=sc.Model(sc.actual_components(129,Q(21,50),Q(1001,1000),args.row_parity),inner.actual(),
                   int(Q(args.distance)*sc.N),args.minimum_groups,tilt=Q(args.base_tilt),inner=inner,
                   variance_shuffle=args.variance_shuffle,variance_bins=args.variance_bins)
    rows=[]
    for x in map(Q,args.points):
        row=probe(model,x);rows.append(row);print(json.dumps(row),flush=True)
    if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,row_parity=args.row_parity,
        distance=args.distance,minimum_groups=args.minimum_groups,variance_shuffle=args.variance_shuffle,
        variance_bins=args.variance_bins,base_tilt=args.base_tilt,return_sources='all nonzero coordinates',rows=rows),indent=2)+'\n')


if __name__=='__main__':main()
