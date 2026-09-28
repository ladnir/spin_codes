"""Selected composition bounds for the complete BCH row-mixture envelope.

All mixture counts carry their multinomial location multiplicity. A
uniformly shuffled heterogeneous input is dominated using a three-category
tilt and the explicit loss in heterogeneous.py. Selected compositions do
not certify the remaining compositions or the full code.
"""
import argparse
from fractions import Fraction as Q
from math import factorial,log

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
from flint import arb,ctx

import model
import iid_kernel
import density_tangent
from bch_joint_support import authenticated_caps
from row_mixture import envelope,pair_components
from heterogeneous import density_loss,tilted_reference,capped_density_loss
from probe import aq


def prepare(components,data):
    coefficients=np.array([log(c.numerator)-log(c.denominator) for _,c,*_ in components])
    probabilities=np.array([[float(x) for x in row[2:5]] for row in components])
    loss=density_loss(4096)
    return coefficients,probabilities,(log(loss.numerator)-log(loss.denominator)),data


def posterior_density_anchors(prepared,counts,witness,anchors):
    """Propose integer anchors from the tilted envelope, never proof inputs.

    Logarithmic derivatives estimate category counts per region under the
    output tilt. Any returned bounded integer anchors remain valid witnesses;
    neither these numerical derivatives nor a posterior interpretation are
    used by outward verification.
    """
    _,probabilities,_,data=prepared
    lam,t1,t2=(float(Q(x)) for x in witness);z=np.array([1.,t1,t2])
    weights=(np.asarray(counts,dtype=float)/4096/(probabilities@z))@probabilities
    _,logtau=density_tangent.floating(4096,anchors)
    weights*=np.exp(logtau)
    def moment(ws):
        scale=float(np.sum(ws));v=float((ws[1]+ws[2])/scale)
        r=float(ws[2]/(ws[1]+ws[2])) if v else 0.
        return 4096*256*log(scale)+model.log_power_moment(
            iid_kernel.floating(data,v,r,lam),16384,np.ones(2))
    result=[];step=1e-5
    for j in range(3):
        if weights[j]==0:result.append(0);continue
        plus=weights.copy();minus=weights.copy()
        plus[j]*=np.exp(step);minus[j]*=np.exp(-step)
        estimate=(moment(plus)-moment(minus))/(2*step*256)
        if not np.isfinite(estimate):raise ArithmeticError('nonfinite density anchor proposal')
        result.append(max(0,min(4096,int(round(estimate)))))
    return tuple(result)


def propose(prepared,counts,*,threshold=model.THRESHOLD,density_anchors=None,capped=False):
    coeff,probs,logloss,data=prepared;n=np.asarray(counts,dtype=float);f=n/4096
    if np.any(n<0) or np.sum(n)!=4096:raise ValueError('complete nonnegative composition required')
    logtau=np.zeros(3)
    if density_anchors is not None:
        if len(density_anchors)!=3:raise ValueError('three density anchors required')
        logloss,logtau=density_tangent.floating(4096,density_anchors);logtau=np.array(logtau)
    elif capped:
        caps=tuple(sum(int(c) for c,p in zip(counts,probs[:,j]) if p>0) for j in range(3))
        loss=capped_density_loss(4096,caps);logloss=log(loss.numerator)-log(loss.denominator)
    multiplicity=gammaln(4097)-np.sum(gammaln(n+1))
    outer=multiplicity+float(n@coeff)+256*logloss
    def objective(w):
        lam=w[0];tilt=np.exp(np.r_[0.,w[1:]])
        Z=probs@tilt;weights=((f/Z)@probs)*np.exp(logtau);scale=np.sum(weights)
        v=(weights[1]+weights[2])/scale;r=weights[2]/(weights[1]+weights[2]) if v else 0.
        inner=model.log_power_moment(iid_kernel.floating(data,v,r,lam),16384,np.ones(2))
        return outer+lam*threshold+256*(float(n@np.log(Z))+4096*log(scale))+inner
    best=None
    for lam0,t0 in ((.064,-1),(.256,-2),(1.,-1),(2.,0.)):
        fit=minimize(objective,[lam0,t0,t0],method='L-BFGS-B',bounds=[(.000001,5),(-8,8),(-8,8)],
                     options={'maxiter':160,'ftol':1e-12})
        if best is None or fit.fun<best.fun:best=fit
    witness=[Q(round(best.x[0]*10**9),10**9)]+[Q(round(np.exp(x)*10**9),10**9) for x in best.x[1:]]
    w=[float(witness[0])]+[log(float(x)) for x in witness[1:]]
    value=objective(w)/log(2)
    if not np.isfinite(value):raise ArithmeticError('nonfinite composition proposal')
    return value,witness


def outward(components,data,counts,witness,*,threshold=model.THRESHOLD,density_anchors=None,capped=False):
    if (len(counts)!=len(components) or any(type(c) is not int or c<0 for c in counts)
            or sum(counts)!=4096 or len(witness)!=3 or any(Q(x)<=0 for x in witness)):
        raise ValueError('valid complete composition and positive witness required')
    lam,t1,t2=map(Q,witness)
    if type(threshold) is not int or not 0<=threshold<2097152:raise ValueError('bounded integer output cutoff required')
    if density_anchors is not None and len(density_anchors)!=3:raise ValueError('three density anchors required')
    probabilities=[row[2:5] for row in components]
    if density_anchors is not None or capped:
        normalizers=[sum(p*t for p,t in zip(row,[1,t1,t2])) for row in probabilities]
        weights=[sum(Q(c,4096)*row[j]/z for c,row,z in zip(counts,probabilities,normalizers)) for j in range(3)]
        if density_anchors is not None:
            logloss,logtau=density_tangent.outward(4096,density_anchors)
            # Preserve monotonicity of the positive input polynomial by
            # rounding each unnormalized, tilted category weight upward.
            rounded=[]
            for w,tau in zip(weights,logtau):
                mantissa,exponent=map(int,model.up(aq(w)*tau.exp()).man_exp())
                rounded.append(Q(mantissa)*Q(2)**exponent)
            weights=rounded
        else:
            caps=tuple(sum(c for c,row in zip(counts,probabilities) if row[j]>0) for j in range(3))
            logloss=aq(capped_density_loss(4096,caps)).log()
        scale=sum(weights);v=(weights[1]+weights[2])/scale;r=weights[2]/(weights[1]+weights[2]) if v else Q(0)
        matrix=iid_kernel.outward(data,v,r,lam)**16384;moment=matrix[0,0]+matrix[0,1]
        if not moment>0:raise ArithmeticError('positive composition moment required')
        multiplicity=factorial(4096)
        for count in counts:multiplicity//=factorial(count)
        outer=arb(multiplicity).log()+sum((n*(aq(row[1]).log()+256*aq(z).log())
                   for row,n,z in zip(components,counts,normalizers) if n),arb(0))
        return model.up((outer+256*logloss+4096*256*aq(scale).log()+moment.log()+aq(lam)*threshold).exp())
    reference,factor=tilted_reference(probabilities,counts,[1,t1,t2])
    p0,p1,p2=reference;v=p1+p2;r=p2/v if v else Q(0)
    matrix=iid_kernel.outward(data,v,r,lam)**16384
    value=model.up((matrix[0,0]+matrix[0,1])*(aq(lam)*threshold).exp()*aq(factor)**256)
    multiplicity=factorial(4096)
    for count in counts:multiplicity//=factorial(count)
    value=model.up(value*multiplicity)
    for row,count in zip(components,counts):value=model.up(value*aq(row[1])**count)
    return value


def counts_from_items(components,items):
    names={row[0]:i for i,row in enumerate(components)};counts=[0]*len(components);seen=set()
    for item in items:
        name,raw=item.split('=');count=int(raw)
        if name not in names or name in seen or count<0:raise ValueError('distinct known component names and nonnegative counts required')
        counts[names[name]]=count;seen.add(name)
    if sum(counts)!=4096:raise ValueError('component counts must sum to 4096')
    return counts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,nargs='+')
    parser.add_argument('--types',nargs='+',default=None,help='Unordered row-component names, e.g. 2,2 for both central')
    parser.add_argument('--central-bits',type=int,default=130)
    parser.add_argument('--theta',default='1/4')
    parser.add_argument('--threshold',type=int,default=model.THRESHOLD)
    parser.add_argument('--updates',type=int,default=2)
    parser.add_argument('--composition',nargs='+',help='Complete component counts, e.g. 0,0=3695 0,3=400 4,4=1')
    parser.add_argument('--witness',nargs=3,help='Fixed positive output tilt and two input tilts instead of optimization')
    loss=parser.add_mutually_exclusive_group()
    loss.add_argument('--density-anchors',nargs=3,type=int)
    loss.add_argument('--capped',action='store_true')
    parser.add_argument('--outward',action='store_true')
    args=parser.parse_args()
    if args.composition and (args.groups or args.types):parser.error('explicit composition cannot be combined with groups/types')
    if args.witness and not args.composition:parser.error('a fixed witness requires one explicit composition')
    args.groups=args.groups or [129,256,512,1024,2048,4096]
    if any(not 1<=q<=4096 for q in args.groups) or not 0<=args.threshold<2097152 or not 1<=args.updates<=32:
        parser.error('valid active-label occupancy and output cutoff required')
    caps=authenticated_caps();ctx.prec=192
    rows=envelope(caps,1<<args.central_bits,Q(args.theta));components=pair_components(rows)
    images,columns,_=model.maps();data=iid_kernel.prepare(images,columns,19,updates=args.updates)
    prepared=prepare(components,data)
    print('COMPONENTS',[(name,(log(c.numerator)-log(c.denominator))/log(2),tuple(map(str,row)))
                        for name,c,*row in components],flush=True)
    print('Per-region heterogeneous density loss log2',prepared[2]/log(2),flush=True)
    selected=[i for i,row in enumerate(components) if row[-1] and (args.types is None or row[0] in args.types)]
    if args.types and len(selected)!=len(set(args.types)):parser.error('unknown active component name')
    if args.composition:
        try:cases=[('explicit',counts_from_items(components,args.composition))]
        except ValueError as error:parser.error(str(error))
    else:
        cases=[]
        for q in args.groups:
            for i in selected:
                counts=[0]*len(components);counts[0]=4096-q;counts[i]=q
                cases.append((f'q={q},type={components[i][0]}',counts))
    for label,counts in cases:
        if args.witness:witness=list(map(Q,args.witness))
        else:
            value,witness=propose(prepared,counts,threshold=args.threshold,density_anchors=args.density_anchors,capped=args.capped)
            print('MIXTURE COMPOSITION',label,'binary64 log2',value,'witness',list(map(str,witness)),flush=True)
        if args.outward:
            value=outward(components,data,counts,witness,threshold=args.threshold,density_anchors=args.density_anchors,capped=args.capped)
            print('MIXTURE COMPOSITION OUTWARD',label,'updates',args.updates,'log2 upper',model.up(value.log()/arb(2).log()),flush=True)
    print('Selected mixture compositions only; no full-code certificate.',flush=True)


if __name__=='__main__':main()
