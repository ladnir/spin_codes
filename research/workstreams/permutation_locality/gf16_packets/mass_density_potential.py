"""Positive-function bounds retaining a shared lazy mass budget.

Each vector bounds future cost at zero, its nonzero supremum, its
nonzero uniform average, and its supremum on each expansion class.
This is a nonlinear backward bound, NOT a comparison matrix for measures.
The command below is an iid diagnostic, not a regional certificate.
"""
import argparse
from fractions import Fraction as Q
import json
from math import comb
from pathlib import Path
import numpy as np
from flint import arb,ctx
import fiber_density as fiber

aq,up=fiber.aq,fiber.up


def shared_budget(mass,caps,values):
    """Outward maximum of sum x_i v_i, 0<=x_i<=caps_i, sum x_i<=mass."""
    if len(caps)!=len(values) or not mass>=0 or any(not x>=0 for x in (*caps,*values)):
        raise ValueError('matching nonnegative budgets and costs required')
    caps=list(map(up,caps));values=list(map(up,values))
    order=sorted(range(len(values)),key=lambda i:values[i],reverse=True)
    left=up(mass);result=arb(0)
    for i in order:
        take=min(left,up(caps[i]));result=up(result+take*values[i])
        left=max(arb(0),up(left-take))
    return result


class Bound:
    def __init__(self,data,source,z,alternatives=()):
        """Source must precede allocations of its lazy M branch to U/F."""
        import lazy_density
        self.W=data['windows'];self.n=3+len(data['birth_class_levels'])
        self.linear=[source,*alternatives]
        if any(len(table)!=self.W+1 or any(m.nrows()!=self.n or m.ncols()!=self.n for m in table)
               for table in self.linear):raise ValueError('matching birth-class operators required')
        sizes=list(map(int,data['birth_class_counts']));L=(1<<data['bits'])-1
        if any(s<=0 for s in sizes) or sum(sizes)!=L:raise ValueError('complete positive class census required')
        profiles=fiber.profile_caps(data,z);g=lazy_density.density_caps(data,z)
        levels=list(map(int,data['birth_class_levels']))
        selected=[[i for i,w in enumerate(data['image_histogram_weights']) if w==l] for l in levels]
        alpha=arb(2)**-data['updates'];self.mass=[];self.density=[];self.caps=[];self.fixed=[]
        for j,m in enumerate(source):
            maxima=[max(profiles[j]),
                up(sum((int(count)*cap for count,cap in zip(data['histogram_multiplicities'],profiles[j])),arb(0))/L),
                *(max(profiles[j][i] for i in ids) for ids in selected)]
            D=[arb(0),*(up(alpha*L*k) for k in maxima)]
            D[2]=min(D[2],up(alpha*g[j]))
            H=[arb(0),*(m[i,1] if j else arb(0) for i in range(1,self.n))]
            fixed=m*arb(1)
            if j:
                for i in range(1,self.n):fixed[i,1]=arb(0)
            self.mass.append(H);self.density.append(D)
            self.caps.append([[up(d*s/L) for s in sizes] for d in D]);self.fixed.append(fixed)
        self.float_linear=np.array([[[[float(m[i,k]) for k in range(self.n)] for i in range(self.n)] for m in table]
                                    for table in self.linear])
        self.float_fixed=np.array([[[float(m[i,k]) for k in range(self.n)] for i in range(self.n)] for m in self.fixed])
        self.float_mass=np.array(self.mass,dtype=float);self.float_density=np.array(self.density,dtype=float)
        self.float_caps=np.array(self.caps,dtype=float)

    def outward(self,j,values):
        if type(j) is not int or not 0<=j<=self.W or len(values)!=self.n or any(not x>=0 for x in values):
            raise ValueError('valid occupancy and nonnegative future-cost vector required')
        v=list(map(up,values));result=[]
        for i in range(self.n):
            previous=min(up(sum((m[j][i,k]*v[k] for k in range(self.n)),arb(0))) for m in self.linear)
            lazy=min(up(self.mass[j][i]*v[1]),up(self.density[j][i]*v[2]),
                shared_budget(self.mass[j][i],self.caps[j][i],v[3:]))
            value=up(sum((self.fixed[j][i,k]*v[k] for k in range(self.n)),arb(0))+lazy)
            result.append(min(previous,value))
        return result

    def floating(self,j,values):
        """Batch proposals. Every retained result requires outward replay."""
        v=np.asarray(values,dtype=float);single=v.ndim==1
        if single:v=v[None,:]
        if (v.ndim!=2 or v.shape[1]!=self.n or not 0<=j<=self.W
                or not np.isfinite(v).all() or (v<0).any()):raise ValueError('valid nonnegative vectors required')
        previous=np.minimum.reduce([v@m.T for m in self.float_linear[:,j]])
        order=np.argsort(-v[:,3:],axis=1);left=np.broadcast_to(self.float_mass[j],(len(v),self.n)).copy()
        budget=np.zeros_like(left)
        for k in range(self.n-3):
            index=order[:,k];take=np.minimum(left,self.float_caps[j][:,index].T)
            budget+=take*v[np.arange(len(v)),index+3,None];left-=take
        lazy=np.minimum(np.minimum(self.float_mass[j]*v[:,1,None],self.float_density[j]*v[:,2,None]),budget)
        result=np.minimum(previous,v@self.float_fixed[j].T+lazy)
        return result[0] if single else result


def iid_proposal(bound,p,steps=16384,iterations=200):
    weights=np.array([comb(bound.W,j)*p**j*(1-p)**(bound.W-j) for j in range(bound.W+1)])
    v=np.ones(bound.n);best=None
    for iteration in range(iterations):
        image=sum((w*bound.floating(j,v) for j,w in enumerate(weights)),np.zeros(bound.n))
        rho=float(np.max(image/v))
        score=(steps*np.log(rho)+np.log(v[0]/v.min()))/np.log(2)
        if best is None or score<best['log2']:best=dict(log2=float(score),potential=v.tolist(),rho=rho,iteration=iteration)
        new=image/image.max()
        if not np.isfinite(new).all() or new.min()<=0:raise ArithmeticError('positive potential iteration failed')
        if np.max(np.abs(np.log(new/v)))<1e-12:break
        v=new
    return best


def main():
    import birth_classes
    import scalar_cover as sc
    import regional_count as regional
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('--activities',nargs='+',default=['.35','.4'])
    parser.add_argument('--precision',type=int,default=256);parser.add_argument('--output',type=Path)
    args=parser.parse_args();source=json.loads(args.input.read_text())
    if source.get('schema')!='gf16-regional-count-screen-1' or not source.get('diagnostic_only'):
        parser.error('regional floating screen required')
    if args.precision<128:parser.error('precision at least 128 required')
    data=fiber.attach(birth_classes.actual(source['updates']));ctx.prec=args.precision
    model=sc.Model(sc.actual_components(129,Q(21,50),Q(1001,1000),True),data,source['threshold'],
        source['minimum_groups'],tilt=Q(3,16),inner=birth_classes,variance_shuffle=True,variance_bins=16,regional_count=True)
    rows=[]
    for p in map(Q,args.activities):
        matching=[row for row in source['rows'] if Q(row['activity'])==p]
        if not matching:parser.error('requested activity missing')
        row=min(matching,key=lambda r:r['proposal']);witness=row['witness'];lam=Q(witness['parameters'][0])
        baseline=regional.local_operators(model,witness)
        raw={k:v for k,v in witness.items() if k not in (
            'regional_lazy_density_through','regional_feedback_classes_from','regional_feedback_classes_through',
            'regional_feedback_uniform_classes','regional_feedback_uniform_replace')}
        original=regional.local_operators(model,raw);bound=Bound(data,original,(-aq(lam)).exp(),[baseline])
        weights=np.array([comb(bound.W,j)*float(p)**j*(1-float(p))**(bound.W-j) for j in range(bound.W+1)])
        old=sc.log_power(np.tensordot(weights,bound.float_linear[1],axes=1))/np.log(2)
        proposed=iid_proposal(bound,float(p));result=dict(activity=str(p),tilt=str(lam),baseline_inner_log2=float(old),
            potential=proposed,gain_bits=float(old-proposed['log2']))
        rows.append(result);print('SHARED MASS IID',json.dumps(result),flush=True)
        if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,scope='iid inner only',
            precision=args.precision,updates=source['updates'],rows=rows),indent=2)+'\n')
    print('IID diagnostic only; no regional placement, outer count, or complete certificate.',flush=True)


if __name__=='__main__':main()
