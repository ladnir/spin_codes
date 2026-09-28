"""Three-coordinate mixture cover retaining the central-row fraction.

No mixture compositions are dropped. Floating solvers propose affine and
exponential witnesses; rational inequalities and Arb verify their bounds.
"""
import argparse
import json
from fractions import Fraction as Q
from itertools import product
from math import log
from pathlib import Path

import numpy as np
from scipy.optimize import linprog,minimize,minimize_scalar
from scipy.special import logsumexp
from flint import arb,ctx

import model
import iid_kernel
from mixture_cover import CoverModel,cover,replay,children,G,L,PACKETS,resolve_threshold
from row_mixture import envelope,pair_components
from bch_joint_support import authenticated_caps
from probe import aq
from integer_projection import integer_projection_empty,boundary_composition_empty

SCHEMA='two-bit-lifted-mixture-cover-1'


def support(cell,vector,maximum):
    choose=max if maximum else min
    return sum((choose(e*lo,e*hi) for e,lo,hi in zip(vector,cell[::2],cell[1::2])),Q(0))


def affine_upper(values,features,active,rho,cell,dual):
    *eta,mu=map(Q,dual)
    if len(eta)*2!=len(cell) or mu<0:raise ValueError('compatible affine witness required')
    intercept=max(v-sum(e*x for e,x in zip(eta,f))+mu*b for v,f,b in zip(values,features,active))
    return intercept+support(cell,eta,True)-mu*rho


def outer_log(coefficients,features,active,n,q_min,cell,dual):
    *eta,mu=map(Q,dual)
    if len(eta)*2!=len(cell) or mu<0:raise ValueError('compatible exponential witness required')
    partition=sum((aq(c)*aq(sum(e*x for e,x in zip(eta,f))+mu*b).exp()
                   for c,f,b in zip(coefficients,features,active)),arb(0))
    return model.up(n*partition.log()-aq(n*support(cell,eta,False)+mu*q_min))


class LiftedModel:
    root=(Q(0),Q(1))*3

    def __init__(self,components,data,q_min,input_tilt,*,threshold=model.THRESHOLD):
        self.base=CoverModel(components,data,q_min,input_tilt,threshold=threshold)
        self.threshold=self.base.threshold
        self.components=components;self.data=data;self.q_min=q_min;self.active=self.base.active
        central=[Q(row[0].split(',').count('2'),2) for row in components]
        self.features=[(f[1],f[2],c) for f,c in zip(self.base.features,central)]
        self.array=np.array([[*map(float,f),b] for f,b in zip(self.features,self.active)])
        self.base_cache={};self.empty_cache={}
        self.minimum_positive=[min(f[j] for f in self.features if f[j]>0)/G for j in range(2)]

    def integer_cell(self,cell):
        """Contract a rectangle while retaining every integer composition.

        The central coordinate is an integer multiple of 1/(2G). A
        positive category coordinate needs at least one capable component.
        The split tree still describes the original, uncontracted cube.
        """
        if len(cell)!=6 or any(not 0<=lo<=hi<=1 for lo,hi in zip(cell[::2],cell[1::2])):
            raise ValueError('bounded three-dimensional cell required')
        result=list(cell)
        lo,hi=2*G*cell[4],2*G*cell[5]
        result[4]=Q(-(-lo.numerator//lo.denominator),2*G)
        result[5]=Q(hi.numerator//hi.denominator,2*G)
        if result[4]>result[5]:return None
        for j,minimum in enumerate(self.minimum_positive):
            if result[2*j+1]<minimum:
                if result[2*j]>0:return None
                result[2*j+1]=Q(0)
        return tuple(result)

    @staticmethod
    def children(cell):
        # Preserve the economical two-dimensional cover until its rectangles
        # become small, then split the added coordinate and all three axes.
        if max(cell[1]-cell[0],cell[3]-cell[2])>Q(1,16):
            return [c+tuple(cell[4:]) for c in children(cell[:4])]
        return children(cell)

    def constraints(self,cell):
        rows=[];caps=[]
        for j,(lo,hi) in enumerate(zip(cell[::2],cell[1::2])):
            rows.extend([self.array[:,j],-self.array[:,j]]);caps.extend([float(hi),-float(lo)])
        rows.append(-self.array[:,-1]);caps.append(-self.q_min/G)
        return np.array(rows),np.array(caps)

    @staticmethod
    def rational_dual(marginals):
        values=[marginals[2*j+1]-marginals[2*j] for j in range((len(marginals)-1)//2)]
        values.append(max(0.,-marginals[-1]))
        return [Q(round(float(x)*10**9),10**9) for x in values]

    def empty(self,cell):
        cell=self.integer_cell(cell)
        if cell is None:return True
        if integer_projection_empty(self.features,cell,G):return True
        if boundary_composition_empty(self.features,self.active,cell,G,self.q_min):return True
        if self.base.empty(cell[:4]):return True
        if cell in self.empty_cache:return self.empty_cache[cell]
        rows,caps=self.constraints(cell);count=len(self.components)
        # Minimize a common slack in the coordinate inequalities. Occupancy
        # is not relaxed. A positive float optimum is never itself evidence.
        slack=-np.ones(len(rows));slack[-1]=0
        fit=linprog(np.r_[np.zeros(count),1],A_ub=np.column_stack((rows,slack)),b_ub=caps,
                    A_eq=np.array([[*np.ones(count),0]]),b_eq=[1],bounds=(0,None),method='highs')
        result=False
        if fit.success and fit.fun>1e-10:
            dual=self.rational_dual(fit.ineqlin.marginals)
            # Every feasible composition would imply 0 <= this exact upper.
            result=affine_upper([Q(0)]*count,self.features,self.active,Q(self.q_min,G),cell,dual)<0
        self.empty_cache[cell]=result
        return result

    def alternate(self,cell,tilt,duals=None):
        coefficients,all_values,_=self.base.family(tilt)
        rows,caps=self.constraints(cell);witnesses=[];weights=[]
        if duals is not None and (len(duals)!=3 or any(len(d)!=4 for d in duals)):
            raise ValueError('three four-coordinate majorants required')
        for j in range(3):
            values=all_values[j]
            if duals is None:
                fit=linprog(-np.array(list(map(float,values))),A_ub=rows,b_ub=caps,
                            A_eq=np.ones((1,len(values))),b_eq=[1],bounds=(0,None),method='highs')
                witness=self.rational_dual(fit.ineqlin.marginals) if fit.success else [Q(0)]*4
            else:witness=list(map(Q,duals[j]))
            weights.append(affine_upper(values,self.features,self.active,Q(self.q_min,G),cell,witness))
            witnesses.append(witness)
        if min(weights)<0 or sum(weights)<=0:raise ArithmeticError('invalid comparison mass in retained cell')
        return weights,coefficients,witnesses

    def proposal_with(self,cell,weights,coefficients,logcoeff=None):
        if logcoeff is None:logcoeff=np.array([log(c.numerator)-log(c.denominator) for c in coefficients])
        lows=np.array(list(map(float,cell[::2])));highs=np.array(list(map(float,cell[1::2])))
        best=None;seen=set()
        # The cell support function is smooth in each orthant. This avoids
        # chasing an infeasible center and paying its excessive radius cost.
        prior=np.exp(logcoeff-logsumexp(logcoeff))@self.array
        queue=[tuple(prior[:3]<(lows+highs)/2)]
        while queue:
            signs=queue.pop()
            if signs in seen:continue
            seen.add(signs)
            endpoint=np.where(signs,lows,highs);target=np.r_[endpoint,self.q_min/G]
            def objective(w):
                logs=logcoeff+self.array@w;norm=logsumexp(logs)
                return norm-target@w,np.exp(logs-norm)@self.array-target
            bounds=[(0,20000) if sign else (-20000,0) for sign in signs]+[(0,20000)]
            fit=minimize(objective,np.zeros(4),method='L-BFGS-B',jac=True,bounds=bounds,
                         options={'maxiter':250,'ftol':1e-14,'gtol':1e-9})
            dual=[Q(round(float(x)*10**8),10**8) for x in fit.x]
            value=G*objective(np.array(list(map(float,dual))))[0]
            if best is None or value<best[0]:best=value,dual
            logs=logcoeff+self.array@fit.x;mean=np.exp(logs-logsumexp(logs))@self.array
            crossing=[j for j in range(3) if abs(fit.x[j])<1e-6
                      and (mean[j]<lows[j]-1e-8 or mean[j]>highs[j]+1e-8)]
            if not crossing:break
            for j in crossing:
                neighbor=list(signs);neighbor[j]=not neighbor[j];queue.append(tuple(neighbor))
            # Optimality of the proposal is not trusted by the verifier.
            # If numerical iteration cycles, finish the small sign search.
            if not any(s not in seen for s in queue):
                queue.extend(s for s in product((False,True),repeat=3) if s not in seen)
        outer,dual=best
        weights=np.array(list(map(float,weights)));scale=weights.sum()
        v=(weights[1]+weights[2])/scale;r=weights[2]/(weights[1]+weights[2]) if v else 0.
        def moment(lam):
            return model.log_power_moment(iid_kernel.floating(self.data,v,r,lam),16384,np.ones(2))+lam*self.threshold
        grid=[.000001,.004,.016,.064,.256,1.,2.,3.,5.]
        index=min(range(len(grid)),key=lambda i:moment(grid[i]))
        fit=minimize_scalar(moment,bounds=(grid[max(0,index-1)],grid[min(len(grid)-1,index+1)]),method='bounded')
        lam=Q(round(float(fit.x)*10**9),10**9);loss=self.base.cell_loss(cell[:4])
        value=(outer+moment(float(lam))+PACKETS*np.log(scale)+L*(log(loss.numerator)-log(loss.denominator)))/np.log(2)
        if not np.isfinite(value):raise ArithmeticError('nonfinite lifted proposal')
        return value,[lam,*dual]

    def proposal(self,cell):
        cell=self.integer_cell(cell)
        if cell is None:raise ValueError('empty integer-composition cell')
        small=cell[:4]
        if small not in self.base_cache:self.base_cache[small]=self.base.proposal(small)
        value,witness=self.base_cache[small]
        best=value,{'base':witness if isinstance(witness,dict) else list(map(str,witness))}
        if value<-100:return best
        for tilt in ((1,Q(1,4),Q(1,4)),(1,Q(1,4),Q(1,2)),(1,Q(1,4),Q(3,4)),(1,Q(1,4),1)):
            weights,coefficients,duals=self.alternate(cell,tilt)
            value,parameters=self.proposal_with(cell,weights,coefficients,self.base.family(tilt)[2])
            if value<best[0]:
                best=value,dict(input_tilt=list(map(str,tilt)),weights_dual=[list(map(str,d)) for d in duals],
                                parameters=list(map(str,parameters)))
            if best[0]<-100:break
        return best

    def outward(self,cell,witness):
        cell=self.integer_cell(cell)
        if cell is None:raise ValueError('empty integer-composition cell')
        if 'base' in witness:return self.base.outward(cell[:4],witness['base'])
        weights,coefficients,_=self.alternate(cell,witness['input_tilt'],witness['weights_dual'])
        lam,*dual=map(Q,witness['parameters'])
        if lam<=0 or len(dual)!=4 or dual[-1]<0:raise ValueError('valid lifted witness required')
        scale=sum(weights);v=(weights[1]+weights[2])/scale;r=weights[2]/(weights[1]+weights[2]) if v else Q(0)
        kernel=iid_kernel.outward(self.data,v,r,lam)**16384;moment=kernel[0,0]+kernel[0,1]
        if not moment>0:raise ArithmeticError('positive inner enclosure required')
        bound=(moment.log()+PACKETS*aq(scale).log()+aq(lam)*self.threshold+L*aq(self.base.cell_loss(cell[:4])).log()
               +outer_log(coefficients,self.features,self.active,G,self.q_min,cell,dual))
        return model.up(bound.exp())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--minimum-groups',type=int,default=1024)
    parser.add_argument('--theta',default='1/3');parser.add_argument('--central-bits',type=int,default=130)
    parser.add_argument('--input-tilt',nargs=2,default=['1/4','1/4'])
    parser.add_argument('--max-cells',type=int,default=5000);parser.add_argument('--max-depth',type=int,default=48)
    parser.add_argument('--max-unresolved',type=int,default=1);parser.add_argument('--target-bits',type=int,default=70)
    parser.add_argument('--threshold',type=int,help='Bad output weight cutoff; default 209715 or the saved cutoff')
    parser.add_argument('--updates',type=int,help='Independent transvections per step; default two or the saved value')
    parser.add_argument('--screen-only',action='store_true');parser.add_argument('--output',type=Path)
    saved=parser.add_mutually_exclusive_group()
    saved.add_argument('--replay',type=Path)
    saved.add_argument('--resume',type=Path,help='Reverify accepted cells and continue the saved unresolved partition')
    saved.add_argument('--retarget',type=Path,help='Change the cutoff, recheck every old leaf, and requeue failures')
    args=parser.parse_args();record=None
    if args.replay or args.resume or args.retarget:
        record=json.loads((args.replay or args.resume or args.retarget).read_text())
        if record.get('schema')!=SCHEMA:parser.error('unsupported lifted schema')
        for name in ('minimum_groups','theta','central_bits','input_tilt'):setattr(args,name,record[name])
    try:
        args.threshold=resolve_threshold(args.threshold,record,resume=bool(args.resume),retarget=bool(args.retarget))
        args.updates=model.resolve_updates(args.updates,record,retarget=bool(args.retarget))
    except ValueError as error:parser.error(str(error))
    if not 1<=args.minimum_groups<=G or min(args.max_cells,args.max_depth,args.max_unresolved,args.target_bits)<1:
        parser.error('bounded occupancy and positive cover limits required')
    caps=authenticated_caps();ctx.prec=192
    components=pair_components(envelope(caps,1<<args.central_bits,Q(args.theta)))
    images,columns,_=model.maps();data=iid_kernel.prepare(images,columns,19,updates=args.updates)
    instance=LiftedModel(components,data,args.minimum_groups,[1,*map(Q,args.input_tilt)],threshold=args.threshold)
    print('OUTPUT WEIGHT THRESHOLD',args.threshold,'of',2*PACKETS,flush=True)
    print('INDEPENDENT UPDATES PER STEP',args.updates,flush=True)
    if args.replay:
        bound=replay(instance,record)
        if not 0<bound<arb(2)**-40:raise ArithmeticError('lifted replay misses dense budget')
        print('REPLAYED COMPLETE LIFTED MIXTURE COVER q>=',args.minimum_groups,'margin',-bound.log()/arb(2).log(),flush=True)
        return
    result=cover(instance,target_bits=args.target_bits,max_cells=args.max_cells,max_depth=args.max_depth,
                 max_unresolved=args.max_unresolved,screen_only=args.screen_only,resume=record,retarget=bool(args.retarget))
    result.update(schema=SCHEMA,minimum_groups=args.minimum_groups,theta=args.theta,central_bits=args.central_bits,input_tilt=args.input_tilt,threshold=args.threshold,updates=args.updates)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print('No complete dense certificate.' if result['unresolved'] or args.screen_only else 'Complete dense range; complementary occupancies still required.',flush=True)


if __name__=='__main__':main()
