"""One-dimensional whole-domain proof attempt for GF(16)-randomized packets.

This is a new ensemble. It does not import the lane-permutation sparse
certificates. All retained bounds use outward arithmetic on the actual maps.
"""
import argparse
from fractions import Fraction as Q
import heapq
import hashlib
import json
from math import comb
from pathlib import Path
import sys
from time import monotonic

import numpy as np
from scipy.optimize import brentq,minimize_scalar,linprog
from scipy.special import logsumexp
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'independent_rows'/'dense_closure'))
import kernel
from mixture import group_components,G,N,REGIONS,PACKETS,EPOCHS,logq,capped_density_loss
from row_mixture import envelope
from bch_joint_support import authenticated_caps
from mixture_cover import partition

SCHEMA='gf16-packet-scalar-cover-1'


def retarget_claim(record,updates=None,distance=None):
    """Select a new claim only; saved bounds must still be recomputed."""
    threshold=record['threshold'];count=record['updates']
    if type(threshold) is not int or not 0<=threshold<N or type(count) is not int or count<1:
        raise ValueError('valid saved threshold and update count required')
    if updates is not None:
        if type(updates) is not int or updates not in (2,3,4):
            raise ValueError('retargeting supports two through four updates')
        count=updates
    if distance is not None:
        delta=Q(distance)
        if not 0<delta<1:raise ValueError('retargeted relative distance must lie in (0,1)')
        threshold=int(delta*N)
    return threshold,count


def actual_components(central_bits=130,theta=Q(2,5),central_scale=Q(1),row_parity=False):
    if type(central_bits) is not int or not 1<=central_bits<=256 or Q(central_scale)<=0:
        raise ValueError('positive bounded central exponent and rational scale required')
    if type(row_parity) is not bool:raise ValueError('boolean row-parity option required')
    if row_parity:
        import parity_mixture
        return parity_mixture.actual_components(central_bits,theta,central_scale)
    return group_components(envelope(authenticated_caps(),Q(central_scale)*(1<<central_bits),Q(theta)))


def log_power(matrix,n=EPOCHS):
    """Rescaled nonnegative matrix power; proposals only."""
    v=np.zeros(len(matrix));v[0]=1.;logv=0.;a=matrix.copy();loga=0.
    while n:
        if n&1:
            v=v@a;s=v.max()
            if s<=0:return -np.inf
            v/=s;logv+=loga+np.log(s)
        n>>=1
        if n:
            a=a@a;s=a.max()
            if s<=0:return -np.inf
            a/=s;loga=2*loga+np.log(s)
    return logv+np.log(v.sum())


def multiply(a,b):
    if not 0<=a<16 or not 0<=b<16:raise ValueError('four-bit field elements required')
    result=0
    for _ in range(4):
        if b&1:result^=a
        b>>=1;a<<=1
        if a&16:a^=0x13
    return result


def probabilities(activity):
    return [1-activity,*(activity*Q(comb(4,j),15) for j in range(1,5))]


def components_after_mixing(components):
    # Equal Bernoulli activity laws with the same active label can be
    # combined by adding their positive outer coefficients.
    merged={}
    for _,coefficient,law,active in components:
        key=1-law[0],active
        merged[key]=merged.get(key,Q(0))+coefficient
    return [(c,p,a) for (p,a),c in sorted(merged.items())]


def outer_dual(logs,features,active,target,occupancy,positive):
    """Eliminate the binary occupancy multiplier before solving the mean dual.

    Floating proposals only. Every rounded dual is checked by outward().
    The remaining convex derivative is monotone; a root bracket avoids
    false two-variable L-BFGS convergence near the occupancy boundary.
    """
    active=np.asarray(active,dtype=bool);features=np.asarray(features)
    lo,hi=(0.,20000.) if positive else (-20000.,0.)
    def evaluate(eta):
        scores=logs+eta*features;mu=0.
        if not active.all():
            inactive_log=logsumexp(scores[~active]);active_log=logsumexp(scores[active])
            if occupancy>=1:mu=20000.
            elif occupancy>0:
                mu=np.clip(np.log(occupancy)-np.log1p(-occupancy)+inactive_log-active_log,0.,20000.)
        scores=scores+mu*active;norm=logsumexp(scores)
        gradient=float(np.exp(scores-norm)@features-target)
        return norm-eta*target-mu*occupancy,gradient,mu
    if evaluate(lo)[1]>=0:eta=lo
    elif evaluate(hi)[1]<=0:eta=hi
    else:eta=brentq(lambda e:evaluate(e)[1],lo,hi,xtol=1e-10)
    return eta,float(evaluate(eta)[2])


class Model:
    def __init__(self,components,data,threshold,q_min=1,tilt=Q(1,32),inner=kernel,variance_shuffle=False,variance_bins=0,regional_count=False):
        if type(q_min) is not int or not 1<=q_min<=G or type(threshold) is not int or not 0<=threshold<N or tilt<=0:
            raise ValueError('valid occupancy, output cutoff, and positive activity tilt required')
        self.components=components_after_mixing(components);self.data=data;self.inner=inner
        self.threshold=threshold;self.q_min=q_min;self.tilt=Q(tilt)
        self.active=[a for _,_,a in self.components]
        self.features=[p*self.tilt/(1-p+p*self.tilt) for _,p,_ in self.components]
        self.minimum=Q(q_min,G)*min(f for f,a in zip(self.features,self.active) if a)
        self.root=(self.minimum,Q(1));self.array=np.array([[float(f),a] for f,a in zip(self.features,self.active)])
        self.min_positive=(min(1-f for f in self.features if f<1),min(f for f in self.features if f>0))
        self.cache={}
        if type(variance_shuffle) is not bool:raise ValueError('boolean variance-shuffle option required')
        self.variance_shuffle=variance_shuffle;self.variance_cache={};self.loss_cache={}
        if type(variance_bins) is not int or not 0<=variance_bins<=64 or (variance_bins and not variance_shuffle):
            raise ValueError('variance bins require the variance-shuffle option; choose 0 through 64')
        self.variance_bins=variance_bins
        if (type(regional_count) is not bool or (regional_count and
                (not variance_bins or inner.__name__!='birth_classes' or data.get('windows')!=32))):
            raise ValueError('regional counts require birth classes, 32-packet steps, and variance partitions')
        self.regional_count=regional_count

    def family(self,tilt):
        tilt=Q(tilt)
        if tilt<=0:raise ValueError('positive activity tilt required')
        if tilt not in self.cache:
            normalizers=[1-p+p*tilt for _,p,_ in self.components]
            cs=[c*z**REGIONS for (c,_,_),z in zip(self.components,normalizers)]
            values=[[((1-p) if j==0 else p)/z for (_,p,_),z in zip(self.components,normalizers)] for j in (0,1)]
            self.cache[tilt]=(cs,values,np.array([logq(c) for c in cs]))
        return self.cache[tilt]

    def loss(self,cell):
        lo,hi=cell
        caps=tuple(max(0,min(G,int(G*x/p))) for x,p in zip((1-lo,hi),self.min_positive))
        return capped_density_loss(G,caps)

    def shuffle_dual(self,cell):
        import shuffle_variance
        if cell not in self.variance_cache:
            self.variance_cache[cell]=shuffle_variance.variance_dual(self.features,self.active,cell,Q(self.q_min,G))[1]
        return self.variance_cache[cell]

    def comparison_loss(self,cell,tilt,dual=None):
        original=self.loss(cell)
        # The scalar cell controls the mean under its base tilt only.
        # Alternate activity tilts retain the existing universal bound.
        if not self.variance_shuffle or Q(tilt)!=self.tilt or not 0<cell[0]<=cell[1]<1:return original
        import shuffle_variance
        dual=tuple(map(Q,dual)) if dual is not None else self.shuffle_dual(cell)
        key=cell,dual,ctx.prec
        if key not in self.loss_cache:
            lower,_=shuffle_variance.variance_dual(self.features,self.active,cell,Q(self.q_min,G),dual)
            if not lower:self.loss_cache[key]=original
            else:
                value=shuffle_variance.density_interval_upper(G,[G*x for x in cell],G*lower)
                m,e=value.upper().man_exp()
                self.loss_cache[key]=min(original,Q(int(m))*Q(2)**int(e))
        return self.loss_cache[key]

    def empty(self,cell):return cell[1]<self.minimum or cell[0]>1

    def weights(self,cell,tilt,duals=None):
        if Q(tilt)==self.tilt and duals is None:return [1-cell[0],cell[1]/self.tilt],None
        _,values,_=self.family(tilt)
        constraints=np.array([self.array[:,0],-self.array[:,0],-self.array[:,1]])
        caps=[float(cell[1]),-float(cell[0]),-self.q_min/G]
        if duals is not None and (len(duals)!=2 or any(len(d)!=2 for d in duals)):
            raise ValueError('two exact affine weight duals required')
        weights=[];witnesses=[]
        for j,column in enumerate(values):
            if duals is None:
                fit=linprog(-np.array(list(map(float,column))),A_ub=constraints,b_ub=caps,
                            A_eq=np.ones((1,len(column))),b_eq=[1],bounds=(0,None),method='highs')
                if fit.success:
                    d=fit.ineqlin.marginals
                    dual=[Q(round(float(d[1]-d[0])*10**9),10**9),Q(round(max(0.,-float(d[2]))*10**9),10**9)]
                else:dual=[Q(0),Q(0)]
            else:dual=list(map(Q,duals[j]))
            eta,mu=dual
            if mu<0:raise ValueError('nonnegative occupancy dual required')
            intercept=max(v-eta*f+mu*a for v,f,a in zip(column,self.features,self.active))
            weights.append(max(Q(0),intercept+max(eta*cell[0],eta*cell[1])-mu*Q(self.q_min,G)))
            witnesses.append(dual)
        return weights,witnesses

    def propose_with(self,cell,tilt):
        _,_,logs=self.family(tilt);weights,weight_duals=self.weights(cell,tilt)
        best=None
        for sign in (False,True):
            target=np.array([float(cell[0] if sign else cell[1]),self.q_min/G])
            def objective(w):
                scores=logs+self.array@w;norm=logsumexp(scores)
                return norm-target@w,np.exp(scores-norm)@self.array-target
            proposed=outer_dual(logs,self.array[:,0],self.active,*target,sign)
            dual=[Q(round(float(x)*10**8),10**8) for x in proposed]
            outside=G*objective(np.array(list(map(float,dual))))[0]
            if best is None or outside<best[0]:best=outside,dual
        outside,dual=best;scale=float(sum(weights));activity=float(weights[1]/sum(weights))
        probs=np.array(list(map(float,probabilities(Q(activity)))))
        def inner(lam):return log_power(self.inner.floating(self.data,probs,lam))+lam*self.threshold
        grid=[.000001,.0001,.001,.004,.016,.064,.256,1.,4.]
        scores=[inner(lam) for lam in grid];i=min(range(len(grid)),key=scores.__getitem__)
        fit=minimize_scalar(inner,bounds=(grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]),method='bounded')
        lam=Q(round(float(fit.x)*10**10),10**10)
        value=(outside+inner(float(lam))+PACKETS*np.log(scale)+REGIONS*logq(self.comparison_loss(cell,tilt)))/np.log(2)
        witness=dict(tilt=str(tilt),parameters=list(map(str,[lam,*dual])))
        if weight_duals is not None:witness['weights_dual']=[list(map(str,d)) for d in weight_duals]
        if self.variance_shuffle and Q(tilt)==self.tilt:witness['variance_dual']=list(map(str,self.shuffle_dual(cell)))
        return value,witness

    def proposal(self,cell):
        stop=proposal_cutoff(self)
        best=self.propose_with(cell,self.tilt)
        if best[0]>stop and cell[1]-cell[0]<=Q(1,16):
            if self.variance_bins and 0<cell[0]<=cell[1]<1:
                import variance_partition
                candidate=variance_partition.propose(self,cell,best)
                if candidate[0]<best[0]:best=candidate
                if best[0]<stop:return best
                if self.regional_count and cell[1]-cell[0]<=Q(1,1024):
                    import regional_count
                    candidate=regional_count.propose(self,cell,dict(candidate[1],
                        regional_tilted_atom=True,regional_fine_tilts=True,regional_tilted_variance=True,
                        regional_direct_counts=True,regional_exact_zero=True))
                    if candidate[0]<best[0]:best=candidate
                    if best[0]<stop:return best
                    refined=regional_count.propose(self,cell,dict(candidate[1],regional_lazy_density_through=6))
                    if refined[0]<best[0]:best=refined
                    if best[0]<stop:return best
                    combined=regional_count.propose(self,cell,dict(refined[1],regional_joint_return_through=3))
                    if combined[0]<best[0]:best=combined
                    if best[0]<stop:return best
                    classified=regional_count.propose(self,cell,dict(combined[1],
                        regional_feedback_classes_from=3,regional_feedback_classes_through=32))
                    if classified[0]<best[0]:best=classified
                    if best[0]<stop:return best
                    uniform=regional_count.propose(self,cell,dict(classified[1],
                        regional_feedback_uniform_classes=True,regional_feedback_uniform_replace=True))
                    if uniform[0]<best[0]:best=uniform
                    if best[0]<stop:return best
                    # The iid tilt is only a starting proposal. Conditional
                    # regional placement can prefer a smaller output tilt.
                    # Keep every previous candidate and retain exact rational
                    # parameters for fresh outward verification.
                    lam=Q(uniform[1]['parameters'][0])
                    for factor in (Q(9,10),Q(4,5),Q(39,40)):
                        tuned=regional_count.propose(self,cell,dict(uniform[1],
                            parameters=[str(lam*factor),*uniform[1]['parameters'][1:]]))
                        if tuned[0]<best[0]:best=tuned
                        if best[0]<stop:return best
            for tilt in (Q(1,128),Q(1,8),Q(1,2),Q(1),Q(2)):
                candidate=self.propose_with(cell,tilt)
                if candidate[0]<best[0]:best=candidate
                if best[0]<stop:break
        return best

    def outward(self,cell,witness):
        if len(witness['parameters'])!=3:raise ValueError('output tilt and two outer duals required')
        lam,eta,mu=map(Q,witness['parameters']);tilt=Q(witness['tilt'])
        if lam<=0 or mu<0:raise ValueError('positive tilt and nonnegative occupancy dual required')
        cs,_,_=self.family(tilt)
        # An alternate tilt must provide its checked weight majorants.
        if tilt!=self.tilt and 'weights_dual' not in witness:raise ValueError('alternate tilt missing weight duals')
        if self.variance_shuffle and tilt==self.tilt and 'variance_dual' not in witness:raise ValueError('base tilt missing variance dual')
        if 'regional_count_parts' in witness:
            if not self.regional_count:raise ValueError('regional count witness requires enabled model')
            import regional_count
            return regional_count.outward(self,cell,witness)[0]
        weights,_=self.weights(cell,tilt,witness.get('weights_dual'));scale=sum(weights)
        matrix=self.inner.outward(self.data,probabilities(weights[1]/scale),lam)**EPOCHS
        inner=sum((matrix[0,j] for j in range(matrix.ncols())),arb(0))
        if 'variance_partition' in witness:
            import variance_partition
            outside=variance_partition.outward_outer(self,cell,witness)
            value=inner.log()+PACKETS*kernel.aq(scale).log()+kernel.aq(lam)*self.threshold+outside.log()
            return kernel.up(value.exp())
        total=sum((kernel.aq(c)*kernel.aq(eta*f+mu*a).exp() for c,f,a in zip(cs,self.features,self.active)),arb(0))
        value=(inner.log()+PACKETS*kernel.aq(scale).log()+kernel.aq(lam)*self.threshold
               +REGIONS*kernel.aq(self.comparison_loss(cell,tilt,witness.get('variance_dual'))).log()+G*total.log()
               -kernel.aq(G*min(eta*cell[0],eta*cell[1])+mu*self.q_min))
        return kernel.up(value.exp())


def proposal_cutoff(model):
    """Search effort only; acceptance always requires a separate outward bound."""
    bits=getattr(model,'proposal_stop_bits',90)
    if type(bits) is not int or bits<=0:raise ValueError('positive integer proposal-stop bits required')
    return -bits


def run(model,max_cells,max_depth,target_bits,checkpoint=None,resume=None,coalesce=True,repropose=False):
    if type(coalesce) is not bool:raise ValueError('boolean resume-coalescing option required')
    if type(repropose) is not bool:raise ValueError('boolean reproposal option required')
    if repropose and resume is None:raise ValueError('reproposing a partition requires a saved partition')
    pending=[(0,'',model.root)];leaves={};unresolved={};visited=0
    previous=0
    if resume is not None:
        cells=partition(model,resume['leaves'],resume['unresolved']);pending=[]
        previous=resume.get('visited',0)
        if type(previous) is not int or previous<0:raise ValueError('invalid prior work count')
        for path,cell in cells.items():
            old=resume['leaves'].get(path)
            # Reconstruct coordinates and recompute every retained bound.
            # A checkpoint's accepted labels are not proof inputs.
            if old is not None and not repropose:
                upper=model.outward(cell,old['witness'])
                if not upper>0:raise ArithmeticError('positive resumed cell bound required')
                if upper<arb(2)**-target_bits:
                    leaves[path]=dict(old,proposal=float(upper.log()/arb(2).log()));continue
            heapq.heappush(pending,(len(path),path,cell))
        # Old unsuccessful subdivisions can be coalesced without removing
        # any accepted cell. Reconstruct the full partition afterwards.
        retry={path for _,path,_ in pending}
        if coalesce:
            for path in sorted(tuple(retry),key=len,reverse=True):
                current=path
                while current and current in retry:
                    sibling=current[:-1]+('1' if current[-1]=='0' else '0')
                    if sibling not in retry:break
                    retry.remove(current);retry.remove(sibling);current=current[:-1];retry.add(current)
        cells=partition(model,leaves,{path:{} for path in retry})
        pending=[(len(path),path,cells[path]) for path in retry];heapq.heapify(pending)
        print('GF16 RESUME rechecked',0 if repropose else len(resume['leaves']),
              'retained',len(leaves),'pending',len(pending),'repropose all',repropose,flush=True)
    while pending and visited<max_cells:
        depth,path,cell=heapq.heappop(pending);visited+=1
        score,witness=model.proposal(cell)
        if score<-(target_bits+2):
            upper=model.outward(cell,witness)
            if not 0<upper<arb(2)**-target_bits:raise ArithmeticError('outward cell misses target')
            leaves[path]=dict(witness=witness,proposal=float(score))
        elif depth>=max_depth:unresolved[path]=dict(cell=list(map(str,cell)),proposal=float(score),witness=witness)
        else:
            mid=sum(cell)/2
            heapq.heappush(pending,(depth+1,path+'0',(cell[0],mid)))
            heapq.heappush(pending,(depth+1,path+'1',(mid,cell[1])))
        if visited%20==0:
            print('GF16 COVER visited',previous+visited,'accepted',len(leaves),'pending',len(pending),'failed',len(unresolved),'last',score,flush=True)
            if checkpoint:checkpoint(dict(leaves=leaves,unresolved={**unresolved,**{p:dict(cell=list(map(str,c))) for _,p,c in pending}},visited=previous+visited))
    unresolved.update({p:dict(cell=list(map(str,c))) for _,p,c in pending})
    return dict(leaves=leaves,unresolved=unresolved,visited=previous+visited)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distance',default='1/20')
    parser.add_argument('--minimum-groups',type=int,default=1)
    parser.add_argument('--updates',type=int,default=2)
    parser.add_argument('--kernel',choices=('two-state','refresh','feedback-refresh','gf-refresh','weighted-return','profile-return','shape-return','conditioned-return','trimmed-return','rank-return','birth-refresh','birth-classes'),default='two-state')
    parser.add_argument('--base-tilt',default='1/32')
    parser.add_argument('--central-bits',type=int,default=130)
    parser.add_argument('--central-scale',default='1')
    parser.add_argument('--row-bias',default='2/5')
    parser.add_argument('--row-parity',action='store_true',help='Retain even BCH row weights after GF randomization')
    parser.add_argument('--variance-shuffle',action='store_true',help='Use a certified variance-sensitive shuffle comparison at the base tilt')
    parser.add_argument('--variance-bins',type=int,default=0,help='Partition variance when the base-tilt cell needs sharpening; requires --variance-shuffle')
    parser.add_argument('--regional-count',action='store_true',help='Retain conditional regional counts in difficult birth-class cells')
    parser.add_argument('--max-cells',type=int,default=300,help='Proposal budget for this invocation, additional when resuming')
    parser.add_argument('--max-depth',type=int,default=20)
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--target-bits',type=int,default=70)
    parser.add_argument('--output',type=Path)
    prior=parser.add_mutually_exclusive_group()
    prior.add_argument('--replay',type=Path)
    prior.add_argument('--resume',type=Path,help='Reload model metadata, recheck retained cells, and continue the full partition')
    parser.add_argument('--keep-unresolved-partition',action='store_true',
        help='On resume, retain previous subdivisions instead of retrying coalesced unresolved cells')
    parser.add_argument('--repropose-all',action='store_true',
        help='Resume only the partition: recompute proposals for every cell, including previously accepted cells')
    parser.add_argument('--retarget-updates',type=int,choices=(2,3,4),
        help='Resume rational witnesses against a changed inner; freshly recheck every retained leaf')
    parser.add_argument('--retarget-distance',
        help='Resume at an explicitly changed relative distance; freshly recheck every retained leaf')
    args=parser.parse_args()
    if args.precision<128 or args.max_cells<1 or args.max_depth<1 or args.target_bits<40:parser.error('valid work limits and precision required')
    if args.keep_unresolved_partition and not args.resume:parser.error('--keep-unresolved-partition requires --resume')
    if args.repropose_all and (not args.resume or not args.keep_unresolved_partition):
        parser.error('--repropose-all requires --resume and --keep-unresolved-partition')
    retarget=args.retarget_updates is not None or args.retarget_distance is not None
    if retarget and (not args.resume or not args.output or args.resume.resolve()==args.output.resolve()):
        parser.error('retargeting requires --resume and a separate --output; replay cannot change its claim')
    source=args.replay or args.resume
    raw=source.read_bytes() if source else None
    saved=json.loads(raw) if raw is not None else None
    if saved is not None and saved.get('schema')!=SCHEMA:parser.error('wrong ensemble/schema')
    provenance=None if saved is None else saved.get('retargeted_from')
    if saved is not None:
        threshold=saved['threshold'];q_min=saved['minimum_groups'];updates=saved['updates']
        args.kernel=saved.get('kernel','two-state')
        args.base_tilt=saved.get('base_tilt','1/32')
        args.central_bits=saved.get('central_bits',130)
        args.central_scale=saved.get('central_scale','1')
        args.row_bias=saved.get('row_bias','2/5')
        args.row_parity=saved.get('row_parity',False)
        args.variance_shuffle=saved.get('variance_shuffle',False)
        args.variance_bins=saved.get('variance_bins',0)
        args.regional_count=saved.get('regional_count',False)
        if args.kernel not in ('two-state','refresh','feedback-refresh','gf-refresh','weighted-return','profile-return','shape-return','conditioned-return','trimmed-return','rank-return','birth-refresh','birth-classes'):parser.error('unknown inner kernel')
        if retarget:
            try:threshold,updates=retarget_claim(saved,args.retarget_updates,args.retarget_distance)
            except ValueError as error:parser.error(str(error))
            provenance=dict(sha256=hashlib.sha256(raw).hexdigest(),updates=saved['updates'],threshold=saved['threshold'])
            print('GF16 RETARGET updates',saved['updates'],'to',updates,'cutoff',saved['threshold'],'to',threshold,
                  '; every saved leaf will be recomputed',flush=True)
    else:
        distance=Q(args.distance)
        if not 0<distance<1:parser.error('distance in (0,1) required')
        threshold=int(distance*N);q_min=args.minimum_groups;updates=args.updates
    if args.kernel=='birth-classes':
        import birth_classes
        inner=birth_classes
    elif args.kernel=='birth-refresh':
        import birth_refresh
        inner=birth_refresh
    elif args.kernel=='rank-return':
        import rank_return
        inner=rank_return
    elif args.kernel=='trimmed-return':
        import trimmed_return
        inner=trimmed_return
    elif args.kernel=='conditioned-return':
        import conditioned_return
        inner=conditioned_return
    elif args.kernel=='shape-return':
        import shape_return
        inner=shape_return
    elif args.kernel=='profile-return':
        import profile_return
        inner=profile_return
    elif args.kernel=='weighted-return':
        import weighted_return
        inner=weighted_return
    elif args.kernel=='gf-refresh':
        import gf_refresh
        inner=gf_refresh
    elif args.kernel=='feedback-refresh':
        import feedback_refresh
        inner=feedback_refresh
    elif args.kernel=='refresh':
        import refresh_kernel
        inner=refresh_kernel
    else:inner=kernel
    if (type(args.central_bits) is not int or not 1<=args.central_bits<=256
            or Q(args.central_scale)<=0 or not 0<Q(args.row_bias)<Q(1,2)):
        parser.error('positive bounded central exponent and row bias in (0,1/2) required')
    components=actual_components(args.central_bits,Q(args.row_bias),Q(args.central_scale),args.row_parity);data=inner.actual(updates);ctx.prec=args.precision
    model=Model(components,data,threshold,q_min,tilt=Q(args.base_tilt),inner=inner,
                variance_shuffle=args.variance_shuffle,variance_bins=args.variance_bins,regional_count=args.regional_count)
    print('GF16 RANDOMIZER DOMAIN',q_min,'through',G,'cutoff',threshold,'activity types',len(model.components),flush=True)
    if args.replay:
        cells=partition(model,saved['leaves'],saved['unresolved']);total=arb(0)
        for path,row in saved['leaves'].items():total=kernel.up(total+model.outward(cells[path],row['witness']))
        print('REPLAY accepted aggregate log2',total.log()/arb(2).log(),'unresolved',len(saved['unresolved']),flush=True)
        if saved['unresolved']:print('INCOMPLETE: no full certificate.',flush=True)
        elif not 0<total<arb(2)**-40:raise ArithmeticError('complete cover misses 40 bits')
        return
    def save(result):
        result.update(schema=SCHEMA,threshold=threshold,minimum_groups=q_min,updates=updates,precision=args.precision,
                      kernel=args.kernel,base_tilt=str(model.tilt),central_bits=args.central_bits,
                      central_scale=str(Q(args.central_scale)),row_bias=str(Q(args.row_bias)),row_parity=args.row_parity,
                      variance_shuffle=args.variance_shuffle,variance_bins=args.variance_bins,regional_count=args.regional_count)
        if provenance is not None:result['retargeted_from']=provenance
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    start=monotonic();result=run(model,args.max_cells,args.max_depth,args.target_bits,save,
        saved if args.resume else None,coalesce=not args.keep_unresolved_partition,
        repropose=args.repropose_all);save(result)
    partition(model,result['leaves'],result['unresolved'])
    print('GF16 COVER FINISHED accepted',len(result['leaves']),'unresolved',len(result['unresolved']),'seconds',monotonic()-start,flush=True)
    if result['unresolved']:print('INCOMPLETE: no full certificate.',flush=True)
    else:print('All composition cells covered; fresh replay and aggregate required.',flush=True)


if __name__=='__main__':main()
