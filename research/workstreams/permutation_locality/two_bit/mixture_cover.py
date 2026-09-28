"""Two-dimensional cover of every row-mixture composition above q_min.

The cover sums whole cells using a multinomial generating function. It
does not enumerate or sample the 15-component count simplex. Every cell
uses one fixed input tilt and one rational Chernoff witness. Saved bounds
are not trusted by the replay verifier.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path
from math import log

import numpy as np
from scipy.optimize import minimize,minimize_scalar,linprog
from scipy.special import logsumexp
from flint import arb,ctx

import model
import iid_kernel
from row_mixture import envelope,pair_components
from heterogeneous import density_loss,capped_density_loss
from bch_joint_support import authenticated_caps
from probe import aq

G=4096
L=256
PACKETS=G*L
SCHEMA='two-bit-mixture-cover-1'


def checked_threshold(value):
    if type(value) is not int or not 0<=value<2*PACKETS:
        raise ValueError('integer output-weight threshold in [0,N) required')
    return value


def resolve_threshold(requested,record=None,*,resume=False,retarget=False):
    """Old records used 209715; replay never silently changes their claim."""
    saved=checked_threshold(model.THRESHOLD if record is None else record.get('threshold',model.THRESHOLD))
    result=saved if requested is None else checked_threshold(requested)
    if retarget and (record is None or requested is None):
        raise ValueError('retarget requires a saved partition and an explicit threshold')
    if record is not None and not retarget and (result>saved or (not resume and result!=saved)):
        raise ValueError('replay must preserve the threshold; resume may only lower it')
    return result


def orientation(a,b,p):return (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])


def convex_hull(points):
    points=sorted(set(points))
    if len(points)<3:raise ValueError('a two-dimensional comparison domain is required')
    def half(sequence):
        result=[]
        for point in sequence:
            while len(result)>1 and orientation(result[-2],result[-1],point)<=0:result.pop()
            result.append(point)
        return result
    return half(points)[:-1]+half(reversed(points))[:-1]


def separated(cell,polygon):
    x0,x1,y0,y1=cell;corners=[(x,y) for x in (x0,x1) for y in (y0,y1)]
    if max(p[0] for p in polygon)<x0 or min(p[0] for p in polygon)>x1:return True
    if max(p[1] for p in polygon)<y0 or min(p[1] for p in polygon)>y1:return True
    return any(max(orientation(a,b,p) for p in corners)<0 for a,b in zip(polygon,polygon[1:]+polygon[:1]))


def outer_log_bound(coefficients,features,active,n,q_min,cell,dual):
    e1,e2,mu=map(Q,dual)
    if mu<0:raise ValueError('nonnegative occupancy dual required')
    partition=sum((aq(c)*aq(e1*f[1]+e2*f[2]+mu*a).exp()
                   for c,f,a in zip(coefficients,features,active)),arb(0))
    x0,x1,y0,y1=cell
    minimum=min(e1*x0,e1*x1)+min(e2*y0,e2*y1)
    return model.up(n*partition.log()-aq(n*minimum+mu*q_min))


def linear_upper(values,features,active,rho,cell,dual):
    """Exact affine majorant, valid for every composition in the cell."""
    e1,e2,mu=map(Q,dual)
    if mu<0:raise ValueError('nonnegative occupancy dual required')
    intercept=max(v-e1*f[1]-e2*f[2]+mu*b for v,f,b in zip(values,features,active))
    x0,x1,y0,y1=cell
    return intercept+max(e1*x0,e1*x1)+max(e2*y0,e2*y1)-mu*rho


class CoverModel:
    def __init__(self,components,data,q_min,input_tilt,*,threshold=model.THRESHOLD):
        self.components=components;self.data=data;self.q_min=q_min
        self.threshold=checked_threshold(threshold)
        self.tilt=tuple(map(Q,input_tilt))
        if not 1<=q_min<=G or len(self.tilt)!=3 or self.tilt[0]!=1 or min(self.tilt)<=0:
            raise ValueError('bounded active-label threshold and positive normalized input tilt required')
        self.normalizers=[sum(p*t for p,t in zip(row[2:5],self.tilt)) for row in components]
        self.features=[tuple(p*t/z for p,t in zip(row[2:5],self.tilt))
                       for row,z in zip(components,self.normalizers)]
        self.coefficients=[row[1]*z**L for row,z in zip(components,self.normalizers)]
        self.active=[row[-1] for row in components]
        self.minimum_density=Q(q_min,G)*min(f[1]+f[2] for f,a in zip(self.features,self.active) if a)
        rho=Q(q_min,G)
        self.domain=convex_hull([(scale*f[1],scale*f[2]) for f,a in zip(self.features,self.active) if a for scale in (rho,Q(1))])
        self.loss=density_loss(G)
        self.logcoeff=np.array([log(c.numerator)-log(c.denominator) for c in self.coefficients])
        self.array=np.array([[float(f[1]),float(f[2]),a] for f,a in zip(self.features,self.active)])
        self.logloss=float(aq(self.loss).log())
        self.family_cache={}

    def family(self,tilt):
        tilt=tuple(map(Q,tilt))
        if len(tilt)!=3 or tilt[0]!=1 or min(tilt)<=0:raise ValueError('positive normalized input tilt required')
        if tilt not in self.family_cache:
            z=[sum(p*t for p,t in zip(row[2:5],tilt)) for row in self.components]
            coefficients=[row[1]*zi**L for row,zi in zip(self.components,z)]
            values=[[row[2+j]/zi for row,zi in zip(self.components,z)] for j in range(3)]
            self.family_cache[tilt]=(coefficients,values,np.array([log(c.numerator)-log(c.denominator) for c in coefficients]))
        return self.family_cache[tilt]

    def alternate(self,cell,tilt,duals=None):
        coefficients,all_values,_=self.family(tilt)
        x0,x1,y0,y1=map(float,cell)
        constraints=np.array([self.array[:,0],-self.array[:,0],self.array[:,1],-self.array[:,1],-self.array[:,2]])
        caps=[x1,-x0,y1,-y0,-self.q_min/G]
        witnesses=[];weights=[]
        for j in range(3):
            values=all_values[j]
            if duals is None:
                fit=linprog(-np.array(list(map(float,values))),A_ub=constraints,b_ub=caps,
                            A_eq=np.ones((1,len(values))),b_eq=[1],bounds=(0,None),method='highs')
                if fit.success:
                    d=fit.ineqlin.marginals
                    witness=[Q(round(v*10**9),10**9) for v in (d[1]-d[0],d[3]-d[2],max(0.,-d[4]))]
                else:witness=[Q(0)]*3
            else:witness=list(map(Q,duals[j]))
            weights.append(linear_upper(values,self.features,self.active,Q(self.q_min,G),cell,witness))
            witnesses.append(witness)
        if min(weights)<0 or sum(weights)<=0:raise ArithmeticError('invalid comparison mass in retained cell')
        return tuple(weights),coefficients,witnesses

    def empty(self,cell):
        x0,x1,y0,y1=cell
        return x0+y0>1 or x1+y1<self.minimum_density or separated(cell,self.domain) or sum(self.category_caps(cell))<G

    def category_caps(self,cell):
        x0,x1,y0,y1=cell;maxima=(1-x0-y0,x1,y1)
        # A type with a positive category probability contributes at least
        # this minimum to its mean feature. Hence only this many slots can
        # be capable of producing that category, under any positive tilt.
        return tuple(max(0,min(G,int(G*m/min(f[j] for f in self.features if f[j]>0))))
                     for j,m in enumerate(maxima))

    def cell_loss(self,cell):
        caps=self.category_caps(cell)
        if min(caps)>= (G+2)//3:return self.loss
        return capped_density_loss(G,tuple(sorted(caps)))

    def upper_weights(self,cell):
        x0,x1,y0,y1=cell
        if self.empty(cell):raise ValueError('empty comparison cell')
        return (1-x0-y0,x1/self.tilt[1],y1/self.tilt[2])

    def proposal_with(self,cell,weights,logcoeff):
        x0,x1,y0,y1=map(float,cell);center=np.array([(x0+x1)/2,(y0+y1)/2])
        radius=np.array([(x1-x0)/2,(y1-y0)/2])
        # The center dual is smooth. The final bound pays the entire cell's
        # support function, so a center outside the feasible polygon is safe.
        def dual(w):
            logs=logcoeff+self.array@w;norm=logsumexp(logs)
            value=G*(norm-center@w[:2])-self.q_min*w[2]
            gradient=G*(np.exp(logs-norm)@self.array-np.r_[center,self.q_min/G])
            return value,gradient
        fit=minimize(dual,np.zeros(3),method='L-BFGS-B',jac=True,
                     bounds=[(-20000,20000),(-20000,20000),(0,20000)],options={'maxiter':200,'ftol':1e-13,'gtol':1e-7})
        eta=[Q(round(float(x)*10**8),10**8) for x in fit.x]
        e=np.array(list(map(float,eta)))
        outer=G*(logsumexp(logcoeff+self.array@e)-center@e[:2]+radius@np.abs(e[:2]))-self.q_min*e[2]
        weights=np.array(list(map(float,weights)));scale=weights.sum()
        v=(weights[1]+weights[2])/scale;r=weights[2]/(weights[1]+weights[2]) if v else 0.
        def inner(lam):
            matrix=iid_kernel.floating(self.data,v,r,lam)
            return model.log_power_moment(matrix,16384,np.ones(2))+lam*self.threshold
        grid=[.000001,.004,.016,.064,.256,1.,2.,3.,5.]
        scores=[inner(t) for t in grid];best=min(range(len(scores)),key=scores.__getitem__)
        low,high=grid[max(0,best-1)],grid[min(len(grid)-1,best+1)]
        opt=minimize_scalar(inner,bounds=(low,high),method='bounded')
        lam=Q(round(float(opt.x)*10**9),10**9)
        loss=self.cell_loss(cell)
        value=(outer+inner(float(lam))+PACKETS*np.log(scale)+L*(log(loss.numerator)-log(loss.denominator)))/np.log(2)
        if not np.isfinite(value):raise ArithmeticError('nonfinite cover proposal')
        return value,[lam,*eta]

    def proposal(self,cell):
        best=self.proposal_with(cell,self.upper_weights(cell),self.logcoeff)
        x0,x1,y0,y1=cell
        if best[0]>-100 and max(x1-x0,y1-y0)<=Q(1,16):
            for tilt in ((1,Q(1,4),Q(1,2)),(1,Q(1,4),Q(3,8)),(1,Q(1,4),Q(5,8)),
                         (1,Q(1,4),Q(3,4)),(1,Q(1,4),Q(1)),(1,Q(3,10),Q(1,2)),
                         (1,Q(2,5),Q(4,5)),(1,Q(1,5),Q(1,4))):
                weights,coefficients,duals=self.alternate(cell,tilt)
                logcoeff=self.family(tilt)[2]
                value,parameters=self.proposal_with(cell,weights,logcoeff)
                if value<best[0]:
                    best=value,dict(input_tilt=list(map(str,tilt)),weights_dual=[list(map(str,d)) for d in duals],
                                    parameters=list(map(str,parameters)))
                if best[0]<-100:break
        return best

    def outward(self,cell,witness):
        if isinstance(witness,dict):
            if len(witness['weights_dual'])!=3 or any(len(d)!=3 for d in witness['weights_dual']):
                raise ValueError('three affine weight majorants required')
            weights,coefficients,_=self.alternate(cell,witness['input_tilt'],witness['weights_dual'])
            witness=witness['parameters']
        else:weights,coefficients=self.upper_weights(cell),self.coefficients
        lam,e1,e2,mu=map(Q,witness)
        if lam<=0 or mu<0:raise ValueError('positive output tilt and nonnegative occupancy dual required')
        scale=sum(weights);v=(weights[1]+weights[2])/scale
        r=weights[2]/(weights[1]+weights[2]) if v else Q(0)
        matrix=iid_kernel.outward(self.data,v,r,lam)**16384
        inner=matrix[0,0]+matrix[0,1]
        if not inner>0:raise ArithmeticError('positive inner enclosure required')
        # The categorical weights are per two-bit packet, not per bit.
        logbound=(inner.log()+PACKETS*aq(scale).log()+aq(lam)*self.threshold+L*aq(self.cell_loss(cell)).log()
                  +outer_log_bound(coefficients,self.features,self.active,G,self.q_min,cell,(e1,e2,mu)))
        return model.up(logbound.exp())


def children(cell):
    if not cell or len(cell)%2:raise ValueError('nonempty rectangular cell required')
    axis=max(range(len(cell)//2),key=lambda j:cell[2*j+1]-cell[2*j])
    mid=(cell[2*axis]+cell[2*axis+1])/2
    left,right=list(cell),list(cell);left[2*axis+1]=mid;right[2*axis]=mid
    return [tuple(left),tuple(right)]


def partition(model_,leaves,unresolved):
    """Reconstruct all terminal cells; saved cell coordinates are not trusted."""
    if set(leaves)&set(unresolved):raise ValueError('cell is both accepted and unresolved')
    terminals=set(leaves)|set(unresolved)
    if not terminals or any(not isinstance(p,str) or set(p)-{'0','1'} for p in terminals):
        raise ValueError('nonempty binary split-tree partition required')
    root=getattr(model_,'root',(Q(0),Q(1),Q(0),Q(1)));split=getattr(model_,'children',children)
    pending=[('',root)];visited=0;cells={}
    while pending:
        path,cell=pending.pop();visited+=1
        if visited>2*len(terminals)-1:raise ValueError('incomplete atlas partition')
        if path in terminals:cells[path]=cell
        else:
            for digit,child in reversed(list(enumerate(split(cell)))):pending.append((path+str(digit),child))
    if set(cells)!=terminals:raise ValueError('overlapping or unreachable atlas cells')
    return cells


def replay_leaves(model_,leaves,cells):
    total=arb(0)
    for path,leaf in leaves.items():
        cell=cells[path]
        # A stronger exact separation may now prune an old accepted cell.
        if model_.empty(cell):continue
        if leaf.get('empty'):raise ValueError('incorrectly pruned cell')
        else:
            bound=model_.outward(cell,leaf['witness'])
            if not bound>0:raise ArithmeticError('positive outward cell bound required')
            total=model.up(total+bound)
    return total


def cover(model_,*,target_bits=70,max_cells=20000,max_depth=28,max_unresolved=32,screen_only=False,resume=None,retarget=False,reuse_witnesses=False):
    root=getattr(model_,'root',(Q(0),Q(1),Q(0),Q(1)))
    split=getattr(model_,'children',children)
    pending=[('',root)];leaves={};unresolved={};total=arb(0);visited=0;seeds={};reused=0
    if reuse_witnesses and screen_only:raise ValueError('witness reuse requires outward verification')
    if retarget and resume is None:raise ValueError('retarget requires a saved partition')
    if resume is not None:
        if resume.get('screen_only') or screen_only:raise ValueError('resume requires outward-verified cells')
        cells=partition(model_,resume['leaves'],resume['unresolved'])
        retry=set(resume['unresolved'])
        if reuse_witnesses:
            seeds={p:leaf['witness'] for p,leaf in resume['unresolved'].items() if 'witness' in leaf}
        if retarget:
            # Raising the cutoff can invalidate every old accepted witness.
            # Recompute each bound and requeue failures; keep the whole tree.
            for index,(path,leaf) in enumerate(resume['leaves'].items(),1):
                cell=cells[path]
                if model_.empty(cell):leaves[path]={'empty':True}
                else:
                    if leaf.get('empty'):raise ValueError('incorrectly pruned cell')
                    bound=model_.outward(cell,leaf['witness'])
                    if not bound>0:raise ArithmeticError('positive outward cell bound required')
                    if bound<arb(2)**-target_bits:
                        leaves[path]=leaf;total=model.up(total+bound)
                    else:
                        retry.add(path)
                        if reuse_witnesses:seeds[path]=leaf['witness']
                if index%100==0:print('RETARGET RECHECKED',index,'retained',len(leaves),'retry',len(retry),flush=True)
        else:
            leaves=dict(resume['leaves'])
            total=replay_leaves(model_,leaves,cells)
        pending=[(path,cells[path]) for path in reversed(cells) if path in retry]
        print('REPLAYED RESUME accepted',len(leaves),'pending',len(pending),flush=True)
    while pending and visited<max_cells and len(unresolved)<max_unresolved:
        path,cell=pending.pop();visited+=1
        if model_.empty(cell):leaves[path]={'empty':True};continue
        if reuse_witnesses and path in seeds:
            witness=seeds.pop(path);bound=model_.outward(cell,witness)
            if not bound>0:raise ArithmeticError('positive reused witness bound required')
            if bound<arb(2)**-target_bits:
                total=model.up(total+bound);leaves[path]={'witness':witness};reused+=1
                continue
        value,witness=model_.proposal(cell)
        if value<-(target_bits+2):
            bound=None if screen_only else model_.outward(cell,witness)
            if bound is not None and not 0<bound<arb(2)**-target_bits:
                raise ArithmeticError('outward cell misses its budget')
            if bound is not None:total=model.up(total+bound)
            leaves[path]={'witness':witness if isinstance(witness,dict) else list(map(str,witness)),'proposal':value}
        elif len(path)>=max_depth:
            unresolved[path]={'cell':list(map(str,cell)),'proposal':value,
                              'witness':witness if isinstance(witness,dict) else list(map(str,witness))}
        else:
            for digit,child in reversed(list(enumerate(split(cell)))):
                next_path=path+str(digit);pending.append((next_path,child))
                if reuse_witnesses:seeds[next_path]=witness
        if visited%100==0:
            print('MIXTURE ATLAS cells',visited,'passed',sum('witness' in l for l in leaves.values()),
                  'pending',len(pending),'unresolved',len(unresolved),'last log2',value,
                  'root volume remaining',sum(2.**-len(p) for p,_ in pending)+sum(2.**-len(p) for p in unresolved),flush=True)
    for path,cell in pending:
        unresolved[path]={'cell':list(map(str,cell))}
        if path in seeds:unresolved[path]['witness']=seeds[path]
    print('MIXTURE ATLAS FINISHED visited',visited,'leaves',len(leaves),'unresolved',len(unresolved),flush=True)
    if reuse_witnesses:print('CHILD CELLS VERIFIED WITH REUSED WITNESSES',reused,flush=True)
    if not screen_only and total>0:print('OUTWARD COVERED-CELL AGGREGATE margin',-total.log()/arb(2).log(),flush=True)
    return dict(leaves=leaves,unresolved=unresolved,screen_only=screen_only,
                upper_dyadic=None if screen_only else [int(x) for x in total.upper().man_exp()])


def replay(model_,record):
    if record.get('screen_only') or record.get('unresolved'):raise ValueError('only complete outward atlases may be replayed as complete')
    cells=partition(model_,record['leaves'],{})
    return replay_leaves(model_,record['leaves'],cells)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--minimum-groups',type=int,default=1024)
    parser.add_argument('--theta',default='1/3')
    parser.add_argument('--central-bits',type=int,default=130)
    parser.add_argument('--input-tilt',nargs=2,default=['1/4','1/4'])
    parser.add_argument('--target-bits',type=int,default=70)
    parser.add_argument('--threshold',type=int,help='Bad output weight cutoff; default 209715 or the saved cutoff')
    parser.add_argument('--updates',type=int,help='Independent transvections per step; default two or the saved value')
    parser.add_argument('--max-cells',type=int,default=20000)
    parser.add_argument('--max-depth',type=int,default=28)
    parser.add_argument('--max-unresolved',type=int,default=32,help='Stop after this many depth-limited failures for targeted diagnosis')
    parser.add_argument('--screen-only',action='store_true')
    parser.add_argument('--output',type=Path)
    saved=parser.add_mutually_exclusive_group()
    saved.add_argument('--replay',type=Path)
    saved.add_argument('--resume',type=Path,help='Reverify accepted cells and continue the saved unresolved partition')
    saved.add_argument('--retarget',type=Path,help='Change the cutoff, recheck every old leaf, and requeue failures')
    args=parser.parse_args()
    record=None
    if args.replay or args.resume or args.retarget:
        record=json.loads((args.replay or args.resume or args.retarget).read_text())
        if record.get('schema')!=SCHEMA:parser.error('unsupported atlas schema')
        args.minimum_groups=record['minimum_groups'];args.theta=record['theta'];args.central_bits=record['central_bits'];args.input_tilt=record['input_tilt']
    try:
        args.threshold=resolve_threshold(args.threshold,record,resume=bool(args.resume),retarget=bool(args.retarget))
        args.updates=model.resolve_updates(args.updates,record,retarget=bool(args.retarget))
    except ValueError as error:parser.error(str(error))
    if not 1<=args.minimum_groups<=G or args.target_bits<1 or args.max_cells<1 or args.max_depth<1 or args.max_unresolved<1:
        parser.error('valid occupancy threshold and positive cover limits required')
    caps=authenticated_caps();ctx.prec=192
    rows=envelope(caps,1<<args.central_bits,Q(args.theta));components=pair_components(rows)
    images,columns,_=model.maps();data=iid_kernel.prepare(images,columns,19,updates=args.updates)
    model_=CoverModel(components,data,args.minimum_groups,[1,*map(Q,args.input_tilt)],threshold=args.threshold)
    print('OUTPUT WEIGHT THRESHOLD',args.threshold,'of',2*PACKETS,flush=True)
    print('INDEPENDENT UPDATES PER STEP',args.updates,flush=True)
    if args.replay:
        bound=replay(model_,record)
        if not 0<bound<arb(2)**-40:raise ArithmeticError('replayed dense atlas does not close its target')
        print('REPLAYED COMPLETE DENSE MIXTURE ATLAS q>=',args.minimum_groups,'margin',-bound.log()/arb(2).log(),flush=True)
        return
    result=cover(model_,target_bits=args.target_bits,max_cells=args.max_cells,max_depth=args.max_depth,
                 max_unresolved=args.max_unresolved,screen_only=args.screen_only,resume=record,retarget=bool(args.retarget))
    result.update(schema=SCHEMA,minimum_groups=args.minimum_groups,theta=args.theta,central_bits=args.central_bits,input_tilt=args.input_tilt,threshold=args.threshold,updates=args.updates)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    if result['unresolved'] or args.screen_only:print('No complete dense-range certificate from this run.',flush=True)
    else:print('Complete dense-range atlas; full code still requires the complementary sparse range.',flush=True)


if __name__=='__main__':main()
