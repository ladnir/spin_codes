"""Couple the inner moment and outer count through a verified affine bound.

Convexity belongs to the exact nonnegative input polynomial, not to the
two-state envelope. The envelope is used only at the cell vertices.
"""
import argparse
import json
from fractions import Fraction as Q
from itertools import product
from math import log
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
from flint import arb,ctx

import model,iid_kernel,density_tangent
from lifted_mixture import LiftedModel,support,SCHEMA as LIFTED_SCHEMA
from mixture_cover import cover,replay,G,L,PACKETS,resolve_threshold
from row_mixture import envelope,pair_components
from bch_joint_support import authenticated_caps
from probe import aq
from polytope import ClippedBox
from exact_composition import boundary_compositions
import mixture_probe

SCHEMA='two-bit-affine-mixture-cover-1'


def resolve_minimum_groups(requested,record=None,*,retarget=False):
    """A retarget may restrict the domain; old empty cells stay empty."""
    saved=512 if record is None else record['minimum_groups']
    result=saved if requested is None else requested
    if any(type(x) is not int or not 1<=x<=G for x in (saved,result)):
        raise ValueError('integer occupancy threshold in 1..4096 required')
    if record is not None and result!=saved and (not retarget or result<saved):
        raise ValueError('only retargeting may raise the occupancy threshold; lowering it needs a new cover')
    return result


def sparse_anchor_candidates(counts):
    # A rare category can disappear under the low-output tilt. Its zero
    # anchor can be substantially better than rounding its unconditioned mean.
    return tuple(product(*[(0,n) if 0<n<=4 else (n,) for n in counts]))


def vertices(cell):return list(dict.fromkeys(product(*zip(cell[::2],cell[1::2]))))


def affine_maps(values,features,active,rho,duals):
    """Componentwise upper weights as affine functions of the cell mean."""
    result=[]
    if len(values)!=len(duals):raise ValueError('one dual per category required')
    for column,dual in zip(values,duals):
        *eta,mu=map(Q,dual)
        if mu<0 or any(len(f)!=len(eta) for f in features):raise ValueError('valid affine weight dual required')
        h=max(v-sum(e*x for e,x in zip(eta,f))+mu*b for v,f,b in zip(column,features,active))
        result.append((h-mu*rho,tuple(eta)))
    return result


def weight_at(affine,point):
    offset,gradient=affine
    return offset+sum(e*x for e,x in zip(gradient,point))


def tangent_anchors(maps,cell):
    center=[(lo+hi)/2 for lo,hi in zip(cell[::2],cell[1::2])]
    corners=vertices(cell)
    # An infeasible center can give a nonpositive affine weight. Anchoring
    # only at a tiny positive floor then produces enormous exponentials.
    return [max(Q(1,10**8),weight_at(f,center),max(weight_at(f,p) for p in corners)/2) for f in maps]


def tangent_weights(maps,anchors,point):
    """Exp of a tangent upper bound on log of each affine weight."""
    if len(maps)!=len(anchors) or min(anchors)<=0:raise ValueError('positive tangent anchors required')
    return [aq(a)*aq((weight_at(f,point)-a)/a).exp() for f,a in zip(maps,anchors)]


def dyadic_upper(value):
    mantissa,exponent=map(int,model.up(value).man_exp())
    return Q(mantissa)*Q(2)**exponent


def coupled_outer_log(coefficients,features,active,n,q_min,cell,slope,dual):
    *eta,mu=map(Q,dual);slope=list(map(Q,slope))
    if mu<0 or len(slope)!=len(eta) or len(cell)!=2*len(eta):raise ValueError('compatible affine outer witness required')
    partition=sum((aq(c)*aq(sum((s/n+e)*f for s,e,f in zip(slope,eta,feature))+mu*b).exp()
                   for c,feature,b in zip(coefficients,features,active)),arb(0))
    return model.up(n*partition.log()-aq(n*support(cell,eta,False)+mu*q_min))


def outer_proposal(instance,cell,logs):
    lows=np.array(list(map(float,cell[::2])));highs=np.array(list(map(float,cell[1::2])))
    best=None
    for signs in product((False,True),repeat=3):
        target=np.r_[np.where(signs,lows,highs),instance.q_min/G]
        def objective(w):
            terms=logs+instance.array@w;norm=logsumexp(terms)
            return norm-target@w,np.exp(terms-norm)@instance.array-target
        bounds=[(0,20000) if sign else (-20000,0) for sign in signs]+[(0,20000)]
        fit=minimize(objective,np.zeros(4),method='L-BFGS-B',jac=True,bounds=bounds,
                     options={'maxiter':180,'ftol':1e-13,'gtol':1e-8})
        dual=[Q(round(float(x)*10**8),10**8) for x in fit.x]
        value=G*objective(np.array(list(map(float,dual))))[0]
        if best is None or value<best[0]:best=value,dual
    return best


class AffineModel(LiftedModel):
    clip_moments=False
    polish_tilt=False
    exact_boundary=False

    def composition_rows(self,cell,witness):
        complete=boundary_compositions(self.features,self.active,cell,G,self.q_min)
        rows=witness['compositions']
        if complete is None or not complete:raise ValueError('complete nonempty integer composition list required')
        recorded=[tuple(row['counts']) for row in rows]
        if len(set(recorded))!=len(recorded) or not set(complete)<=set(recorded):
            raise ValueError('composition witness has gaps or duplicates')
        # A parent's complete list may be reused on a child. Reconstruct the
        # child's full list and retain exactly those entries, never a subset.
        lookup={tuple(row['counts']):row for row in rows}
        return [lookup[counts] for counts in complete]

    def composition_proposal(self,cell):
        counts=boundary_compositions(self.features,self.active,cell,G,self.q_min)
        if not counts:return None
        if not hasattr(self,'composition_cache'):
            self.composition_cache={};self.composition_data=mixture_probe.prepare(self.components,self.data)
        rows=[];values=[]
        for composition in counts:
            if composition not in self.composition_cache:
                value,parameters=mixture_probe.propose(self.composition_data,composition,threshold=self.threshold,capped=True)
                best=value,dict(counts=list(composition),parameters=list(map(str,parameters)),capped=True)
                _,probabilities,_,_=self.composition_data
                z=np.array([1,float(parameters[1]),float(parameters[2])]);normalizers=probabilities@z
                weights=(np.array(composition)/G/normalizers)@probabilities
                seen=set()
                for bias in (float(parameters[0]),0.):
                    p=weights*np.exp(-bias*np.arange(3));p/=p.sum()
                    rounded=tuple(min(G,max(0,int(round(G*x)))) for x in p)
                    for anchors in sparse_anchor_candidates(rounded):
                        if anchors in seen:continue
                        seen.add(anchors)
                        value,witness=mixture_probe.propose(self.composition_data,composition,threshold=self.threshold,density_anchors=anchors)
                        if value<best[0]:best=value,dict(counts=list(composition),parameters=list(map(str,witness)),density_anchors=list(anchors))
                        if best[0]<-120:break
                    if best[0]<-120:break
                # Retune against the low-output tilted envelope. Raw input
                # means can put the single-packet anchor far above the
                # profile that dominates this bound. Only integer witnesses
                # survive; replay never trusts the numerical derivative.
                for _ in range(3):
                    if best[0]<-120 or 'density_anchors' not in best[1]:break
                    row=best[1]
                    try:
                        rounded=mixture_probe.posterior_density_anchors(
                            self.composition_data,composition,row['parameters'],row['density_anchors'])
                    except ArithmeticError:break
                    candidates=[ks for ks in sparse_anchor_candidates(rounded) if ks not in seen]
                    if not candidates:break
                    for anchors in candidates:
                        seen.add(anchors)
                        value,witness=mixture_probe.propose(self.composition_data,composition,threshold=self.threshold,density_anchors=anchors)
                        if value<best[0]:best=value,dict(counts=list(composition),parameters=list(map(str,witness)),density_anchors=list(anchors))
                self.composition_cache[composition]=best
            value,row=self.composition_cache[composition];values.append(value);rows.append(row)
        return float(logsumexp(np.array(values)*log(2))/log(2)),{'compositions':rows}

    def moment_vertices(self,cell,clipped):
        if not clipped:return vertices(cell)
        if not hasattr(self,'clipped_domain'):
            rho=Q(self.q_min,G)
            points=[tuple(scale*x for x in f) for f,b in zip(self.features,self.active) if b for scale in (rho,Q(1))]
            self.clipped_domain=ClippedBox(points)
        result=self.clipped_domain.vertices(cell)
        if not result:raise ArithmeticError('no feasible moment vertices')
        return result

    def setup_plane(self,cell,tilt,duals=None,anchors=None):
        _,coefficients,duals=self.alternate(cell,tilt,duals)
        maps=affine_maps(self.base.family(tilt)[1],self.features,self.active,Q(self.q_min,G),duals)
        if anchors is None:anchors=tangent_anchors(maps,cell)
        return coefficients,maps,duals,list(map(Q,anchors))

    def plane_proposal(self,cell,tilt,lam,density_anchors=None,*,clipped=True):
        coefficients,maps,duals,anchors=self.setup_plane(cell,tilt)
        corners=self.moment_vertices(cell,clipped)
        if clipped:
            center=tuple(sum(p[j] for p in corners)/len(corners) for j in range(3))
            anchors=[max(Q(1,10**8),weight_at(f,center),max(weight_at(f,p) for p in corners)/2) for f in maps]
        if density_anchors is None:
            loss=self.base.cell_loss(cell[:4]);logloss=log(loss.numerator)-log(loss.denominator);logtau=np.zeros(3)
        else:
            if len(density_anchors)!=3:raise ValueError('three density anchors required')
            logloss,logtau=density_tangent.floating(G,density_anchors);logtau=np.array(logtau)
        logs=[]
        for corner in corners:
            logweights=np.array([log(a.numerator)-log(a.denominator)+float((weight_at(f,corner)-a)/a) for f,a in zip(maps,anchors)])
            if not np.all(np.isfinite(logweights)) or np.max(np.abs(logweights))>10000:
                raise ArithmeticError('numerically unusable tangent; split the cell instead')
            logweights+=logtau
            logscale=logsumexp(logweights);probabilities=np.exp(logweights-logscale)
            v=probabilities[1]+probabilities[2];r=probabilities[2]/v if v else 0.
            logs.append(PACKETS*logscale+model.log_power_moment(iid_kernel.floating(self.data,v,r,float(lam)),16384,np.ones(2)))
        if not np.all(np.isfinite(logs)):raise ArithmeticError('nonfinite tangent vertex proposal')
        if clipped:
            # Any slope is valid; fit only proposes one. The verifier takes
            # the maximum residual at every exact vertex of the enclosure.
            design=np.column_stack((np.ones(len(corners)),np.array(corners,dtype=float)))
            fitted=np.linalg.lstsq(design,np.array(logs),rcond=None)[0][1:]
            slope=[Q(round(float(value)*10**6),10**6) for value in fitted]
        else:
            slope=[]
            for j,(lo,hi) in enumerate(zip(cell[::2],cell[1::2])):
                value=0. if lo==hi else (np.mean([v for p,v in zip(corners,logs) if p[j]==hi])
                                        -np.mean([v for p,v in zip(corners,logs) if p[j]==lo]))/float(hi-lo)
                slope.append(Q(round(float(value)*10**6),10**6))
        offset=max(value-float(sum(s*x for s,x in zip(slope,p))) for p,value in zip(corners,logs))
        coefficient_logs=self.base.family(tilt)[2]+self.array[:,:3]@np.array(list(map(float,slope)))/G
        outer,outer_dual=outer_proposal(self,cell,coefficient_logs)
        value=(offset+outer+float(lam)*self.threshold+L*logloss)/log(2)
        if not np.isfinite(value):raise ArithmeticError('nonfinite affine moment proposal')
        witness=dict(input_tilt=list(map(str,tilt)),weights_dual=[list(map(str,d)) for d in duals],
                     anchors=list(map(str,anchors)),slope=list(map(str,slope)),tilt=str(lam),outer_dual=list(map(str,outer_dual)))
        if density_anchors is not None:witness['density_anchors']=list(density_anchors)
        if clipped:witness['clipped']=True
        return value,{'plane':witness}

    def density_anchors(self,cell,tilt,lam):
        _,maps,_,anchors=self.setup_plane(cell,tilt)
        center=[(lo+hi)/2 for lo,hi in zip(cell[::2],cell[1::2])]
        weights=np.array([max(0.,float(weight_at(f,center))) for f in maps])
        if weights.sum()==0:weights=np.array(list(map(float,anchors)))
        result=[]
        for bias in (0.,float(lam)):
            tilted=weights*np.exp(-bias*np.arange(3));probabilities=tilted/tilted.sum()
            counts=tuple(min(G,max(0,int(round(G*p)))) for p in probabilities)
            if counts not in result:result.append(counts)
        return result

    def proposal(self,cell):
        cell=self.integer_cell(cell)
        if cell is None:raise ValueError('empty integer-composition cell')
        exact=self.composition_proposal(cell) if self.exact_boundary else None
        if exact is not None and exact[0]<-100:return exact
        value,witness=super().proposal(cell);best=value,{'old':witness}
        if exact is not None and exact[0]<best[0]:best=exact
        if value<-100 or max(cell[1]-cell[0],cell[3]-cell[2])>Q(1,16):return best
        for tilt in ((1,Q(1,4),Q(1,4)),(1,Q(1,4),Q(1,2)),(1,Q(1,4),1)):
            weights,coefficients,_=self.alternate(cell,tilt)
            _,parameters=self.proposal_with(cell,weights,coefficients,self.base.family(tilt)[2])
            for factor in (Q(1),Q(4,5),Q(6,5)):
                lam=parameters[0]*factor
                for density_anchors in [None,*self.density_anchors(cell,tilt,lam)]:
                    try:candidate=self.plane_proposal(cell,tilt,lam,density_anchors,clipped=self.clip_moments)
                    except ArithmeticError:continue
                    if candidate[0]<best[0]:best=candidate
                    if best[0]<-100:return best
        if self.polish_tilt and best[0]>=-100 and 'plane' in best[1]:
            best=self.polish_proposal(cell,best)
        return best

    def polish_proposal(self,cell,best):
        """Optimize the actual coupled witness, not the earlier uncoupled bound.

        This is proposal work only. Every final vertex and dual is replayed
        outward. Keeping the option off preserves the established search.
        """
        from scipy.optimize import minimize_scalar
        w=best[1]['plane'];tilt=list(map(Q,w['input_tilt']))
        lam=Q(w['tilt']);density=w.get('density_anchors')
        def objective(value):
            nonlocal best
            candidate_lam=Q(round(float(value)*10**9),10**9)
            if candidate_lam<=0:return float('inf')
            try:candidate=self.plane_proposal(cell,tilt,candidate_lam,density,clipped=self.clip_moments)
            except ArithmeticError:return float('inf')
            if candidate[0]<best[0]:best=candidate
            return candidate[0]
        minimize_scalar(objective,bounds=(float(lam)/3,3*float(lam)),method='bounded',
                        options={'maxiter':24,'xatol':1e-6})
        return best

    def outward(self,cell,witness):
        cell=self.integer_cell(cell)
        if cell is None:raise ValueError('empty integer-composition cell')
        if 'compositions' in witness:
            total=arb(0)
            for row in self.composition_rows(cell,witness):
                bound=mixture_probe.outward(self.components,self.data,row['counts'],row['parameters'],threshold=self.threshold,
                                          density_anchors=row.get('density_anchors'),capped=row.get('capped',False))
                total=model.up(total+bound)
            return total
        if 'old' in witness:return super().outward(cell,witness['old'])
        w=witness['plane'];lam=Q(w['tilt']);slope=list(map(Q,w['slope']));dual=list(map(Q,w['outer_dual']))
        if lam<=0 or len(slope)!=3 or len(dual)!=4 or dual[-1]<0:raise ValueError('valid moment plane and outer dual required')
        coefficients,maps,_,anchors=self.setup_plane(cell,w['input_tilt'],w['weights_dual'],w['anchors'])
        if len(anchors)!=3 or min(anchors)<=0:raise ValueError('three positive anchors required')
        if 'density_anchors' in w:
            if len(w['density_anchors'])!=3:raise ValueError('three density anchors required')
            logloss,logtau=density_tangent.outward(G,w['density_anchors'])
            tau=[s.exp() for s in logtau]
        else:
            logloss=aq(self.base.cell_loss(cell[:4])).log();tau=[arb(1)]*3
        offsets=[]
        for corner in self.moment_vertices(cell,w.get('clipped',False)):
            weights=[dyadic_upper(value*t) for value,t in zip(tangent_weights(maps,anchors,corner),tau)]
            scale=sum(weights);v=(weights[1]+weights[2])/scale;r=weights[2]/(weights[1]+weights[2])
            kernel=iid_kernel.outward(self.data,v,r,lam)**16384;moment=kernel[0,0]+kernel[0,1]
            if not moment>0:raise ArithmeticError('positive vertex moment required')
            offsets.append(model.up(moment.log()+PACKETS*aq(scale).log()-aq(sum(s*x for s,x in zip(slope,corner)))))
        outer=coupled_outer_log(coefficients,self.features,self.active,G,self.q_min,cell,slope,dual)
        result=max(offsets)+outer+aq(lam)*self.threshold+L*logloss
        return model.up(result.exp())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--minimum-groups',type=int,help='Minimum active-pair count; default 512 or the saved value');parser.add_argument('--theta',default='2/5')
    parser.add_argument('--central-bits',type=int,default=130);parser.add_argument('--input-tilt',nargs=2,default=['1/4','1/4'])
    parser.add_argument('--probe-cell',nargs=6,action='append',help='Selected rectangles only; never a full-domain certificate')
    parser.add_argument('--max-cells',type=int,default=2000);parser.add_argument('--max-depth',type=int,default=48)
    parser.add_argument('--max-unresolved',type=int,default=1);parser.add_argument('--target-bits',type=int,default=70)
    parser.add_argument('--threshold',type=int,help='Bad output weight cutoff; default 209715 or the saved cutoff')
    parser.add_argument('--updates',type=int,help='Independent transvections per step; default two or the saved value')
    parser.add_argument('--clip-domain',action='store_true',help='Experimental: exclude impossible moment vertices using exact supporting halfspaces')
    parser.add_argument('--polish-tilt',action='store_true',help='Optimize the final coupled witness on difficult cells')
    parser.add_argument('--reuse-witnesses',action='store_true',help='Recompute parent witnesses on child cells before proposing new ones')
    parser.add_argument('--exact-boundary',action='store_true',help='Sum complete small integer composition lists near the sparse boundary')
    parser.add_argument('--output',type=Path)
    saved=parser.add_mutually_exclusive_group();saved.add_argument('--replay',type=Path);saved.add_argument('--resume',type=Path)
    saved.add_argument('--retarget',type=Path,help='Change the cutoff, recheck every old leaf, and requeue failures')
    args=parser.parse_args();record=None
    if args.replay or args.resume or args.retarget:
        record=json.loads((args.replay or args.resume or args.retarget).read_text())
        if (args.resume or args.retarget) and record.get('schema')==LIFTED_SCHEMA:
            for leaf in record['leaves'].values():
                if not leaf.get('empty'):leaf['witness']={'old':leaf['witness']}
            record['schema']=SCHEMA
        if record.get('schema')!=SCHEMA:parser.error('unsupported affine schema')
        for key in ('theta','central_bits','input_tilt'):setattr(args,key,record[key])
    try:
        args.minimum_groups=resolve_minimum_groups(args.minimum_groups,record,retarget=bool(args.retarget))
        args.threshold=resolve_threshold(args.threshold,record,resume=bool(args.resume),retarget=bool(args.retarget))
        args.updates=model.resolve_updates(args.updates,record,retarget=bool(args.retarget))
    except ValueError as error:parser.error(str(error))
    if not 1<=args.minimum_groups<=G or min(args.max_cells,args.max_depth,args.max_unresolved,args.target_bits)<1:
        parser.error('bounded occupancy and positive cover limits required')
    if args.probe_cell and (record is not None or args.output):parser.error('selected probes do not save or replay full covers')
    caps=authenticated_caps();ctx.prec=192
    components=pair_components(envelope(caps,1<<args.central_bits,Q(args.theta)))
    images,columns,_=model.maps();data=iid_kernel.prepare(images,columns,19,updates=args.updates)
    instance=AffineModel(components,data,args.minimum_groups,[1,*map(Q,args.input_tilt)],threshold=args.threshold)
    instance.clip_moments=args.clip_domain
    instance.polish_tilt=args.polish_tilt
    instance.exact_boundary=args.exact_boundary
    print('OUTPUT WEIGHT THRESHOLD',args.threshold,'of',2*PACKETS,flush=True)
    print('MINIMUM ACTIVE-PAIR COUNT',args.minimum_groups,flush=True)
    print('INDEPENDENT UPDATES PER STEP',args.updates,flush=True)
    if args.probe_cell:
        for raw in args.probe_cell:
            cell=tuple(map(Q,raw))
            if instance.empty(cell):print('EMPTY SELECTED CELL',raw,flush=True);continue
            value,witness=instance.proposal(cell);bound=instance.outward(cell,witness)
            print('SELECTED AFFINE CELL',raw,'proposal',value,'log2 upper',bound.log()/arb(2).log(),flush=True)
        print('Selected rectangles only; no complete occupancy-range certificate.',flush=True);return
    if args.replay:
        bound=replay(instance,record)
        if not 0<bound<arb(2)**-40:raise ArithmeticError('replay misses dense budget')
        print('REPLAYED COMPLETE AFFINE COVER q>=',args.minimum_groups,'margin',-bound.log()/arb(2).log(),flush=True);return
    result=cover(instance,target_bits=args.target_bits,max_cells=args.max_cells,max_depth=args.max_depth,
                 max_unresolved=args.max_unresolved,resume=record,retarget=bool(args.retarget),reuse_witnesses=args.reuse_witnesses)
    result.update(schema=SCHEMA,minimum_groups=args.minimum_groups,theta=args.theta,central_bits=args.central_bits,input_tilt=args.input_tilt,threshold=args.threshold,updates=args.updates)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print('No complete dense certificate.' if result['unresolved'] else 'Complete dense range; complementary occupancies still required.',flush=True)


if __name__=='__main__':main()
