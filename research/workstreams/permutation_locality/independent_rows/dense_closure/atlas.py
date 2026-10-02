"""Four-dimensional, complete composition cover for four-bit row mixtures.

Floating point chooses witnesses; Arb replays every retained cell. An
unfinished split tree is deliberately not a minimum-distance certificate.
"""
import argparse
from fractions import Fraction as Q
import json
from itertools import product
from pathlib import Path

import numpy as np
from flint import arb,ctx
from scipy.optimize import minimize,minimize_scalar,linprog
from scipy.special import logsumexp
from scipy.spatial import ConvexHull

import kernel
import boundary
import clipping
import selected as composition_probe
from mixture import (actual_components,G,REGIONS,PACKETS,EPOCHS,N,logq,log_power,
                     density_loss,capped_density_loss)
from mixture_cover import cover,replay,partition
import density_tangent
import shuffle_comparison

SCHEMA='four-bit-dense-composition-cover-1'


def resolve_minimum_groups(requested,record=None,*,retarget=False):
    """Retargeting may restrict the domain, never silently enlarge it."""
    saved=59 if record is None else record['minimum_groups']
    result=saved if requested is None else requested
    if any(type(x) is not int or not 1<=x<=G for x in (saved,result)):
        raise ValueError('integer occupancy threshold in 1..2048 required')
    if record is not None and result!=saved and (not retarget or result<saved):
        raise ValueError('only retargeting may raise the occupancy threshold; lowering it needs a new cover')
    return result


class Model:
    root=tuple([Q(0),Q(1)]*4)
    coupled_alternates=False
    posterior_shuffle=False

    def __init__(self,components,data,q_min,threshold,input_tilt):
        if (type(q_min) is not int or not 1<=q_min<=G or type(threshold) is not int
                or not 0<=threshold<N or len(input_tilt)!=5 or input_tilt[0]!=1
                or min(input_tilt)<=0):
            raise ValueError('valid geometry, threshold, and positive normalized tilt required')
        self.components=components;self.data=data;self.q_min=q_min;self.threshold=threshold
        self.tilt=tuple(map(Q,input_tilt))
        z=[sum(p*t for p,t in zip(row[2],self.tilt)) for row in components]
        self.features=[tuple(p*t/zi for p,t in zip(row[2],self.tilt)) for row,zi in zip(components,z)]
        self.coefficients=[row[1]*zi**REGIONS for row,zi in zip(components,z)]
        self.active=[row[3] for row in components]
        self.array=np.array([[*map(float,f[1:]),a] for f,a in zip(self.features,self.active)])
        self.logcoeff=np.array([logq(c) for c in self.coefficients])
        self.minimum_density=Q(q_min,G)*min(sum(f[1:]) for f,a in zip(self.features,self.active) if a)
        self.min_category=[min(f[j] for f in self.features if f[j]>0) for j in range(5)]
        self.loss=density_loss(G,5)
        # Float hull normals only propose separating directions. Each
        # intercept is recomputed on every exact rational vertex, so the
        # prune does not trust Qhull tolerances or its facet completeness.
        points=sorted({tuple(scale*x for x in f[1:]) for f,a in zip(self.features,self.active)
                       if a for scale in (Q(q_min,G),Q(1))})
        hull=ConvexHull(np.array([list(map(float,p)) for p in points]))
        normals=sorted({tuple(Q(round(float(x)*10**8),10**8) for x in row[:4]) for row in hull.equations})
        self.facets=[(e,max(sum(x*y for x,y in zip(e,p)) for p in points)) for e in normals]
        self.facet_array=np.array([[*map(float,e),float(b)] for e,b in self.facets])
        self.boundary_cache={};self.composition_cache={};self.vertex_cache={};self.family_cache={};self.seeds={}
        self.affine_cache={}

    def family(self,tilt):
        tilt=tuple(map(Q,tilt))
        if len(tilt)!=5 or tilt[0]!=1 or min(tilt)<=0:raise ValueError('positive normalized category tilt required')
        if tilt not in self.family_cache:
            z=[sum(p*t for p,t in zip(row[2],tilt)) for row in self.components]
            coefficients=[row[1]*zi**REGIONS for row,zi in zip(self.components,z)]
            values=[[row[2][j]/zi for row,zi in zip(self.components,z)] for j in range(5)]
            self.family_cache[tilt]=(coefficients,values,np.array([logq(c) for c in coefficients]))
        return self.family_cache[tilt]

    def alternate(self,cell,tilt,duals=None):
        coefficients,values,logs=self.family(tilt)
        constraints=np.array([sign*self.array[:,j] for j in range(4) for sign in (1,-1)]+[-self.array[:,4]])
        caps=[float(v) for lo,hi in zip(cell[::2],cell[1::2]) for v in (hi,-lo)]+[-self.q_min/G]
        if duals is not None and (len(duals)!=5 or any(len(d)!=5 for d in duals)):
            raise ValueError('five affine weight majorants required')
        weights=[];witnesses=[]
        for j,column in enumerate(values):
            if duals is None:
                fit=linprog(-np.array(list(map(float,column))),A_ub=constraints,b_ub=caps,
                            A_eq=np.ones((1,len(column))),b_eq=[1],bounds=(0,None),method='highs')
                if fit.success:
                    d=fit.ineqlin.marginals
                    dual=[Q(round(float(x)*10**9),10**9) for x in
                          [*(d[2*k+1]-d[2*k] for k in range(4)),max(0.,-d[8])]]
                else:dual=[Q(0)]*5
            else:dual=list(map(Q,duals[j]))
            *eta,mu=dual
            if mu<0:raise ValueError('nonnegative occupancy dual required')
            intercept=max(v-sum(e*x for e,x in zip(eta,f[1:]))+mu*a
                          for v,f,a in zip(column,self.features,self.active))
            upper=intercept+sum(max(e*lo,e*hi) for e,lo,hi in zip(eta,cell[::2],cell[1::2]))-mu*Q(self.q_min,G)
            weights.append(max(Q(0),upper));witnesses.append(dual)
        return weights,coefficients,logs,witnesses

    def moment_vertices(self,cell):
        if cell not in self.vertex_cache:
            points=clipping.vertices(cell,self.facets)
            # A proved-empty box needs no bound. Retaining the full box
            # here is also safe and avoids a special witness convention.
            self.vertex_cache[cell]=points or tuple(product(*zip(cell[::2],cell[1::2])))
        return self.vertex_cache[cell]

    def boundary(self,cell):
        if cell not in self.boundary_cache:
            self.boundary_cache[cell]=boundary.compositions([f[1:] for f in self.features],self.active,cell,G,self.q_min)
        return self.boundary_cache[cell]

    def category_caps(self,cell):
        maxima=(1-sum(cell[::2]),*cell[1::2])
        return tuple(max(0,min(G,int(G*x/p))) for x,p in zip(maxima,self.min_category))

    def empty(self,cell):
        if (sum(cell[::2])>1 or sum(cell[1::2])<self.minimum_density
                or sum(self.category_caps(cell))<G):return True
        e=self.facet_array[:,:4];b=self.facet_array[:,4]
        lo=np.array(list(map(float,cell[::2])));hi=np.array(list(map(float,cell[1::2])))
        candidates=np.flatnonzero(np.sum(np.where(e>=0,e*lo,e*hi),axis=1)>b-1e-12)
        for j in candidates:
            eta,bound=self.facets[j]
            if sum(min(x*l,x*h) for x,l,h in zip(eta,cell[::2],cell[1::2]))>bound:return True
        return self.boundary(cell)==()

    def cell_loss(self,cell):
        caps=self.category_caps(cell)
        return self.loss if min(caps)>=(G+4)//5 else capped_density_loss(G,tuple(sorted(caps)))

    def upper_weights(self,cell):
        return (1-sum(cell[::2]),*(x/t for x,t in zip(cell[1::2],self.tilt[1:])))

    def outer_proposal(self,cell,logcoeff=None):
        if logcoeff is None:logcoeff=self.logcoeff
        lows=np.array(list(map(float,cell[::2])));highs=np.array(list(map(float,cell[1::2])))
        best=None
        # Optimize the actual box support function, not a center-point dual
        # with its box penalty added afterward. Each orthant is smooth.
        for signs in product((False,True),repeat=4):
            target=np.r_[np.where(signs,lows,highs),self.q_min/G]
            def objective(w):
                logs=logcoeff+self.array@w;norm=logsumexp(logs)
                return norm-target@w,np.exp(logs-norm)@self.array-target
            fit=minimize(objective,np.zeros(5),jac=True,method='L-BFGS-B',
                         bounds=[(0,20000) if sign else (-20000,0) for sign in signs]+[(0,20000)],
                         options={'maxiter':160,'ftol':1e-13,'gtol':1e-8})
            witness=[Q(round(float(x)*10**8),10**8) for x in fit.x]
            value=G*objective(np.array(list(map(float,witness))))[0]
            if best is None or value<best[0]:best=value,witness
        return best

    def weights_at(self,point):
        return (1-sum(point),*(x/t for x,t in zip(point,self.tilt[1:])))

    def affine_weights(self,cell,tilt,duals=None):
        """Keep each alternate-weight majorant affine until the moment bound.

        The exact intercept is checked against all component types. Only
        the nonnegative occupancy multiplier is replaced by q_min/G.
        """
        key=cell,tuple(map(Q,tilt))
        proposed=duals is None
        if proposed and key in self.affine_cache:return self.affine_cache[key]
        _,coefficients,logs,duals=self.alternate(cell,tilt,duals)
        _,values,_=self.family(tilt)
        maps=[]
        for column,dual in zip(values,duals):
            *eta,mu=map(Q,dual)
            intercept=max(v-sum(e*x for e,x in zip(eta,f[1:]))+mu*a
                          for v,f,a in zip(column,self.features,self.active))
            maps.append((intercept-mu*Q(self.q_min,G),tuple(eta)))
        result=maps,coefficients,logs,duals
        if proposed:self.affine_cache[key]=result
        return result

    @staticmethod
    def mapped_weights(maps,point):
        return tuple(h+sum(e*x for e,x in zip(eta,point)) for h,eta in maps)

    def plane_proposal(self,cell,lam,density_anchors=None,*,input_tilt=None):
        coefficient_logs=self.logcoeff
        weights_at=self.weights_at
        if input_tilt is not None:
            maps,_,coefficient_logs,duals=self.affine_weights(cell,input_tilt)
            weights_at=lambda p:self.mapped_weights(maps,p)
        corners=self.moment_vertices(cell)
        center=[sum(p[j] for p in corners)/len(corners) for j in range(4)]
        center_weights=weights_at(center)
        upper_weights=[max(weights_at(p)[j] for p in corners) for j in range(5)]
        anchors=[max(Q(1,10**8),x,u/2) for x,u in zip(center_weights,upper_weights)]
        if density_anchors is None:logloss=logq(self.cell_loss(cell));logtau=np.zeros(5)
        else:logloss,logtau=density_tangent.floating(G,density_anchors)
        logs=[]
        for point in corners:
            logweights=np.array([logq(a)+float((w-a)/a) for w,a in zip(weights_at(point),anchors)])+logtau
            if not np.all(np.isfinite(logweights)) or np.max(np.abs(logweights))>10000:
                raise ArithmeticError('numerically unusable tangent; split the cell instead')
            logscale=logsumexp(logweights);probabilities=np.exp(logweights-logscale)
            logs.append(PACKETS*logscale+log_power(kernel.floating(self.data,probabilities,float(lam))))
        design=np.column_stack((np.ones(len(corners)),np.array(corners,dtype=float)))
        fit=np.linalg.lstsq(design,np.array(logs),rcond=None)[0][1:]
        slope=[Q(round(float(x)*10**6),10**6) for x in fit]
        offset=max(value-float(sum(s*x for s,x in zip(slope,p))) for p,value in zip(corners,logs))
        outer,dual=self.outer_proposal(cell,coefficient_logs+self.array[:,:4]@np.array(list(map(float,slope)))/G)
        witness=dict(tilt=str(lam),anchors=list(map(str,anchors)),slope=list(map(str,slope)),dual=list(map(str,dual)),clipped=True)
        if input_tilt is not None:
            witness.update(input_tilt=list(map(str,input_tilt)),weights_dual=[list(map(str,d)) for d in duals])
        if density_anchors is not None:witness['density_anchors']=list(density_anchors)
        return (offset+outer+float(lam)*self.threshold+REGIONS*logloss)/np.log(2),dict(plane=witness)

    def moment_proposal(self,weights):
        weights=np.array(list(map(float,weights)));scale=weights.sum()
        probabilities=weights/scale
        def objective(lam):
            return log_power(kernel.floating(self.data,probabilities,lam))+lam*self.threshold
        grid=[.000001,.004,.016,.064,.256,1.,3.,6.]
        scores=[objective(t) for t in grid];best=min(range(len(grid)),key=scores.__getitem__)
        lo,hi=grid[max(0,best-1)],grid[min(len(grid)-1,best+1)]
        opt=minimize_scalar(objective,bounds=(lo,hi),method='bounded',options={'xatol':1e-8})
        lam=Q(round(float(opt.x)*10**9),10**9)
        return objective(float(lam))+PACKETS*np.log(scale),lam

    def posterior_anchors(self,weights,lam,anchors):
        _,tau=density_tangent.floating(G,anchors)
        weights=np.array(list(map(float,weights)))*np.exp(tau)
        def moment(w):
            scale=w.sum()
            return PACKETS*np.log(scale)+log_power(kernel.floating(self.data,w/scale,float(lam)))
        result=[];step=1e-5
        for j in range(5):
            if weights[j]==0:result.append(0);continue
            plus=weights.copy();minus=weights.copy()
            plus[j]*=np.exp(step);minus[j]*=np.exp(-step)
            value=(moment(plus)-moment(minus))/(2*step*REGIONS)
            if not np.isfinite(value):raise ArithmeticError('nonfinite anchor proposal')
            result.append(max(0,min(G,int(round(value)))))
        return tuple(result)

    def retune_seed(self,cell,witness):
        if isinstance(witness,dict) and ('plane' in witness or 'compositions' in witness):return None
        weights=self.upper_weights(cell);logs=self.logcoeff;anchors=None
        if isinstance(witness,dict):
            parameters=witness['parameters'];anchors=witness.get('density_anchors')
            if 'input_tilt' in witness:
                weights,_,logs,_=self.alternate(cell,witness['input_tilt'],witness['weights_dual'])
        else:parameters=witness
        _,*dual=map(Q,parameters);eta=np.array(list(map(float,dual)))
        minimum=sum(min(e*lo,e*hi) for e,lo,hi in zip(dual[:4],cell[::2],cell[1::2]))
        outside=G*(logsumexp(logs+self.array@eta)-float(minimum))-self.q_min*eta[4]
        if anchors is None:logloss=logq(self.cell_loss(cell))
        else:
            logloss,logtau=density_tangent.floating(G,anchors)
            weights=np.array(list(map(float,weights)))*np.exp(logtau)
        inside,lam=self.moment_proposal(weights)
        parameters=list(map(str,[lam,*dual]))
        updated=dict(witness,parameters=parameters) if isinstance(witness,dict) else parameters
        return (outside+inside+REGIONS*logloss)/np.log(2),updated

    def composition_proposal(self,counts):
        key=counts,self.posterior_shuffle
        if key not in self.composition_cache:
            options=[{}]
            if self.posterior_shuffle:
                options.append(dict(posterior_shuffle=True))
                try:shuffle_comparison.structure(self.components,counts)
                except ValueError:pass
                else:options.append(dict(exact_shuffle=True))
            best=None
            for mode in options:
                result=composition_probe.selected(self.data,self.components,counts,self.threshold,**mode)
                raw=result['tilt']
                parameters=[Q(round(raw[0]*10**9),10**9),
                            *(Q(round(float(np.exp(x))*10**12),10**12) for x in raw[1:])]
                row=dict(counts=list(counts),parameters=list(map(str,parameters)),**mode)
                if best is None or result['log2_upper']<best[0]:best=result['log2_upper'],row
            self.composition_cache[key]=best
        return self.composition_cache[key]

    def composition_cell_proposal(self,cell):
        compositions=self.boundary(cell)
        if not compositions or len(compositions)>(32 if self.posterior_shuffle else 2):return None
        values=[self.composition_proposal(counts) for counts in compositions]
        return float(logsumexp([v*np.log(2) for v,_ in values])/np.log(2)),dict(compositions=[row for _,row in values])

    def proposal(self,cell):
        if cell in self.seeds:
            candidate=self.retune_seed(cell,self.seeds[cell])
            if candidate is not None and candidate[0]<-100:return candidate
        candidate=self.composition_cell_proposal(cell)
        if candidate is not None and candidate[0]<-100:return candidate
        outer,dual=self.outer_proposal(cell)
        weights=self.upper_weights(cell)
        moment,lam=self.moment_proposal(weights)
        value=(outer+moment+REGIONS*logq(self.cell_loss(cell)))/np.log(2)
        if not np.isfinite(value):raise ArithmeticError('nonfinite proposal')
        best=value,[lam,*dual]
        alternate_candidates=[]
        if value>-100:
            for tilt in ((1,Q(1,32),Q(1,8),Q(1,4),1),
                         (1,Q(1,16),Q(1,4),Q(1,2),1),
                         (1,Q(1,16),Q(1,16),Q(1,16),1),
                         (1,Q(1,16),Q(1,16),Q(1,16),Q(1,16)),
                         (1,Q(1,8),Q(1,4),Q(1,2),1),(1,1,1,1,1)):
                ws,_,logs,duals=self.alternate(cell,tilt)
                outside,eta=self.outer_proposal(cell,logs)
                inside,lambda_=self.moment_proposal(ws)
                candidate=(outside+inside+REGIONS*logq(self.cell_loss(cell)))/np.log(2)
                alternate_candidates.append((candidate,tilt,lambda_,ws))
                if candidate<-100:
                    return candidate,dict(input_tilt=list(map(str,tilt)),weights_dual=[list(map(str,d)) for d in duals],
                                          parameters=list(map(str,[lambda_,*eta])))
        if value>-100 and max(hi-lo for lo,hi in zip(cell[::2],cell[1::2]))<=Q(1,32):
            p=np.array(list(map(float,weights)));p/=p.sum();seen=set()
            for bias in (0.,float(lam)):
                tilted=p*np.exp(-bias*np.arange(5));tilted/=tilted.sum()
                rounded=tuple(int(round(G*x)) for x in tilted)
                for anchors in (rounded,tuple(0 if x<=4 else x for x in rounded)):
                    if anchors in seen:continue
                    seen.add(anchors)
                    logloss,logtau=density_tangent.floating(G,anchors)
                    w=np.array(list(map(float,weights)))*np.exp(logtau)
                    moment,lambda_=self.moment_proposal(w)
                    candidate=(outer+moment+REGIONS*logloss)/np.log(2)
                    if candidate<best[0]:
                        best=candidate,dict(parameters=list(map(str,[lambda_,*dual])),density_anchors=list(anchors))
                    if best[0]<-100:return best
            # Conditioning on low output can shift category frequencies far
            # from the raw reference mean. These derivatives choose integer
            # witnesses only; the proof does not rely on their accuracy.
            current=best[1]
            anchors=current.get('density_anchors',tuple(int(round(G*x)) for x in p)) if isinstance(current,dict) else tuple(int(round(G*x)) for x in p)
            lambda_=Q(current['parameters'][0]) if isinstance(current,dict) else current[0]
            for _ in range(4):
                anchors=self.posterior_anchors(weights,lambda_,anchors)
                if anchors in seen:break
                seen.add(anchors)
                logloss,logtau=density_tangent.floating(G,anchors)
                moment,lambda_=self.moment_proposal(np.array(list(map(float,weights)))*np.exp(logtau))
                candidate=(outer+moment+REGIONS*logloss)/np.log(2)
                if candidate<best[0]:best=candidate,dict(parameters=list(map(str,[lambda_,*dual])),density_anchors=list(anchors))
                if best[0]<-100:return best
        if best[0]>-100 and max(hi-lo for lo,hi in zip(cell[::2],cell[1::2]))<=Q(1,16):
            witness=best[1]
            density=witness.get('density_anchors') if isinstance(witness,dict) else None
            lambda_=Q(witness['parameters'][0]) if isinstance(witness,dict) else witness[0]
            for factor in (Q(1),Q(4,5),Q(6,5)):
                candidate=self.plane_proposal(cell,lambda_*factor,density)
                if candidate[0]<best[0]:best=candidate
                if best[0]<-100:break
        if (self.coupled_alternates and best[0]>-100
                and max(hi-lo for lo,hi in zip(cell[::2],cell[1::2]))<=Q(1,16)):
            for _,tilt,lambda_,ws in sorted(alternate_candidates)[:3]:
                probabilities=np.array(list(map(float,ws)));probabilities/=probabilities.sum()
                counts=tuple(int(round(G*p)) for p in probabilities)
                anchor_options=[None,counts,tuple(0 if x<=4 else x for x in counts)]
                for density in dict.fromkeys(anchor_options):
                    for factor in (Q(1),Q(4,5),Q(6,5)):
                        try:candidate=self.plane_proposal(cell,lambda_*factor,density,input_tilt=tilt)
                        except ArithmeticError:continue
                        if candidate[0]<best[0]:best=candidate
                        if best[0]<-100:return best
        return best

    def outward(self,cell,witness):
        if isinstance(witness,dict) and 'plane' in witness:return self.outward_plane(cell,witness['plane'])
        if isinstance(witness,dict) and 'compositions' in witness:
            complete=self.boundary(cell)
            rows=witness['compositions']
            if complete is None or sorted(tuple(row['counts']) for row in rows)!=list(complete):
                raise ValueError('saved composition list is not the complete cell')
            return kernel.up(sum((composition_probe.outward(self.data,self.components,row['counts'],
                               self.threshold,row['parameters'],exact_shuffle=row.get('exact_shuffle',False),
                               posterior_shuffle=row.get('posterior_shuffle',False)) for row in rows),arb(0)))
        anchors=None
        weights=self.upper_weights(cell);coefficients=self.coefficients
        if isinstance(witness,dict):
            if 'input_tilt' in witness:
                weights,coefficients,_,_=self.alternate(cell,witness['input_tilt'],witness['weights_dual'])
            anchors=witness.get('density_anchors');witness=witness['parameters']
            if anchors is not None and len(anchors)!=5:raise ValueError('five density anchors required')
        if len(witness)!=6:raise ValueError('output tilt and five outer dual coordinates required')
        lam,*eta,mu=map(Q,witness)
        if lam<=0 or mu<0:raise ValueError('positive output tilt and nonnegative occupancy dual required')
        if anchors is None:logloss=kernel.aq(self.cell_loss(cell)).log()
        else:
            logloss,logtau=density_tangent.outward(G,anchors)
            uppers=[kernel.up(kernel.aq(w)*t.exp()).man_exp() for w,t in zip(weights,logtau)]
            weights=[Q(int(m))*Q(2)**int(e) for m,e in uppers]
        scale=sum(weights)
        probabilities=[w/scale for w in weights]
        matrix=kernel.outward(self.data,probabilities,lam)**EPOCHS
        inner=matrix[0,0]+matrix[0,1]
        minimum=sum(min(e*lo,e*hi) for e,lo,hi in zip(eta,cell[::2],cell[1::2]))
        partition=sum((kernel.aq(c)*kernel.aq(sum(e*x for e,x in zip(eta,f[1:]))+mu*a).exp()
                       for c,f,a in zip(coefficients,self.features,self.active)),arb(0))
        value=(inner.log()+PACKETS*kernel.aq(scale).log()+kernel.aq(lam)*self.threshold
               +REGIONS*logloss+G*partition.log()
               -kernel.aq(G*minimum+mu*self.q_min))
        return kernel.up(value.exp())

    def outward_plane(self,cell,witness):
        lam=Q(witness['tilt']);anchors=list(map(Q,witness['anchors']))
        slope=list(map(Q,witness['slope']));*eta,mu=map(Q,witness['dual'])
        if lam<=0 or len(anchors)!=5 or min(anchors)<=0 or len(slope)!=4 or len(eta)!=4 or mu<0:
            raise ValueError('valid affine moment and outer witnesses required')
        weights_at=self.weights_at;coefficients=self.coefficients
        if 'input_tilt' in witness:
            maps,coefficients,_,_=self.affine_weights(cell,witness['input_tilt'],witness['weights_dual'])
            weights_at=lambda p:self.mapped_weights(maps,p)
        if 'density_anchors' in witness:
            if len(witness['density_anchors'])!=5:raise ValueError('five density anchors required')
            logloss,logtau=density_tangent.outward(G,witness['density_anchors'])
        else:logloss=kernel.aq(self.cell_loss(cell)).log();logtau=[arb(0)]*5
        offsets=[]
        points=self.moment_vertices(cell) if witness.get('clipped') else product(*zip(cell[::2],cell[1::2]))
        for point in points:
            uppers=[kernel.up(kernel.aq(a)*(kernel.aq((w-a)/a)+t).exp()).man_exp()
                    for w,a,t in zip(weights_at(point),anchors,logtau)]
            weights=[Q(int(m))*Q(2)**int(e) for m,e in uppers];scale=sum(weights)
            matrix=kernel.outward(self.data,[w/scale for w in weights],lam)**EPOCHS
            moment=matrix[0,0]+matrix[0,1]
            offsets.append(kernel.up(moment.log()+PACKETS*kernel.aq(scale).log()
                                     -kernel.aq(sum(s*x for s,x in zip(slope,point)))))
        minimum=sum(min(e*lo,e*hi) for e,lo,hi in zip(eta,cell[::2],cell[1::2]))
        partition=sum((kernel.aq(c)*kernel.aq(sum((s/G+e)*x for s,e,x in zip(slope,eta,f[1:]))+mu*a).exp()
                       for c,f,a in zip(coefficients,self.features,self.active)),arb(0))
        value=(max(offsets)+kernel.aq(lam)*self.threshold+REGIONS*logloss+G*partition.log()
               -kernel.aq(G*minimum+mu*self.q_min))
        return kernel.up(value.exp())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--minimum-groups',type=int,help='Minimum active-group count; default 59 or the saved value')
    parser.add_argument('--threshold',type=int,help='Integer bad-weight cutoff; default 20971 or the saved value')
    parser.add_argument('--input-tilt',nargs=4,default=['1/4','1/16','1/64','1/256'])
    parser.add_argument('--theta',default='2/5')
    parser.add_argument('--central-bits',type=int,default=130)
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--target-bits',type=int,default=70)
    parser.add_argument('--max-cells',type=int,default=10000)
    parser.add_argument('--max-depth',type=int,default=56)
    parser.add_argument('--max-unresolved',type=int,default=8)
    parser.add_argument('--screen-only',action='store_true')
    parser.add_argument('--coupled-alternates',action='store_true',help='Couple alternate input tilts to outer counts with affine moment bounds')
    parser.add_argument('--posterior-shuffle',action='store_true',help='Use sharper shuffle bounds on completely enumerated boundary cells')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--replay',type=Path)
    parser.add_argument('--resume',type=Path)
    parser.add_argument('--retarget',type=Path,help='Change the cutoff, recheck accepted cells, and retry failures')
    args=parser.parse_args();ctx.prec=args.precision
    if sum(bool(x) for x in (args.replay,args.resume,args.retarget))>1:parser.error('choose replay, resume, or retarget')
    record=None
    if args.replay or args.resume or args.retarget:
        record=json.loads((args.replay or args.resume or args.retarget).read_text())
        if record['schema']!=SCHEMA:parser.error('unsupported schema')
        if args.replay and (record.get('unresolved') or record.get('screen_only')):
            parser.error('replay requires a complete outward atlas')
        if args.retarget:
            if args.threshold is None:parser.error('retarget requires an explicit threshold')
        elif args.threshold is not None and args.threshold!=record['threshold']:
            parser.error('changing the saved cutoff requires retarget')
        else:args.threshold=record['threshold']
        for key in ('input_tilt','theta','central_bits'):
            setattr(args,key,record[key])
    elif args.threshold is None:args.threshold=20971
    try:args.minimum_groups=resolve_minimum_groups(args.minimum_groups,record,retarget=bool(args.retarget))
    except ValueError as error:parser.error(str(error))
    components=actual_components(args.central_bits,Q(args.theta));data=kernel.actual()
    # Authentication uses its own precision. Restore the requested replay
    # precision after loading the exact integer/rational inputs.
    ctx.prec=args.precision
    model=Model(components,data,args.minimum_groups,args.threshold,[Q(1),*map(Q,args.input_tilt)])
    model.coupled_alternates=args.coupled_alternates
    model.posterior_shuffle=args.posterior_shuffle
    if record is not None:
        cells=partition(model,record['leaves'],record['unresolved'])
        for path,leaf in {**record['leaves'],**record['unresolved']}.items():
            if 'witness' in leaf:model.seeds[cells[path]]=leaf['witness']
    print('FOUR-BIT DENSE DOMAIN q>=',args.minimum_groups,'threshold',args.threshold,'two updates',flush=True)
    print('OUTWARD PRECISION',ctx.prec,flush=True)
    if args.replay:
        bound=replay(model,record)
        if not 0<bound<arb(2)**-40:raise ArithmeticError('complete replay misses the 40-bit target')
        print('REPLAYED COMPLETE DENSE DOMAIN margin',-bound.log()/arb(2).log(),flush=True)
        return
    result=cover(model,target_bits=args.target_bits,max_cells=args.max_cells,max_depth=args.max_depth,
                 max_unresolved=args.max_unresolved,screen_only=args.screen_only,resume=record,retarget=bool(args.retarget))
    result.update(schema=SCHEMA,minimum_groups=args.minimum_groups,threshold=args.threshold,
                  input_tilt=args.input_tilt,theta=args.theta,central_bits=args.central_bits)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    if result['unresolved'] or args.screen_only:print('No complete certificate from this run.',flush=True)
    else:print('Dense domain closed; combine only with the same-setup complementary sparse proof.',flush=True)


if __name__=='__main__':main()
