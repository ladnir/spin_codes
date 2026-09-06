"""Dense counting with fixed row-type counts and convex simplex-box bounds.

Zero rows are category 0. Each other category is a deterministic outer
weight-band counting measure dominated by Gamma_g Bernoulli(p_g). A
single proposal and tilt are held fixed over every vertex of each box.
All numerical results are nearest binary64 diagnostics.
"""
import heapq
import itertools
import math

import numpy as np
from scipy.special import gammaln,logsumexp
from scipy.optimize import minimize_scalar,minimize

import activation_occupation as general
import dense_occupation_ranges as dense
from composition_occupation import terminal_logs
from activation_q1 import region_logs as marked_power


def weight_categories(counts,block):
    bands=[[],[],[],[]]
    for w in sorted(counts):
        if w==block:bands[3].append(w)
        elif 4*w<block:bands[0].append(w)
        elif 4*w<=3*block:bands[1].append(w)
        else:bands[2].append(w)
    return [band for band in bands if band]


def category_parameters(counts,block,bands=None):
    bands=weight_categories(counts,block) if bands is None else bands
    if sorted(w for band in bands for w in band)!=sorted(counts):raise ValueError('bands must partition the spectrum')
    probabilities=[0.];costs=[0.]
    for band in bands:
        if band==[block]:
            probabilities.append(1.);costs.append(math.log(counts[block]));continue
        weights=np.array(band)
        offsets=np.array([math.log(counts[w])-math.log(math.comb(block,w)) for w in band])
        def objective(eta):return float(block*np.logaddexp(0.,eta)+max(offsets-weights*eta))
        optimum=minimize_scalar(objective,bounds=(-16.,16.),method='bounded',options={'xatol':1e-10})
        eta=float(optimum.x)
        probabilities.append(math.exp(-float(np.logaddexp(0.,-eta))))
        costs.append(objective(eta))
    return bands,np.array(probabilities),np.array(costs)


def vertices(lower,upper,total):
    lower=np.asarray(lower,dtype=np.int64);upper=np.asarray(upper,dtype=np.int64)
    if lower.shape!=upper.shape or lower.ndim!=1 or len(lower)<2:
        raise ValueError('invalid count box')
    if np.any(lower<0) or np.any(lower>upper) or lower.sum()>total or upper.sum()<total:
        return np.empty((0,len(lower)),dtype=np.int64)
    result=[];categories=len(lower)
    for free in range(categories):
        fixed=[g for g in range(categories) if g!=free]
        for choices in itertools.product((0,1),repeat=categories-1):
            row=lower.copy()
            for g,choice in zip(fixed,choices):row[g]=upper[g] if choice else lower[g]
            row[free]=total-int(row[fixed].sum())
            if lower[free]<=row[free]<=upper[free]:result.append(row)
    return np.unique(np.array(result,dtype=np.int64),axis=0) if result else np.empty((0,categories),dtype=np.int64)


def lattice_log_count(lower,upper):
    widths=np.asarray(upper,dtype=np.int64)-np.asarray(lower,dtype=np.int64)+1
    if np.any(widths<1):raise ValueError('empty count box')
    logs=np.log(widths.astype(float))
    # Any choice of all but one coordinate fixes the last coordinate.
    return float(logs.sum()-logs.max())


def point_logs(count_vectors,total,block,log_gammas,proposal,moment_log,cutoff,lam):
    counts=np.asarray(count_vectors,dtype=np.int64)
    proposal=np.asarray(proposal,dtype=float)
    if (counts.ndim!=2 or counts.shape[1]!=len(proposal) or np.any(counts<0)
            or not np.all(counts.sum(axis=1)==total) or np.any(proposal<=0)
            or abs(float(proposal.sum())-1.)>1e-12):
        raise ValueError('invalid type counts or positive proposal')
    multinomial=gammaln(total+1)-gammaln(counts+1).sum(axis=1)
    return (cutoff*lam+moment_log+counts@(np.asarray(log_gammas)-block*np.log(proposal))
            -(block-1)*multinomial)


def proposal_for(corners,total):
    proposal=np.maximum(corners.mean(axis=0)/total,1e-12)
    return proposal/proposal.sum()


def split_box(lower,upper,corners):
    # Tighten to the exact coordinate extrema of the simplex intersection.
    lower=corners.min(axis=0);upper=corners.max(axis=0)
    coordinate=int(np.argmax(upper-lower))
    if upper[coordinate]==lower[coordinate]:return []
    midpoint=(int(lower[coordinate])+int(upper[coordinate]))//2
    left_upper=upper.copy();left_upper[coordinate]=midpoint
    right_lower=lower.copy();right_lower[coordinate]=midpoint+1
    return [(lower.copy(),left_upper),(right_lower,upper.copy())]


def moment_and_density(epoch,t,total_bits,theta):
    """Moment and tilted input density, from a positive marked matrix power."""
    if not 0<theta<1 or total_bits%t:raise ValueError('invalid marked moment geometry')
    matrix=np.full((1,3,3),-np.inf);marked=np.full_like(matrix,-np.inf)
    for j in range(t+1):
        term=(math.log(math.comb(t,j))+j*math.log(theta)+(t-j)*math.log1p(-theta))+epoch[j]
        np.logaddexp(matrix[0],term,out=matrix[0])
        if j:np.logaddexp(marked[0],term+math.log(j/t),out=marked[0])
    ordinary,average_mark=marked_power(matrix,marked,total_bits//t)
    value=float(np.logaddexp.reduce(ordinary[0,0]))
    mark=float(np.logaddexp.reduce(average_mark[0,0]))
    return value,min(1.,max(0.,math.exp(mark-value)))


def proposal_objective(logits,corners,total,block,ps,log_gammas,epoch,t):
    """Smooth vertex upper and analytic gradient, normalized by output size."""
    log_proposal=np.asarray(logits)-logsumexp(logits)
    proposal=np.exp(log_proposal);theta=float(proposal@ps)
    value,density=moment_and_density(epoch,t,block*total,theta)
    multinomial=gammaln(total+1)-gammaln(corners+1).sum(axis=1)
    terms=corners@(log_gammas-block*log_proposal)-(block-1)*multinomial
    aggregate=float(logsumexp(terms))
    mean_counts=np.exp(terms-aggregate)@corners
    derivative=(density-theta)/(theta*(1-theta))
    gradient=proposal*(1+derivative*(ps-theta))-mean_counts/total
    return (value+aggregate)/(block*total),gradient


def optimized_proposal(corners,total,block,ps,costs,epoch,t,initial):
    anchor=int(np.argmax(initial));free=np.array([j for j in range(len(initial)) if j!=anchor])
    initial_logits=np.log(initial/initial[anchor])
    def objective(coordinates):
        logits=np.zeros(len(initial));logits[free]=coordinates
        value,gradient=proposal_objective(logits,corners,total,block,ps,costs,epoch,t)
        return value,gradient[free]
    solution=minimize(objective,initial_logits[free],jac=True,method='L-BFGS-B',
                      bounds=[(-28.,28.)]*len(free),options={'maxiter':40,'ftol':1e-11,'gtol':1e-7})
    logits=np.zeros(len(initial));logits[free]=solution.x
    return np.exp(logits-logsumexp(logits))


class TypedDense:
    def __init__(self,counts,block,t,s,a_counts,kernel,length,tilts,probability_scales=(.5,.75,1.),optimize_proposals=True,bands=None):
        if length%t or length<t:raise ValueError('non-native typed geometry')
        self.block=block;self.t=t;self.length=length;self.cutoff=block*length//10
        self.bands,self.ps,self.log_gammas=category_parameters(counts,block,bands)
        self.probability_banks=[]
        for scale in probability_scales:
            probabilities=self.ps.copy();costs=self.log_gammas.copy()
            for g,band in enumerate(self.bands,start=1):
                if band==[block]:continue
                eta=scale*(math.log(self.ps[g])-math.log1p(-self.ps[g]))
                log_p=-float(np.logaddexp(0.,-eta));log_not_p=-float(np.logaddexp(0.,eta))
                probabilities[g]=math.exp(log_p)
                costs[g]=max(math.log(counts[w])-math.log(math.comb(block,w))-w*log_p-(block-w)*log_not_p for w in band)
            self.probability_banks.append((float(scale),probabilities,costs))
        self.tilts=list(tilts)
        self.optimize_proposals=optimize_proposals
        self.epochs=[general.epoch_logs(t,s,a_counts,kernel,math.exp(z),t) for z in tilts]

    def evaluate(self,lower,upper):
        corners=vertices(lower,upper,self.length)
        if not len(corners):return None
        proposal=proposal_for(corners,self.length)
        best=math.inf;witness=None;candidates=[]
        penalty=lattice_log_count(lower,upper)
        for tilt,epoch in zip(self.tilts,self.epochs):
            thetas=np.array([proposal@p for scale,p,costs in self.probability_banks])
            matrices=dense.epoch_mixture_logs(epoch,self.t,thetas)
            moments=terminal_logs(matrices,self.block*self.length//self.t)
            for bank,((scale,p,costs),moment) in enumerate(zip(self.probability_banks,moments)):
                values=point_logs(corners,self.length,self.block,costs,proposal,float(moment),self.cutoff,math.exp(tilt))
                bound=float(max(values))+penalty
                candidates.append((bound,tilt,bank))
                if bound<best:
                    best=bound;witness=dict(log_surprisal=tilt,proposal=proposal.tolist(),theta=float(thetas[bank]),
                                           probability_scale=scale,probability_bank=bank,
                                           maximum_vertex=corners[int(np.argmax(values))].tolist())
        if self.optimize_proposals:
            # One or two cold optimizer calls per box; all hot moment arithmetic
            # has an analytic derivative and uses fixed three-state contractions.
            for _,tilt,bank in sorted(candidates)[:2]:
                scale,p,costs=self.probability_banks[bank];epoch=self.epochs[self.tilts.index(tilt)]
                improved=optimized_proposal(corners,self.length,self.block,p,costs,epoch,self.t,proposal)
                theta=float(improved@p)
                moment,_=moment_and_density(epoch,self.t,self.block*self.length,theta)
                values=point_logs(corners,self.length,self.block,costs,improved,moment,self.cutoff,math.exp(tilt))
                bound=float(max(values))+penalty
                if bound<best:
                    best=bound;witness=dict(log_surprisal=tilt,proposal=improved.tolist(),theta=theta,
                        probability_scale=scale,probability_bank=bank,proposal_optimized=True,
                        maximum_vertex=corners[int(np.argmax(values))].tolist())
        return dict(lower=lower.tolist(),upper=upper.tolist(),corners=corners,
                    own_log_bound=best,best_log_bound=best,witness=witness,parent=None,children=[])

    def search(self,minimum=17,maximum_nodes=1023,target_bits=80):
        if not 1<=minimum<=self.length or maximum_nodes<1:raise ValueError('invalid typed search range')
        lower=np.zeros(len(self.ps),dtype=np.int64);upper=np.full(len(self.ps),self.length,dtype=np.int64)
        upper[0]=self.length-minimum
        root=self.evaluate(lower,upper)
        if root is None:raise ValueError('empty typed domain')
        nodes=[root];heap=[(-root['own_log_bound'],0)]
        while heap and len(nodes)+2<=maximum_nodes and nodes[0]['best_log_bound']>-target_bits*math.log(2):
            _,index=heapq.heappop(heap);node=nodes[index]
            children=split_box(node['lower'],node['upper'],node['corners'])
            for lo,hi in children:
                child=self.evaluate(lo,hi)
                if child is None:continue
                child['parent']=index;child_index=len(nodes);nodes.append(child)
                node['children'].append(child_index)
                heapq.heappush(heap,(-child['own_log_bound'],child_index))
            if not node['children']:continue
            changed=index
            while changed is not None:
                ancestor=nodes[changed]
                if ancestor['children']:
                    combined=float(np.logaddexp.reduce([nodes[j]['best_log_bound'] for j in ancestor['children']]))
                    ancestor['best_log_bound']=min(ancestor['own_log_bound'],combined)
                changed=ancestor['parent']
        selected=[]
        def extract(index):
            node=nodes[index]
            if not node['children'] or node['own_log_bound']<=node['best_log_bound']:
                selected.append(index)
            else:
                for child in node['children']:extract(child)
        extract(0)
        return dict(log_union_upper=nodes[0]['best_log_bound'],nodes_evaluated=len(nodes),selected_boxes=[
            {k:v for k,v in nodes[i].items() if k in ('lower','upper','own_log_bound','witness')} for i in selected],
            bands=self.bands,probabilities=self.ps.tolist(),log_density_costs=self.log_gammas.tolist(),
            probability_banks=[dict(scale=scale,probabilities=p.tolist(),log_density_costs=c.tolist()) for scale,p,c in self.probability_banks],
            occupation_min=minimum,occupation_max=self.length,arithmetic='nearest binary64 diagnostic')
