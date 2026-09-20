"""All-occupancy diagnostic adapter for the IMT parameter-slice maps.

The counting and partition rules are reused; both sparse and dense transfers
are rebuilt from the separate IMT maps. No exact-refresh transfer is used.
"""
import argparse
from collections import Counter,defaultdict
from functools import lru_cache
import itertools
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
import parameter_q1 as grid
import typed_dense_boxes as typed
import balanced_occupation as balanced

search = grid.ladder.model.independent.search
g = search.g
import screen_dense as direct


@lru_cache(maxsize=40)
def prepare(t,s):
    record = grid.inner(t,s)
    a,b = record['expansion_columns'],record['feedback_columns']
    exact_b = grid.maps.enumerate_spectrum(grid.maps.generator_words(b,s),s,t)
    kernel = grid.maps.macwilliams_kernel(exact_b,s)
    b_spectrum = {w:n for w,n in enumerate(exact_b) if w and n}
    caps = g.fiber_caps(t,s,b_spectrum,kernel)
    images = [sum(((c&q).bit_count() & 1) << j for j,c in enumerate(a)) for q in b]
    low = {}
    for degree in (1,2):
        by_syndrome,by_weight = defaultdict(Counter),defaultdict(Counter)
        for support in itertools.combinations(range(t),degree):
            syndrome,image,mask = 0,0,0
            for j in support:
                syndrome ^= b[j]
                image ^= images[j]
                mask |= 1 << j
            if syndrome:
                weight = (image ^ mask).bit_count()
                by_syndrome[syndrome][weight] += 1
                by_weight[image.bit_count()][weight] += 1
        low[degree] = dict(patterns=sorted({tuple(sorted(h.items())) for h in by_syndrome.values()}),
                           by_weight=dict(by_weight),nonzero_input_count=sum(sum(h.values()) for h in by_syndrome.values()))
    assert kernel[:3] == [1,0,0]
    return record,b_spectrum,kernel,caps,low


def sparse(t,s,b,exponent,maximum):
    record,bs,kernel,caps,low = prepare(t,s)
    counts,_ = grid.outer(b)
    rows = (1 << exponent)//(b//2)
    cutoff = b*rows//10
    best = np.full(maximum-1,np.inf)
    witnesses = [None]*(maximum-1)
    for tilt in np.arange(-12.,-1.99,0.5):
        lam = math.exp(tilt)
        epochs = g.epochs(record['spectrum'],kernel,record['feedback_columns'],caps,low,lam,maximum=min(t,maximum))
        regions = g.regions(epochs,t,rows,maximum)
        for scale in (0.,0.5,1.):
            bands,roots,active,inactive = balanced.density_roots(counts,b,scale)
            envelope = balanced.Envelope(roots,active,inactive)
            current = regions.copy()
            for q in range(1,maximum+1):
                current = envelope.apply(current[:-1],current[1:])
                if q == 1:
                    continue
                bound = (math.log(math.comb(rows,q))+q*math.log(len(bands))+cutoff*lam
                         +g.terminal(current[0],b))
                if bound < best[q-2]:
                    best[q-2] = bound
                    witnesses[q-2] = dict(tilt=float(tilt),scale=scale)
        print('sparse tilt',tilt,'best union',-float(np.logaddexp.reduce(best))/math.log(2),flush=True)
    return dict(first=2,last=maximum,log_bounds=best.tolist(),witnesses=witnesses,
                margin_bits=-float(np.logaddexp.reduce(best))/math.log(2))


class Dense(typed.TypedDense):
    def __init__(self,t,s,b,exponent,tilts):
        self.record,self.bs,self.kernel,_,_ = prepare(t,s)
        counts,_ = grid.outer(b)
        rows = (1 << exponent)//(b//2)
        # Empty tilt list: the historical constructor builds only counting
        # banks. It never evaluates its exact-refresh epoch implementation.
        super().__init__(counts,b,t,s,self.record['spectrum'],self.kernel,rows,[],optimize_proposals=False)
        self.tilts = list(tilts)
        self.evaluations = 0

    def moment(self,theta,tilt):
        matrix = direct.bernoulli(self.record['spectrum'],self.bs,self.kernel,float(theta),math.exp(tilt))
        return g.terminal(matrix,self.block*self.length//self.t)

    def evaluate(self,lower,upper):
        corners = typed.vertices(lower,upper,self.length)
        if not len(corners):
            return None
        initial = typed.proposal_for(corners,self.length)
        penalty = typed.lattice_log_count(lower,upper)
        ranked = []
        def evaluate(proposal,tilt,bank):
            scale,p,costs = self.probability_banks[bank]
            moment = self.moment(proposal@p,tilt)
            values = typed.point_logs(corners,self.length,self.block,costs,proposal,moment,self.cutoff,math.exp(tilt))
            return float(max(values))+penalty
        for index,tilt in enumerate(self.tilts):
            for bank in range(len(self.probability_banks)):
                ranked.append((evaluate(initial,tilt,bank),index,bank,initial))
        best,index,bank,proposal = min(ranked,key=lambda row:row[0])
        # Direct transfer evaluations, not interpolated moments, determine
        # every accepted diagnostic bound. One cold proposal refinement.
        anchor = int(np.argmax(initial))
        free = [j for j in range(len(initial)) if j != anchor]
        def unpack(coordinates):
            logits = np.zeros(len(initial)); logits[free] = coordinates
            return np.exp(logits-logsumexp(logits))
        optimum = minimize(lambda x:evaluate(unpack(x),self.tilts[index],bank)/(self.block*self.length),
            np.log(initial/initial[anchor])[free],method='Nelder-Mead',
            options=dict(maxiter=60,xatol=1e-4,fatol=1e-10))
        updated = unpack(optimum.x)
        value = evaluate(updated,self.tilts[index],bank)
        if value < best:
            best,proposal = value,updated
        self.evaluations += 1
        if self.evaluations % 32 == 0:
            print('IMT typed boxes',self.evaluations,flush=True)
        return dict(lower=lower.tolist(),upper=upper.tolist(),corners=corners,
            own_log_bound=best,best_log_bound=best,parent=None,children=[],
            witness=dict(tilt=self.tilts[index],probability_bank=bank,proposal=proposal.tolist()))


def run(output,b,t,s,exponent,maximum,nodes):
    assert not output.exists()
    record,bs,kernel,caps,low = prepare(t,s)
    q1 = grid.screen(record,b,exponent,np.arange(-180,1,dtype=float)/10)
    sparse_result = sparse(t,s,b,exponent,maximum)
    dense = Dense(t,s,b,exponent,np.arange(-7.,0.01,0.5)).search(maximum+1,nodes,60)
    union = float(np.logaddexp.reduce([-q1['q1_margin_bits']*math.log(2),
        *sparse_result['log_bounds'],dense['log_union_upper']]))
    model = grid.ladder.model
    sources = grid.ladder.candidate.sources()
    _,path = grid.outer(b)
    sources[path.relative_to(grid.ROOT).as_posix()] = grid.maps.sha(path)
    model.base.write_new(output,dict(status='BINARY64_IMT_FULL_PARAMETER_DIAGNOSTIC',
        geometry=[b,t,s,exponent],inner=record,q1=q1,sparse=sparse_result,dense=dense,
        full_margin_bits=-union/math.log(2),full_distance_proved=False,
        source_sha256=sources))
    print('IMT full diagnostic margin',-union/math.log(2),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--b',type=int,choices=(64,128),default=128)
    p.add_argument('--t',type=int,choices=(64,128,256),default=64)
    p.add_argument('--s',type=int,default=20)
    p.add_argument('--m',type=int,default=20)
    p.add_argument('--maximum-sparse',type=int,default=8)
    p.add_argument('--nodes',type=int,default=127)
    a = p.parse_args()
    assert (a.b,a.t,a.s,a.m) in grid.geometries() and 2 <= a.maximum_sparse
    run(a.output.resolve(),a.b,a.t,a.s,a.m,a.maximum_sparse,a.nodes)
