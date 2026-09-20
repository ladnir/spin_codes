"""Continuously refine retained three-state dense bounds; never alter the code."""
import argparse
import copy
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
import evaluate_smaller_margins as base


class Epochs:
    """Prepared hypergeometric laws for the existing three-state transfer."""
    def __init__(self,a,kernel):
        t,s = 128,19
        base.fixed.general.epoch_logs(t,s,a,kernel,1.,0)
        self.weights = np.arange(t+1)
        support = sorted(a)
        self.mass = np.array([math.log(a[w])-math.log((1<<s)-1) for w in support])
        self.hyper = np.full((t+1,len(support),t+1),-np.inf)
        self.beta,self.nonkernel = np.full(t+1,-np.inf),np.full(t+1,-np.inf)
        for j in range(t+1):
            total = math.comb(t,j)
            if kernel[j]: self.beta[j] = math.log(kernel[j])-math.log(total)
            if kernel[j]<total: self.nonkernel[j] = math.log(total-kernel[j])-math.log(total)
            for i,w in enumerate(support):
                for v in range(max(0,j-t+w),min(w,j)+1):
                    self.hyper[j,i,w+j-2*v] = math.log(math.comb(w,v)*math.comb(t-w,j-v))-math.log(total)
        self.cache = {}

    def at(self,z):
        if z not in self.cache:
            lam = math.exp(z)
            moments = logsumexp(self.hyper-lam*self.weights,axis=2)
            arbitrary = np.max(moments,axis=1)
            uniform = logsumexp(moments+self.mass,axis=1)
            log_m = math.log((1<<19)-1); log_kappa = math.log1p(1/((1<<19)-2))
            out = np.full((129,3,3),-np.inf)
            out[:,0,0] = self.beta-lam*self.weights
            out[:,0,1] = self.nonkernel-lam*self.weights
            out[:,1,0] = np.minimum(self.nonkernel,arbitrary)-log_m
            out[:,1,2] = arbitrary
            out[:,2,0] = log_kappa+np.minimum(self.nonkernel,uniform)-log_m
            out[:,2,2] = log_kappa+uniform
            if len(self.cache)>=128: self.cache.clear()
            self.cache[z] = out
        return self.cache[z]


class Refiner:
    def __init__(self,row,a,kernel):
        self.row = row
        self.length = row['outer_rows']
        self.epochs = Epochs(a,kernel)
        self.banks = row['dense']['probability_banks']

    def direct(self,box,witness):
        bank = self.banks[witness['probability_bank']]
        proposal = np.array(witness['proposal'])
        theta = float(proposal@bank['probabilities'])
        matrix = base.typed.dense.epoch_mixture_logs(self.epochs.at(witness['log_surprisal']),128,np.array([theta]))
        moment = float(base.fixed.composition.terminal_logs(matrix,self.row['output_bits']//128)[0])
        corners = base.typed.vertices(box['lower'],box['upper'],self.length)
        values = base.typed.point_logs(corners,self.length,128,np.array(bank['log_density_costs']),proposal,moment,
                                      self.row['bad_weight'],math.exp(witness['log_surprisal']))
        penalty = base.typed.lattice_log_count(box['lower'],box['upper'])
        return float(max(values))+penalty,float(logsumexp(values))+penalty

    def refine(self,box):
        original = self.direct(box,box['witness'])[0]
        base.fixed.outer.require(abs(original-box['own_log_bound'])<2e-6,'input witness mismatch')
        current = copy.deepcopy(box)
        initial = np.array(box['witness']['proposal'])
        anchor = int(np.argmax(initial)); free = [j for j in range(3) if j!=anchor]
        z = box['witness']['log_surprisal']
        start = np.r_[z,np.log(initial/initial[anchor])[free]]
        def objective(coordinates):
            logits = np.zeros(3); logits[free] = coordinates[1:]
            proposal = np.exp(logits-logsumexp(logits))
            witness = dict(box['witness'],log_surprisal=float(coordinates[0]),proposal=proposal.tolist(),
                theta=float(proposal@self.banks[box['witness']['probability_bank']]['probabilities']),continuous_refinement=True)
            bound,smooth = self.direct(box,witness)
            if bound<current['own_log_bound']:
                current.update(own_log_bound=bound,witness=witness)
            return smooth/self.row['output_bits']
        result = minimize(objective,start,method='L-BFGS-B',
            bounds=[(max(-18.,z-.75),min(2.,z+.75))]+[(-28.,28.)]*2,
            options={'maxiter':80,'ftol':1e-14,'gtol':1e-10,'maxls':30})
        objective(result.x)
        return current

    def root(self,box,target,budget):
        remaining = budget
        def recurse(box,local_target):
            nonlocal remaining
            remaining -= 1
            current = self.refine(box)
            if current['own_log_bound']<=local_target or remaining<2:
                return current['own_log_bound'],[current]
            corners = base.typed.vertices(current['lower'],current['upper'],self.length)
            children = base.typed.split_box(current['lower'],current['upper'],corners)
            if not children: return current['own_log_bound'],[current]
            bounds,leaves = [],[]
            for lo,hi in children:
                child = dict(lower=lo.tolist(),upper=hi.tolist(),witness=current['witness'])
                child['own_log_bound'] = self.direct(child,child['witness'])[0]
                if remaining:
                    value,selected = recurse(child,local_target-math.log(len(children)))
                else:
                    value,selected = child['own_log_bound'],[child]
                bounds.append(value); leaves.extend(selected)
            combined = float(np.logaddexp.reduce(bounds))
            if combined<current['own_log_bound']: return combined,leaves
            return current['own_log_bound'],[current]
        return recurse(box,target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=base.fixed.HERE/'SMALLER_MARGIN_SCREEN.json')
    parser.add_argument('--output',type=Path,default=base.fixed.HERE/'SMALLER_MARGIN_REFINED.json')
    parser.add_argument('--nodes-per-root',type=int,default=7)
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    for name,digest in payload['source_sha256'].items():
        base.fixed.outer.require(base.fixed.sha(base.fixed.ROOT/name)==digest,'changed dependency: '+name)
    a,kernel,_ = base.fixed.load_inner()
    counts = {w:n for w,n in enumerate(base.smaller_outer.spectrum()) if w and n}
    # Check the prepared evaluator against the original transfer before use.
    prepared = Epochs(a,kernel)
    for z in (-10.,-5.,-1.,.5):
        np.testing.assert_allclose(prepared.at(z),base.fixed.general.epoch_logs(128,19,a,kernel,math.exp(z),128),atol=2e-12,rtol=2e-13)
    for row in payload['results']:
        if row['dense_union_margin_bits']>=60: continue
        model = Refiner(row,a,kernel)
        boxes = row['dense']['selected_boxes']
        target = -60*math.log(2)-math.log(len(boxes))
        groups = [[box] for box in boxes]
        for index in sorted(range(len(boxes)),key=lambda j:boxes[j]['own_log_bound'],reverse=True):
            if boxes[index]['own_log_bound']<=target: continue
            value,groups[index] = model.root(boxes[index],target,args.nodes_per_root)
            combined = float(np.logaddexp.reduce([b['own_log_bound'] for group in groups for b in group]))
            print(dict(root=index,root_margin=-value/math.log(2),dense_margin=-combined/math.log(2)),flush=True)
            if combined<=-60*math.log(2): break
        result = row['dense']
        result['selected_boxes'] = [b for group in groups for b in group]
        result['log_union_upper'] = float(np.logaddexp.reduce([b['own_log_bound'] for b in result['selected_boxes']]))
        dense = base.replay_dense(result,counts,a,kernel,row)
        row['dense_union_margin_bits'] = -dense/math.log(2)
        row['combined_margin_bits'] = -float(np.logaddexp(dense,-row['sparse_union_margin_bits']*math.log(2)))/math.log(2)
        row['dense_refinement'] = dict(nodes_per_root=args.nodes_per_root,input_sha256=base.fixed.sha(args.input))
    payload['source_sha256'][Path(__file__).relative_to(base.fixed.ROOT).as_posix()] = base.fixed.sha(Path(__file__))
    payload['source_sha256'][args.input.resolve().relative_to(base.fixed.ROOT).as_posix()] = base.fixed.sha(args.input)
    args.output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('FINAL',[(r['distance_target'],r['combined_margin_bits']) for r in payload['results']],flush=True)


if __name__ == '__main__':
    main()
