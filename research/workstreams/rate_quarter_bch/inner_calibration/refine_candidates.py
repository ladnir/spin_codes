"""Bounded continuous refinement of candidate-specific dense witnesses."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
import check_candidates as check
import refine_smaller_dense as original


class Epochs:
    def __init__(self,record):
        self.t,self.s = t,s = record['step_bits'],record['state_bits']
        a = {w:n for w,n in enumerate(record['a_counts']) if w and n}
        kernel = record['kernel_counts']
        check.fixed.general.epoch_logs(t,s,a,kernel,1.,0)
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
            denominator = (1<<self.s)-1
            log_m,log_kappa = math.log(denominator),math.log1p(1/(denominator-1))
            out = np.full((self.t+1,3,3),-np.inf)
            out[:,0,0] = self.beta-lam*self.weights
            out[:,0,1] = self.nonkernel-lam*self.weights
            out[:,1,0] = np.minimum(self.nonkernel,arbitrary)-log_m
            out[:,1,2] = arbitrary
            out[:,2,0] = log_kappa+np.minimum(self.nonkernel,uniform)-log_m
            out[:,2,2] = log_kappa+uniform
            if len(self.cache)>=128: self.cache.clear()
            self.cache[z] = out
        return self.cache[z]


class Refiner(original.Refiner):
    def __init__(self,row,record,epochs):
        self.row = dict(row,outer_rows=check.L,output_bits=check.N)
        delta = check.Fraction(row['distance_target'])
        self.row['bad_weight'] = check.N*delta.numerator//delta.denominator
        self.length,self.t = check.L,record['step_bits']
        self.epochs = epochs
        self.banks = row['dense']['probability_banks']

    def direct(self,box,witness):
        bank = self.banks[witness['probability_bank']]
        proposal = np.array(witness['proposal'])
        theta = float(proposal@bank['probabilities'])
        matrix = check.typed.dense.epoch_mixture_logs(self.epochs.at(witness['log_surprisal']),self.t,np.array([theta]))
        moment = float(check.fixed.composition.terminal_logs(matrix,check.N//self.t)[0])
        corners = check.typed.vertices(box['lower'],box['upper'],self.length)
        values = check.typed.point_logs(corners,self.length,check.B,np.array(bank['log_density_costs']),proposal,moment,
                                        self.row['bad_weight'],math.exp(witness['log_surprisal']))
        penalty = check.typed.lattice_log_count(box['lower'],box['upper'])
        return float(max(values))+penalty,float(logsumexp(values))+penalty


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag',default='t128_s18_nested')
    parser.add_argument('--nodes-per-root',type=int,default=3)
    args = parser.parse_args()
    folder = check.calibration.HERE
    source = folder/f'{args.tag}_dense.json'
    payload = json.loads(source.read_text())
    for name,digest in payload['source_sha256'].items():
        assert check.fixed.sha(check.fixed.ROOT/name)==digest,name
    map_path = folder/'maps'/f'{args.tag}.json'
    record = json.loads(map_path.read_text())
    epochs = Epochs(record)
    a = {w:n for w,n in enumerate(record['a_counts']) if w and n}
    for z in (-10.,-5.,-1.,.5):
        np.testing.assert_allclose(epochs.at(z),check.fixed.general.epoch_logs(record['step_bits'],record['state_bits'],a,
                                   record['kernel_counts'],math.exp(z),record['step_bits']),atol=2e-11,rtol=2e-13)
    for row in payload['results']:
        if row['dense_union_margin_bits']>=60: continue
        model = Refiner(row,record,epochs)
        boxes = row['dense']['selected_boxes']
        groups = [[box] for box in boxes]
        target = -60*math.log(2)-math.log(len(boxes))
        for i in sorted(range(len(boxes)),key=lambda j:-boxes[j]['own_log_bound']):
            if boxes[i]['own_log_bound']<=target: continue
            _,groups[i] = model.root(boxes[i],target,args.nodes_per_root)
            total = float(np.logaddexp.reduce([b['own_log_bound'] for g in groups for b in g]))
            print(args.tag,row['distance_target'],'root',i,'dense margin',-total/math.log(2),flush=True)
            if total<=-60*math.log(2): break
        dense = row['dense']
        dense['selected_boxes'] = [b for g in groups for b in g]
        check.check_coverage(dense['selected_boxes'],check.L,dense['occupation_min'])
        for box in dense['selected_boxes']:
            assert abs(model.direct(box,box['witness'])[0]-box['own_log_bound'])<2e-6
        dense['log_union_upper'] = float(np.logaddexp.reduce([b['own_log_bound'] for b in dense['selected_boxes']]))
        row['dense_union_margin_bits'] = -dense['log_union_upper']/math.log(2)
    for path in (Path(__file__),Path(original.__file__),source,map_path):
        payload['source_sha256'][path.relative_to(check.fixed.ROOT).as_posix()] = check.fixed.sha(path)
    payload['refinement_nodes_per_root'] = args.nodes_per_root
    (folder/f'{args.tag}_dense_refined.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('FINAL',[(r['distance_target'],r['dense_union_margin_bits']) for r in payload['results']],flush=True)


if __name__=='__main__': main()
