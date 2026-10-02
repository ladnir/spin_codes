from copy import deepcopy
from fractions import Fraction as Q
from itertools import product
from math import comb,factorial
import unittest
from unittest.mock import patch
import numpy as np
from flint import arb,ctx
import scalar_cover as sc
import variance_partition as vp
from test_scalar_cover import endpoint


class VariancePartitionTests(unittest.TestCase):
    def test_full_interval_outer_bound_against_exact_compositions(self):
        ctx.prec=192;n=4;regions=2
        probabilities=[Q(0),Q(1,4),Q(3,4),Q(1)]
        components=[(str(i),Q(i+1),sc.probabilities(p),int(i!=0)) for i,p in enumerate(probabilities)]
        with patch.object(sc,'G',n),patch.object(sc,'REGIONS',regions):
            model=sc.Model(components,{},0,q_min=2,tilt=Q(1),variance_shuffle=True,variance_bins=4)
            cell=(Q(1,4),Q(3,4));cs,_,logs=model.family(model.tilt)
            parts=[]
            for j in range(4):
                interval=(Q(j,16),Q(j+1,16))
                _,dual=vp.outer_witness(logs,model.features,model.active,cell,Q(1,2),interval)
                parts.append(dict(interval=list(map(str,interval)),dual=list(map(str,dual))))
            witness=dict(tilt='1',variance_dual=list(map(str,model.shuffle_dual(cell))),variance_partition=parts)
            bound=endpoint(vp.outward_outer(model,cell,witness));expected=Q(0)
            for counts in product(range(n+1),repeat=4):
                if sum(counts)!=n or sum(a*c for a,c in zip(model.active,counts))<2:continue
                mean=sum(f*c for f,c in zip(model.features,counts))
                if not n*cell[0]<=mean<=n*cell[1]:continue
                masses=[Q(1)];coefficient=Q(factorial(n))
                for c,weight,p in zip(counts,cs,model.features):
                    coefficient*=weight**c/factorial(c)
                    for _ in range(c):
                        updated=[Q(0)]*(len(masses)+1)
                        for k,m in enumerate(masses):updated[k]+=m*(1-p);updated[k+1]+=m*p
                        masses=updated
                p=mean/n
                ratio=max(m/(comb(n,k)*p**k*(1-p)**(n-k)) for k,m in enumerate(masses))
                expected+=coefficient*ratio**regions
            self.assertGreaterEqual(bound,expected)
            for changed in (parts[:-1],parts[1:],list(reversed(parts)),parts+parts):
                with self.assertRaises(ValueError):vp.outward_outer(model,cell,dict(witness,variance_partition=changed))
            changed=deepcopy(parts);changed[0]['dual'][1]='-1'
            with self.assertRaises(ValueError):vp.outward_outer(model,cell,dict(witness,variance_partition=changed))
            changed=deepcopy(parts);changed[1]['interval'][0]='0'
            with self.assertRaises(ValueError):vp.outward_outer(model,cell,dict(witness,variance_partition=changed))

    def test_real_size_proposal_and_outward_agree(self):
        ctx.prec=192
        components=sc.group_components([(Q(1),p) for p in (Q(0),Q(1),Q(1,2),Q(2,5),Q(3,5))])
        data=sc.kernel.prepare(list(range(8)),[1,2,4,3]*32,3)
        model=sc.Model(components,data,10485,q_min=129,tilt=Q(1,8),variance_shuffle=True,variance_bins=4)
        cell=(Q(1,25),Q(401,10000))
        base=model.propose_with(cell,model.tilt)
        score,witness=vp.propose(model,cell,base)
        bound=model.outward(cell,witness)
        self.assertAlmostEqual(score,float(bound.log()/arb(2).log()),places=5)
        ctx.prec=256
        precise=model.outward(cell,witness)
        self.assertAlmostEqual(score,float(precise.log()/arb(2).log()),places=5)
        for changed in (dict(witness,tilt='1',weights_dual=[['0','0'],['0','0']]),):
            with self.assertRaises(ValueError):model.outward(cell,changed)
        model.variance_bins=0
        with self.assertRaises(ValueError):model.outward(cell,witness)

    def test_invalid_modes(self):
        components=[('0',Q(1),sc.probabilities(Q(0)),0),('1',Q(1),sc.probabilities(Q(1,2)),1)]
        for value in (-1,65,True,1.5):
            with self.assertRaises(ValueError):sc.Model(components,{},0,variance_shuffle=True,variance_bins=value)
        with self.assertRaises(ValueError):sc.Model(components,{},0,variance_bins=2)


if __name__=='__main__':unittest.main()
