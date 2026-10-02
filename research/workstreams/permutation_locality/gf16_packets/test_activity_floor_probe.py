from fractions import Fraction as Q
from itertools import product
from math import factorial,exp,log
import unittest
from unittest.mock import patch
from pathlib import Path
import sys
import numpy as np
import activity_floor_probe as floor_probe
import scalar_cover as sc
import variance_partition as variance

# The proof modules add their dense-reference directories to sys.path.
# Keep this suite first so unittest discovery does not import a sibling
# test_assemble.py with the same module name.
sys.path.insert(0,str(Path(__file__).resolve().parent))


class ActivityFloorTests(unittest.TestCase):
    def test_reoptimized_counts_cover_marked_composition_sum(self):
        n=4;ps=list(map(Q,('0','1/4','1/2','3/4')))
        components=[(str(i),Q(i+1),sc.probabilities(p),int(i!=0)) for i,p in enumerate(ps)]
        with patch.object(sc,'G',n),patch.object(sc,'REGIONS',2):
            model=sc.Model(components,{},0,q_min=2,tilt=Q(1))
            cell=(Q(1,4),Q(3,4));interval=(Q(0),Q(1,4))
            parts=[(interval,(Q(0),Q(0),Q(0)))];cs,_,logs=model.family(model.tilt)
            for gamma in map(Q,('0','1/8','1/2')):
                actual=floor_probe.marked_counts(model,cell,parts,gamma,floor=2)[0]
                total=0.
                for counts in product(range(n+1),repeat=len(ps)):
                    q=sum(a*c for a,c in zip(model.active,counts))
                    if sum(counts)!=n or q<model.q_min:continue
                    mean=sum(f*c for f,c in zip(model.features,counts))/n
                    if not cell[0]<=mean<=cell[1]:continue
                    weight=float(factorial(n))
                    for c,coefficient in zip(counts,cs):weight*=float(coefficient)**c/factorial(c)
                    total+=weight*exp(-2*float(gamma)*q)
                self.assertGreaterEqual(actual+1e-10,log(total))
                if gamma==0:
                    expected=n*variance.outer_witness(logs,np.array(list(map(float,model.features))),
                        np.array(model.active),cell,model.q_min/n,interval)[0]
                    self.assertAlmostEqual(actual,expected,places=10)

    def test_negative_mark_or_floor_rejected(self):
        for gamma,floor in ((Q(-1),2),(Q(1),-1)):
            with self.assertRaises(ValueError):floor_probe.marked_counts(None,None,None,gamma,floor)


if __name__=='__main__':unittest.main()
