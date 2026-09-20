import itertools
import math
import unittest
import csv
import tempfile
from pathlib import Path

import numpy as np

import activation_occupation as general
import composition_boxes as boxes
import composition_occupation as sparse


class CompositionBoxesTest(unittest.TestCase):
    def test_multinomial_bound_covers_every_feasible_composition(self):
        for total in range(1,9):
            points=[c for c in itertools.product(range(total+1),repeat=3) if sum(c)==total]
            for lower,upper in (([0,0,0],[total]*3),([0,0,0],[total,1,2]),([0,1,0],[total,total,total])):
                feasible=[c for c in points if all(a<=x<=b for a,x,b in zip(lower,c,upper))]
                if not feasible:continue
                values=[math.factorial(total)//math.prod(math.factorial(x) for x in c) for c in feasible]
                mode=boxes.balanced_counts(lower,upper,total)
                self.assertEqual(math.factorial(total)//math.prod(math.factorial(int(x)) for x in mode),max(values))
                self.assertGreaterEqual(boxes.log_assignment_count_upper(lower,upper,total)+1e-12,math.log(sum(values)))

    def test_box_matrix_dominates_each_fixed_type_mixture(self):
        epoch=general.epoch_logs(8,4,{4:14,8:1},[1,0,0,0,14,0,0,0,1],.7,8)
        regions=general.region_logs(epoch[:7],8,16,6)
        counts={2:2,4:7,8:1};bands=[[2],[4],[8]];ps=np.array([.25,.5,.85])
        model=boxes.CompositionBoxes(counts,8,bands,ps,6)
        for lower,total in (([0,0,0],5),([1,1,1],6),([2,2,2],6)):
            bound=model.matrix(regions,lower,total)
            for c in itertools.product(range(total+1),repeat=3):
                if sum(c)!=total or any(x<a for x,a in zip(c,lower)):continue
                indices=np.array([[g for g,n in enumerate(c) for _ in range(n)]])
                p=sparse.probability_logs(indices,np.log(ps),np.log1p(-ps))
                actual=sparse.mixture_logs(regions,p)[0]+np.array(c)@model.log_gamma/8
                self.assertTrue(np.all(actual<=bound+2e-12))
                if sum(lower)==total:np.testing.assert_allclose(actual,bound,atol=2e-12)

    def test_selected_boxes_partition_all_counts(self):
        epoch=general.epoch_logs(8,4,{4:14,8:1},[1,0,0,0,14,0,0,0,1],.7,8)
        regions=general.region_logs(epoch[:6],8,16,5)
        model=boxes.CompositionBoxes({2:2,4:7,8:1},8,[[2],[4],[8]],[.25,.5,.85],5)
        result=boxes.search([(model,regions,.7,'test')],5,16,12,maximum_nodes=31,target_bits=1000)
        for c in itertools.product(range(6),repeat=3):
            if sum(c)!=5:continue
            self.assertEqual(sum(all(a<=x<=b for a,x,b in zip(row['lower'],c,row['upper'])) for row in result['boxes']),1)
        self.assertAlmostEqual(float(np.logaddexp.reduce([r['log_bound'] for r in result['boxes']])),result['log_union_upper'],places=10)

    def test_all_one_singleton_is_exact(self):
        epoch=general.epoch_logs(8,4,{4:14,8:1},[1,0,0,0,14,0,0,0,1],.7,4)
        regions=general.region_logs(epoch,8,16,4)
        model=boxes.CompositionBoxes({4:14,8:1},8,[[4],[8]],[.5,1.],4)
        np.testing.assert_allclose(model.matrix(regions,[0,4],4),regions[4],atol=1e-12)

    def test_producer_authenticates_each_refined_occupation(self):
        import run_threeband_grid as producer
        _,observations,_=producer.coverage.snapshot()
        row=next(r for r in observations.values() if int(r['block_bits'])==512 and int(r['step_bits'])==64
                 and int(r['state_bits'])==18 and int(r['message_exponent'])==16 and r['outer_model']!='random-ensemble-average')
        arguments=dict(maximum_occupation=3,low_probabilities=[.25],high_probability=.8677722630069483,
            log_surprisals=[-4.],maximum_nodes=3,refine_minimum=3,target_bits=1000.,witness_window=.8,setup_failure_bits=60)
        with tempfile.TemporaryDirectory(dir=producer.grid.HERE) as temporary:
            directory=Path(temporary)/'t64_s18'
            self.assertEqual(producer.run_batch(directory,[row],arguments,{}),2)
            self.assertEqual(producer.verify(directory,arguments)['row_count'],2)
            with (directory/'occupations.csv').open() as handle:rows=list(csv.DictReader(handle))
            self.assertTrue(rows[-1]['range_witness_source'])


if __name__=='__main__':unittest.main()
