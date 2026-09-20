import itertools
import math
import unittest
import csv
import tempfile
from pathlib import Path

import numpy as np

import activation_occupation as general
import batched_typed_ranges as batched
import dense_occupation_ranges as dense
import typed_dense_boxes as typed
from composition_occupation import terminal_logs


class BatchedTypedTest(unittest.TestCase):
    def test_partition_covers_every_feasible_dense_type_once(self):
        for depth in (0,1,3,7):
            parts=batched.partition(8,4,3,depth)
            for c in itertools.product(range(9),repeat=4):
                if sum(c)!=8 or c[0]>5:continue
                self.assertEqual(sum(all(a<=x<=b for a,x,b in zip(lo,c,hi)) for lo,hi in parts),1)

    def test_batched_bounds_match_individual_vertex_evaluations(self):
        counts={1:1,4:6,7:1,8:1}
        engine=batched.TypedPartition(counts,8,16,minimum=3,depth=3,scales=(.5,1.),eligibility_shifts=(0.,.5))
        epoch=general.epoch_logs(8,4,{4:14,8:1},[1,0,0,0,14,0,0,0,1],.7,8)
        engine.observe(epoch,8,.7,math.log(.7))
        for j,(lo,hi) in enumerate(engine.boxes):
            corners=typed.vertices(lo,hi,16);values=[]
            for proposal,p,cost in zip(engine.proposals[:,j],engine.probabilities,engine.costs):
                theta=float(proposal@p)
                moment=terminal_logs(dense.epoch_mixture_logs(epoch,8,np.array([theta])),16)[0]
                bound=max(typed.point_logs(corners,16,8,cost,proposal,moment,12,.7))+typed.lattice_log_count(lo,hi)
                values.append(bound)
            self.assertAlmostEqual(engine.best[j],min(values),places=10)
        result=engine.result()
        self.assertAlmostEqual(result['log_union_upper'],float(np.logaddexp.reduce(engine.best)),places=10)

    def test_producer_records_full_range_and_shared_random_event(self):
        import run_typed_range_grid as producer
        _,observations,_=producer.coverage.snapshot()
        rows=[r for r in observations.values() if int(r['block_bits'])==512 and int(r['step_bits'])==64
              and int(r['state_bits'])==12 and int(r['message_exponent'])==16]
        self.assertEqual(len(rows),2)
        arguments=dict(minimum_occupation=65,partition_depth=0,probability_scales=[.75],
                       eligibility_shifts=[0.],log_surprisals=[-.5],setup_failure_bits=60)
        with tempfile.TemporaryDirectory(dir=producer.grid.HERE) as temporary:
            directory=Path(temporary)/'t64_s12'
            self.assertEqual(producer.run_batch(directory,rows,arguments,{}),2)
            receipt=producer.common.verify(directory,arguments)
            self.assertEqual(receipt['row_count'],2)
            with (directory/'occupations.csv').open() as handle:output=list(csv.DictReader(handle))
            self.assertTrue(all((r['occupation_min'],r['occupation_max'])==('65','256') for r in output))
            self.assertEqual(sum(bool(r['setup_event_id']) for r in output),1)
            self.assertEqual(producer.run_batch(directory,rows,arguments,{}),2)


if __name__=='__main__':unittest.main()
