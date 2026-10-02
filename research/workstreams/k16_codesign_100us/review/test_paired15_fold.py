"""Independent symbolic interface audit of the final s15 folded kernel."""
import random
import re
from types import SimpleNamespace
import unittest

import test_paired15_interface as model

PERM=(0,1,2,7,3,4,5,6,8,15,10,12,11,13,14,9)


def permuted(value):
    return sum(((value>>PERM[j])&1)<<j for j in range(16))


class Paired15FoldTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _,cls.record=model.restriction.prepare((10,))
        cls.rows=[int(x,16) for x in cls.record['expansion_rows_hex']]
        cls.columns=cls.record['feedback_columns']
        cls.source=(model.KERNEL/'K16Paired15Fold.cpp').read_text()
        cls.setup=(model.KERNEL/'K16Paired15ShuffleSetup.cpp').read_text()
        cls.header=(model.KERNEL/'T64Paired15ShuffleMap.h').read_text()
        cls.expansion=model.body(cls.header,'streamHigh')
        cls.feedback=model.body(cls.source,'finishWide')
        cls.total_body=model.body(cls.source,'totals4')

    @classmethod
    def folded_totals(cls,a,b,c,d):
        prefix,returned=cls.total_body.split('return',1)
        env=dict(a=a,b=b,c=c,d=d)
        model.execute(prefix,env)
        return model.evaluate(returned.strip().removesuffix(';'),env)

    def test_actual_folded_totals(self):
        inputs=[tuple(1<<(4*i+j) for j in range(4)) for i in range(4)]
        self.assertEqual(self.folded_totals(*inputs),tuple(sum(v) for v in inputs))
        original=(model.KERNEL/'K16Paired15.cpp').read_text()
        self.assertEqual(model.body(self.source,'transpose4'),model.body(original,'transpose4'))

    def test_actual_expansion_feedback_and_dummy_slot9(self):
        state=[1<<(64+PERM[j]) for j in range(16)]
        for raw_feedback in (False,True):
            env=dict(state=SimpleNamespace(v=[tuple(state[i:i+4]) for i in range(0,16,4)]),
                raw=model.Raw([1<<p for p in range(64)]),packetBase=0,
                RawFeedback=raw_feedback,high={},emitted={},
                totals4=self.folded_totals,
                _mm512_maskz_shuffle_i32x4=lambda k,a,b,i:model.mask((0,)*4,k,model.shuffle(a,b,i),4))
            model.execute(self.expansion,env)
            emitted=sum((env['emitted'][j] for j in range(16)),())
            self.assertEqual(emitted,tuple((1<<p)^(column<<64) for p,column in enumerate(self.columns)))
            env['feedback']=SimpleNamespace(v={})
            model.execute(self.feedback,env)
            actual=sum((env['feedback'].v[j] for j in range(4)),())
            self.assertEqual(actual,tuple((self.rows+[0])[j] for j in PERM))
            self.assertEqual(actual[9],0)

    def test_noninvolutory_setup_permutation_and_GL15_embedding(self):
        text=re.search(r'permutation\[16\]=\{([^}]+)\}',self.setup).group(1)
        self.assertEqual(tuple(map(int,text.split(','))),PERM)
        self.assertEqual(sorted(PERM),list(range(16)))
        self.assertEqual(PERM[9],15)
        self.assertIn('customizePaired15(plan,seed,tables,optimized);',self.setup)
        self.assertIn('plan.reverseMatrices[epoch][permutation[r]]>>permutation[c]',self.setup)
        rng=random.Random(52015)
        for _ in range(16):
            while True:
                rows=[rng.getrandbits(15) for _ in range(15)]
                if model.rank(rows,15)==15:break
            rows.append(1<<15)
            fast=[permuted(rows[PERM[j]]) for j in range(16)]
            self.assertEqual(fast[9],1<<9)
            self.assertTrue(all(not(row>>9&1) for i,row in enumerate(fast) if i!=9))
            for j in range(16):
                self.assertEqual(model.apply(fast,permuted(1<<j)),permuted(model.apply(rows,1<<j)))

    def test_peeled_epoch_boundaries_and_update_indices(self):
        peeled=model.body(self.source,'reversePeeled')
        self.assertIn('const auto first=epochs-1',peeled)
        self.assertIn('if(!first)return',peeled)
        self.assertIn('for(std::size_t epoch=first;--epoch;)',peeled)
        self.assertIn('tables.updates.data()+4*epoch',peeled)
        self.assertIn('paired15shuffle::streamHigh<true>(state,input,0,emit,high)',peeled)
        rng=random.Random(52123)
        for epochs in (1,2,3,4,17):
            raw=[rng.getrandbits(64) for _ in range(epochs)]
            matrices=[]
            for _ in range(epochs):
                while True:
                    rows=[rng.getrandbits(15) for _ in range(15)]
                    if model.rank(rows,15)==15:break
                matrices.append(rows)
            expand=lambda s:sum(((s&c).bit_count()%2)<<p for p,c in enumerate(self.columns))
            feedback=lambda x:model.apply(self.rows,x)
            direct=[0]*epochs;state=0
            for epoch in reversed(range(epochs)):
                direct[epoch]=raw[epoch]^expand(state)
                if epoch:state=model.apply(matrices[epoch],state)^feedback(raw[epoch])
            peeled_output=[0]*epochs
            peeled_output[-1]=raw[-1]
            if epochs>1:
                state=feedback(raw[-1])
                for epoch in range(epochs-2,0,-1):
                    peeled_output[epoch]=raw[epoch]^expand(state)
                    state=model.apply(matrices[epoch],state)^feedback(raw[epoch])
                peeled_output[0]=raw[0]^expand(state)
            self.assertEqual(peeled_output,direct)


if __name__=='__main__':
    unittest.main()
