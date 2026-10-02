from fractions import Fraction as Q
import unittest
from unittest.mock import patch
from flint import arb,ctx
import birth_classes as birth
import scalar_cover as sc
from test_scalar_cover import endpoint


class BirthClassOutwardTests(unittest.TestCase):
    def test_class_and_tail_envelopes_share_only_fresh_same_call_fourier_parts(self):
        ctx.prec=192
        data=birth.prepare(list(range(8)),[1,2,4,3,5,7,6,1],3)
        p=Q(1,5);z=arb(3)/4;probabilities=sc.probabilities(p)
        original=birth.rank_return.outward_at_z(data,probabilities,z)
        expected=birth.base.outward_candidate(data,p,z,2,original)[0,2]
        actual_parts=birth.base.fourier_parts
        with patch.object(birth.base,'fourier_parts',wraps=actual_parts) as parts:
            matrix=birth.outward_candidate(data,probabilities,z,'classes_below_2',original)
            self.assertEqual(matrix[0,2],expected)
            parts.assert_called_once_with(data,p,z,1)
            again=birth.outward_candidate(data,probabilities,z,'classes_below_2',original)
            self.assertEqual(parts.call_count,2)
            self.assertEqual(matrix,again)

    def test_exact_class_masses_and_repeated_transitions(self):
        ctx.prec=192;z=Q(3,4)
        images=[sum(((a>>i)&1)*v for i,v in enumerate((1,6,120))) for a in range(8)]
        for columns in ([1,2,4,3,5,7,6,1],[0]*8,[1]*8):
            data=birth.prepare(images,columns,3)
            levels=list(map(int,data['birth_class_levels']))
            for p in (Q(0),Q(1,5),Q(15,16),Q(1)):
                exact=[[Q(0)]*8 for _ in range(8)]
                masses={k:[Q(0)]*len(levels) for k in (None,1,2,3)}
                for x in range(256):
                    probability=(p/15 if x&15 else 1-p)*(p/15 if x>>4 else 1-p)
                    feedback=0
                    for bit,c in enumerate(columns):
                        if x>>bit&1:feedback^=c
                    occupancy=int(bool(x&15))+int(bool(x>>4))
                    if feedback:
                        index=levels.index(images[feedback].bit_count())
                        for cutoff,values in masses.items():
                            if cutoff is None or occupancy<cutoff:values[index]+=probability*z**x.bit_count()
                    for state in range(8):
                        weight=probability*z**(images[state]^x).bit_count()
                        if not state:exact[state][feedback]+=weight
                        else:
                            exact[state][state^feedback]+=weight/4
                            for refreshed in range(1,8):exact[state][refreshed^feedback]+=weight*Q(3,28)
                for cutoff,expected in masses.items():
                    actual=birth.weighted_classes(data,p,arb(3)/4,cutoff)
                    for value,truth in zip(actual,expected):
                        self.assertGreaterEqual(endpoint(value),truth)
                        self.assertLessEqual(endpoint(value)-truth,Q(1,10**40))
                    name='all_classes' if cutoff is None else f'classes_below_{cutoff}'
                    matrix=birth.outward_candidate(data,sc.probabilities(p),arb(3)/4,name)
                    rows=[[endpoint(matrix[i,j]) for j in range(matrix.ncols())] for i in range(matrix.nrows())]
                    values=[Q(1)]*8;upper=[Q(1)]*len(rows)
                    for _ in range(3):
                        values=[sum(a*v for a,v in zip(row,values)) for row in exact]
                        upper=[sum(a*v for a,v in zip(row,upper)) for row in rows]
                        self.assertGreaterEqual(upper[0],values[0])
                        self.assertGreaterEqual(upper[1],max(values[1:]))
                        self.assertGreaterEqual(upper[2],sum(values[1:])/7)
                        for i,level in enumerate(levels,3):
                            self.assertGreaterEqual(upper[i],max(values[a] for a in range(1,8) if images[a].bit_count()==level))

    def test_selection_and_invalid_cutoff(self):
        ctx.prec=192;data=birth.prepare(list(range(8)),[1,2,4,3,5,7,6,1],3)
        probabilities=sc.probabilities(Q(1,5))
        for cutoff in (0,4,1.5):
            with self.assertRaises(ValueError):birth.weighted_classes(data,Q(1,5),arb(3)/4,cutoff)
        for z in (arb(0),arb(2)):
            with self.assertRaises(ValueError):birth.outward_at_z(data,probabilities,z)
        with self.assertRaises(ValueError):birth.outward_candidate(data,probabilities,arb(3)/4,'unknown')
        matrix=birth.outward(data,probabilities,Q(1,8))
        self.assertTrue(all(matrix[i,j]>=0 for i in range(matrix.nrows()) for j in range(matrix.ncols())))


if __name__=='__main__':unittest.main()
