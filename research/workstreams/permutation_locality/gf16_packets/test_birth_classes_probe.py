from fractions import Fraction as Q
import unittest
import numpy as np
import birth_classes_probe as probe
import scalar_cover as sc


class BirthClassesTests(unittest.TestCase):
    def test_integer_walsh(self):
        values=np.array([1,2,0,-1,3,1,2,0])
        expected=[sum(int(v)*(-1)**((i&j).bit_count()%2) for j,v in enumerate(values)) for i in range(8)]
        np.testing.assert_array_equal(probe.walsh(values),expected)
        np.testing.assert_array_equal(probe.walsh(probe.walsh(values)),8*values)
        with self.assertRaises(ValueError):probe.walsh([1,2,3])

    def test_class_masses_and_multistep_against_all_states(self):
        images=[sum(((s>>i)&1)*v for i,v in enumerate((1,6,120))) for s in range(8)]
        for columns in ([1,2,4,3,5,7,6,1],[0]*8,[1]*8):
            for updates in (1,2,3):
                data=probe.prepare(images,columns,3,updates)
                for p in (Q(0),Q(1,5),Q(15,16),Q(1)):
                    z=Q(3,4);masses={int(level):Q(0) for level in data['birth_class_levels']}
                    exact=[[Q(0) for _ in range(8)] for _ in range(8)]
                    for x in range(256):
                        probability=(p/15 if x&15 else 1-p)*(p/15 if x>>4 else 1-p)
                        feedback=0
                        for bit,c in enumerate(columns):
                            if x>>bit&1:feedback^=c
                        if feedback:masses[images[feedback].bit_count()]+=probability*z**x.bit_count()
                        for state in range(8):
                            weight=probability*z**(images[state]^x).bit_count()
                            if not state:exact[state][feedback]+=weight
                            else:
                                exact[state][state^feedback]+=weight/Q(2**updates)
                                for refreshed in range(1,8):
                                    exact[state][refreshed^feedback]+=weight*(1-Q(1,2**updates))/7
                    np.testing.assert_allclose(probe.class_masses(data,float(p),float(z)),list(map(float,masses.values())),atol=1e-14)
                    exact=np.array(exact,dtype=float)
                    options=[matrix for refine in (False,True)
                        for matrix in probe.candidates(data,sc.probabilities(p),-np.log(float(z)),class_returns=refine).values()]
                    for matrix in options:
                        values=np.ones(8);upper=np.ones(len(matrix))
                        for step in range(8):
                            values=exact@values;upper=matrix@upper
                            self.assertGreaterEqual(upper[0]+1e-13,values[0])
                            self.assertGreaterEqual(upper[1]+1e-13,max(values[1:]))
                            self.assertGreaterEqual(upper[2]+1e-13,np.mean(values[1:]))
                            for i,level in enumerate(data['birth_class_levels'],3):
                                if i<len(matrix):self.assertGreaterEqual(upper[i]+1e-13,max(values[s] for s in range(1,8) if images[s].bit_count()==level))


if __name__=='__main__':unittest.main()
