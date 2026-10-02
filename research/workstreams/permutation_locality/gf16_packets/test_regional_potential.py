from types import SimpleNamespace
from fractions import Fraction as Q
from math import comb
import unittest
import numpy as np
import regional_count_probe as linear
import regional_potential as nonlinear


class RegionalPotentialTests(unittest.TestCase):
    def test_nonlinear_conditional_region_bounds_exact_small_code(self):
        from flint import ctx
        import birth_classes
        import occupancy_birth_classes
        import fiber_density
        import mass_density_potential as potential
        ctx.prec=192;z=Q(3,4)
        images=[0x93*(a&1)^0x65*((a>>1)&1)^0x3c*((a>>2)&1) for a in range(8)]
        columns=[1,2,4,3,6,5,7,1]
        data=fiber_density.attach(birth_classes.prepare(images,columns,3,2))
        source=occupancy_birth_classes.outward_at_z(data,potential.aq(z))
        bound=potential.Bound(data,source,potential.aq(z));levels=list(data['birth_class_levels'])
        exact=[[[Q(0)]*8 for _ in range(8)] for _ in range(3)]
        for x in range(256):
            j=int(bool(x&15))+int(bool(x>>4));chance=Q(1,comb(2,j)*15**j);feedback=0
            for bit,c in enumerate(columns):
                if x>>bit&1:feedback^=c
            for a in range(8):
                weight=chance*z**(images[a]^x).bit_count()
                if not a:exact[j][a][feedback]+=weight
                else:
                    exact[j][a][a^feedback]+=weight/4
                    for b in range(1,8):exact[j][a][b^feedback]+=weight*Q(3,28)
        terminal=[Q(1+a*a,7) for a in range(8)]
        def norms(v):return [v[0],max(v[1:]),sum(v[1:])/7,
            *(max(v[a] for a in range(1,8) if images[a].bit_count()==l) for l in levels)]
        values,logs=nonlinear.conditional(bound,np.array(list(map(float,norms(terminal)))),2)
        for k in range(5):
            expected=[Q(0)]*8
            for j in range(max(0,k-2),min(2,k)+1):
                chance=Q(comb(2,j)*comb(2,k-j),comb(4,k))
                suffix=[sum(c*v for c,v in zip(row,terminal)) for row in exact[k-j]]
                result=[sum(c*v for c,v in zip(row,suffix)) for row in exact[j]]
                expected=[a+chance*b for a,b in zip(expected,result)]
            self.assertTrue((values[k]*np.exp(logs[k])+1e-12>=np.array(list(map(float,norms(expected))))).all())

    def test_linear_special_case_matches_all_conditional_products(self):
        matrices=np.array([[[1.,0],[0,.8]],[[.5,.2],[.3,.5]],[[.2,.3],[.1,.4]]])
        v=np.array([1.,.4])
        for scale in (1.,1e-80):
            local=matrices*scale
            bound=SimpleNamespace(W=2,n=2,floating=lambda j,values:values@local[j].T)
            for epochs in (1,2,3,4):
                region,logs=nonlinear.conditional(bound,v,epochs)
                expected,expected_logs=linear.scaled_placement(local,epochs)
                truth=expected@v;peaks=truth.max(axis=1)
                np.testing.assert_allclose(logs,expected_logs+np.log(peaks),rtol=1e-12,atol=1e-10)
                np.testing.assert_allclose(region,truth/peaks[:,None],rtol=1e-12,atol=1e-12)
                weights=np.linspace(-5,0,len(region))
                image,shift=nonlinear.weighted_image(region,logs,weights)
                matrix,matrix_shift=linear.weighted_region(expected,expected_logs+weights)
                reference=matrix@v;peak=reference.max()
                self.assertAlmostEqual(shift,matrix_shift+np.log(peak),places=9)
                np.testing.assert_allclose(image,reference/peak,rtol=1e-12,atol=1e-12)

    def test_positive_potential_gives_an_upper_bound_on_finite_linear_power(self):
        local=np.array([[[1.,0],[0,.8]],[[.3,.4],[.2,.5]]])
        bound=SimpleNamespace(W=1,n=2,floating=lambda j,values:values@local[j].T)
        epochs=3;regions=4;weights=np.log([.1,.2,.3,.4]);offset=.5;count=.7
        score,parts=nonlinear.propose(bound,[(weights,count)],offset,[1.,1.],regions,epochs,4)
        matrices,logs=linear.scaled_placement(local,epochs)
        matrix,shift=linear.weighted_region(matrices,logs+weights)
        exact=(count+offset+regions*shift+np.log(np.linalg.matrix_power(matrix,regions)[0].sum()))/np.log(2)
        self.assertGreaterEqual(score+1e-12,exact)
        self.assertTrue(parts[0]['potential'])


if __name__=='__main__':unittest.main()
