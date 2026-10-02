from fractions import Fraction as Q
from math import comb
import unittest
from unittest.mock import patch
from types import SimpleNamespace
from flint import arb,arb_mat,ctx
import occupancy_birth_classes as birth
import birth_classes
from test_scalar_cover import endpoint


class FixedBirthClassTests(unittest.TestCase):
    def test_refined_builder_regenerates_census_and_applies_returns_before_density(self):
        data={'windows':32};local=[arb_mat([[1]])];returned=[arb_mat([[2]])];refined=[arb_mat([[3]])]
        checked=[{'checked':'fresh'}];order=[]
        args=SimpleNamespace(exact_feedback=True,updates=2,precision=192,groups=1,
            tilts=['.1','.2'],joint_return_through=3,lazy_density_through=6)
        def returns(d,m,c,z):
            self.assertIs(d,data);self.assertIs(m,local);self.assertIs(c,checked)
            order.append('return');return returned
        def density(d,m,z,through):
            self.assertIs(d,data);self.assertIs(m,returned);self.assertEqual(through,6)
            order.append('density');return refined
        with patch.object(birth,'actual',return_value=data),patch.object(birth,'outward',return_value=local), \
             patch('return_moment.actual_census',return_value=checked) as census, \
             patch('return_moment.refine_class_returns',side_effect=returns), \
             patch('lazy_density.refine_actual',side_effect=density), \
             patch('occupancy_model.placement',side_effect=lambda m,**kw:m) as placement:
            result=birth.build_operators(args)
            self.assertEqual(order,['return','density']*2)
            self.assertEqual(census.call_count,1)
            self.assertEqual(placement.call_count,2)
            self.assertIs(result['.1','1'][0],refined)
            birth.build_operators(args)
            self.assertEqual(census.call_count,2)
        for field,limit in (('joint_return_through',4),('lazy_density_through',32)):
            for invalid in (-1,limit+1,True,'1'):
                with self.assertRaises(ValueError):birth.build_operators(SimpleNamespace(**dict(vars(args),**{field:invalid})))

    def test_exact_masses_and_ordered_transition_products(self):
        ctx.prec=192;z=Q(3,4)
        images=[sum(((a>>i)&1)*v for i,v in enumerate((1,6,120))) for a in range(8)]
        for columns in ([1,2,4,3,5,7,6,1],[0]*8,[1]*8):
            for updates in (1,2,3,4):
                data=birth.prepare(images,columns,3,updates);levels=list(map(int,data['birth_class_levels']))
                classes=birth.class_masses(data,arb(3)/4);matrices=birth.outward_at_z(data,arb(3)/4)
                complete=birth.class_masses(data,arb(3)/4,include_zero=True)
                refined=birth.refine_zero(data,matrices,arb(3)/4)
                exact=[[[Q(0)]*8 for _ in range(8)] for _ in range(3)]
                masses=[[Q(0)]*len(levels) for _ in range(3)]
                for x in range(256):
                    j=int(bool(x&15))+int(bool(x>>4));chance=Q(1,comb(2,j)*15**j)
                    feedback=0
                    for bit,c in enumerate(columns):
                        if x>>bit&1:feedback^=c
                    if feedback:masses[j][levels.index(images[feedback].bit_count())]+=chance*z**x.bit_count()
                    for state in range(8):
                        weight=chance*z**(images[state]^x).bit_count()
                        if not state:exact[j][state][feedback]+=weight
                        else:
                            exact[j][state][state^feedback]+=weight/Q(2**updates)
                            for refreshed in range(1,8):
                                exact[j][state][refreshed^feedback]+=weight*(1-Q(1,2**updates))/7
                for row,truths in zip(classes,masses):
                    for value,truth in zip(row,truths):
                        self.assertGreaterEqual(endpoint(value),truth)
                        self.assertLessEqual(endpoint(value)-truth,Q(1,10**40))
                for j,(old,new,row) in enumerate(zip(matrices,refined,complete)):
                    self.assertGreaterEqual(endpoint(row[0]),exact[j][0][0])
                    self.assertLessEqual(endpoint(row[0])-exact[j][0][0],Q(1,10**40))
                    moment=(((1+z)**4-1)/15)**j
                    self.assertGreaterEqual(sum(map(endpoint,row)),moment)
                    self.assertLessEqual(sum(map(endpoint,row))-moment,Q(1,10**40))
                    self.assertGreaterEqual(endpoint(new[0,0]),exact[j][0][0])
                    self.assertLessEqual(endpoint(new[0,0]),endpoint(old[0,0]))
                    self.assertLessEqual(endpoint(new[0,0])-exact[j][0][0],Q(1,10**40))
                    for a in range(old.nrows()):
                        for b in range(old.ncols()):
                            if a or b:self.assertEqual(new[a,b],old[a,b])
                matrices=refined
                bounds=[[[endpoint(matrix[i,k]) for k in range(matrix.ncols())]
                         for i in range(matrix.nrows())] for matrix in matrices]
                for sequence in ((0,0,0),(1,1,1),(2,2,2),(0,1,2),(2,0,1),(1,0,0)):
                    values=[Q(1)]*8;upper=[Q(1)]*len(bounds[0])
                    for j in sequence:
                        values=[sum(a*v for a,v in zip(row,values)) for row in exact[j]]
                        upper=[sum(a*v for a,v in zip(row,upper)) for row in bounds[j]]
                        self.assertGreaterEqual(upper[0],values[0])
                        self.assertGreaterEqual(upper[1],max(values[1:]))
                        self.assertGreaterEqual(upper[2],sum(values[1:])/7)
                        for i,level in enumerate(levels,3):
                            self.assertGreaterEqual(upper[i],max(values[a] for a in range(1,8) if images[a].bit_count()==level))

    def test_invalid_output_weight(self):
        data=birth.prepare(list(range(8)),[1,2,4,3,5,7,6,1],3)
        for z in (arb(0),arb(2)):
            with self.assertRaises(ValueError):birth.class_masses(data,z)
        for tilt in (0,-1):
            with self.assertRaises(ValueError):birth.outward(data,Q(tilt))
        with self.assertRaises(ValueError):birth.class_masses(data,arb(1),include_zero=1)
        with self.assertRaises(ValueError):birth.refine_zero(data,[],arb(1))

    def test_fixed_occupancy_mixture_matches_iid_births(self):
        ctx.prec=192
        images=[sum(((a>>i)&1)*v for i,v in enumerate((1,6,120))) for a in range(8)]
        data=birth.prepare(images,[1,2,4,3,5,7,6,1],3)
        W=data['windows']
        for z in (Q(1,2),Q(3,4),Q(1)):
            conditioned=birth.class_masses(data,birth.aq(z))
            for p in (Q(0),Q(1,5),Q(1,2),Q(1)):
                weights=[comb(W,j)*p**j*(1-p)**(W-j) for j in range(W+1)]
                for cutoff in (None,1,2,3):
                    iid=birth_classes.weighted_classes(data,p,birth.aq(z),cutoff)
                    terms=W+1 if cutoff is None else cutoff
                    for i,value in enumerate(iid):
                        mixture=sum((birth.aq(weights[j])*conditioned[j][i] for j in range(terms)),arb(0))
                        self.assertLess(abs(endpoint(value)-endpoint(mixture)),Q(1,10**40))


if __name__=='__main__':unittest.main()
