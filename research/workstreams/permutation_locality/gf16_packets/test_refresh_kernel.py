from fractions import Fraction as Q
from collections import Counter
import unittest
import numpy as np
from flint import arb,arb_mat,ctx
import scalar_cover as sc
import refresh_kernel as rk
import feedback_refresh
import gf_refresh
import weighted_return
import profile_return
import shape_return
import conditioned_return
import trimmed_return
import rank_return
import birth_refresh
from test_scalar_cover import endpoint


class RefreshKernelTests(unittest.TestCase):
    kernel=rk
    def setUp(self):
        ctx.prec=192
        self.images=[sum(((s>>i)&1)*v for i,v in enumerate((1,6,120))) for s in range(8)]
        self.columns=[1,2,4,3,5,7,6,1]

    def test_fast_histograms_match_existing_and_counts(self):
        for W in (1,2,32):
            columns=(self.columns*16)[:4*W]
            images=list(range(8)) if W==1 else self.images
            data=self.kernel.prepare(images,columns,3)
            old=sc.kernel.prepare(images,columns,3)
            for name in ('records','multiplicities','histograms'):np.testing.assert_array_equal(data[name],old[name])
            expected=Counter(tuple(sum(((x>>(4*w))&15).bit_count()==j for w in range(W)) for j in range(5)) for x in images[1:])
            self.assertEqual(dict(zip(map(tuple,data['histograms']),map(int,data['histogram_multiplicities']))),expected)

    def test_every_start_state_and_uniform_density_through_eight_steps(self):
        for updates in (1,2,3):
            data=self.kernel.prepare(self.images,self.columns,3,updates)
            for activity in (Q(0),Q(1,5),Q(15,16),Q(1)):
                z=Q(3,4);probabilities=sc.probabilities(activity)
                matrix=self.kernel.outward_at_z(data,probabilities,arb(3)/4)
                exact=[[Q(0) for _ in range(8)] for _ in range(8)]
                lazy=[[Q(0) for _ in range(8)] for _ in range(8)]
                fresh=[[Q(0) for _ in range(8)] for _ in range(8)]
                for state in range(8):
                    for x in range(256):
                        probability=Q(1)
                        for offset in (0,4):probability*=activity/15 if (x>>offset)&15 else 1-activity
                        feedback=0
                        for b,c in enumerate(self.columns):
                            if x>>b&1:feedback^=c
                        weight=probability*z**(self.images[state]^x).bit_count()
                        if state==0:exact[state][feedback]+=weight
                        else:
                            lazy[state][state^feedback]+=weight/Q(2**updates)
                            for refreshed in range(1,8):fresh[state][refreshed^feedback]+=weight*(1-Q(1,2**updates))/7
                    if state:
                        exact[state]=[a+b for a,b in zip(lazy[state],fresh[state])]
                        self.assertGreaterEqual(endpoint(matrix[1,0]),exact[state][0])
                        self.assertGreaterEqual(endpoint(matrix[1,1]),sum(lazy[state][1:]))
                        for target in range(1,8):self.assertGreaterEqual(endpoint(matrix[1,2])/7,fresh[state][target])
                self.assertGreaterEqual(endpoint(matrix[0,0]),exact[0][0])
                # New births may be represented by arbitrary mass plus a
                # uniform-density envelope. With U=0 this is the original
                # arbitrary-mass check, unchanged for the earlier kernels.
                residual=sum(max(Q(0),value-endpoint(matrix[0,2])/7) for value in exact[0][1:])
                self.assertGreaterEqual(endpoint(matrix[0,1]),residual)
                self.assertGreaterEqual(endpoint(matrix[2,0]),sum(row[0] for row in exact[1:])/7)
                self.assertGreaterEqual(endpoint(matrix[2,1]),sum(sum(row[1:]) for row in lazy[1:])/7)
                for target in range(1,8):self.assertGreaterEqual(endpoint(matrix[2,2])/7,sum(row[target] for row in fresh[1:])/7)
                values=[Q(1)]*8
                for steps in range(1,9):
                    values=[sum(a*b for a,b in zip(row,values)) for row in exact]
                    upper=matrix**steps
                    sums=[endpoint(sum((upper[i,j] for j in range(3)),arb(0))) for i in range(3)]
                    self.assertGreaterEqual(sums[0],values[0])
                    for value in values[1:]:self.assertGreaterEqual(sums[1],value)
                    self.assertGreaterEqual(sums[2],sum(values[1:])/7)

    def test_float_and_outward_agree(self):
        data=self.kernel.prepare(self.images,self.columns,3)
        for activity in (Q(1,100),Q(1,2),Q(15,16),Q(1)):
            p=sc.probabilities(activity)
            for tilt in (Q(1,100),Q(1,4),Q(1)):
                proposed=self.kernel.floating(data,p,float(tilt));exact=self.kernel.outward(data,p,tilt)
                np.testing.assert_allclose(proposed,np.array([[float(exact[i,j]) for j in range(3)] for i in range(3)]),rtol=2e-12,atol=1e-15)


class FeedbackRefreshTests(RefreshKernelTests):
    kernel=feedback_refresh

    def test_empty_input_cannot_cancel_a_nonzero_state(self):
        data=self.kernel.prepare(self.images,self.columns,3)
        matrix=self.kernel.outward_at_z(data,sc.probabilities(Q(0)),arb(3)/4)
        self.assertEqual(matrix[1,0],0)
        self.assertEqual(matrix[2,0],0)


class GFRefreshTests(FeedbackRefreshTests):
    kernel=gf_refresh

    def test_wrong_packet_law_is_rejected(self):
        data=self.kernel.prepare(self.images,self.columns,3)
        with self.assertRaises(ValueError):self.kernel.outward_at_z(data,[0,1,0,0,0],arb(3)/4)


class WeightedReturnTests(GFRefreshTests):
    kernel=weighted_return

    def test_tilted_feedback_moments_directly(self):
        data=self.kernel.prepare(self.images,self.columns,3)
        for p in (Q(0),Q(1,5),Q(1)):
            for h in (Q(1),Q(4,3),Q(10)):
                exact=[Q(0)]*8
                for x in range(256):
                    probability=Q(1);syndrome=0
                    for offset in (0,4):probability*=p/15 if (x>>offset)&15 else 1-p
                    for bit,c in enumerate(self.columns):
                        if x>>bit&1:syndrome^=c
                    exact[syndrome]+=probability*h**x.bit_count()
                atom,nonzero=self.kernel.feedback_moments(data,sc.probabilities(p),arb(h.numerator)/h.denominator)
                self.assertGreaterEqual(endpoint(atom),max(exact[1:]))
                self.assertGreaterEqual(endpoint(nonzero),sum(exact[1:]))
                self.assertLess(endpoint(nonzero)-sum(exact[1:]),Q(1,2)**120)


class ProfileReturnTests(GFRefreshTests):
    kernel=profile_return

    def test_raw_lazy_bound_for_every_nonzero_state(self):
        data=self.kernel.prepare(self.images,self.columns,3);z=Q(3,4)
        for p in (Q(0),Q(1,5),Q(15,16)):
            bound=endpoint(self.kernel.lazy_bound(data,p,arb(3)/4))
            for state in range(1,8):
                exact=Q(0)
                for x in range(256):
                    probability=Q(1);syndrome=0
                    for offset in (0,4):probability*=p/15 if (x>>offset)&15 else 1-p
                    for bit,c in enumerate(self.columns):
                        if x>>bit&1:syndrome^=c
                    if syndrome==state:exact+=probability*z**(self.images[state]^x).bit_count()
                self.assertGreaterEqual(bound,exact)

    def test_canonical_profile_dominates_every_small_weight_assignment(self):
        from itertools import product
        p=Q(1,5);z=Q(3,4);a=1-16*p/15;b=p/15
        for W in (1,2,3):
            characters=list(product(range(5),repeat=W))
            records=np.array([[r.count(j) for j in range(5)] for r in characters])
            factors=[[a*z**w+b*(1+z)**(4-r)*(1-z)**r for r in range(5)] for w in range(5)]
            patterns=list(product(range(5),repeat=W))
            for distance in range(4*W+1):
                weights,profiles=self.kernel.profiles(records,distance)
                bound=[]
                for row in profiles:
                    value=Q(1)
                    for n,f in zip(row,[factors[w][r] for w in weights for r in range(5)]):value*=f**int(n)
                    bound.append(value-a**W*z**distance)
                for pattern in patterns:
                    if sum(pattern)<distance:continue
                    for index,character in enumerate(characters):
                        value=Q(1)
                        for w,r in zip(pattern,character):value*=factors[w][r]
                        self.assertGreaterEqual(bound[index],value-a**W*z**sum(pattern))


class ShapeReturnTests(GFRefreshTests):
    kernel=shape_return

    def test_exceptional_states_and_average_against_full_enumeration(self):
        images=[0,15,240,255];columns=[1,2,3,1,3,2,1,3];z=Q(3,4)
        for cols in (columns,[0]*8):
            data=self.kernel.prepare(images,cols,2,exact_exceptions=True)
            self.assertEqual([row['state'] for row in data['exceptions']],[1,2])
            self.assertEqual(sum(row['multiplicity'] for row in data['shapes']),1)
            for p in (Q(0),Q(1,5),Q(15,16)):
                maximum,mean=self.kernel.bounds(data,p,arb(3)/4)
                exact=[Q(0)]*4
                for x in range(256):
                    probability=Q(1);syndrome=0
                    for offset in (0,4):probability*=p/15 if (x>>offset)&15 else 1-p
                    for bit,c in enumerate(cols):
                        if x>>bit&1:syndrome^=c
                    if syndrome:exact[syndrome]+=probability*z**(images[syndrome]^x).bit_count()
                self.assertGreaterEqual(endpoint(maximum),max(exact[1:]))
                self.assertGreaterEqual(endpoint(mean),sum(exact[1:])/3)
                c=1-16*p/15;b=p/15
                factors=[arb(c.numerator)/c.denominator*(arb(3)/4)**w
                         +arb(b.numerator)/b.denominator*((-1)**r if w else 1)*(arb(7)/4)**(4-r)*(arb(1)/4)**r
                         for w in (0,4) for r in range(5)]
                for row in data['exceptions']:
                    values=self.kernel.arithmetic.arb_products(factors,row['tuples'],2)
                    result=sum((int(n)*value for n,value in zip(row['counts'],values)),arb(0))/4
                    self.assertGreaterEqual(endpoint(result),exact[row['state']])
                    self.assertGreaterEqual(endpoint(-result),-exact[row['state']])
                    moments=[]
                    for i,marginal in enumerate(row['marginals']):
                        squares=[v*v for v in factors[5*i:5*i+5]]
                        values=self.kernel.arithmetic.arb_products(squares,marginal['tuples'],2)
                        moments.append(sum((int(n)*value for n,value in zip(marginal['counts'],values)),arb(0))/4)
                    self.assertGreaterEqual(endpoint((moments[0]*moments[1]).sqrt()),exact[row['state']])


class ConditionedReturnTests(GFRefreshTests):
    kernel=conditioned_return

    def test_fixed_occupancy_moments_and_raw_return(self):
        data=self.kernel.prepare(self.images,self.columns,3)
        z=Q(3,4)
        for p in (Q(0),Q(1,5),Q(15,16),Q(1)):
            maximum,mean=self.kernel.bounds(data,p,arb(3)/4)
            exact=[Q(0)]*8
            for x in range(256):
                probability=Q(1);syndrome=0
                for offset in (0,4):probability*=p/15 if (x>>offset)&15 else 1-p
                for bit,c in enumerate(self.columns):
                    if x>>bit&1:syndrome^=c
                if syndrome:exact[syndrome]+=probability*z**(self.images[syndrome]^x).bit_count()
            self.assertGreaterEqual(endpoint(maximum),max(exact[1:]))
            self.assertGreaterEqual(endpoint(mean),sum(exact[1:])/7)
        for power in (1,*self.kernel.POWERS):
            fast=self.kernel.float_moments(data,float(z**power))
            exact=self.kernel.fixed.moments(data,(arb(3)/4)**power)
            for a,b in zip(fast,exact):np.testing.assert_allclose(a,list(map(float,b)),rtol=2e-13)


class TrimmedReturnTests(GFRefreshTests):
    kernel=trimmed_return

    def test_integer_census_matches_every_small_input(self):
        for y in (0,1,15,86,255):
            histogram=[sum(((y>>offset)&15).bit_count()==w for offset in (0,4)) for w in range(5)]
            counts=self.kernel.output_counts(histogram)
            exact=np.zeros((3,9),dtype=object)
            for x in range(256):exact[int(bool(x&15))+int(bool(x>>4)),(y^x).bit_count()]+=1
            np.testing.assert_array_equal(counts,exact)
        big=self.kernel.output_counts([32,0,0,0,0])
        self.assertEqual(sum(map(sum,big)),1<<128)

    def test_lightest_bound_dominates_every_event(self):
        from itertools import combinations
        counts=[1,0,2,1];weights=[0,2,2,3];z=Q(3,4)
        for size in range(5):
            trimmed=self.kernel.lightest(counts,size)
            bound=sum(n*z**w for w,n in enumerate(trimmed))
            self.assertEqual(sum(trimmed),size)
            for selected in combinations(range(4),size):self.assertGreaterEqual(bound,sum(z**weights[i] for i in selected))
        for size in (-1,5):
            with self.assertRaises(ValueError):self.kernel.lightest(counts,size)

    def test_each_fixed_occupancy_return_and_uniform_event(self):
        data=self.kernel.prepare(self.images,self.columns,3);z=Q(3,4)
        exact=[[Q(0)]*8 for _ in range(3)]
        for x in range(256):
            syndrome=0;j=int(bool(x&15))+int(bool(x>>4))
            for bit,c in enumerate(self.columns):
                if x>>bit&1:syndrome^=c
            if syndrome:exact[j][syndrome]+=z**(self.images[syndrome]^x).bit_count()
        hist_to_index={tuple(h):i for i,h in enumerate(data['histograms'])}
        for j in range(3):
            total=0
            for state,image in enumerate(self.images[1:],1):
                h=tuple(sum(((image>>offset)&15).bit_count()==w for offset in (0,4)) for w in range(5))
                row=data['trimmed'][hist_to_index[h]][j]
                bound=sum(n*z**w for w,n in enumerate(row));self.assertGreaterEqual(bound,exact[j][state]);total+=bound
            uniform=sum(n*z**w for w,n in enumerate(data['trimmed_uniform'][j]))
            self.assertGreaterEqual(uniform,sum(exact[j]))
            self.assertGreaterEqual(total,sum(exact[j]))


class RankReturnTests(GFRefreshTests):
    kernel=rank_return

    def test_raw_fixed_occupancy_rank_bound(self):
        z=Q(3,4)
        for columns in (self.columns,[0]*8,[1]*8):
            data=self.kernel.prepare(self.images,columns,3)
            hist_to_index={tuple(h):i for i,h in enumerate(data['histograms'])}
            for state,image in enumerate(self.images[1:],1):
                actual=[Q(0)]*3
                for x in range(256):
                    syndrome=0;j=int(bool(x&15))+int(bool(x>>4))
                    for bit,c in enumerate(columns):
                        if x>>bit&1:syndrome^=c
                    if syndrome==state:actual[j]+=z**(image^x).bit_count()
                h=tuple(sum(((image>>offset)&15).bit_count()==w for offset in (0,4)) for w in range(5))
                for j,row in enumerate(data['rank_counts'][hist_to_index[h]]):
                    bound=(1+z)**(4*j-3)*sum(int(n)*z**w for w,n in enumerate(row))
                    self.assertGreaterEqual(bound,actual[j])


class BirthRefreshTests(RefreshKernelTests):
    kernel=birth_refresh


if __name__=='__main__':unittest.main()
