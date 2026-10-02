"""Independent small-code checks of GF16 shell / MacWilliams constraints."""
import itertools
import copy
import unittest
from unittest.mock import patch
from fractions import Fraction as Q
from math import comb
from shared_relaxed_alternative import krawtchouk,shell_caps,Constraints
from shared_relaxed_alternative_mixture import candidates
from shared_mixture import verify
from shared_relaxed_alternative_counts import build_components


def enumerator(words,n,g):
    result=[0]*(n+1)
    for rows in itertools.product(words,repeat=g):
        union=0
        for row in rows:union|=row
        result[union.bit_count()]+=1
    return result


class Check(unittest.TestCase):
    def test_recurrence(self):
        for q in (2,4,16):
            for n in range(1,9):
                for u in range(n+1):
                    actual=krawtchouk(n,u,n,q)
                    expected=[sum((-1)**l*(q-1)**(j-l)*comb(u,l)*comb(n-u,j-l)
                                  for l in range(max(0,j-(n-u)),min(u,j)+1)) for j in range(n+1)]
                    self.assertEqual(actual,expected)

    def test_scalar_extension_constraints(self):
        for n,generators in ((3,[7]),(3,[3,5]),(4,[3,12])):
            words=[0]
            for row in generators:words += [v^row for v in words]
            dual=[x for x in range(1<<n) if all((x&r).bit_count()%2==0 for r in generators)]
            k=len(generators)
            for g in (1,2,4):
                q=1<<g
                actual=[enumerator(w,n,g) for w in (words,dual)]
                cdfs=[];caps=[]
                for side,rows in enumerate(actual):
                    cdf=[];total=0
                    for v in rows:total+=v;cdf.append(total)
                    cdfs.append(cdf)
                    other=dual if side==0 else words
                    d=min((v.bit_count() for v in other if v),default=n+1)
                    cap=shell_caps(n,k if side==0 else n-k,d,q)
                    cap[0]=1
                    self.assertTrue(all(a<=b for a,b in zip(rows,cap)))
                    caps.append(cap)
                model=Constraints(cdfs,caps,k,q)
                model.check_values(actual)
                for u in range(1,n+1):
                    row=model.solve(u,.5)
                    self.assertGreaterEqual(row['cap'],actual[0][u])
                    self.assertGreaterEqual(model.transform_cap(u),actual[0][u])

    def test_bad_multipliers(self):
        model=Constraints([[1,1,4],[1,1,4]],[[1,0,3],[1,0,3]],1,4)
        with self.assertRaises(ValueError):model.verify(2,[[0,'-1']])
        with self.assertRaises(ValueError):model.verify(2,[[0,'1'],[0,'2']])

    def test_pruned_method_selection(self):
        caps=[0,3,9,7]
        rows=list(candidates(caps,zero_bits=None,methods=['pruned-1/4']))
        self.assertEqual([row[0] for row in rows],['pruned-1/4'])
        verify(caps,rows[0][1])
        with self.assertRaises(ValueError):list(candidates(caps,zero_bits=None,methods=['unknown']))

    def test_saved_mixture_rechecks_caps(self):
        caps=[0]+[int(u>=128) for u in range(1,257)]
        source=dict(witnesses=[],iterations=0,stages=[])
        # p=99/128 is on the actual step-eight grid starting at38/256.
        saved=[dict(mass=str(1<<300),activity='99/128')]
        with patch('shared_relaxed_alternative_counts.actual',return_value=(caps,source)):
            components,replayed,mixture,metadata=build_components([],saved_mixture=saved)
            self.assertEqual(caps,replayed)
            self.assertEqual(metadata['mixture_mode'],'saved-rationals-freshly-verified')
            self.assertEqual(len(metadata['cap_sha256']),64)
            self.assertEqual(len(components),2)
            verify(caps,mixture)
            with self.assertRaises(ArithmeticError):
                build_components([],saved_mixture=[dict(mass='1',activity='99/128')])
            with self.assertRaises(ValueError):
                build_components([],saved_mixture=[dict(mass=str(1<<300),activity='3/4')])

    def test_pruning_saved_start_is_componentwise(self):
        caps=[0]+[int(u>=128) for u in range(1,257)]
        source=dict(witnesses=[],iterations=0,stages=[])
        start=[dict(mass=str(1<<300),activity='99/128')]
        with patch('shared_relaxed_alternative_counts.actual',return_value=(caps,source)):
            _,_,mixture,metadata=build_components([],starting_mixture=start)
            verify(caps,mixture)
            self.assertEqual(metadata['mixture_mode'],'saved-start-componentwise-pruned')
            self.assertTrue(all(p==Q(99,128) and c<=1<<300 for c,p in mixture))
            with self.assertRaises(ValueError):
                build_components([],starting_mixture=start,saved_mixture=start)
            with patch('positive_prune.prune',return_value=([(Q(1<<301),Q(99,128))],{})):
                with self.assertRaises(ArithmeticError):
                    build_components([],starting_mixture=start)

    def test_joint_complement_constraints_when_both_codes_contain_ones(self):
        from constraint_counts import Constraints as JointConstraints
        from dual_moments import exact_tuple_cdfs
        from bch_joint_support import rank_total
        from test_constraint_counts import dimensions
        for n,generators in ((4,[3,12]),(6,[3,12,48])):
            words=[0]
            for row in generators:words += [v^row for v in words]
            dual=[x for x in range(1<<n) if all((x&r).bit_count()%2==0 for r in generators)]
            self.assertIn((1<<n)-1,words)
            self.assertIn((1<<n)-1,dual)
            cdfs=[exact_tuple_cdfs(code,n,4) for code in (words,dual)]
            inflated=[[[2*x for x in row] for row in side] for side in cdfs]
            spectra=[[sum(w.bit_count()==u for w in code) for u in range(n+1)]
                     for code in (words,dual)]
            for last in (n//2,n):
                system=JointConstraints(*inflated,dimensions(words,n),dimensions(dual,n),
                                        len(generators),last,spectra=spectra,ones=(True,True))
                values=[Q(cdfs[side][h-1][u]//rank_total(h,4,h),system.caps[side][h-1][u])
                        for side,h,u in system.keys]
                for row,rhs in system.equalities:
                    self.assertEqual(sum(a*values[i] for i,a in row.items()),rhs)
                for row,rhs in system.inequalities:
                    self.assertLessEqual(sum(a*values[i] for i,a in row.items()),rhs)

    def test_exact_composition_of_coordinate_duals(self):
        from constraint_counts import normalize,verify_dual
        from shared_relaxed_alternative_dual_repair import compose
        rows=[normalize({0:1,1:1},1),normalize({0:4},1)]
        y,z,upper,trace=compose(rows,[],{0:Q(1)},[Q(1),Q(0)],[],
                                {0:([Q(0),Q(2)],[])},2,4)
        self.assertEqual(upper,Q(5,8))
        self.assertEqual(verify_dual(rows,[],{0:Q(1)},y,z,2)[0],upper)
        self.assertEqual(len(trace),2)
        for x in map(Q,(0,.125,.25)):
            for yy in map(Q,(0,.125,.5,.75)):
                if x+yy<=1:self.assertLessEqual(x,upper)
        with self.assertRaises(ValueError):
            compose(rows,[],{0:Q(1)},[Q(1),Q(0)],[],{0:([Q(0),Q(-2)],[])},2)

    def test_saved_probe_retarget_clears_bounds_and_preserves_source(self):
        from shared_relaxed_alternative_probe_saved import prepare_record
        source=dict(schema='shared-gf16-relaxed-dense-1',updates=2,K=1<<20,N=1<<21,
            minimum_groups=33,zero_bits=64,variance_bins=16,base_tilt='3/16',cost_tilt='1/4',
            pruned=True,mixture=[dict(mass='10',activity='99/128')],cap_sha256='b'*64,
            results=[dict(distance='7/100',threshold=146800,root=['0','1'],
                probes=[dict(mean='1/32',upper=[1,-999])],cover=dict(leaves={'':{}},unresolved={}),
                rechecked_precision=384)])
        before=copy.deepcopy(source)
        for distance,cutoff in ((None,146800),('3/40',157286),('2/25',167772)):
            result=prepare_record(source,'a'*64,distance)
            row=result['results'][0]
            self.assertEqual(source,before)
            self.assertEqual(result['mixture'],source['mixture'])
            self.assertEqual(row['threshold'],cutoff)
            self.assertEqual(row['probes'],[])
            self.assertNotIn('cover',row)
            self.assertNotIn('rechecked_precision',row)
            self.assertTrue(result['screen_only'])
            if distance is not None:
                self.assertEqual(row['retargeted_from']['source_sha256'],'a'*64)
                self.assertEqual(row['retargeted_from']['prior_probes'],before['results'][0]['probes'])
        with self.assertRaises(ValueError):prepare_record(source,'a'*64,'1/2')


if __name__=='__main__':unittest.main()
