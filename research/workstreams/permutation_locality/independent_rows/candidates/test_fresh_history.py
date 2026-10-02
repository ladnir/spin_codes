from collections import Counter,defaultdict
from fractions import Fraction as Q
from math import comb
from itertools import combinations,product
import unittest

from flint import arb,arb_mat,ctx
import fresh_history as history
from group_moment import maps
import cancellation_joint
import window_histogram


def endpoint(value):
    m,e=map(int,value.upper().man_exp())
    return Q(m)*Q(2)**e


class FreshHistoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec=192
        cls.images,columns,cls.spectrum=maps()
        cls.packets=defaultdict(list);seen=set()
        for window in range(32):
            for bits in range(1,16):
                state=0
                for bit in range(4):
                    if bits>>bit&1:state^=columns[4*window+bit]
                assert state and state not in seen;seen.add(state)
                cls.packets[bits.bit_count()].append((bits<<(4*window),state))
        cls.tilt=(arb(8)/7).log()
        joint=cancellation_joint.census(1,by_fresh_class=True)
        histogram=window_histogram.census()
        cls.details={'histograms':histogram,'moments':window_histogram.moments(histogram,cls.tilt,1),
                     'cancellations':joint[:2],'fresh_classes':joint[2],'spectrum':cls.spectrum}
        cls.feedback={}
        for a,packets in cls.packets.items():
            cls.feedback[a,]=(0,1,len(packets),dict(Counter(cls.images[s].bit_count() for _,s in packets)))

    def test_class_partition(self):
        for split in (False,True):
            records=history.classes(self.details,split)
            for a in range(1,5):
                self.assertEqual(sum(r['n'] for r in records if r['a']==a),32*comb(4,a))

    def test_every_single_packet_history_component(self):
        records=history.classes(self.details,True)
        single=history.single_density(records,self.tilt)
        powers=[7**w*8**(128-w) for w in range(129)]
        for b in range(1,5):
            old=arb_mat([[1]*11 for _ in range(11)])
            old[history.Z,history.F]=(-self.tilt*b).exp()
            matrix=history.family_matrix(old,1,(b,),self.details,self.feedback,records,
                                         self.tilt,'1',2,single)
            for index,r in enumerate(records,11):
                source=[s for _,s in self.packets[r['a']] if self.images[s].bit_count()==r['v']]
                self.assertEqual(len(source),r['n'])
                weighted=Counter()
                for s in source:
                    for x,feedback in self.packets[b]:
                        weighted[s^feedback]+=powers[(self.images[s]^x).bit_count()]
                denominator=len(source)*len(self.packets[b])*8**128
                total=Q(sum(weighted.values()),denominator);zero=Q(weighted[0],denominator)
                self.assertGreaterEqual(endpoint(matrix[index,history.Z]),zero/4+3*total/(4*((1<<19)-1)))
                self.assertGreaterEqual(endpoint(matrix[index,history.M]),(total-zero)/4)
                self.assertGreaterEqual(endpoint(matrix[index,history.C]),
                                        Q(max(n for s,n in weighted.items() if s),4*denominator))
                for target,cut in ((9,48),(10,56)):
                    exact=Q(sum(n for s,n in weighted.items() if s and self.images[s].bit_count()<=cut),4*denominator)
                    self.assertGreaterEqual(endpoint(matrix[index,target]),exact)
                for ell,v in enumerate(sorted(self.spectrum)):
                    self.assertGreaterEqual(endpoint(matrix[index,history.U+ell]),
                                            3*total*self.spectrum[v]/(4*((1<<19)-1)))
            self.assertEqual(matrix[history.Z,history.F],0)
            activation=sum((matrix[history.Z,i] for i in range(11,matrix.ncols())),arb(0))
            self.assertTrue(activation>=(-self.tilt*b).exp() or activation.overlaps((-self.tilt*b).exp()))

    def test_empty_epochs_preserve_fresh_class(self):
        records=history.classes(self.details,True)
        for rounds in (2,4,None):
            a=Q(0) if rounds is None else Q(1,2**rounds)
            matrix=history.family_matrix(arb_mat(11,11),0,None,self.details,self.feedback,
                                         records,self.tilt,'1',rounds)
            for i,r in enumerate(records,11):
                f=Q(7,8)**r['v']
                self.assertGreaterEqual(endpoint(matrix[i,i]),a*f)
                self.assertEqual(matrix[i,history.Z],0)
                for ell,v in enumerate(sorted(self.spectrum)):
                    self.assertGreaterEqual(endpoint(matrix[i,history.U+ell]),
                                            (1-a)*f*self.spectrum[v]/((1<<19)-1))

    def test_fresh_class_counts_preserve_existing_marginals(self):
        for shape,(_,prior,_) in self.details['cancellations'][0].items():
            records=self.details['fresh_classes'][shape]
            for a in range(1,5):
                self.assertEqual([sum(row[w] for (b,v),row in records.items() if b==a) for w in range(129)],prior[a-1])

    def test_two_window_fresh_census_against_independent_scalar_enumeration(self):
        data,_,refined=cancellation_joint.census(2,check_old=False,by_fresh_class=True)
        lookup={s:(a,self.images[s].bit_count()) for a,packets in self.packets.items() for _,s in packets}
        windows=defaultdict(list)
        for a,packets in self.packets.items():
            for word,s in packets:windows[(word.bit_length()-1)//4].append((a,word,s))
        expected=defaultdict(Counter)
        for left,right in combinations(range(32),2):
            for (a,x,s),(b,y,t) in product(windows[left],windows[right]):
                state=s^t
                if state in lookup:
                    incoming,v=lookup[state]
                    expected[tuple(sorted((a,b)))][incoming,v,(self.images[state]^x^y).bit_count()]+=1
        for shape,records in refined.items():
            if len(shape)!=2:continue
            for (a,v),row in records.items():
                self.assertEqual(row,[expected[shape][a,v,w] for w in range(129)])


if __name__=='__main__':unittest.main()
