from fractions import Fraction as Q
from math import comb
import unittest
import numpy as np
from flint import arb,ctx
import refresh_kernel
import occupancy_kernel as kernel
import density_kernel
import single_packet
import occupancy_rank
from test_scalar_cover import endpoint


class FixedOccupancyTests(unittest.TestCase):
    kernel=kernel
    refined=False
    classes=False
    exact_feedback=False
    def test_feedback_moments_and_repeated_transition(self):
        ctx.prec=192
        images=[sum(((s>>i)&1)*v for i,v in enumerate((1,6,120))) for s in range(8)]
        columns=[1,2,4,3,5,7,6,1]
        z=Q(3,4)
        for updates in (1,2,3):
            data=kernel.prepare(refresh_kernel.prepare(images,columns,3,updates),self.exact_feedback)
            if self.refined:
                single=single_packet.census(images,columns,arb(3)/4)
                operators=self.kernel.outward_at_z(data,arb(3)/4,single,classes=self.classes)
            else:operators=self.kernel.outward_at_z(data,arb(3)/4)
            moments=[kernel.polynomial(h,arb(3)/4) for h in data['histograms']]
            for q in range(3):
                inputs=[x for x in range(256) if sum(bool((x>>offset)&15) for offset in (0,4))==q]
                self.assertEqual(len(inputs),comb(2,q)*15**q)
                feedback=[]
                for x in inputs:
                    syndrome=0
                    for b,c in enumerate(columns):
                        if x>>b&1:syndrome^=c
                    feedback.append(syndrome)
                self.assertEqual(data['zero_probabilities'][q],Q(feedback.count(0),len(inputs)))
                for target in range(1,8):
                    self.assertGreaterEqual(data['nonzero_atom_caps'][q],Q(feedback.count(target),len(inputs)))
                exact=[[Q(0)]*8 for _ in range(8)]
                for state in range(8):
                    total=Q(0)
                    for x,syndrome in zip(inputs,feedback):
                        weight=z**(images[state]^x).bit_count()/len(inputs);total+=weight
                        if state==0:exact[state][syndrome]+=weight
                        else:
                            exact[state][state^syndrome]+=weight/Q(2**updates)
                            for refreshed in range(1,8):
                                exact[state][refreshed^syndrome]+=weight*(1-Q(1,2**updates))/7
                    if state:
                        histogram=tuple(sum(((images[state]>>offset)&15).bit_count()==w for offset in (0,4)) for w in range(5))
                        index=list(map(tuple,data['histograms'])).index(histogram)
                        self.assertLessEqual(-endpoint(-moments[index][q]),total)
                        self.assertGreaterEqual(endpoint(moments[index][q]),total)
                self.assertGreaterEqual(endpoint(operators[q][0,0]),exact[0][0])
                self.assertGreaterEqual(endpoint(operators[q][0,1]),sum(exact[0][1:]))
                if len(self.kernel.TERMINAL)==3:
                    for state in range(1,8):self.assertGreaterEqual(endpoint(operators[q][1,0]),exact[state][0])
                    self.assertGreaterEqual(endpoint(operators[q][2,0]),sum(row[0] for row in exact[1:])/7)
                values=[Q(1)]*8
                for steps in range(1,9):
                    values=[sum(a*b for a,b in zip(row,values)) for row in exact]
                    upper=operators[q]**steps
                    terminal=self.kernel.TERMINAL
                    if self.classes:
                        terminal=np.ones(upper.nrows());terminal[2]=0
                    sums=[endpoint(sum((upper[i,j] for j,t in enumerate(terminal) if t),arb(0))) for i in range(len(terminal))]
                    self.assertGreaterEqual(sums[0],values[0])
                    arbitrary=sums[1] if len(terminal)==3 else sums[1]+sums[2]
                    self.assertGreaterEqual(arbitrary,max(values[1:]))
                    if self.classes:
                        for i,v in enumerate(sorted(set(x.bit_count() for x in images[1:]))):
                            selected=[values[a] for a in range(1,8) if images[a].bit_count()==v]
                            self.assertGreaterEqual(sums[3+i],sum(selected)/len(selected))
                    else:self.assertGreaterEqual(sums[-1],sum(values[1:])/7)
            self.assertEqual(operators[0][1,0],0)
            self.assertEqual(operators[0][2,0],0)


class DensityOccupancyTests(FixedOccupancyTests):
    kernel=density_kernel


class RefinedDensityTests(DensityOccupancyTests):
    refined=True


class ExpansionClassTests(RefinedDensityTests):
    classes=True


class ExactFeedbackOccupancyTests(FixedOccupancyTests):
    exact_feedback=True


class ExactFeedbackDensityTests(ExpansionClassTests):
    exact_feedback=True


class RankOccupancyTests(FixedOccupancyTests):
    kernel=occupancy_rank
    exact_feedback=True


if __name__=='__main__':unittest.main()
