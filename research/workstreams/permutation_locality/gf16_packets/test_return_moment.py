from fractions import Fraction as Q
import unittest
from unittest.mock import patch
from flint import arb,ctx
import return_moment
import fixed_joint_return_probe as joint
from test_scalar_cover import endpoint


class JointReturnTests(unittest.TestCase):
    def test_actual_census_checks_map_identity(self):
        images=list(range(8));columns=[1,2,4,3,5,7,6,1];data={'columns':columns,'bits':3}
        with patch('group_moment.maps',return_value=(images,columns,None)):
            rows=return_moment.actual_census(data,2)
            self.assertEqual([row['denominator'] for row in rows],[1,30,225])
        with patch('group_moment.maps',return_value=(images,[0]*8,None)):
            with self.assertRaises(ArithmeticError):return_moment.actual_census(data,2)

    def test_four_packet_joint_census(self):
        images=[sum(((s>>i)&1)*v for i,v in enumerate((1,6,65528))) for s in range(8)]
        columns=[1,2,4,3,5,7,6,1,2,3,4,5,6,7,1,4]
        row=return_moment.census(images,columns,3,4)[4]
        expected=[[0]*17 for _ in images]
        for x in range(1<<16):
            if any(not ((x>>shift)&15) for shift in (0,4,8,12)):continue
            feedback=0
            for bit,c in enumerate(columns):
                if x>>bit&1:feedback^=c
            expected[feedback][(images[feedback]^x).bit_count()]+=1
        actual=[[0]*17 for _ in images]
        for i,target in enumerate(row['targets']):
            end=row['starts'][i+1] if i+1<len(row['starts']) else len(row['counts'])
            for k in range(row['starts'][i],end):actual[int(target)][int(row['weights'][k])]=int(row['counts'][k])
        self.assertEqual(actual,expected)
        self.assertEqual(sum(map(sum,actual)),15**4)

    def test_refined_class_return_entries_dominate_exact_transitions(self):
        ctx.prec=192;z=Q(3,4)
        images=[sum(((s>>i)&1)*v for i,v in enumerate((1,6,120))) for s in range(8)]
        for columns in ([1,2,4,3,5,7,6,1],[0]*8,[1]*8):
            census=return_moment.census(images,columns,3,2)
            for updates in (2,3,4):
                data=joint.base.prepare(images,columns,3,updates)
                local=joint.base.outward_at_z(data,arb(3)/4)
                refined=joint.refine(data,local,census,arb(3)/4)
                truth=[[Q(0)]*8 for _ in range(3)]
                for x in range(256):
                    j=int(bool(x&15))+int(bool(x>>4));feedback=0
                    for bit,c in enumerate(columns):
                        if x>>bit&1:feedback^=c
                    chance=Q(1,census[j]['denominator'])
                    for a in range(8):
                        mass=chance*z**(images[a]^x).bit_count()
                        probability=(Q(int(feedback==0)) if not a else
                                     Q(int(feedback==a),2**updates)+(1-Q(1,2**updates))*Q(int(feedback!=0),7))
                        truth[j][a]+=mass*probability
                for j,matrix in enumerate(refined):
                    targets=[truth[j][0],max(truth[j][1:]),sum(truth[j][1:])/7,
                             *(max(truth[j][a] for a in range(1,8) if images[a].bit_count()==level)
                               for level in data['birth_class_levels'])]
                    for i,target in enumerate(targets):
                        self.assertGreaterEqual(endpoint(matrix[i,0]),target)
                        self.assertLessEqual(endpoint(matrix[i,0]),endpoint(local[j][i,0]))
                        for k in range(1,matrix.ncols()):self.assertEqual(matrix[i,k],local[j][i,k])

    def test_every_small_joint_count_and_moment(self):
        ctx.prec=192;images=[sum(((s>>i)&1)*v for i,v in enumerate((1,6,120))) for s in range(8)]
        columns=[1,2,4,3,5,7,6,1];rows=return_moment.census(images,columns,3,2)
        for j,row in enumerate(rows):
            hist=[[0]*9 for _ in images]
            for x in range(256):
                if bool(x&15)+bool(x>>4)!=j:continue
                s=0
                for b,c in enumerate(columns):
                    if x>>b&1:s^=c
                hist[s][(images[s]^x).bit_count()]+=1
            reconstructed=[[0]*9 for _ in images]
            for i,s in enumerate(row['targets']):
                end=row['starts'][i+1] if i+1<len(row['starts']) else len(row['counts'])
                for k in range(row['starts'][i],end):reconstructed[int(s)][int(row['weights'][k])]=int(row['counts'][k])
            self.assertEqual(reconstructed,hist)
            for z in (Q(1),Q(3,4),Q(1,10)):
                exact=[sum(Q(n)*z**w for w,n in enumerate(h))/row['denominator'] for h in hist]
                bounded=return_moment.outward(row,arb(z.numerator)/z.denominator)
                for label,target in (('zero',exact[0]),('maximum',max(exact[1:])),('total',sum(exact[1:]))):
                    self.assertGreaterEqual(endpoint(bounded[label]),target)
                    self.assertLess(endpoint(bounded[label])-target,Q(1,2)**40)
                for level,value in bounded['classes'].items():
                    target=sum(exact[s] for s in range(1,8) if images[s].bit_count()==level)
                    self.assertGreaterEqual(endpoint(value),target)
                    maximum=max(exact[s] for s in range(1,8) if images[s].bit_count()==level)
                    self.assertGreaterEqual(endpoint(bounded['class_maxima'][level]),maximum)
                    self.assertLess(endpoint(bounded['class_maxima'][level])-maximum,Q(1,2)**40)


if __name__=='__main__':unittest.main()
