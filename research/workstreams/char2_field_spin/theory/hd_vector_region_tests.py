from fractions import Fraction as F
from itertools import combinations,product
from math import comb
import unittest

from hd_vector_region import local,mm,mv,IDENTITY,exact_region,exact_profile,exact_active_envelope


def mul(a,b):
    out=0
    while b:
        if b&1:out^=a
        b>>=1;a<<=1
        if a&4:a^=7
    return out


def dot(a,b):
    out=0
    for x,y in zip(a,b):out^=mul(x,y)
    return out


def rank(matrix):
    matrix=[list(row) for row in matrix]
    if not matrix:return 0
    r=0
    for j in range(len(matrix[0])):
        pivot=next((i for i in range(r,len(matrix)) if matrix[i][j]),None)
        if pivot is None:continue
        matrix[r],matrix[pivot]=matrix[pivot],matrix[r]
        inv=next(v for v in (1,2,3) if mul(v,matrix[r][j])==1)
        matrix[r]=[mul(inv,x) for x in matrix[r]]
        for i in range(len(matrix)):
            if i!=r and matrix[i][j]:
                a=matrix[i][j]
                matrix[i]=[x^mul(a,y) for x,y in zip(matrix[i],matrix[r])]
        r+=1
        if r==len(matrix):break
    return r


class VectorRegionTests(unittest.TestCase):
    def test_exact_local_law_all_small_states(self):
        y=F(2,5)
        for r in (1,2):
            vectors=list(product(range(4),repeat=r))
            for old in vectors:
                for x in (0,1):
                    actual=[F(0),F(0)]
                    for d,f,c in product(vectors,vectors,range(1,4)):
                        output=dot(f,old)^mul(c,x)
                        following=tuple(s^mul(a,x) for s,a in zip(old,d))
                        actual[int(any(following))]+=y**bool(output)
                    denominator=len(vectors)**2*3
                    self.assertEqual(tuple(v/denominator for v in actual),
                                     local(y,4,r)[bool(x)][bool(any(old))])
            for tau in (F(1),F(2),F(4)):
                self.assertEqual(local(y,4,r,tau)[1][0],local(y,4,r)[1][0])

    def test_region_keeps_chronology_and_binomial_normalization(self):
        y=F(1,3);t=local(y,4,2)
        for rows in range(1,6):
            for h in range(rows+1):
                total=[[F(0),F(0)],[F(0),F(0)]]
                for positions in combinations(range(rows),h):
                    m=IDENTITY
                    for i in range(rows):m=mm(m,t[i in positions])
                    for i,j in product(range(2),repeat=2):total[i][j]+=m[i][j]
                expected=tuple(tuple(v/comb(rows,h) for v in row) for row in total)
                self.assertEqual(exact_region(rows,h,y,4,2),expected)
        self.assertNotEqual(exact_profile([0,1],3,y,4,2),exact_profile([1,0],3,y,4,2))

    def test_whole_state_rank_charge_and_joint_return(self):
        states=(((1,0),(0,1)),((1,0),(1,0)),((1,0,0),(0,1,0)),
                ((1,0,1),(0,1,1)),((1,1,0),(2,2,0)))
        for s in states:
            width=len(s[0])
            inputs=[tuple(int(i==j) for i in range(width)) for j in range(width)]
            for x in inputs:
                for tau in (F(1),F(2),F(4)):
                    g=1/tau+(1-1/tau)/16
                    output_total=F(0);joint_total=F(0)
                    for f,c in product(product(range(4),repeat=2),range(1,4)):
                        eq=tuple(mul(c,v)^dot(f,column) for v,column in zip(x,zip(*s)))
                        r0=rank((eq,))
                        live=rank((eq,)+s)>r0 and rank((eq,x))>r0
                        if live:output_total+=tau**(-r0)
                        for d in product(range(4),repeat=2):
                            equations=(eq,)+tuple(tuple(a^mul(b,v) for a,v in zip(row,x))
                                                  for row,b in zip(s,d))
                            rr=rank(equations)
                            if rank(equations+s)>rr and rank(equations+(x,))>rr:
                                joint_total+=tau**(-rr)
                    self.assertLessEqual(output_total/48,g)
                    self.assertLessEqual(joint_total/(48*16),g/tau**2)
                    self.assertLessEqual(joint_total/(48*16),1/(4*tau**2))
        # Whole-state liveness is strictly weaker than a fixed live scalar pivot.
        self.assertGreater(F(1,16)+F(15,16)/2,F(1,2))

    def test_exact_active_majorant_for_repetition_rows(self):
        for rows in (1,2,3):
            for a in range(1,rows+1):
                y=F(3,7)
                actual=comb(rows,a)*3**(a-1)*exact_profile([a,a],rows,y,4,2)
                for lam in (F(1),F(3,2)):
                    bound=exact_active_envelope(a,2,1,rows,y,4,2,lam=lam)
                    self.assertLessEqual(actual,bound)

    def test_complete_positive_operator_with_live_history(self):
        z=(0,0,0);u=(1,0,0);v=(0,1,0);w=(0,0,1)
        cases=(((z,z),z,()),((z,z),u,((v,),)),
               ((u,z),z,()),((u,z),v,()),((u,u),u,((v,),)),
               ((u,v),z,()),((u,v),u,()),((u,v),w,((u,v),)))
        y=F(2,5)
        for s,x,past in cases:
            active=any(x);old=rank(s)>0
            required=list(past)+([(x,)] if active else [])+([s] if old else [])
            totals=[[F(0),F(0)] for _ in (1,2,4)]
            for f,d,c in product(product(range(4),repeat=2),product(range(4),repeat=2),range(1,4)):
                output=tuple(mul(c,a)^dot(f,column) for a,column in zip(x,zip(*s)))
                new=tuple(tuple(a^mul(b,t) for a,t in zip(row,x)) for row,b in zip(s,d))
                for select_output,weight in ((False,y),(True,1-y)):
                    base=(output,) if select_output else ()
                    for new_live in (False,True):
                        equations=base if new_live else base+new
                        rr=rank(equations)
                        histories=required+([new] if new_live else [])
                        if all(rank(equations+tuple(m))>rr for m in histories):
                            for i,tau in enumerate((1,2,4)):
                                totals[i][new_live]+=weight*F(1,tau**rr)
            for tau,actual in zip((F(1),F(2),F(4)),totals):
                bound=local(y,4,2,tau)[active][old]
                for observed,limit in zip(actual,bound):
                    self.assertLessEqual(observed/768,limit)


if __name__=='__main__':
    unittest.main()
