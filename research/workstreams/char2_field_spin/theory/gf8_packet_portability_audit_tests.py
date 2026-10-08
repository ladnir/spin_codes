"""Independent audit controls for the width-eight packet transfer.

These checks support the algebraic portability argument. They are not a
finite distance certificate and do not change the screened construction.
"""
from fractions import Fraction as F
from itertools import combinations, product
import unittest

from diagonal_refresh_rank import rank_potential_bound
import gf8_packet_maps as maps
from gf8_packet_width_screen import joint_terms, local_matrix


def scalar_rank(matrix):
    """Gauss-Jordan reference, separate from the candidate's rank helper."""
    if not matrix:return 0
    work=[list(row) for row in matrix]
    found=0
    for column in range(len(work[0])):
        pivots=[i for i in range(found,len(work)) if work[i][column]]
        if not pivots:continue
        pivot=pivots[0]
        work[pivot],work[found]=work[found],work[pivot]
        inverse=maps.INV[work[found][column]]
        work[found]=[maps.MUL[inverse][x] for x in work[found]]
        for i in range(len(work)):
            if i==found:continue
            coefficient=work[i][column]
            work[i]=[x^maps.MUL[coefficient][y] for x,y in zip(work[i],work[found])]
        found+=1
        if found==len(work):break
    return found


class PacketPortabilityAudit(unittest.TestCase):
    def test_width8_all_return_support_pairs_and_joint_patterns(self):
        a,c=maps.matrices(8,4)
        # An invertible nonidentity diagonal state change gives A',C' that
        # remain superregular and C'A'=0, without requiring C'=(A')^T.
        scale=(1,2,3,4)
        a=tuple(tuple(maps.MUL[value][maps.INV[scale[j]]] for j,value in enumerate(row)) for row in a)
        c=tuple(tuple(maps.MUL[scale[i]][value] for value in row) for i,row in enumerate(c))
        self.assertNotEqual(c,tuple(zip(*a)))
        self.assertTrue(all(x==0 for row in maps.matmul(c,a) for x in row))
        count=0
        for state_mask in range(1,16):
            hs=[j for j in range(4) if state_mask>>j&1]
            h=len(hs)
            for input_mask in range(1,256):
                js=[j for j in range(8) if input_mask>>j&1]
                w=len(js)
                if h+w<=4:continue
                inactive=[j for j in range(8) if not input_mask>>j&1]
                for active_zeros in ((),tuple(js[::2]),tuple(js)):
                    for z in range(min(h-1,len(inactive))+1):
                        zero_rows=list(active_zeros)+inactive[:z]
                        return_rows=[[int(v==j) for j in hs]+[c[v][j] for j in js] for v in range(4)]
                        output_rows=[[maps.MUL[2][a[v][j]] for j in hs]+[int(v==j) for j in js] for v in zero_rows]
                        value=scalar_rank(return_rows+output_rows)
                        i=len(active_zeros)
                        full=max(4,i+z,h+w-(8-i-z))
                        relaxed=max(4,i+z,w-8+i+2*z+1)
                        self.assertGreaterEqual(value,full)
                        self.assertGreaterEqual(full,relaxed)
                        self.assertLessEqual(4-h,min(w-1,3-z))
                        count+=1
        self.assertGreater(count,20000)

    def test_width8_all_joint_multiplicities(self):
        minima=maps.proved_return_minima(8,4)
        for weight in range(1,9):
            active=set(range(weight))
            observed={}
            for mask in range(256):
                i=sum(bool(mask>>j&1) for j in active)
                z=sum(bool(mask>>j&1) for j in range(weight,8))
                if z<=min(3,8-weight) and 8-i-z>=minima[weight-1]:
                    observed[i,z]=observed.get((i,z),0)+1
            analytic={(i,z):count for i,z,count,_ in joint_terms(weight,8,4,minima)}
            self.assertEqual(observed,analytic)

    def test_small_field_fresh_diagonal_rank_moment(self):
        def multiply(a,b):
            result=0
            while b:
                if b&1:result^=a
                a<<=1
                if a&4:a^=7
                b>>=1
            return result
        inverse=(0,1,3,2)
        # One Cauchy row at node2 against state nodes0,1.
        row=(inverse[2],inverse[3])
        spaces=[((1,t),) for t in range(4)]+[((0,1),),((1,0),(0,1))]
        cases=0
        for basis in spaces:
            rank=len(basis)
            for tau in (F(1),F(2),F(3)):
                total=F(0)
                for diagonal in product(range(1,4),repeat=2):
                    equations=[multiply(multiply(row[0],diagonal[0]),v[0])
                               ^multiply(multiply(row[1],diagonal[1]),v[1]) for v in basis]
                    successor=rank-int(any(equations))
                    if successor:total+=tau**(successor-rank)
                self.assertLessEqual(total/9,rank_potential_bound(rank,1,tau,4,2))
                cases+=1
        self.assertEqual(cases,18)


if __name__=='__main__':unittest.main()
