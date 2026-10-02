from fractions import Fraction as Q
from itertools import combinations,product
import unittest
import clipping


def solve(rows,rhs):
    a=[list(map(Q,row))+[Q(y)] for row,y in zip(rows,rhs)];n=len(rows)
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        if pivot is None:return None
        a[j],a[pivot]=a[pivot],a[j];v=a[j][j];a[j]=[x/v for x in a[j]]
        for i in range(n):
            if i!=j:
                v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[j])]
    return tuple(row[-1] for row in a)


class ClippingTests(unittest.TestCase):
    def test_simplex(self):
        box=(Q(0),Q(1))*4
        points=clipping.vertices(box,[((1,1,1,1),Q(1))])
        self.assertEqual(set(points),{(0,0,0,0),*(tuple(int(i==j) for j in range(4)) for i in range(4))})

    def test_clipping_against_all_constraint_intersections(self):
        box=(Q(0),Q(1))*4
        axes=[(tuple(sign*int(j==i) for j in range(4)),int(sign==1)) for i in range(4) for sign in (-1,1)]
        families=[[(tuple([1]*4),Q(1)),((1,0,0,0),Q(1,2)),((-1,1,1,0),Q(1,3))],
                  [((1,2,-1,0),Q(1,3)),((0,1,1,2),Q(3,2)),((-2,0,1,1),Q(1,4))],
                  [((1,1,1,1),Q(0))],
                  [((1,1,1,1),Q(-1))]]
        for extra in families:
            constraints=axes+extra;expected=set()
            for chosen in combinations(constraints,4):
                p=solve([a for a,b in chosen],[b for a,b in chosen])
                if p is not None and all(clipping.dot(a,p)<=b for a,b in constraints):expected.add(p)
            self.assertEqual(set(clipping.vertices(box,extra)),expected)


if __name__=='__main__':unittest.main()
