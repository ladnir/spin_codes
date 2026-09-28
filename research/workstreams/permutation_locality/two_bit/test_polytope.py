from fractions import Fraction as Q
from itertools import product
import unittest

from polytope import ClippedBox,dot


class ClippedBoxTests(unittest.TestCase):
    def setUp(self):
        self.points=[(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
        self.domain=ClippedBox(self.points)

    def test_tetrahedron_and_empty_intersection(self):
        self.assertEqual(set(self.domain.vertices((0,1)*3)),set(self.points))
        self.assertEqual(self.domain.vertices((Q(1,2),1)*3),())
        self.assertEqual(len(self.domain.vertices((0,Q(1,2))*3)),7)
        for p in self.points:
            for plane in self.domain.planes:self.assertLessEqual(dot(plane[:3],p),plane[3])

    def test_degenerate_cells_and_boundary(self):
        p=(Q(1,2),Q(1,3),Q(1,6));cell=tuple(y for x in p for y in (x,x))
        self.assertEqual(self.domain.vertices(cell),(p,))
        self.assertEqual(set(self.domain.vertices((0,0,0,0,0,1))),{(0,0,0),(0,0,1)})
        self.assertEqual(self.domain.vertices((Q(1,2),Q(1,2))*3),())

    def test_small_feasible_grid_is_in_every_vertex_support_bound(self):
        cell=(Q(1,8),Q(5,8),Q(0),Q(3,4),Q(0),Q(1,2))
        vertices=self.domain.vertices(cell)
        for p in product([Q(k,8) for k in range(9)],repeat=3):
            if sum(p)>1 or any(not lo<=x<=hi for x,lo,hi in zip(p,cell[::2],cell[1::2])):continue
            for direction in ((1,1,1),(3,-4,7),(-2,5,-8)):
                self.assertLessEqual(dot(direction,p),max(dot(direction,v) for v in vertices))


if __name__=='__main__':unittest.main()
