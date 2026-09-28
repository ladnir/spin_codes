"""Exact vertices of a box clipped by verified supporting halfspaces.

Qhull only proposes facets. Every retained inequality is checked against
all source points with rational arithmetic. Missing facets merely weaken
the enclosure. Vertex enumeration and feasibility checks are exact.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import gcd,lcm

import numpy as np
from scipy.spatial import ConvexHull


def dot(a,b):return sum(x*y for x,y in zip(a,b))


def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def normalize(normal,bound):
    values=(*normal,bound);den=lcm(*(v.denominator for v in values))
    integers=[int(v*den) for v in values];divisor=gcd(*integers)
    return tuple(v//divisor for v in integers)


def supporting_planes(points):
    points=sorted(set(tuple(map(Q,p)) for p in points))
    if len(points)<4 or any(len(p)!=3 for p in points):raise ValueError('three-dimensional source points required')
    hull=ConvexHull(np.array(points,dtype=float));planes=set()
    for ids in hull.simplices:
        a,b,c=(points[i] for i in ids)
        normal=cross(tuple(y-x for x,y in zip(a,b)),tuple(y-x for x,y in zip(a,c)))
        if not any(normal):continue
        bound=dot(normal,a);sides=[dot(normal,p)-bound for p in points]
        if max(sides)<=0:planes.add(normalize(normal,bound))
        elif min(sides)>=0:planes.add(normalize(tuple(-v for v in normal),-bound))
        # A floating facet cutting through the source hull is not trusted.
    return tuple(sorted(planes))


class ClippedBox:
    def __init__(self,points):
        self.planes=supporting_planes(points)
        axes=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
        normals=axes+tuple(p[:3] for p in self.planes)
        self.bases=[];self.cache={}
        for indices in combinations(range(len(normals)),3):
            a,b,c=(normals[i] for i in indices);det=dot(a,cross(b,c))
            if not det:continue
            inverse=tuple(tuple(Q(x,det) for x in column) for column in (cross(b,c),cross(c,a),cross(a,b)))
            self.bases.append((indices,inverse))

    def vertices(self,cell):
        cell=tuple(map(Q,cell))
        if len(cell)!=6 or any(lo>hi for lo,hi in zip(cell[::2],cell[1::2])):
            raise ValueError('nonempty three-dimensional rectangle required')
        if cell in self.cache:return self.cache[cell]
        bounds=(cell[1],-cell[0],cell[3],-cell[2],cell[5],-cell[4],*(Q(p[3]) for p in self.planes))
        result=set()
        for indices,columns in self.bases:
            rhs=[bounds[i] for i in indices]
            point=tuple(sum(columns[j][k]*rhs[j] for j in range(3)) for k in range(3))
            if any(not lo<=x<=hi for x,lo,hi in zip(point,cell[::2],cell[1::2])):continue
            if all(dot(p[:3],point)<=p[3] for p in self.planes):result.add(point)
        result=tuple(sorted(result))
        if len(self.cache)>=4096:self.cache.clear()
        self.cache[cell]=result
        return result
