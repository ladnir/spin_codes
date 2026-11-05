"""Exact box-clipped profile cells for the adaptive all-pair coverage lemma.

Interior coordinates are (r,a,b); even/odd face coordinates are (u,k).
Ranges are inclusive integer intervals, but vertices describe the complete
real polytope. Splits partition its integer lattice, not its real points.
Only the fixed domain inequalities and the supplied box are imposed.
This module has no FLOAT geometry, probability bounds, or coverage driver.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from functools import cached_property
from itertools import combinations


KINDS = ('interior', 'even_face', 'odd_face')


def _exact_point(point, dimension):
    try:
        result = tuple(point)
    except TypeError as error:
        raise ValueError('Expected an exact coordinate tuple') from error
    if (len(result) != dimension
            or any(type(value) not in (int, F) for value in result)):
        raise ValueError('Coordinates must be integers or Fractions; FLOAT is not accepted')
    return result


def _intersection(rows, dimension):
    """Return integer Cramer numerators and a positive denominator, or None.

    Dimension is fixed at two or three. Integer feasibility checks precede
    Fraction construction, avoiding rational elimination in each candidate.
    """
    if dimension == 2:
        (a, b), x = rows[0]
        (c, d), y = rows[1]
        denominator = a*d-b*c
        if not denominator:
            return None
        numerators = (x*d-b*y, a*y-x*c)
    else:
        (a, b, c), x = rows[0]
        (d, e, f), y = rows[1]
        (g, h, i), z = rows[2]
        bc = (e*i-f*h, f*g-d*i, d*h-e*g)
        denominator = a*bc[0]+b*bc[1]+c*bc[2]
        if not denominator:
            return None
        ca = (h*c-i*b, i*a-g*c, g*b-h*a)
        ab = (b*f-c*e, c*d-a*f, a*e-b*d)
        numerators = (x*bc[0]+y*ca[0]+z*ab[0],
                      x*bc[1]+y*ca[1]+z*ab[1],
                      x*bc[2]+y*ca[2]+z*ab[2])
    if denominator < 0:
        return tuple(-value for value in numerators), -denominator
    return numerators, denominator


@dataclass(frozen=True)
class Cell:
    """A fixed profile-domain intersection with an inclusive integer box.

    ``constraints`` is a tuple of (coefficients,rhs), each meaning
    dot(coefficients,point)<=rhs. ``contains`` tests integer cell coordinates,
    whereas ``contains_real`` accepts exact rational coordinates. Face points
    must be supplied as (u,k), not as their three-coordinate profile images.
    """

    kind: str
    ranges: tuple
    pairs: int = 4096

    def __post_init__(self):
        if self.kind not in KINDS or type(self.pairs) is not int or self.pairs < 1:
            raise ValueError('Expected a known cell kind and a positive integer pair count')
        try:
            ranges = tuple(tuple(row) for row in self.ranges)
        except TypeError as error:
            raise ValueError('Expected inclusive integer coordinate ranges') from error
        if (len(ranges) != self.dimension
                or any(len(row) != 2 or any(type(value) is not int for value in row)
                       or row[0] > row[1] for row in ranges)):
            raise ValueError('Expected one ordered inclusive integer range per coordinate')
        object.__setattr__(self, 'ranges', ranges)

    @property
    def dimension(self):
        return 3 if self.kind == 'interior' else 2

    @property
    def coordinate_names(self):
        return ('r', 'a', 'b') if self.kind == 'interior' else ('u', 'k')

    @cached_property
    def constraints(self):
        p = self.pairs
        if self.kind == 'interior':
            rows = [((1, 0, 0), 2*p-1), ((-1, 0, 0), -1),
                    ((0, 1, 0), p), ((0, -1, 0), 0),
                    ((0, 0, 1), p), ((0, 0, -1), 0),
                    ((-1, 2, 0), -1), ((-1, 0, 2), -1),
                    ((1, -1, 0), p), ((1, 0, -1), p),
                    ((-1, 1, 1), -1)]
        else:
            rows = [((1, 0), p-1), ((-1, 0), -1), ((0, -1), -1),
                    ((-1, 1), 0), ((1, 1), p)]
        for axis, (low, high) in enumerate(self.ranges):
            upper = tuple(int(index == axis) for index in range(self.dimension))
            rows.append((upper, high))
            rows.append((tuple(-value for value in upper), -low))
        # Parallel inequalities with the same normal need only their tightest rhs.
        tight = {}
        for coefficients, rhs in rows:
            tight[coefficients] = min(rhs, tight.get(coefficients, rhs))
        return tuple(sorted(tight.items()))

    def contains_real(self, point):
        point = _exact_point(point, self.dimension)
        return all(sum(a*x for a, x in zip(coefficients, point)) <= rhs
                   for coefficients, rhs in self.constraints)

    def contains(self, point):
        point = _exact_point(point, self.dimension)
        return all(F(value).denominator == 1 for value in point) and self.contains_real(point)

    def profile(self, point):
        """Affine map to (r,a,b), without asserting that point lies in this cell."""
        point = _exact_point(point, self.dimension)
        if self.kind == 'interior':
            return point
        u, k = point
        return (2*u, u, u-k) if self.kind == 'even_face' else (2*u, u-k, u)

    @cached_property
    def _vertices(self):
        vertices = set()
        dimension, constraints = self.dimension, self.constraints
        for rows in combinations(constraints, dimension):
            intersection = _intersection(rows, dimension)
            if intersection is None:
                continue
            numerators, denominator = intersection
            if all(sum(a*x for a, x in zip(coefficients, numerators)) <= rhs*denominator
                   for coefficients, rhs in constraints):
                vertices.add(tuple(F(value, denominator) for value in numerators))
        return tuple(sorted(vertices))

    def vertices(self):
        """Complete exact vertices, including segments/points; () means real-empty.

        The finite box makes the polytope bounded. At each vertex, active
        constraint normals span the ambient space, even in lower dimension.
        Enumerating every independent ambient-size subset therefore finds all
        vertices. A nonempty bounded polytope has a vertex.
        """
        return self._vertices

    def real_empty(self):
        return not self.vertices()

    def integer_point(self):
        """Return an integer witness or None, using exact native-domain bounds."""
        p = self.pairs
        if self.kind == 'interior':
            (rlow, rhigh), (alow, ahigh), (blow, bhigh) = self.ranges
            alow, ahigh = max(0, alow), min(p, ahigh)
            blow, bhigh = max(0, blow), min(p, bhigh)
            rlow = max(1, rlow, 2*alow+1, 2*blow+1)
            rhigh = min(2*p-1, rhigh, p+ahigh, p+bhigh)
            if alow > ahigh or blow > bhigh or rlow > rhigh:
                return None
            return rlow, max(alow, rlow-p), max(blow, rlow-p)
        (ulow, uhigh), (klow, khigh) = self.ranges
        klow = max(1, klow)
        ulow, uhigh = max(1, ulow, klow), min(p-1, uhigh, p-klow)
        return None if klow > khigh or ulow > uhigh else (ulow, klow)

    def integer_empty(self):
        return self.integer_point() is None

    def split(self, axis, integer_cut):
        """Two disjoint integer children: coordinate<=cut and coordinate>=cut+1.

        The cut must lie strictly before the axis range's last integer. Either
        child may be empty after intersecting the native profile domain.
        """
        if (type(axis) is not int or not 0 <= axis < self.dimension
                or type(integer_cut) is not int):
            raise ValueError('Expected a valid axis index and an integer cut')
        low, high = self.ranges[axis]
        if not low <= integer_cut < high:
            raise ValueError('Cut must give two nonempty integer box ranges')
        left, right = list(self.ranges), list(self.ranges)
        left[axis], right[axis] = (low, integer_cut), (integer_cut+1, high)
        return Cell(self.kind, tuple(left), self.pairs), Cell(self.kind, tuple(right), self.pairs)

    def histogram_bounds(self, parity, *, integer=True):
        """Bounds ((single_min,single_max),(full_min,full_max)), or None if empty.

        Integer bounds are ceil(minimum),floor(maximum) of exact real extrema,
        enclosing every attainable integer histogram. They need not be jointly
        attainable. With integer=False the endpoints remain Fractions.
        """
        if type(parity) is not int or parity not in (0, 1) or type(integer) is not bool:
            raise ValueError('Expected parity0/1 and a boolean integer-envelope flag')
        vertices = self.vertices()
        if not vertices:
            return None
        profiles = tuple(self.profile(point) for point in vertices)
        singles = tuple(profile[0]-2*profile[parity+1] for profile in profiles)
        full = tuple(profile[parity+1] for profile in profiles)
        bounds = ((min(singles), max(singles)), (min(full), max(full)))
        if not integer:
            return bounds
        return tuple((-((-low.numerator)//low.denominator), high.numerator//high.denominator)
                     for low, high in bounds)

    def descriptor(self):
        """JSON-safe constructor data, with no numerical or coverage claim."""
        return dict(kind=self.kind, ranges=[list(row) for row in self.ranges], pairs=self.pairs)
