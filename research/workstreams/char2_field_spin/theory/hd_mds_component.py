"""Exact MDS-aware component transfers, scoped to bounded source classes.

Rows are EXACT active source rows, unlike permitted rows in the older
Bernoulli expansion. No scalable all-support certificate is claimed.
"""

from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import product

from hd_affine_proof import row_probability
from hd_component_rank import exact_kernels


def scalar(value, reference):
    """Embed a rational in the kernel's number type, including directed Arb."""
    value = F(value)
    if isinstance(reference, (int,F)):
        return value
    return (reference*0+value.numerator)/value.denominator


def cut_choices(old,new,q):
    """(nonzero difference, multiplicity) for one new free component value."""
    if old not in (0,1) or new not in (0,1) or q < 2:
        raise ValueError("invalid live flags or alphabet")
    if not old and not new:
        return ((0,1),)
    if not old:
        return ((1,q-1),)
    if not new:
        return ((1,1),)
    return ((0,1),(1,q-2))


def region_polynomial(cuts,rows,q,kernel):
    """Positive polynomial in one difference-degree variable per source row.

    cuts lists (row label, position); labels must be exactly range(a).
    Coefficient mask increments those labelled row degrees once. Its matrix
    acts on the carried component live flag and post-state support class.
    """
    cuts = tuple(sorted(cuts,key=lambda item:item[1]))
    active = len(cuts)
    if (not cuts or sorted(r for r,_ in cuts) != list(range(active))
            or len({j for _,j in cuts}) != active
            or not 0 <= cuts[0][1] <= cuts[-1][1] < rows):
        raise ValueError("distinct valid cuts with consecutive row labels required")
    initial = kernel(0)
    post,zero = len(initial),initial[0][0]*0
    result = {}
    for incoming in (0,1):
        # current flag, prefix occupancy, difference mask -> multiplicity
        frontier = {(incoming,cuts[0][1]*incoming,0):1}
        for index,(row,start) in enumerate(cuts):
            end = cuts[index+1][1] if index+1 < active else rows
            following = defaultdict(int)
            for (old,occupancy,mask),multiplicity in frontier.items():
                for new in (0,1):
                    for delta,count in cut_choices(old,new,q):
                        following[new,occupancy+(end-start)*new,mask|(delta<<row)] += multiplicity*count
            frontier = following
        for (outgoing,occupancy,mask),multiplicity in frontier.items():
            if not multiplicity:
                continue
            matrix = result.setdefault(mask,[[zero]*(2*post) for _ in range(2*post)])
            local = kernel(occupancy)
            for old in range(post):
                for new in range(post):
                    matrix[incoming*post+old][outgoing*post+new] += multiplicity*local[old][new]
    return result


def degree_distribution(cuts_by_region,rows,q,kernel,cap=None):
    """Exact coefficient vectors; optional saturation preserves threshold tests."""
    active = len(cuts_by_region[0])
    initial = kernel(0)
    post,zero = len(initial),initial[0][0]*0
    start = (zero+1,)+(zero,)*(2*post-1)
    values = {(0,)*active:start}
    for cuts in cuts_by_region:
        polynomial = region_polynomial(cuts,rows,q,kernel)
        following = {}
        for degrees,vector in values.items():
            for mask,matrix in polynomial.items():
                degree = tuple(d+((mask>>r)&1) for r,d in enumerate(degrees))
                if cap is not None:
                    degree = tuple(min(cap,d) for d in degree)
                target = following.setdefault(degree,[zero]*(2*post))
                for old,weight in enumerate(vector):
                    if weight:
                        for new in range(2*post):
                            target[new] += weight*matrix[old][new]
        values = following
    return {degrees:sum(vector) for degrees,vector in values.items()}


def source_moment(cuts_by_region,rows,q,dimension,kernel,feasible_only=False):
    """Nonprojective exact-active source moment (or its feasible z envelope)."""
    length,active = len(cuts_by_region),len(cuts_by_region[0])
    if not 1 <= dimension < length:
        raise ValueError("invalid outer dimension")
    minimum = length-dimension+1
    distribution = degree_distribution(cuts_by_region,rows,q,kernel,
                                       cap=minimum if feasible_only else None)
    zero = kernel(0)[0][0]*0
    if feasible_only:
        total = sum((mass for degrees,mass in distribution.items()
                     if all(d==minimum for d in degrees)),zero)
        return scalar(F(1,(q-1)**((length-dimension)*active)),zero)*total
    result = zero
    for degrees,mass in distribution.items():
        if any(d < minimum for d in degrees):
            continue
        for degree in degrees:
            mass *= scalar(row_probability(length,dimension,degree,q),zero)
        result += mass
    return result


def marked_upper(cuts_by_region,rows,q,dimension,kernel,degree_mark):
    """Positive total minimum-degree upper, not a per-row exact filter.

    Drop the all-zero assignment before multiplying by lambda^(-a*d).
    Exact-active MDS feasibility implies total difference degree >= a*d.
    """
    degree_mark = F(degree_mark)
    if degree_mark < 1:
        raise ValueError("minimum-degree mark must be at least one")
    length,active = len(cuts_by_region),len(cuts_by_region[0])
    minimum = length-dimension+1
    initial = kernel(0)
    size,zero = 2*len(initial),initial[0][0]*0
    vector = [zero+1]+[zero]*(size-1)
    for cuts in cuts_by_region:
        matrix = [[zero]*size for _ in range(size)]
        for mask,coefficient in region_polynomial(cuts,rows,q,kernel).items():
            factor = scalar(degree_mark**mask.bit_count(),zero)
            for i in range(size):
                for j in range(size):
                    matrix[i][j] += factor*coefficient[i][j]
        vector = [sum(vector[i]*matrix[i][j] for i in range(size)) for j in range(size)]
    prefactor = F(1,(q-1)**((length-dimension)*active))/degree_mark**(active*minimum)
    return (sum(vector)-1)*scalar(prefactor,zero)


def direct_value_moment(cuts_by_region,rows,q,dimension,kernels):
    """GF4 exhaustive component-value oracle; XOR is field addition."""
    if q != 4:
        raise ValueError("this independent bounded value oracle is GF4 only")
    length,active = len(cuts_by_region),len(cuts_by_region[0])
    @lru_cache(None)
    def post_moment(occupancies):
        vector = [F(1)]+[F(0)]*(len(kernels[0])-1)
        for occupancy in occupancies:
            matrix = kernels[occupancy]
            vector = [sum(vector[i]*matrix[i][j] for i in range(len(vector)))
                      for j in range(len(vector))]
        return sum(vector)
    result = F(0)
    for components in product(range(q),repeat=length*active):
        degrees = [0]*active
        carried,index = 0,0
        occupancies = []
        for cuts in cuts_by_region:
            cuts = sorted(cuts,key=lambda item:item[1])
            occupancy = cuts[0][1]*bool(carried)
            for j,(row,position) in enumerate(cuts):
                new = components[index]
                index += 1
                degrees[row] += (new != carried)
                end = cuts[j+1][1] if j+1 < active else rows
                occupancy += (end-position)*bool(new)
                carried = new
            occupancies.append(occupancy)
        if not all(degrees):
            continue
        probability = F(1)
        for degree in degrees:
            probability *= row_probability(length,dimension,degree,q)
        if probability:
            result += probability*post_moment(tuple(occupancies))
    return result


def bounded_receipts():
    from hd_component_rank import component_moment
    cases = (
        ("one row; common cuts",4,2,[[(0,0)]]*4),
        ("one row; moving cuts",4,2,[[(0,j)] for j in (0,1,2,3)]),
        ("two rows; reversed middle order",3,1,
         [[(0,0),(1,1)],[(1,0),(0,2)],[(0,1),(1,2)]]),
    )
    for name,rows,dimension,cuts in cases:
        q,mark = 4,F(1,4)
        kernels = exact_kernels(rows,q,mark)
        exact = source_moment(cuts,rows,q,dimension,kernels.__getitem__)
        feasible = source_moment(cuts,rows,q,dimension,kernels.__getitem__,True)
        active,length = len(cuts[0]),len(cuts)
        positions = [[j for _,j in sorted(c,key=lambda item:item[1])] for c in cuts]
        old = component_moment(positions,rows,q,kernels.__getitem__)/(q-1)**((length-dimension)*active)
        marks = {str(lam):str(marked_upper(cuts,rows,q,dimension,kernels.__getitem__,lam))
                 for lam in (1,2,4)}
        yield dict(name=name,q=q,mark=str(mark),exact_source_moment=str(exact),
                   feasible_source_upper=str(feasible),unfiltered_source_upper=str(old),
                   feasibility_ratio=float(feasible/old),exact_ratio=float(exact/old),
                   degree_mark_uppers=marks,whole_certificate=False)


if __name__ == '__main__':
    import json
    for receipt in bounded_receipts():
        print(json.dumps(receipt),flush=True)
