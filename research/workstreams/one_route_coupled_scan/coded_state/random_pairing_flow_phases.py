"""Exact population cells for the final matching-exposure phase.

Research proof geometry only. No files, no native benchmark. Closed real
cells overlap on boundaries; integer counts assign a boundary to its first
phase. The proof covers every integer population without double-counting.
"""
import math
import shutil
import sys
import time
from fractions import Fraction as F
from itertools import combinations, product
import random_pairing_flow_envelope_verify as v


def phase_inequalities(bounds,priority,phase,n=512):
    assert sorted(priority)==list(range(5)) and 0<=phase<5
    inequalities=[]
    for i,(lo,hi) in enumerate(bounds):
        unit=tuple(int(i==j) for j in range(5))
        inequalities.extend(((unit,hi),(tuple(-a for a in unit),-lo)))
    cap=n//2-1;assert n>=2 and n%2==0
    before=tuple(int(i in priority[:phase]) for i in range(5))
    after=tuple(int(i in priority[:phase+1]) for i in range(5))
    inequalities.extend(((before,cap),(tuple(-a for a in after),-cap)))
    return tuple(inequalities)


def solve_vertex(rows):
    matrix=[[F(x) for x in row] for row in rows]
    for col in range(5):
        pivot=next((i for i in range(col,5) if matrix[i][col]),None)
        if pivot is None:return None
        matrix[col],matrix[pivot]=matrix[pivot],matrix[col]
        scale=matrix[col][col];matrix[col]=[x/scale for x in matrix[col]]
        for i in range(5):
            if i!=col:
                scale=matrix[i][col]
                if scale:matrix[i]=[x-scale*y for x,y in zip(matrix[i],matrix[col])]
    return tuple(row[-1] for row in matrix)


def phase_vertices(bounds,priority,phase,n=512):
    inequalities=phase_inequalities(bounds,priority,phase,n);vertices=set()
    for active in combinations(inequalities,4):
        point=solve_vertex([(1,1,1,1,1,n)]+[tuple(a)+(b,) for a,b in active])
        if point is None or any(sum(x*y for x,y in zip(a,point))>b for a,b in inequalities):continue
        # Box plus nested prefix constraints is integral. Verify this
        # property explicitly rather than silently rounding a vertex.
        assert all(x.denominator==1 for x in point),'nonintegral phase vertex'
        vertices.add(tuple(int(x) for x in point))
    return tuple(sorted(vertices))


def phase_reference(vertices,inequalities,n=512):
    assert vertices
    total=[sum(x[i] for x in vertices) for i in range(5)];size=len(vertices)
    base=[x//size for x in total];remainder=n-sum(base)
    for raised in combinations(range(5),remainder):
        point=tuple(x+int(i in raised) for i,x in enumerate(base))
        if all(sum(x*y for x,y in zip(a,point))<=b for a,b in inequalities):return point
    # A vertex is always a valid integral fallback, never a rounded
    # infeasible center. No certificate depends on center optimality.
    return vertices[0]


def phase_geometry(bounds,priority,phase,n=512,star=False):
    vertices=phase_vertices(bounds,priority,phase,n)
    if not vertices:return (),(),None
    inequalities=phase_inequalities(bounds,priority,phase,n)
    normals=tuple(a for a,_ in inequalities if any(a))
    simplices=v.pulling_simplices(vertices,normals)
    reference=phase_reference(vertices,inequalities,n)
    if star and reference not in vertices:
        # Cone a relative-interior point over every boundary facet.
        # If rounding puts the reference on a facet, use the pulling cover.
        strict=all(sum(a*x for a,x in zip(normal,reference))<bound
                   for normal,bound in inequalities
                   if any(sum(a*x for a,x in zip(normal,point))<bound for point in vertices))
        if strict:
            dimension=v.exact_rank([[x-y for x,y in zip(point,vertices[0])] for point in vertices[1:]])
            facets=set()
            for normal,bound in inequalities:
                ids=tuple(i for i,point in enumerate(vertices) if sum(a*x for a,x in zip(normal,point))==bound)
                if ids and v.exact_rank([[x-y for x,y in zip(vertices[i],vertices[ids[0]])] for i in ids[1:]])==dimension-1:
                    facets.add(ids)
            assert facets
            simplices=tuple(sorted({tuple(sorted(tuple(ids[i] for i in s)+(len(vertices),)))
                                   for ids in facets for s in v.pulling_simplices(tuple(vertices[i] for i in ids),normals)}))
            vertices=vertices+(reference,)
    return vertices,simplices,reference


def sum_counts(bounds,n):
    values=[1]+[0]*n
    for lo,hi in bounds:
        prefix=[0]
        for value in values:prefix.append(prefix[-1]+value)
        values=[prefix[max(0,s-lo+1)]-prefix[max(0,s-hi)] for s in range(n+1)]
    return values


def phase_count(bounds,priority,phase,n=512):
    """Disjoint integer ownership: the first prefix reaching n/2-1."""
    cap=n//2-1
    if phase and cap==0:return 0
    before=sum_counts([bounds[i] for i in priority[:phase]],n)
    after=sum_counts([bounds[i] for i in priority[phase+1:]],n)
    lo,hi=bounds[priority[phase]];result=0
    for prefix in range(cap if phase else 1):
        for current in range(max(lo,cap-prefix),min(hi,n-prefix)+1):
            result+=before[prefix]*after[n-prefix-current]
    return result


def tests(star=False):
    checks=0
    for n in (2,4,6):
        for bounds in (((0,n),)*5,((0,n//2),)*4+((0,n),)):
            for priority in ((0,1,2,3,4),(4,2,0,3,1)):
                points=[p+(n-sum(p),) for p in product(range(n+1),repeat=4)
                        if sum(p)<=n and all(lo<=x<=hi for x,(lo,hi) in zip(p+(n-sum(p),),bounds))]
                counts=[]
                for phase in range(5):
                    inequalities=phase_inequalities(bounds,priority,phase,n)
                    vertices,simplices,reference=phase_geometry(bounds,priority,phase,n,star)
                    owned=[p for p in points if
                           (not phase or sum(p[i] for i in priority[:phase])<n//2-1)
                           and sum(p[i] for i in priority[:phase+1])>=n//2-1]
                    count=phase_count(bounds,priority,phase,n)
                    assert count==len(owned);counts.append(count);checks+=1
                    if not vertices:assert not owned;continue
                    assert sum(reference)==n and all(sum(a*b for a,b in zip(normal,reference))<=bound for normal,bound in inequalities)
                    for point in vertices:
                        assert point in points and all(sum(a*b for a,b in zip(normal,point))<=bound for normal,bound in inequalities);checks+=1
                    # Verify every exact point has a convex representation
                    # in one returned simplex using independent elimination.
                    for point in points:
                        expected=all(sum(a*b for a,b in zip(normal,point))<=bound for normal,bound in inequalities)
                        actual=any(simplex_contains([vertices[i] for i in s],point) for s in simplices)
                        assert actual==expected;checks+=1
                assert sum(counts)==len(points);checks+=1
    return checks


def simplex_contains(vertices,point):
    size=len(vertices);rows=[[F(v[i]) for v in vertices]+[F(point[i])] for i in range(5)]
    rows.append([F(1)]*(size+1))
    for col in range(size):
        pivot=next((i for i in range(col,len(rows)) if rows[i][col]),None)
        assert pivot is not None,'degenerate simplex'
        rows[col],rows[pivot]=rows[pivot],rows[col]
        scale=rows[col][col];rows[col]=[x/scale for x in rows[col]]
        for i in range(len(rows)):
            if i!=col:
                scale=rows[i][col]
                if scale:rows[i]=[x-scale*y for x,y in zip(rows[i],rows[col])]
    return all(not row[-1] for row in rows[size:]) and all(row[-1]>=0 for row in rows[:size])


def native_quota_tests():
    """Check the production cap255 affine quotas against greedy exposure."""
    box=((0,79),(80,159),(72,179),(0,287),(0,287));checks=0
    for priority in ((0,1,2,3,4),(4,3,2,1,0)):
        owned=0
        for phase in range(5):
            owned+=phase_count(box,priority,phase)
            vertices,_,reference=phase_geometry(box,priority,phase)
            if not vertices:continue
            polynomials,stable=v.phase_quotas(reference,[v.I()]*4,priority,final_phase=phase,vertices=vertices)
            assert stable
            for point in vertices+(reference,):
                delta=[v.I(point[i]-reference[i]) for i in (0,1,2,4)]
                remaining=255
                for color,poly in zip(priority,polynomials):
                    expected=min(remaining,point[color]);remaining-=expected
                    value=v.evaluate_polynomial(poly,delta)
                    assert value.lo<=expected<=value.hi and value.hi-value.lo<v.D('1e-25')
                    checks+=1
                assert remaining==0;checks+=1
        assert owned==sum_counts(box,512)[512]==155373250;checks+=1
    try:
        v.phase_quotas((512,0,0,0,0),[v.I()]*4,(0,1,2,3,4),final_phase=4,vertices=((512,0,0,0,0),))
    except AssertionError:checks+=1
    else:raise AssertionError('invalid exposure cell accepted')
    return checks


def rational_split_tests():
    """Independent barycentric coverage across a nonintegral edge split."""
    vertices=((0,0,0,0,2),(1,0,0,0,1),(0,1,0,0,1))
    split,mesh=v.bisect_simplex_edge(vertices,((0,1,2),),(0,1));checks=0
    assert split[-1]==(F(1,2),0,0,0,F(3,2));checks+=1
    for denominator in (2,3,5):
        for a in range(denominator+1):
            for b in range(denominator-a+1):
                weights=(F(a,denominator),F(b,denominator),F(denominator-a-b,denominator))
                point=tuple(sum(w*p[i] for w,p in zip(weights,vertices)) for i in range(5))
                assert any(simplex_contains([split[i] for i in s],point) for s in mesh);checks+=1
    # Reconstruct a mixed-sign quadratic using the rational mesh controls.
    polynomial={(2,0,0,0):v.I(2),(0,1,1,0):v.I(-3),(0,0,0,0):v.I(7)}
    controls=v.polar_coefficients(polynomial,[p[:4] for p in split],control_pairs=v.triangulation_pairs(mesh))
    for s in mesh:
        weights={i:F(1,len(s)) for i in s}
        point=[sum(weights[i]*split[i][a] for i in s) for a in range(4)]
        value=sum((controls[i,j]*v.I(weights[i]*weights[j]*(1 if i==j else 2)) for i in s for j in s if i<=j),v.I())
        truth=2*point[0]**2-3*point[1]*point[2]+7
        assert F(value.lo)<=truth<=F(value.hi);checks+=1
    return checks


def main():
    assert shutil.disk_usage('.').free>=20*(1<<30),'disk guard'
    probe=next((arg.split('=',1)[1] for arg in sys.argv if arg.startswith('--probe=')),None)
    if probe is not None:
        from random_pairing_flow_quadratic import RETUNED_LOG_WEIGHT_NUMERATORS
        assert probe in ('pulling','nested','star','integrated','split','point')
        precision=80 if '--precision80' in sys.argv else 40
        limit=float(next((arg.split('=',1)[1] for arg in sys.argv if arg.startswith('--deadline-seconds=')),'95'))
        assert 1<=limit<=120
        v.getcontext().prec=precision;start=time.monotonic();v.VERIFY_DEADLINE=start+limit
        kwargs=dict(matrix_mode='star' if probe=='star' else 'triangulated',control_degree=3,
                    log_weights=RETUNED_LOG_WEIGHT_NUMERATORS,bounded_switch=True,
                    exposure_priority=(0,1,2,3,4),exposure_phase=2)
        box=((0,79),(80,159),(72,179),(0,287),(0,287))
        if probe=='nested':kwargs['force_nested_phase']=True
        if probe in ('integrated','split'):kwargs['integrated_degree']=8
        if probe=='split':kwargs['simplex_splits']=((0,15),(10,15))
        if probe=='point':
            box=tuple((x,x) for x in (53,88,114,86,171));kwargs['exposure_phase']=None
            kwargs['center_override']=('96.88028197774','79.139755188407','95.896577195838','132.178071834463')
            if '--recenter' in sys.argv:
                h=tuple(lo for lo,_ in box);u=[x/10**9 for x in RETUNED_LOG_WEIGHT_NUMERATORS]
                _,proposal,_=v.variational(h,u,.1,.2,.08);center,_=v.moments(h,v.tables(proposal))
                kwargs['center_override']=tuple(str(round(x,12)) for x in center[:4])
        result=v.verify_region(box,**kwargs)
        print(dict(probe=probe,seconds=time.monotonic()-start,result=result),flush=True)
        assert time.monotonic()-start<limit
        return
    start=time.monotonic();print(dict(exact_phase_geometry_checks=tests(),star_geometry_checks=tests(True),native_quota_checks=native_quota_tests(),rational_split_checks=rational_split_tests()),flush=True)
    box=((0,79),(80,159),(72,179),(0,287),(0,287));priority=(0,1,2,3,4)
    for phase in range(5):
        count=phase_count(box,priority,phase)
        if count:
            vertices,simplices,reference=phase_geometry(box,priority,phase)
            print(dict(phase=phase,histograms=count,vertices=len(vertices),simplices=len(simplices),reference=reference),flush=True)
    print(dict(seconds=time.monotonic()-start),flush=True)
    assert time.monotonic()-start<55


if __name__=='__main__':main()
