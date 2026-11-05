"""Bounded, exactly accounted coverage diagnostics for the full histogram simplex.

No search tree or raw output is saved. A partial run reports uncovered mass;
neither complete accounting nor selected passing regions mean full coverage.
"""
import math
import shutil
import sys
import time
from decimal import getcontext
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from random_pairing_flow_envelope_verify import (
    I,D,TYPES,dot,tighten_region,region_count,region_reference,verify_region,
    polynomial_coefficients,evaluate_polynomial,local_data,LOG_WEIGHT_NUMERATORS,simplex_polynomial,poly_multiply,
    interval_data,variational,moments,tables,
    coupled_radial_bound,bernstein_polynomial,poly_add,
    bernstein_coefficients,polar_coefficients,bernstein_gaussian_shell,region_vertices,
    weak_compositions,global_covariance_polynomials,OUTPUT,
    polar_square_controls,
    bernstein_coefficients_reference,
    polar_coefficients_reference,
)
from random_pairing_flow_envelope import eigensystem,solve
import random_pairing_flow_envelope_verify as verifier
from random_pairing_flow_quadratic import RETUNED_LOG_WEIGHT_NUMERATORS


def split_region(bounds,n=512):
    bounds=tighten_region(bounds,n);assert bounds is not None
    axis=max(range(5),key=lambda a:(bounds[a][1]-bounds[a][0],-a))
    lo,hi=bounds[axis];assert lo<hi
    middle=(lo+hi)//2;children=[]
    for interval in ((lo,middle),(middle+1,hi)):
        child=list(bounds);child[axis]=interval;child=tighten_region(tuple(child),n)
        if child is not None:children.append(child)
    assert sum(region_count(child,n) for child in children)==region_count(bounds,n)
    return children


def subtract_region(bounds,removed,n=512):
    """Disjoint integer slabs for bounds minus its intersection with removed."""
    bounds=tighten_region(bounds,n);assert bounds is not None
    intersection=tighten_region(tuple((max(a,c),min(b,d)) for (a,b),(c,d) in zip(bounds,removed)),n)
    if intersection is None:return [bounds]
    remaining=list(bounds);pieces=[]
    for axis,(lo,hi) in enumerate(intersection):
        left,right=remaining[axis]
        for interval in ((left,lo-1),(hi+1,right)):
            candidate=remaining.copy();candidate[axis]=interval;candidate=tighten_region(tuple(candidate),n)
            if candidate is not None:pieces.append(candidate)
        remaining[axis]=(lo,hi)
    assert sum(region_count(piece,n) for piece in pieces)+region_count(intersection,n)==region_count(bounds,n)
    return pieces


def interior_partition():
    caps,_=corner_partition('local')
    seeds=[('cap'+str(i),bounds,'radial',2) for i,bounds in enumerate(caps)]
    seeds+=[('middle_band',((80,136),)*5,'polytope',3),
            ('boundary_r4',tighten_region(((181,189),(184,192),(79,87),(0,512),(0,0))),'bernstein',4)]
    pending=[((0,512),)*5]
    for i,(_,bounds,_,_) in enumerate(seeds):
        for _,previous,_,_ in seeds[:i]:
            overlap=tuple((max(a,c),min(b,d)) for (a,b),(c,d) in zip(bounds,previous))
            assert tighten_region(overlap) is None
        pending=[piece for region in pending for piece in subtract_region(region,bounds)]
    covered=sum(region_count(bounds) for _,bounds,_,_ in seeds)
    assert covered==241866501
    assert covered+sum(region_count(region) for region in pending)==math.comb(516,4)
    return seeds,pending


def partition_tests():
    checks=0
    for n in range(6):
        points=[prefix+(n-sum(prefix),) for prefix in product(range(n+1),repeat=4) if sum(prefix)<=n]
        for index in range(12):
            bounds=tuple((0,n if (index+a)%3 else n//2) for a in range(5));bounds=tighten_region(bounds,n)
            if bounds is None:continue
            removed=tuple((min((index+a)%(n+1),(index+2*a)%(n+1)),
                           max((index+a)%(n+1),(index+2*a)%(n+1))) for a in range(5))
            pieces=subtract_region(bounds,removed,n)
            for point in points:
                inside=lambda region:all(a<=x<=b for x,(a,b) in zip(point,region))
                assert sum(inside(piece) for piece in pieces)==int(inside(bounds) and not inside(removed))
            checks+=1
    seeds,pending=interior_partition()
    assert sum(region_count(bounds) for bounds in pending)==2677752444;checks+=1
    return checks


def verify_interior_partition():
    seeds,pending=interior_partition();covered=0
    for name,bounds,mode,degree in seeds:
        assert time.monotonic()<verifier.VERIFY_DEADLINE
        if mode=='radial':result=verify_region(bounds,coarse_only=True,coupled_coarse=True,local_radial='auto')
        else:result=verify_region(bounds,matrix_mode=mode,control_degree=degree)
        assert result['passes'],(name,result['directed_log_upper'])
        covered+=region_count(bounds)
        print(dict(region=name,histograms=region_count(bounds),directed_log_upper=result['directed_log_upper'],passes=True),flush=True)
    outstanding=sum(region_count(bounds) for bounds in pending)
    assert covered+outstanding==math.comb(516,4)
    return dict(certified_histograms=covered,pending_histograms=outstanding,pending_regions=len(pending),
                complete_histogram_cover=False,whole_code_certificate=False)


def retuned_partition():
    """Compact source checkpoint, not a saved search tree or assumed proof."""
    caps,central=corner_partition('retuned-wide')
    seeds=[('wide_cap'+str(i),bounds,'polytope',degree)
           for i,(bounds,degree) in enumerate(zip(caps,(2,2,4,3,3)))]
    seeds.extend(('interior'+str(i),bounds,'triangulated',2) for i,bounds in enumerate((
        ((0,79),(80,159),(0,71),(72,143),(144,215)),
        ((0,79),(80,159),(0,71),(72,143),(216,287)),
        ((0,79),(80,159),(0,71),(0,71),(216,287)),
        ((0,79),(80,159),(0,71),(0,71),(144,215)),
        ((0,79),(80,159),(0,71),(72,143),(60,143)),
        ((0,79),(80,159),(0,71),(216,287),(0,71)),
        ((0,79),(80,159),(0,71),(216,287),(72,143)),
        ((0,79),(80,159),(0,71),(144,215),(72,143)),
        ((0,79),(80,159),(0,71),(144,215),(0,71)),
        ((0,79),(80,159),(0,71),(144,215),(144,215)),
        ((0,72),(80,152),(0,71),(144,215),(216,287)),
        ((0,72),(80,152),(0,71),(216,287),(144,216)),
        ((68,79),(148,159),(60,71),(60,71),(132,143)),
        ((0,79),(80,159),(180,233),(0,71),(72,143)),
    )))
    seeds.extend((name,bounds,'bounded_radial',2) for name,bounds in (
        ('bounded_radial0',((0,79),(80,159),(180,233),(0,71),(144,252))),
        ('bounded_radial1',((0,79),(80,159),(180,233),(0,71),(0,71))),
    ))
    seeds.append(('mixed_degree4',((0,79),(80,159),(180,233),(72,252),(0,180)),
                  'triangulated_ball',4))
    seeds.append(('high_middle_degree3',((0,79),(80,159),(234,287),(0,198),(0,198)),
                  'triangulated_ball',3))
    pending=[central]
    for _,bounds,_,_ in seeds[5:]:
        pending=[piece for region in pending for piece in subtract_region(region,bounds)]
    for i,(_,bounds,_,_) in enumerate(seeds):
        for _,other,_,_ in seeds[:i]:
            assert tighten_region(tuple((max(a,c),min(b,d)) for (a,b),(c,d) in zip(bounds,other))) is None
    covered=sum(region_count(bounds) for _,bounds,_,_ in seeds)
    assert covered==645724890
    assert covered+sum(region_count(bounds) for bounds in pending)==math.comb(516,4)
    return seeds,pending


def verify_retuned_partition(start=0,stop=None):
    seeds,pending=retuned_partition();worst=None;covered=0
    stop=len(seeds) if stop is None else stop
    assert 0<=start<stop<=len(seeds)
    for name,bounds,mode,degree in seeds[start:stop]:
        radial=mode=='bounded_radial';mixed=mode=='triangulated_ball'
        result=verify_region(bounds,matrix_mode='interval' if radial else 'triangulated' if mixed else mode,
                             control_degree=degree,radial_controls=not mixed,
                             coarse_only=radial,bounded_switch=radial or mixed,
                             coupled_coarse=True,local_radial='shells32' if radial else 'auto',
                             log_weights=RETUNED_LOG_WEIGHT_NUMERATORS)
        assert result['passes'],(name,result['directed_log_upper'])
        covered+=region_count(bounds)
        bound=D(result['directed_log_upper']);worst=bound if worst is None else max(worst,bound)
        print(dict(replayed_region=name,histograms=region_count(bounds),log_upper=str(bound)),flush=True)
    return dict(certified_histograms=covered,checkpoint_histograms=645724890,
                replayed_region_range=(start,stop),full_checkpoint_replayed=start==0 and stop==len(seeds),
                pending_histograms=sum(region_count(x) for x in pending),
                pending_regions=len(pending),worst_log_upper=str(worst),
                complete_histogram_cover=False,whole_code_certificate=False)


def pending_cover(max_nodes=32,max_seconds=35):
    """Largest-region-first refinement of the explicitly uncovered domain.

    Existing seeds are excluded, not silently recertified. Every new passing
    region is checked in this call. The caller can replay the old seeds too.
    """
    _,pending=interior_partition();initial=sum(region_count(b) for b in pending)
    start=time.monotonic();covered=0;passed=[];failures=[];nodes=0
    while pending and nodes<max_nodes and time.monotonic()-start<max_seconds:
        pending.sort(key=region_count);bounds=pending.pop();count=region_count(bounds)
        nodes+=1
        result=verify_region(bounds,coarse_only=True,coupled_coarse=True,local_radial='auto')
        if result['passes']:
            covered+=count;passed.append(bounds)
            print(dict(new_region=bounds,histograms=count,log_upper=result['directed_log_upper']),flush=True)
        elif count>1:
            failures.append((count,bounds,result['directed_log_upper']))
            pending.extend(split_region(bounds))
        else:
            pending.append(bounds);break
    outstanding=sum(region_count(b) for b in pending)
    assert covered+outstanding==initial
    print(dict(examined_regions=nodes,newly_certified=covered,previous_seeds_not_replayed=math.comb(516,4)-initial,
               pending_histograms=outstanding,pending_regions=len(pending),seconds=time.monotonic()-start,
               complete_histogram_cover=False,whole_code_certificate=False),flush=True)
    print(dict(largest_failed_regions=sorted(failures,reverse=True)[:3]),flush=True)
    return passed,pending


def retuned_cover_pilot(max_nodes=24,max_seconds=38,hybrid=False,triangulated=False,wide_caps=False,matrix_width=96,resume=False):
    """Exactly accounted pilot from the complete domain, one fixed potential.

    Recheck the five caps; refine their disjoint complement. No old-weight
    seed is imported. This bounded run does not persist a search tree.
    """
    start=time.monotonic();weights=RETUNED_LOG_WEIGHT_NUMERATORS
    hybrid=hybrid or triangulated
    stop_at=min(verifier.VERIFY_DEADLINE-2,start+max_seconds)
    caps,central=corner_partition('retuned-wide' if wide_caps else 'local')
    certified=0;regions=0;attempts=0;worst_upper=None
    previous=0;pending=[central]
    if resume:
        seeds,pending=retuned_partition();previous=sum(region_count(bounds) for _,bounds,_,_ in seeds)
        certified=previous;caps=[]
    for axis,bounds in enumerate(caps):
        degree=(2,2,4,3,3)[axis] if wide_caps else 2
        result=verify_region(bounds,coarse_only=not wide_caps,coupled_coarse=True,local_radial='auto',log_weights=weights,
                             matrix_mode='polytope',control_degree=degree,radial_controls=wide_caps)
        assert result['passes'];certified+=region_count(bounds);regions+=1
        upper=D(result['directed_log_upper']);worst_upper=upper if worst_upper is None else max(worst_upper,upper)
        if wide_caps:
            print(dict(rechecked_wide_cap=axis,histograms=region_count(bounds),degree=degree,
                       log_upper=str(upper),method=result['proof_method']),flush=True)
    seed_count=certified;failed_singletons=[];matrix_attempts=0;inconclusive={};methods={};failed=[];failed_matrices=[]
    while pending and attempts<max_nodes and time.monotonic()<stop_at:
        bounds=pending.pop();count=region_count(bounds);attempts+=1
        degree=3 if hybrid and not triangulated and len(region_vertices(bounds))<=10 else 2
        full_deadline=verifier.VERIFY_DEADLINE;attempt_start=time.monotonic()
        verifier.VERIFY_DEADLINE=min(stop_at,attempt_start+(20 if triangulated else 12)) if hybrid else stop_at
        try:
            result=verify_region(bounds,coarse_only=not hybrid,coupled_coarse=True,local_radial='auto',
                                 log_weights=weights,radial_controls=True,
                                 matrix_mode='triangulated' if triangulated else 'polytope',control_degree=degree,
                                 matrix_gate=(matrix_width,'14') if triangulated else (80,'9.5') if hybrid else None)
            if not result['passes'] and max(hi-lo for lo,hi in bounds)<=32:
                result=verify_region(bounds,coarse_only=True,coupled_coarse=True,local_radial='shells',
                                     log_weights=weights,radial_controls=True)
        except AssertionError as error:
            reason=str(error)
            if reason!='region verification deadline' and not reason.startswith("('nonpositive pivot'"):raise
            reason=reason.split(',')[0];inconclusive[reason]=inconclusive.get(reason,0)+1
            result={'passes':False,'proof_method':'inconclusive'}
            if time.monotonic()>=stop_at:
                pending.append(bounds);break
        finally:verifier.VERIFY_DEADLINE=full_deadline
        matrix_attempts+=int(result['proof_method']=='local_quadratic')
        if result['proof_method']=='local_quadratic' and not result['passes']:
            failed_matrices.append((count,bounds,result['directed_log_upper']))
            failed_matrices=sorted(failed_matrices,reverse=True)[:3]
        method=result['proof_method'];methods[method]=methods.get(method,0)+1
        if result['passes']:
            certified+=count;regions+=1
            bound=D(result['directed_log_upper']);worst_upper=bound if worst_upper is None else max(worst_upper,bound)
            print(dict(new_certified_region=bounds,histograms=count,log_upper=result['directed_log_upper'],
                       method=result['proof_method'],seconds=time.monotonic()-attempt_start),flush=True)
        elif count==1:failed_singletons.append(bounds)
        else:
            if 'directed_log_upper' in result:
                failed.append((count,bounds,result['directed_log_upper'],method))
                failed=sorted(failed,reverse=True)[:3]
            children=split_region(bounds)
            pending.extend(sorted(children,key=region_count))
        assert certified+sum(region_count(x) for x in pending)+len(failed_singletons)==math.comb(516,4)
    assert certified+sum(region_count(x) for x in pending)+len(failed_singletons)==math.comb(516,4)
    print(dict(proof_weights=weights,attempts=attempts,regions_passed=regions,
               certified_this_run=certified-previous,previous_certified_not_replayed=previous,
               total_including_previous_certificates=certified,new_beyond_seeds=certified-seed_count,
               pending_histograms=sum(region_count(x) for x in pending),pending_regions=len(pending),
               hybrid=hybrid,triangulated=triangulated,wide_caps=wide_caps,matrix_width=matrix_width,
               matrix_attempts=matrix_attempts,inconclusive=inconclusive,methods=methods,
               partial_worst_log_upper_this_run=None if worst_upper is None else str(worst_upper),
               failed_singletons=failed_singletons,seconds=time.monotonic()-start,
               complete_histogram_cover=False,whole_code_certificate=False),flush=True)
    print(dict(largest_unproved_bounds=failed,largest_unproved_matrix_bounds=failed_matrices,
               largest_pending=sorted(((region_count(x),x) for x in pending),reverse=True)[:3]),flush=True)
    return pending


def domain_tests():
    checks=0
    for n in range(7):
        root=((0,n),)*5;assert region_count(root,n)==math.comb(n+4,4);checks+=1
        points=[prefix+(n-sum(prefix),) for prefix in product(range(n+1),repeat=4) if sum(prefix)<=n]
        for index in range(25):
            bounds=tuple((min((index+2*a)%(n+1),(3*index+a)%(n+1)),
                          max((index+2*a)%(n+1),(3*index+a)%(n+1))) for a in range(5))
            expected=sum(all(lo<=x<=hi for x,(lo,hi) in zip(point,bounds)) for point in points)
            assert region_count(bounds,n)==expected;checks+=1
            tightened=tighten_region(bounds,n)
            if expected:
                reference=region_reference(tightened,n)
                assert reference in points and all(lo<=x<=hi for x,(lo,hi) in zip(reference,bounds));checks+=1
                if expected>1:
                    children=split_region(tightened,n)
                    for point in points:
                        actual=sum(all(lo<=x<=hi for x,(lo,hi) in zip(point,child)) for child in children)
                        assert actual==int(all(lo<=x<=hi for x,(lo,hi) in zip(point,bounds)))
                    checks+=1
        pending=[root];leaves=set()
        while pending:
            bounds=pending.pop()
            if region_count(bounds,n)==1:
                point=region_reference(bounds,n);assert point not in leaves;leaves.add(point)
            else:pending.extend(split_region(bounds,n))
        assert leaves==set(points);checks+=1
    h=(100,110,120,90,92);scores={key:F(2*key[0]**2-key[1]+3,17) for key in TYPES}
    def exact(v):
        return sum(F(v[a]*(v[a]-1),2*511)*value if a==b else F(v[a]*v[b],511)*value for (a,b),value in scores.items())
    def enclose(x):return I(x.numerator)/x.denominator
    for absorber in range(5):
        axes=tuple(a for a in range(5) if a!=absorber)
        coefficients=polynomial_coefficients(enclose(exact(h)),{key:enclose(x) for key,x in scores.items()},h,axes=axes)
        for delta in product((-1,0,1),repeat=4):
            v=list(h)
            for a,change in zip(axes,delta):v[a]+=change;v[absorber]-=change
            result=evaluate_polynomial(coefficients,[I(x) for x in delta]);truth=exact(v)
            assert F(result.lo)<=truth<=F(result.hi);checks+=1
    linear={tuple(int(i==j) for i in range(4)):I(1) for j in range(4)}
    for polynomial in (linear,poly_multiply(linear,linear)):
        bound=simplex_polynomial(polynomial,[I(0,1)]*4)
        assert bound.lo<=0 and 1<=bound.hi<D('1.000001');checks+=1
    coefficients={powers:F((-1)**sum(powers)*(1+sum((i+1)*k for i,k in enumerate(powers))),13)
                  for powers in product(range(5),repeat=4) if sum(powers)<=4}
    lower=(-2,-1,0,1);box=[I(x,x+2) for x in lower]
    envelope=simplex_polynomial({key:enclose(x) for key,x in coefficients.items()},box)
    for increments in product(range(3),repeat=4):
        if sum(increments)>2:continue
        point=tuple(x+y for x,y in zip(lower,increments))
        truth=sum(value*math.prod(x**k for x,k in zip(point,powers)) for powers,value in coefficients.items())
        assert F(envelope.lo)<=truth<=F(envelope.hi);checks+=1
    caps,central=corner_partition()
    assert sum(region_count(cap) for cap in caps)==16499285
    assert region_count(central)==2903119660;checks+=2
    for profile,expected in (('coupled',111318573),('local',236981157),('retuned-wide',448226845)):
        caps,central=corner_partition(profile)
        assert sum(region_count(cap) for cap in caps)==expected
        assert region_count(central)==math.comb(516,4)-expected;checks+=2
    seeds,pending=retuned_partition()
    assert len(seeds)==23 and sum(region_count(x) for x in pending)==2273894055;checks+=1
    return checks


def radial_tests():
    """Known Gaussian integrals and an exact population-cancellation example."""
    zero=(0,0,0,0);checks=0
    # On this segment x+y=1, while its coordinate box also admits x+y=2.
    g=[{(1,0,0,0):I(1),(0,1,0,0):I(1)}]+[{zero:I()} for _ in range(3)]
    vertices=((1,0,0,0),(0,1,0,0));box=[I(0,1),I(0,1),I(),I()]
    plain=coupled_radial_bound({zero:I()},g,{zero:I()},I('.4'),I(),box,bernstein_polynomial)
    controlled=coupled_radial_bound({zero:I()},g,{zero:I()},I('.4'),I(),box,bernstein_polynomial,vertices=vertices)
    truth=I(1)/(4*I('.4'))-(4*I('.4').sq()).ln()
    assert controlled.hi>=truth.lo and controlled.hi<truth.hi+D('1e-10')
    assert plain.hi>controlled.hi+1;checks+=2
    polynomials=[{(2,0,0,0):I(2),(0,1,0,0):I(-1)},
                 {(1,1,0,0):I(-3),(0,0,0,0):I(2)}]
    vertices=((1,0,0,0),(0,1,0,0),(-1,-1,0,0))
    squares=polar_square_controls(polynomials,vertices)
    # Independent expansion in ordinary barycentric monomials, over Q.
    def multiply_exact(left,right):
        result={}
        for a,x in left.items():
            for b,y in right.items():
                key=tuple(i+j for i,j in zip(a,b));result[key]=result.get(key,F(0))+x*y
        return result
    basis=tuple(tuple(int(i==j) for i in range(3)) for j in range(3))
    total={};unit={key:F(1) for key in basis}
    for polynomial in polynomials:
        expanded={}
        for powers,coefficient in polynomial.items():
            term={(0,0,0):F(coefficient.lo)}
            for axis,power in enumerate(powers):
                linear={key:F(vertices[j][axis]) for j,key in enumerate(basis)}
                for _ in range(power):term=multiply_exact(term,linear)
            for _ in range(2-sum(powers)):term=multiply_exact(term,unit)
            for key,value in term.items():expanded[key]=expanded.get(key,F(0))+value
        for key,value in multiply_exact(expanded,expanded).items():total[key]=total.get(key,F(0))+value
    for alpha,coefficient in squares.items():
        truth=total.get(alpha,F(0))*F(math.prod(math.factorial(x) for x in alpha),24)
        assert F(coefficient.lo)<=truth<=F(coefficient.hi);checks+=1
    for weights in ((F(1),F(0),F(0)),(F(1,3),)*3,(F(1,2),F(1,3),F(1,6))):
        point=[sum(w*v[i] for w,v in zip(weights,vertices)) for i in range(4)]
        truth=sum(sum(F(c.lo)*math.prod(x**p for x,p in zip(point,powers)) for powers,c in poly.items())**2 for poly in polynomials)
        value=I()
        for alpha,c in squares.items():
            weight=F(24,math.prod(math.factorial(x) for x in alpha))*math.prod(w**a for w,a in zip(weights,alpha))
            value+=c*I(weight.numerator)/weight.denominator
        assert F(value.lo)<=truth<=F(value.hi);checks+=1
    for scalar in ('0.1','0.4','0.49'):
        a=I(scalar)
        for vector in ((0,0,0,0),(1,2,-1,3)):
            g=[I(x)/10 for x in vector]
            bound=coupled_radial_bound({zero:I(2)},[{zero:x} for x in g],{zero:I()},a,I(),[I()]*4,bernstein_polynomial)
            exact=I(2)+sum(x.sq() for x in g)/(4*a)-(4*a.sq()).ln()
            assert bound.hi>=exact.lo and bound.hi<exact.hi+D('1e-10');checks+=1
    gradients=[{tuple(int(i==j) for i in range(4)):I(1)} for j in range(4)]
    core={}
    for gradient in gradients:core=poly_add(core,poly_multiply(gradient,gradient),-I(5)/8)
    bound=coupled_radial_bound(core,gradients,{zero:I()},I('0.4'),I(),[I(-1,1)]*4,bernstein_polynomial)
    exact=-(I(4)*I('0.4').sq()).ln()
    assert exact.lo<=bound.hi<exact.hi+D('1e-10');checks+=1
    for lower in (1,2,4,8):
        a=I('0.4')
        bound=coupled_radial_bound({zero:I(2)},[{zero:I()}]*4,{zero:I()},a,I(),[I()]*4,bernstein_polynomial,lower)
        x=a*lower*lower
        exact=I(2)-(4*a.sq()).ln()+((1+x)*(-x).exp()).ln()
        assert exact.lo<=bound.hi<exact.hi+D('1e-10');checks+=1
    return checks


def bernstein_transform_tests():
    """Check every tensor coefficient against an independent rational formula."""
    polynomial={tuple(4*int(i==j) for i in range(4)):F((-1)**j*(j+1),7) for j in range(4)}
    polynomial.update({(1,1,1,1):F(2,3),(2,1,0,0):F(-5,11),(0,0,0,0):F(17,13)})
    lower=(-2,-1,0,1);width=(1,2,3,2);box=[I(x,x+w) for x,w in zip(lower,width)]
    coefficients={key:I(x.numerator)/x.denominator for key,x in polynomial.items()}
    start=time.monotonic();reference=bernstein_coefficients_reference(coefficients,box);reference_time=time.monotonic()-start
    start=time.monotonic();actual=bernstein_coefficients(coefficients,box);actual_time=time.monotonic()-start
    assert actual.keys()==reference.keys();checks=0
    for index,value in actual.items():
        truth=F(0)
        for powers,coefficient in polynomial.items():
            term=coefficient
            for row,power,lo,w in zip(index,powers,lower,width):
                term*=sum(F(math.comb(power,j)*lo**(power-j)*w**j*math.comb(row,j),math.comb(4,j))
                          for j in range(min(power,row)+1))
            truth+=term
        assert F(value.lo)<=truth<=F(value.hi)
        assert F(reference[index].lo)<=truth<=F(reference[index].hi);checks+=1
    print(dict(exact_transform_coefficients=checks,reference_seconds=reference_time,optimized_seconds=actual_time),flush=True)
    quadratic={(2,0,0,0):F(2,7),(0,1,1,0):F(-3,11),(0,0,0,1):F(5,13),(0,0,0,0):F(-7,17)}
    vertices=[(-1,)*4]+[tuple(2 if a==b else -1 for a in range(4)) for b in range(4)]
    coefficients={key:I(x.numerator)/x.denominator for key,x in quadratic.items()}
    actual=polar_coefficients(coefficients,vertices);reference=polar_coefficients_reference(coefficients,vertices)
    def evaluate_exact(point):return sum(c*math.prod(x**p for x,p in zip(point,powers)) for powers,c in quadratic.items())
    for (i,j),value in actual.items():
        left,right=vertices[i],vertices[j];middle=[F(x+y,2) for x,y in zip(left,right)]
        truth=2*evaluate_exact(middle)-(evaluate_exact(left)+evaluate_exact(right))/2
        assert F(value.lo)<=truth<=F(value.hi)
        assert F(reference[i,j].lo)<=truth<=F(reference[i,j].hi);checks+=1
    return checks


def triangulation_tests(star=False,refine=False):
    """Independent exact lattice coverage, slice volume and reconstruction."""
    def determinant(rows):
        a=[[F(x) for x in row] for row in rows];value=F(1)
        for col in range(len(a)):
            pivot=next((i for i in range(col,len(a)) if a[i][col]),None)
            if pivot is None:return F(0)
            if pivot!=col:a[col],a[pivot]=a[pivot],a[col];value=-value
            scale=a[col][col];value*=scale;a[col]=[x/scale for x in a[col]]
            for i in range(col+1,len(a)):
                scale=a[i][col];a[i]=[x-scale*y for x,y in zip(a[i],a[col])]
        return value
    def contains(simplex,point):
        size=len(simplex);rows=[[F(v[a]) for v in simplex]+[F(point[a])] for a in range(5)]
        rows.append([F(1)]*(size+1))
        for col in range(size):
            pivot=next(i for i in range(col,len(rows)) if rows[i][col])
            rows[col],rows[pivot]=rows[pivot],rows[col]
            scale=rows[col][col];rows[col]=[x/scale for x in rows[col]]
            for i in range(len(rows)):
                if i!=col:
                    scale=rows[i][col];rows[i]=[x-scale*y for x,y in zip(rows[i],rows[col])]
        return all(not row[-1] for row in rows[size:]) and all(row[-1]>=0 for row in rows[:size])
    checks=0;star_cases=0;refined_cases=0
    for n in range(6):
        for case in range(4):
            bounds=tuple((0,n if (a+case)%3 else n//2) for a in range(5))
            bounds=tighten_region(bounds,n)
            if bounds is None:continue
            vertices,simplices=verifier.population_triangulation(bounds,n,star)
            star_cases+=int(len(vertices)>len(region_vertices(bounds,n)))
            if refine:
                edges=sorted((i,j) for i,j in verifier.triangulation_pairs(simplices) if i!=j and
                             all((a+b)%2==0 for a,b in zip(vertices[i],vertices[j])))
                if edges:
                    vertices,simplices=verifier.bisect_simplex_edge(vertices,simplices,edges[0]);refined_cases+=1
            pairs=verifier.triangulation_pairs(simplices)
            for prefix in product(range(n+1),repeat=4):
                if sum(prefix)>n:continue
                point=prefix+(n-sum(prefix),)
                expected=all(lo<=x<=hi for x,(lo,hi) in zip(point,bounds))
                actual=any(contains([vertices[i] for i in s],point) for s in simplices)
                assert actual==expected;checks+=1
            if len(simplices[0])==5:
                volume=sum(abs(determinant([[vertices[j][a]-vertices[s[0]][a] for a in range(4)] for j in s[1:]])) for s in simplices)
                offset=n-sum(lo for lo,hi in bounds);width=[hi-lo for lo,hi in bounds]
                exact=sum((-1)**sum(mask)*max(0,offset-sum(w*b for w,b in zip(width,mask)))**4 for mask in product((0,1),repeat=5))
                assert volume==exact;checks+=1
            polynomial={(2,0,0,0):I(2),(0,1,1,0):I(-3),(0,0,0,0):I(7)}
            controls=polar_coefficients(polynomial,[v[:4] for v in vertices],control_pairs=pairs)
            unrestricted=polar_coefficients(polynomial,[v[:4] for v in vertices])
            for key,value in controls.items():
                assert (value.lo,value.hi)==(unrestricted[key].lo,unrestricted[key].hi);checks+=1
            for s in simplices:
                weights={i:F(1,len(s)) for i in s};point=[sum(weights[i]*vertices[i][a] for i in s) for a in range(4)]
                value=sum(controls[i,j]*I((weights[i]*weights[j]*(1 if i==j else 2)).numerator)/
                          (weights[i]*weights[j]*(1 if i==j else 2)).denominator for i in s for j in s if i<=j)
                truth=2*point[0]**2-3*point[1]*point[2]+7
                assert F(value.lo)<=truth<=F(value.hi);checks+=1
    if star:assert star_cases>0
    if refine:assert refined_cases>0
    return checks


def restricted_elevation_tests():
    """Exact blossom formula and reconstruction, independent of elevation weights."""
    from itertools import combinations_with_replacement
    checks=0
    coefficients={(2,0,0,0):F(2,7),(0,1,1,0):F(-3,11),
                  (0,0,0,1):F(5,13),(0,0,0,0):F(-7,17)}
    poly={key:I(c.numerator)/c.denominator for key,c in coefficients.items()}
    def enclosed(interval,truth):
        if truth==0:return interval.lo<=0<=interval.hi
        lower=(interval.lo<=0 and truth>0) or F(interval.lo)<=truth
        upper=(interval.hi>=0 and truth<0) or truth<=F(interval.hi)
        return lower and upper
    for bounds in (((0,2),)*5,((0,1),(0,2),(0,3),(0,3),(0,3))):
        vertices,simplices=verifier.population_triangulation(bounds,4)
        vertices=[v[:4] for v in vertices];pairs=verifier.triangulation_pairs(simplices)
        for degree in (3,4):
            controls=polar_coefficients(poly,vertices,degree,pairs)
            for simplex in simplices:
                for indices in combinations_with_replacement(simplex,degree):
                    alpha=tuple(indices.count(i) for i in range(len(vertices)))
                    assert alpha in controls;checks+=1
            for alpha,value in controls.items():
                points=[vertices[i] for i,a in enumerate(alpha) for _ in range(a)]
                truth=F(0)
                for powers,c in coefficients.items():
                    factors=[i for i,a in enumerate(powers) for _ in range(a)]
                    if not factors:weight=F(1)
                    elif len(factors)==1:weight=sum(F(p[factors[0]]) for p in points)/degree
                    else:
                        a,b=factors
                        weight=sum(F(p[a]*q[b]) for i,p in enumerate(points) for j,q in enumerate(points) if i!=j)/(degree*(degree-1))
                    truth+=c*weight
                assert enclosed(value,truth);checks+=1
            for simplex in simplices:
                for shape in (tuple(range(1,len(simplex)+1)),(1,)+(0,)*(len(simplex)-1)):
                    weights=[F(0)]*len(vertices)
                    for i,w in zip(simplex,shape):weights[i]=F(w,sum(shape))
                    point=[sum(w*v[a] for w,v in zip(weights,vertices)) for a in range(4)]
                    truth=sum(c*math.prod(x**p for x,p in zip(point,key)) for key,c in coefficients.items())
                    result=I();mass=F(0)
                    for alpha,value in controls.items():
                        weight=F(math.factorial(degree),math.prod(math.factorial(a) for a in alpha))
                        weight*=math.prod(w**a for w,a in zip(weights,alpha));mass+=weight
                        if weight:result+=value*I(weight.numerator)/weight.denominator
                    assert mass==1 and enclosed(result,truth);checks+=1
    # A non-clique triple is excluded, but the edge controls remain.
    pairs={(0,0),(1,1),(2,2),(0,1),(1,2)}
    plan=dict(verifier.restricted_elevation_plan(3,3,tuple(sorted(pairs))))
    assert (1,1,1) not in plan and (2,1,0) in plan and (0,1,2) in plan;checks+=1
    return checks


def matrix_control_tests():
    """Exact quadratic reconstruction and analytical Gaussian controls."""
    zero=(0,0,0,0);checks=0
    coefficients={powers:F((-1)**sum(powers)*(1+sum((i+1)*k for i,k in enumerate(powers))),13)
                  for powers in product(range(3),repeat=4) if sum(powers)<=2}
    def enclose(x):return I(x.numerator)/x.denominator
    polynomial={key:enclose(x) for key,x in coefficients.items()}
    def exact(point):return sum(value*math.prod(x**k for x,k in zip(point,powers)) for powers,value in coefficients.items())
    vertices=[(-1,)*4]+[tuple(2 if a==b else -1 for a in range(4)) for b in range(4)]
    controls=polar_coefficients(polynomial,vertices)
    for prefix in product(range(4),repeat=4):
        if sum(prefix)>3:continue
        weights=[F(x,3) for x in prefix+(3-sum(prefix),)]
        point=[sum(weight*v[a] for weight,v in zip(weights,vertices)) for a in range(4)]
        reconstructed=sum(value*enclose(weights[i]*weights[j]*(1 if i==j else 2)) for (i,j),value in controls.items())
        truth=exact(point)
        assert F(reconstructed.lo)<=truth<=F(reconstructed.hi);checks+=1
    for degree in (3,4):
        controls=polar_coefficients(polynomial,vertices,degree)
        for prefix in product(range(4),repeat=4):
            if sum(prefix)>3:continue
            weights=[F(x,3) for x in prefix+(3-sum(prefix),)]
            point=[sum(weight*v[a] for weight,v in zip(weights,vertices)) for a in range(4)]
            reconstructed=I()
            for alpha,value in controls.items():
                weight=F(math.factorial(degree),math.prod(math.factorial(k) for k in alpha))
                weight*=math.prod(x**k for x,k in zip(weights,alpha))
                reconstructed+=value*enclose(weight)
            truth=exact(point)
            assert F(reconstructed.lo)<=truth<=F(reconstructed.hi);checks+=1
    box=[I(-1,1),I(-1,1),I(),I()]
    for degree in (2,4):
        controls=bernstein_coefficients(polynomial,box,(degree,degree,0,0))
        for x,y in ((F(0),F(0)),(F(1,4),F(3,4)),(F(1,2),F(1,2)),(F(1),F(1))):
            reconstructed=I()
            for key,value in controls.items():
                weight=F(1)
                for k,t in zip(key[:2],(x,y)):weight*=math.comb(degree,k)*t**k*(1-t)**(degree-k)
                reconstructed+=value*enclose(weight)
            truth=exact((2*x-1,2*y-1,F(0),F(0)))
            assert F(reconstructed.lo)<=truth<=F(reconstructed.hi);checks+=1
    constant={zero:I()};linear=[{(1,0,0,0):I(1)}]+[{zero:I()}]*3
    precision={(i,j):{zero:I(int(i==j))} for i in range(4) for j in range(4)}
    line_box=[I(-1,1),I(),I(),I()]
    for vertex_set in (None,[(-1,0,0,0),(1,0,0,0)]):
        result=bernstein_gaussian_shell(constant,linear,precision,line_box,None,0,vertices=vertex_set)
        assert D('0.5')<=result.hi<D('0.500000001');checks+=1
    constant={(2,0,0,0):-I(1)/4};precision[0,0][1,0,0,0]=I(1)/5
    result=bernstein_gaussian_shell(constant,linear,precision,line_box,None,0)
    for numerator in range(-10,11):
        x=I(numerator)/10;p=1+x/5
        truth=-x.sq()/4+x.sq()/(2*p)-p.ln()/2
        assert result.hi>=truth.lo;checks+=1
    assert len(region_vertices(((64,128),)*5))==10;checks+=1
    return checks


def global_covariance_tests():
    checks=0
    for size in (1,2,3):
        for counts in weak_compositions(8,size):
            probabilities=[F(x,8) for x in counts]
            for i,p in enumerate(probabilities):
                diagonal=F(1,2)-p+p*p
                assert diagonal>=sum(p*q for j,q in enumerate(probabilities) if i!=j)
            checks+=1
    h=(100,110,120,90,92);fake={};exact_matrices={}
    def enclose(x):return I(x.numerator)/x.denominator
    def pair_count(v,a,b):return F(v[a]*(v[a]-1),2*511) if a==b else F(v[a]*v[b],511)
    for a,b in TYPES:
        vectors=[tuple(F(x) for x in OUTPUT[a,b,j][:4]) for j in range(3)]
        distinct=sorted(set(vectors));mean=[sum(v[i] for v in distinct)/len(distinct) for i in range(4)]
        exact_matrices[a,b]=[[sum((v[i]-mean[i])*(v[j]-mean[j]) for v in distinct)/2 for j in range(4)] for i in range(4)]
        shifted=[([enclose(x+F(1,3)) for x in v],I()) for v in vectors]
        fake[a,b]=(None,None,None,None,enclose(pair_count(h,a,b)),shifted)
    geometry=global_covariance_polynomials(h,fake,(0,1,2,4))
    for delta in ((0,0,0,0),(1,-1,0,0),(-2,1,1,0),(0,0,0,1)):
        v=list(h)
        for axis,amount in zip((0,1,2,4),delta):v[axis]+=amount;v[3]-=amount
        for i in range(4):
            for j in range(4):
                truth=sum(pair_count(v,a,b)*matrix[i][j] for (a,b),matrix in exact_matrices.items())
                result=evaluate_polynomial(geometry[i,j],[I(x) for x in delta])
                assert F(result.lo)<=truth<=F(result.hi);checks+=1
    return checks


def attempt_region(bounds):
    try:
        tightened=tighten_region(bounds);assert tightened is not None
        result=verify_region(tightened,coarse_only=max(hi-lo for lo,hi in tightened)>32,
                             coupled_coarse=True,local_radial='auto')
        return result['passes'],result['directed_log_upper'],None
    except AssertionError as error:
        message=str(error)
        if not (message.startswith("('nonpositive pivot'") or message in
                ('spectral boundary nonpositive','could not certify spectral boundary','tail coefficient nonpositive')):
            raise
        return False,None,message.split(',')[0]


def corner_partition(profile='original'):
    thresholds={'original':(448,416,416,416,416),'coupled':(384,384,352,352,352),
                'local':(352,352,320,320,320),'retuned-wide':(320,320,288,288,288)}[profile]
    caps=[tighten_region(tuple((thresholds[a],512) if a==dominant else (0,512) for a in range(5))) for dominant in range(5)]
    central=tighten_region(tuple((0,x-1) for x in thresholds))
    assert all(thresholds[a]+thresholds[b]>512 for a in range(5) for b in range(a))
    assert sum(region_count(cap) for cap in caps)+region_count(central)==math.comb(516,4)
    return caps,central


def bounded_cover(bounds,max_nodes=6,max_seconds=20,seed_corners=False):
    root=tighten_region(bounds);assert root is not None
    target=region_count(root);pending=[root];certified=0;unclosed_singletons=0
    if seed_corners:
        assert root==((0,512),)*5
        caps,central=corner_partition('local');pending=[central]+list(reversed(caps))
    initial_pending=len(pending)
    attempts=0;regions=0;reasons={};worst=None;started=time.monotonic()
    while pending and attempts<max_nodes and time.monotonic()-started<max_seconds:
        region=pending.pop();size=region_count(region);attempts+=1
        passed,bound,reason=attempt_region(region)
        if passed:
            certified+=size;regions+=1
            if worst is None or D(bound)>D(worst):worst=bound
        elif size==1:unclosed_singletons+=1
        else:pending.extend(reversed(split_region(region)))
        if reason:reasons[reason]=reasons.get(reason,0)+1
        outstanding=sum(region_count(x) for x in pending)
        assert certified+unclosed_singletons+outstanding==target
        assert len(pending)<=max_nodes+initial_pending
    return dict(histograms=target,certified=certified,pending=sum(region_count(x) for x in pending),
                unclosed_singletons=unclosed_singletons,regions_passed=regions,attempts=attempts,
                inconclusive_reasons=reasons,worst_certified_log_bound=worst,
                complete_histogram_cover=certified==target,whole_code_certificate=False,
                seconds=time.monotonic()-started)


def core_hessian_blocks(h,entries,kappa):
    n=sum(h);axes=(0,1,2,4)
    basis=[[int(a==axis)-int(a==3) for a in range(5)] for axis in axes]
    f=[[entries[min(a,b),max(a,b)][0] for b in range(5)] for a in range(5)]
    hh=[[sum(basis[i][a]*basis[j][b]*(f[a][b]/(n-1)-2*kappa/n*int(a==b)) for a in range(5) for b in range(5))
         for j in range(4)] for i in range(4)]
    hz=[[0.]*4 for _ in range(4)]
    for j in range(4):
        g=[[entries[min(a,b),max(a,b)][1][j] for b in range(5)] for a in range(5)]
        gradient=[sum(g[a][b]*h[b] for b in range(5))/(n-1)-g[a][a]/(2*(n-1)) for a in range(5)]
        for i in range(4):hz[i][j]=dot(basis[i],gradient)
    zz=[[sum(value[4]*value[2][i][j] for value in entries.values())-float(i==j) for j in range(4)] for i in range(4)]
    return hh,hz,zz


def joint_curvature_probe(h):
    """Floating test for joint concavity; excludes matching correction."""
    entries,_,_=local_data(h,[x/10**9 for x in LOG_WEIGHT_NUMERATORS],.1,.2,.08)
    hh,hz,zz=core_hessian_blocks(h,entries,.08)
    inverse=[solve(zz,hz[j]) for j in range(4)]
    schur=[[hh[i][j]-dot(hz[i],inverse[j]) for j in range(4)] for i in range(4)]
    eigen,vectors=eigensystem(schur);index=max(range(4),key=lambda i:eigen[i])
    direction=[vectors[i][index] for i in range(4)]
    return dict(histogram=h,core_hh_max_eigenvalue=max(eigensystem(hh)[0]),
                core_joint_schur_max_eigenvalue=eigen[index],histogram_direction=direction,
                fixed_fifth_schur_max_eigenvalue=max(eigensystem([row[:3] for row in schur[:3]])[0]),
                diagnostic_not_certificate=True)


def certify_nonconcavity(h,fixed_fifth=False):
    """A positive rational direction rejects the natural core's concavity.

    This rejects one proof shortcut, not the desired code or its moment bound.
    """
    u=[x/10**9 for x in LOG_WEIGHT_NUMERATORS]
    floating,_,_=local_data(h,u,.1,.2,.08);hh,hz,zz=core_hessian_blocks(h,floating,.08)
    best=None
    for direction in product((-1,0,1),repeat=4):
        if not any(direction) or fixed_fifth and direction[3]:continue
        mixed=[sum(hz[i][j]*direction[i] for i in range(4)) for j in range(4)]
        gaussian=[-x for x in solve(zz,mixed)]
        score=sum(direction[i]*hh[i][j]*direction[j] for i in range(4) for j in range(4))+dot(mixed,gaussian)
        if best is None or score>best[0]:best=(score,direction,gaussian)
    _,direction,gaussian=best
    gaussian=[I(str(round(x,10))) for x in gaussian]
    _,proposal,_=variational(h,u,.1,.2,.08);center,_=moments(h,tables(proposal))
    center=[I(str(round(x,12))) for x in center[:4]];center.append(I(512)-sum(center))
    entries,_,_,_,kappa=interval_data(h,center);hh,hz,zz=core_hessian_blocks(h,entries,kappa)
    curvature=sum(direction[i]*hh[i][j]*direction[j]+2*direction[i]*hz[i][j]*gaussian[j]
                  +gaussian[i]*zz[i][j]*gaussian[j] for i in range(4) for j in range(4))
    assert curvature.lo>0
    return dict(histogram=h,fixed_fifth=fixed_fifth,histogram_direction=direction,
                gaussian_direction=[str(x.lo) for x in gaussian],
                directed_curvature_lower=str(curvature.lo),rejects_natural_core_joint_concavity=True,
                rejects_code=False)


def main():
    assert shutil.disk_usage(Path.cwd()).free>=20*(1<<30),'disk guard'
    assert not ('--precision80' in sys.argv and '--precision40' in sys.argv)
    if '--precision80' in sys.argv:getcontext().prec=80
    elif '--precision40' in sys.argv:getcontext().prec=40
    limit=float(next((x.split('=',1)[1] for x in sys.argv if x.startswith('--deadline-seconds=')),'25'))
    assert 1<=limit<=120
    started=time.monotonic();verifier.VERIFY_DEADLINE=started+limit
    checks=domain_tests()+radial_tests()+matrix_control_tests()+global_covariance_tests()+partition_tests()+bernstein_transform_tests()+triangulation_tests()+triangulation_tests(True)+triangulation_tests(refine=True)+restricted_elevation_tests()
    print(dict(exact_domain_and_coordinate_checks=checks,all_histograms=region_count(((0,512),)*5)),flush=True)
    simplex_probe=next((x.split('=',1)[1] for x in sys.argv if x.startswith('--simplex-probe=')),None)
    if simplex_probe is not None:
        cases={'merged':((0,43),(44,87),(80,159),(0,79),(240,319)),
               'extended':((0,63),(44,87),(80,159),(0,79),(240,319)),
               'balanced':((64,128),)*5,
               'mixed':((0,79),(80,159),(180,233),(72,252),(0,180))}
        bounds=cases[simplex_probe];mode='star' if '--star-controls' in sys.argv else 'triangulated'
        degree=4 if '--degree4' in sys.argv else 3 if '--degree3' in sys.argv else 2
        splits=next((tuple(tuple(int(y) for y in edge.split(',')) for edge in x.split('=',1)[1].split(';'))
                     for x in sys.argv if x.startswith('--simplex-splits=')),())
        result=verify_region(bounds,matrix_mode=mode,log_weights=RETUNED_LOG_WEIGHT_NUMERATORS,
                             coupled_coarse=True,local_radial='auto',radial_controls=degree==2,simplex_splits=splits,
                             control_degree=degree,bounded_switch='--bounded-switch' in sys.argv)
        print(dict(simplex_probe=simplex_probe,mode=mode,bounds=bounds,histograms=region_count(bounds),
                   upper=result['directed_log_upper'],passes=result['passes'],detail=result.get('central_terms')),flush=True)
        if '--control-midpoint' in sys.argv:
            assert degree==2,'edge midpoint diagnostic requires degree two'
            vertices,simplices=verifier.population_triangulation(bounds,star=mode=='star')
            for edge in splits:vertices,simplices=verifier.bisect_simplex_edge(vertices,simplices,edge)
            i,j=result['central_terms']['worst_control']
            assert all((a+b)%2==0 for a,b in zip(vertices[i],vertices[j]))
            midpoint=tuple((a+b)//2 for a,b in zip(vertices[i],vertices[j]))
            point=verifier.verify_histogram(midpoint,log_weights=RETUNED_LOG_WEIGHT_NUMERATORS)
            print(dict(first_shell_worst_edge=(vertices[i],vertices[j]),midpoint=midpoint,
                       point_log_upper=point['directed_log_moment_upper'],passes=point['passes']),flush=True)
    if '--verify-partition' in sys.argv:print(verify_interior_partition(),flush=True)
    if '--verify-retuned-partition' in sys.argv:
        selection=next((tuple(int(y) for y in x.split('=',1)[1].split(':')) for x in sys.argv if x.startswith('--partition-slice=')),(0,None))
        assert len(selection)==2
        print(verify_retuned_partition(*selection),flush=True)
    if '--pending-cover' in sys.argv:pending_cover()
    if '--retuned-cover-pilot' in sys.argv or '--retuned-pending-cover' in sys.argv:
        nodes=next((int(x.split('=',1)[1]) for x in sys.argv if x.startswith('--cover-nodes=')),24)
        assert 1<=nodes<=64
        matrix_width=next((int(x.split('=',1)[1]) for x in sys.argv if x.startswith('--matrix-width=')),96)
        assert 16<=matrix_width<=256
        retuned_cover_pilot(max_nodes=nodes,max_seconds=min(100,limit-8),hybrid='--hybrid-cover' in sys.argv,
                           triangulated='--triangulated-cover' in sys.argv,wide_caps='--wide-caps' in sys.argv,
                           matrix_width=matrix_width,resume='--retuned-pending-cover' in sys.argv)
    if '--retuned-radial-controls' in sys.argv:
        cases=(
            ('low_first_population',((0,79),(0,192),(0,192),(0,192),(0,319))),
            ('balanced_band',((64,128),)*5),
            ('boundary_box',((173,205),(152,184),(99,131),(0,512),(0,16))),
        )
        for name,bounds in cases:
            for controlled in (False,True):
                if controlled and '--degree4' in sys.argv and math.comb(len(region_vertices(tighten_region(bounds)))+3,4)>2000:
                    print(dict(region=name,inconclusive='degree4 polar control count guard'),flush=True);continue
                result=verify_region(bounds,coarse_only=True,coupled_coarse=True,
                                     local_radial='shells' if '--radial-shells' in sys.argv else 'auto',
                                     log_weights=RETUNED_LOG_WEIGHT_NUMERATORS,radial_controls=controlled,
                                     control_degree=4 if '--degree4' in sys.argv else 2)
                print(dict(region=name,controls=controlled,histograms=result['integer_histograms'],
                           log_upper=result['directed_log_upper'],passes=result['passes']),flush=True)
    if '--retuned-caps' in sys.argv:
        for index,bounds in enumerate(corner_partition('local')[0]):
            result=verify_region(bounds,coarse_only=True,coupled_coarse=True,local_radial='auto',
                                 log_weights=RETUNED_LOG_WEIGHT_NUMERATORS)
            print(dict(retuned_cap=index,bounds=bounds,integer_histograms=result['integer_histograms'],
                       directed_log_upper=result['directed_log_upper'],passes=result['passes']),flush=True)
    if '--retuned-probe' in sys.argv:
        weights=RETUNED_LOG_WEIGHT_NUMERATORS
        print(dict(proof_weights=weights,denominator=10**9,allowance_lower=str(verifier.allowance_interval(weights).lo)),flush=True)
        if '--only-retuned-band' not in sys.argv:
            result=verifier.verify_histogram((189,168,115,40,0),log_weights=weights)
            print(dict(retuned_point=result),flush=True)
        cases=(
            ('previously_failing_band',((64,128),)*5),
            ('retuned_boundary_box',((173,205),(152,184),(99,131),(0,512),(0,16))),
        )
        if '--only-retuned-band' in sys.argv:cases=cases[:1]
        for name,bounds in cases:
            result=verify_region(bounds,matrix_mode='polytope',control_degree=3 if '--degree3' in sys.argv else 2,log_weights=weights)
            print(dict(region=name,integer_histograms=result['integer_histograms'],
                       directed_log_upper=result['directed_log_upper'],passes=result['passes']),flush=True)
    if '--probe' in sys.argv:
        probes=(('absent_state3',((124,132),(124,132),(124,132),(0,0),(124,132))),
                ('near_state0_vertex',((508,512),(0,4),(0,4),(0,4),(0,4))),
                ('near_state4_vertex',((0,4),(0,4),(0,4),(0,4),(508,512))))
        for label,bounds in probes:
            assert time.monotonic()-started<24
            passed,bound,reason=attempt_region(bounds)
            print(dict(region=label,histograms=region_count(bounds),passes=passed,log_upper=bound,reason=reason),flush=True)
    if '--root' in sys.argv:
        print(bounded_cover(((0,512),)*5,max_nodes=32,max_seconds=18,seed_corners=True),flush=True)
    if '--matrix-regions' in sys.argv:
        cases=(
            ('boundary_r2',((183,187),(186,190),(81,85),(0,512),(0,0))),
            ('boundary_r4',((181,189),(184,192),(79,87),(0,512),(0,0))),
            ('uniform_box',((88,116),)*5),
        )
        if '--only-boundary-r4' in sys.argv:cases=cases[1:2]
        if '--broad-interior' in sys.argv:cases=(('balanced_band',((64,128),)*5),)
        if '--middle-band' in sys.argv:cases=(('middle_band',((80,136),)*5),)
        if '--narrow-band' in sys.argv:cases=(('narrow_band',((88,124),)*5),)
        if '--dominant-band' in sys.argv:cases=(('dominant0_256_351',((256,351),(0,256),(0,256),(0,256),(0,256))),)
        for label,bounds in cases:
            assert time.monotonic()-started<24
            mode='polytope' if '--polytope-controls' in sys.argv else 'bernstein' if '--integral-controls' in sys.argv else 'common'
            priority=next((tuple(int(c) for c in x.split('=',1)[1]) for x in sys.argv if x.startswith('--exposure-priority=')),None)
            result=verify_region(bounds,matrix_mode=mode,
                                 control_degree=4 if '--degree4' in sys.argv else 3 if '--degree3' in sys.argv else 2,
                                 exposure_priority=priority)
            print(dict(region=label,histograms=result['integer_histograms'],passes=result['passes'],
                       log_upper=result['directed_log_upper'],components=result.get('central_terms')),flush=True)
            if mode=='bernstein' and '--only-boundary-r4' in sys.argv and '--degree4' in sys.argv:assert result['passes']
            if mode=='polytope' and '--narrow-band' in sys.argv:assert result['passes']
    if '--interior' in sys.argv or '--interior-narrow' in sys.argv:
        cases=(
            ('uniform_box',((88,116),)*5),
            ('boundary_box',((177,193),(180,196),(75,91),(0,512),(0,16))),
            ('vertex0_d96',((416,512),(0,96),(0,96),(0,96),(0,96))),
            ('vertex2_d128',((0,128),(0,128),(384,512),(0,128),(0,128))),
        ) if '--interior' in sys.argv else (
            ('uniform_narrow',((96,108),)*5),
            ('boundary_narrow',((183,187),(186,190),(81,85),(0,512),(0,4))),
            ('uniform_singleton',((103,103),(103,103),(102,102),(102,102),(102,102))),
            ('boundary_singleton',((185,185),(188,188),(83,83),(56,56),(0,0))),
        )
        for label,bounds in cases:
            assert time.monotonic()-started<24
            mode='shells' if '--radial-shells' in sys.argv else '--local-radial' in sys.argv
            result=verify_region(bounds,coarse_only=True,coupled_coarse=True,local_radial=mode)
            print(dict(region=label,histograms=result['integer_histograms'],passes=result['passes'],
                       log_upper=result['directed_log_upper'],components=result['coarse_detail']),flush=True)
            if '--interior-narrow' in sys.argv and mode=='shells':
                assert result['passes']==label.startswith('uniform_')
    if '--coupled-ladder' in sys.argv:
        for deficit in (128,160,192):
            for dominant in range(5):
                assert time.monotonic()-started<24
                bounds=tuple((512-deficit,512) if a==dominant else (0,deficit) for a in range(5))
                result=verify_region(bounds,coarse_only=True,coupled_coarse=True)
                print(dict(deficit=deficit,dominant=dominant,histograms=result['integer_histograms'],
                           log_upper=result['directed_log_upper'],passes=result['passes']),flush=True)
    if '--local-cap-controls' in sys.argv:
        for deficit,dominant in ((160,0),(160,1),(192,2),(192,3),(192,4)):
            assert time.monotonic()-started<24
            bounds=tuple((512-deficit,512) if a==dominant else (0,deficit) for a in range(5))
            result=verify_region(bounds,coarse_only=True,coupled_coarse=True,local_radial=True)
            print(dict(deficit=deficit,dominant=dominant,histograms=result['integer_histograms'],
                       log_upper=result['directed_log_upper'],passes=result['passes'],
                       components=result['coarse_detail']),flush=True)
            assert result['passes']
    if '--coarse-vertices' in sys.argv or '--coarse-ladder' in sys.argv:
        for deficit in ((64,96,128) if '--coarse-ladder' in sys.argv else (64,)):
            for dominant in range(5):
                assert time.monotonic()-started<24
                bounds=tuple((512-deficit,512) if a==dominant else (0,deficit) for a in range(5))
                result=verify_region(bounds,coarse_only=True)
                print(dict(deficit=deficit,dominant=dominant,histograms=result['integer_histograms'],log_upper=result['directed_log_upper'],
                           passes=result['passes'],pure_simplex=result['pure_simplex']),flush=True)
    if '--curvature' in sys.argv:
        for h in ((185,188,83,56,0),(185,188,83,55,1),(103,103,102,102,102),
                  (509,1,1,1,0),(1,509,1,1,0),(1,1,509,1,0),(1,1,1,509,0)):
            print(joint_curvature_probe(h),flush=True)
        print(certify_nonconcavity((185,188,83,54,2)),flush=True)
        print(certify_nonconcavity((506,2,2,2,0),True),flush=True)
    assert time.monotonic()-started<limit
    print(dict(seconds=time.monotonic()-started,precision=getcontext().prec,whole_code_certificate=False),flush=True)


if __name__=='__main__':main()
