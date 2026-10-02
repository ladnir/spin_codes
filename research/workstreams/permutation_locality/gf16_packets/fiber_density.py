"""Use GF feedback fibers to retain density from arbitrary entering states.

The fiber target may be zero. These local bounds are not complete
distance certificates; callers must cover and sum the outer domain.
"""
from fractions import Fraction as Q
from math import comb
from flint import arb
import rank_return
import trimmed_return

aq,up=rank_return.aq,rank_return.up


def attach(data):
    """Trim by the largest atom over ALL feedback targets, including zero."""
    W=data['windows'];denoms=[comb(W,j)*15**j for j in range(W+1)]
    budgets=[int(max(Q(a),Q(b))*D) for a,b,D in
             zip(data['nonzero_atom_caps'],data['zero_probabilities'],denoms)]
    counts=[trimmed_return.output_counts(h) for h in data['histograms']]
    trimmed=[[trimmed_return.lightest(row,budgets[j]) for j,row in enumerate(table)] for table in counts]
    return dict(data,all_fiber_trimmed=trimmed)


def profile_caps(data,z):
    """Bound E[z^wt(Aa+X) 1{CX=c}|J=j] for every a!=0 and every c."""
    if not 0<z<=1:raise ValueError('output weight in (0,1] required')
    if 'all_fiber_trimmed' not in data:data=attach(data)
    W=data['windows'];powers=[z**w for w in range(4*W+1)];result=[]
    for j in range(W+1):
        denominator=comb(W,j)*15**j
        factor=(1+z)**(4*j-data['bits'])/denominator
        result.append([min(
            up(factor*sum((int(n)*p for n,p in zip(rank[j],powers) if n),arb(0))),
            up(sum((int(n)*p for n,p in zip(trimmed[j],powers) if n),arb(0))/denominator))
            for rank,trimmed in zip(data['rank_counts'],data['all_fiber_trimmed'])])
    return result


def replace_uniform_classes(data,local,source,z,start,through=None,caps=None):
    """Select complete U rows rebuilt before the lazy density allocation.

    Both local and source must be valid operators for the same step and
    coordinate invariants. Whole-row selection preserves that contract;
    subtracting a rounded lazy-density contribution would not.
    """
    replacement=candidate(data,source,z,start,through,caps,'classes',include_uniform=True)
    through=data['windows'] if through is None else through
    n=3+len(data['birth_class_levels'])
    if len(local)!=len(replacement) or any(m.nrows()!=n or m.ncols()!=n for m in local):
        raise ValueError('matching birth-class operators required')
    result=[]
    for j,matrix in enumerate(local):
        value=matrix*arb(1)
        if start<=j<=through:
            for k in range(n):value[2,k]=replacement[j][2,k]
        result.append(value)
    return result


def candidate(data,local,z,start,through=None,caps=None,allocation='density',include_uniform=False):
    """Allocate each selected lazy branch using density or class mass caps."""
    W=data['windows'];through=W if through is None else through
    if (type(start) is not int or type(through) is not int or not 1<=start<=through<=W
            or len(local)!=W+1):
        raise ValueError('valid positive packet interval and matching operators required')
    if allocation not in ('density','classes'):raise ValueError('density or class allocation required')
    if type(include_uniform) is not bool or (include_uniform and allocation!='classes'):
        raise ValueError('uniform-source option requires class allocation')
    caps=profile_caps(data,z) if caps is None else caps
    if len(caps)!=W+1 or any(len(row)!=len(data['histograms']) for row in caps):
        raise ValueError('one bound per occupancy and expansion profile required')
    levels=list(map(int,data['birth_class_levels']));n=3+len(levels)
    if any(m.nrows()!=n or m.ncols()!=n for m in local):
        raise ValueError('birth-class operators required')
    selected=[[i for i,w in enumerate(data['image_histogram_weights']) if w==level] for level in levels]
    if any(not ids for ids in selected):raise ValueError('empty expansion class')
    L=(1<<data['bits'])-1;alpha=arb(2)**-data['updates'];result=[]
    sizes=list(map(int,data['birth_class_counts']))
    if len(sizes)!=len(levels) or any(n<=0 for n in sizes) or sum(sizes)!=L:
        raise ValueError('complete nonzero-state class census required')
    for j,matrix in enumerate(local):
        value=matrix*arb(1)
        if start<=j<=through:
            bounds=[max(caps[j]),*(max(caps[j][i] for i in ids) for ids in selected)]
            sources=[1,*range(3,n)]
            if include_uniform:
                # An entering measure bounded pointwise by U/L is
                # dominated by the uniform nonzero-state measure.
                bounds.append(up(sum((int(count)*cap for count,cap in
                    zip(data['histogram_multiplicities'],caps[j])),arb(0))/L))
                sources.append(2)
            for source,bound in zip(sources,bounds):
                mass=matrix[source,1]
                # Earlier refinements may already allocate this branch
                # to U or F; leave those allocations untouched.
                if not mass:continue
                value[source,1]=arb(0)
                if allocation=='density':value[source,2]=up(value[source,2]+alpha*L*bound)
                else:
                    # Both inequalities bound the SAME outgoing class
                    # mass: total lazy mass, and size times point cap.
                    # This is not an entrywise min of alternative rows.
                    for target,size in enumerate(sizes,3):
                        cap=min(mass,up(alpha*size*bound))
                        value[source,target]=up(value[source,target]+cap)
        result.append(value)
    return result
