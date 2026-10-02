"""Fixed-occupancy GF16 operators retaining newborn expansion-weight classes.

The packet locations are a uniform subset and their nonzero labels are
independent uniform GF16 values. No iid occupancy approximation is used.
"""
from fractions import Fraction as Q
from math import comb
import numpy as np
from flint import arb,arb_mat,arb_poly,ctx
import birth_classes
import occupancy_kernel
import occupancy_rank
import rank_return

aq,up=occupancy_kernel.aq,occupancy_kernel.up
prepare,actual=birth_classes.prepare,birth_classes.actual


def class_masses(data,z,*,include_zero=False):
    """Weighted class masses, conditional on each exact local occupancy."""
    if not 0<z<=1:raise ValueError('output weight in (0,1] required')
    W=data['windows'];S=1<<data['bits']
    values=[((1+z)**(4-r)*(1-z)**r-1)/15 for r in range(5)]
    # The coefficient of t^j sums over all j-subsets of packet locations.
    # Each a_r already averages the fifteen nonzero packet labels.
    powers=[[arb_poly([1,value])**j for j in range(W+1)] for value in values]
    polynomials=[]
    for profile in data['records']:
        polynomial=arb_poly([1])
        for row,n in zip(powers,profile):polynomial*=row[int(n)]
        polynomials.append([polynomial[j] for j in range(W+1)])
    if type(include_zero) is not bool:raise ValueError('boolean zero-class option required')
    rows=[[int(n) for n in row] for row in data['birth_class_census']]
    # The Walsh transform of the indicator of {0} is one at every
    # character. Its census by character profile is the multiplicity.
    if include_zero:rows.insert(0,list(map(int,data['multiplicities'])))
    census=arb_mat(rows)
    sums=census*arb_mat(polynomials)
    return [[max(arb(0),up(sums[i,j]/(S*comb(W,j)))) for i in range(sums.nrows())]
            for j in range(W+1)]


def refine_zero(data,local,z):
    """Replace only zero-to-zero bounds by their weighted Fourier count."""
    if len(local)!=data['windows']+1 or any(m.nrows()<1 or m.ncols()<1 for m in local):
        raise ValueError('one local matrix per packet occupancy required')
    masses=class_masses(data,z,include_zero=True);result=[]
    for matrix,row in zip(local,masses):
        value=matrix*arb(1);value[0,0]=min(value[0,0],row[0]);result.append(value)
    return result


def outward_at_z(data,z):
    original=occupancy_rank.outward_at_z(data,z)
    births=class_masses(data,z)
    profiles=[occupancy_kernel.polynomial(histogram,z) for histogram in data['histograms']]
    fibers=rank_return.fixed_profile_bounds(data,z)
    levels=list(map(int,data['birth_class_levels']));L=(1<<data['bits'])-1
    alpha=arb(2)**-data['updates'];beta=1-alpha;n=3+len(levels);result=[]
    for j,old in enumerate(original):
        rows=[[arb(0)]*n for _ in range(n)]
        for i in range(3):
            for k in range(3):rows[i][k]=old[i,k]
        rows[0][1]=arb(0);rows[0][2]=arb(0);rows[0][3:]=births[j]
        for i,level in enumerate(levels,3):
            selected=[k for k,w in enumerate(data['image_histogram_weights']) if w==level]
            total=max(up(profiles[k][j]) for k in selected)
            quiet=z**level if j==0 else arb(0)
            lazy_return=max(fibers[j][k] for k in selected)
            rows[i][0]=min(old[1,0],total,up(alpha*lazy_return+beta*total/L))
            rows[i][1]=max(arb(0),up(alpha*(total-quiet)))
            rows[i][2]=up(beta*total);rows[i][i]=up(alpha*quiet)
        matrix=arb_mat(rows)
        if any(matrix[i,k]<0 for i in range(n) for k in range(n)):
            raise ArithmeticError('negative fixed-occupancy class envelope')
        result.append(matrix)
    return result


def outward(data,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,(-aq(tilt)).exp())


def build_operators(args):
    from occupancy_model import placement
    from occupancy_memory import rounded
    if not args.exact_feedback:raise ValueError('exact feedback required for GF birth classes')
    updates=getattr(args,'updates',2)
    if type(updates) is not int or updates not in (2,3,4):raise ValueError('two through four inner updates required')
    joint=getattr(args,'joint_return_through',None)
    density=getattr(args,'lazy_density_through',None)
    for value,limit in ((joint,4),(density,32)):
        if value is not None and (type(value) is not int or not 0<=value<=limit):
            raise ValueError('valid integer local-refinement cutoff required')
    data=actual(updates);ctx.prec=args.precision;result={}
    if joint is not None:
        import return_moment
        checked=return_moment.actual_census(data,joint)
    for tilt in args.tilts:
        local=outward(data,tilt)
        z=(-aq(Q(tilt))).exp()
        # The return bound uses the original U entry as refresh mass.
        # Apply it before adding lazy mass to that coordinate.
        if joint is not None:local=return_moment.refine_class_returns(data,local,checked,z)
        if density is not None:
            import lazy_density
            local=lazy_density.refine_actual(data,local,z,density)
        exact=placement(local,rounding=rounded,maximum_groups=args.groups)
        n=local[0].nrows()
        arrays=[np.array([[float(matrix[i,j]) for j in range(n)] for i in range(n)]) for matrix in exact]
        result[tilt,'1']=exact,arrays
        print('GF16 birth-class fixed-occupancy operators',tilt,'states',n,'degree',args.groups,'updates',updates,
              'joint returns',joint,'lazy density',density,flush=True)
    return result
