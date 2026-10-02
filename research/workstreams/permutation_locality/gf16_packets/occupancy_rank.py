"""Fixed-occupancy operators with trimmed-output and rank return bounds."""
from fractions import Fraction as Q
import numpy as np
from flint import arb,ctx
import occupancy_kernel as base
import trimmed_return
import rank_return

aq,up=base.aq,base.up
TERMINAL=base.TERMINAL


def outward_at_z(data,z):
    if 'rank_counts' not in data:data=rank_return.attach(trimmed_return.attach(data))
    matrices=base.outward_at_z(data,z);maximum,mean=rank_return.fixed_bounds(data,z)
    alpha=arb(2)**-data['updates'];L=(1<<data['bits'])-1
    for j,matrix in enumerate(matrices):
        matrix[1,0]=min(matrix[1,0],up(alpha*maximum[j]+matrix[1,2]/L))
        matrix[2,0]=min(matrix[2,0],up(alpha*mean[j]+matrix[2,2]/L))
    return matrices


def outward(data,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,(-aq(tilt)).exp())


def build_operators(args):
    from occupancy_model import placement
    from occupancy_memory import rounded
    updates=getattr(args,'updates',2)
    if type(updates) is not int or updates not in (2,3,4):raise ValueError('two through four inner updates required')
    data=rank_return.actual(updates);ctx.prec=args.precision;result={}
    for tilt in args.tilts:
        exact=placement(outward(data,tilt),rounding=rounded,maximum_groups=args.groups)
        arrays=[np.array([[float(matrix[i,j]) for j in range(3)] for i in range(3)]) for matrix in exact]
        result[tilt,'1']=exact,arrays
        print('GF16 rank fixed-occupancy operators',tilt,'degree',args.groups,'updates',updates,flush=True)
    return result
