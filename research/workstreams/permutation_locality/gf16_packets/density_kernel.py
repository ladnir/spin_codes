"""Fixed-occupancy GF16 envelope with a persistent pointwise density cap.

Coordinates: zero mass, arbitrary mass, cap per arbitrary state, uniform
envelope mass. The cap is not probability mass and has terminal zero.
"""
from fractions import Fraction as Q
from math import comb
import numpy as np
from flint import arb,arb_mat,ctx
import occupancy_kernel as base

aq,up=base.aq,base.up
prepare=base.prepare
actual=base.actual
TERMINAL=np.array([1.,1.,0.,1.])


def outward_at_z(data,z,single=None,*,classes=False):
    if not 0<z<=1:raise ValueError('output weight in (0,1] required')
    L=(1<<data['bits'])-1;alpha=arb(2)**-data['updates'];beta=1-alpha
    h,mean=base.moments(data,z);h2,mean2=base.moments(data,z*z)
    distance=min(sum(w*int(n) for w,n in enumerate(row)) for row in data['histograms'])
    if classes:
        levels=sorted(set(sum(w*int(n) for w,n in enumerate(row)) for row in data['histograms']))
        sizes=[];class_mean=[];class_mean2=[]
        for level in levels:
            subset=[i for i,row in enumerate(data['histograms']) if sum(w*int(n) for w,n in enumerate(row))==level]
            selected=dict(data,histograms=data['histograms'][subset],
                          histogram_multiplicities=data['histogram_multiplicities'][subset])
            size=sum(int(x) for x in selected['histogram_multiplicities']);sizes.append(size)
            class_mean.append([up(x*L/size) for x in base.moments(selected,z)[1]])
            class_mean2.append([up(x*L/size) for x in base.moments(selected,z*z)[1]])
    else:levels=[None];sizes=[L];class_mean=[mean];class_mean2=[mean2]
    packet=((1+z)**4-1)/15;packet2=((1+z*z)**4-1)/15
    result=[];weights=[arb(1)]
    for q,(zero,atom) in enumerate(zip(data['zero_probabilities'],data['nonzero_atom_caps'])):
        zero=aq(zero);atom=aq(atom);nonzero=1-zero;mass=packet**q
        zz=min(up(mass),up(z**q*zero),up((packet2**q*zero).sqrt()))
        zn=up(mass-z**(4*q)*zero)
        zd=min(up(mass),up(z**q*atom))
        pointwise=z**max(0,distance-4*q)
        g=up(sum((p*z**max(0,distance-w) for w,p in enumerate(weights)),arb(0)))
        cancellation=min(g,up(pointwise*nonzero))
        if single is not None and q:
            if q==1:
                g=min(g,up(arb(single['density'])/single['denominator']))
                cancellation=min(cancellation,up(arb(single['cancellation'])/single['denominator']))
                zd=min(zd,up(arb(single['feedback'])/single['denominator']))
                zz=min(zz,up(arb(single['zero'])/single['denominator']))
            else:
                # Condition on the other q-1 windows. At most q-1 choices
                # are excluded from the one-window average. Shift the target
                # by their feedback and pay their output-weight perturbation.
                factor=arb(data['windows'])/(data['windows']-q+1)*(((1+1/z)**4-1)/15)**(q-1)
                refined=up(factor*single['global_density']/single['denominator'])
                g=min(g,refined);cancellation=min(cancellation,refined)
        fresh=min(h[q],up((h2[q]*nonzero).sqrt()))
        fresh_mean=min(mean[q],up((mean2[q]*nonzero).sqrt()))
        n=3+len(sizes);rows=[[arb(0)]*n for _ in range(n)]
        rows[0][:3]=[zz,zn,zd]
        rows[1][:3]=[up(beta*fresh/L),up(alpha*h[q]),0]
        rows[2][:3]=[up(alpha*cancellation),0,up(alpha*g)]
        for j,size in enumerate(sizes):rows[1][3+j]=up(beta*h[q]*size/L)
        for i,(level,size,averages,averages2) in enumerate(zip(levels,sizes,class_mean,class_mean2)):
            moment=averages[q];moment2=averages2[q]
            if classes and not q:
                # Empty lazy steps preserve the expansion-weight class.
                rows[3+i][3+i]=up(alpha*z**level)
            else:
                if classes:
                    event=min(up(atom),up(nonzero/size))
                    uniform_cancel=min(moment,up(z**max(0,level-4*q)*event),up((moment2*event).sqrt()))
                    class_g=up(sum((p*z**max(0,level-w) for w,p in enumerate(weights)),arb(0)))
                    if single is not None and q:
                        single_g=up(arb(single['density'])/single['denominator']) if q==1 else refined
                        class_g=min(class_g,single_g)
                    uniform_density=up(class_g/size)
                else:uniform_cancel=up(cancellation/L);uniform_density=up(g/L)
                refresh_moment=min(moment,up((moment2*nonzero).sqrt()))
                rows[3+i][:3]=[up(alpha*uniform_cancel+beta*refresh_moment/L),up(alpha*moment),up(alpha*uniform_density)]
            for j,target_size in enumerate(sizes):
                rows[3+i][3+j]=up(rows[3+i][3+j]+beta*moment*target_size/L)
        matrix=arb_mat(rows)
        if any(matrix[i,j]<0 for i in range(n) for j in range(n)):
            raise ArithmeticError('negative density envelope entry')
        result.append(matrix)
        next_weights=[arb(0)]*(len(weights)+4)
        for w,p in enumerate(weights):
            for b in range(1,5):next_weights[w+b]+=p*comb(4,b)/15
        weights=next_weights
    return result


def outward(data,tilt,single=None,*,classes=False):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,(-aq(tilt)).exp(),single,classes=classes)


def build_operators(args,refined=False,classes=False):
    from occupancy_model import placement
    from occupancy_memory import rounded
    data=actual(2,getattr(args,'exact_feedback',False));ctx.prec=args.precision;result={}
    for tilt in args.tilts:
        single=None
        if refined:
            import single_packet
            single=single_packet.actual((-aq(Q(tilt))).exp());ctx.prec=args.precision
            print('GF16 one-packet density census',tilt,single,flush=True)
        local=outward(data,tilt,single,classes=classes)
        exact=placement(local,rounding=rounded,maximum_groups=args.groups)
        n=local[0].nrows()
        arrays=[np.array([[float(matrix[i,j]) for j in range(n)] for i in range(n)]) for matrix in exact]
        result[tilt,'1']=exact,arrays
        print('GF16 density operators',tilt,'degree',args.groups,flush=True)
    return result
