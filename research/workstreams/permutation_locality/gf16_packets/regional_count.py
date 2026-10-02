"""Outward ingredients for count-sensitive regional shuffle comparisons.

Positive polynomial matrix products retain exact regional occupancy.
Finite-feature MGF duals bound Poisson-binomial atoms over mean/variance
intervals. A whole-code certificate still requires complete outer coverage.
"""
from fractions import Fraction as Q
from math import comb
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.special import logsumexp
from flint import arb,arb_mat,arb_poly,ctx
import scalar_cover as sc

aq,up=sc.kernel.aq,sc.kernel.up
TILTS=tuple(map(Q,('-16','-8','-4','-2','-1','-1/2','-1/4','0','1/4','1/2','1','2','4','8','16')))
FINE_TILTS=tuple(sorted(set(TILTS).union(Q(j,8) for j in range(-16,17))))


def polynomial_product(left,right):
    n=len(left);result=[]
    for row in left:
        output=[]
        for j in range(n):
            value=arb_poly([])
            for k in range(n):
                if row[k].degree()>=0 and right[k][j].degree()>=0:value+=row[k]*right[k][j]
            output.append(value)
        result.append(output)
    return result


def placement(local,epochs=64):
    """Coefficients of the ordered local matrix polynomial, via squaring."""
    if type(epochs) is not int or epochs<1 or len(local)<2:raise ValueError('positive placement geometry required')
    W=len(local)-1;n=local[0].nrows()
    if any(m.nrows()!=n or m.ncols()!=n for m in local):raise ValueError('matching square local matrices required')
    if any(m[i,j]<0 for m in local for i in range(n) for j in range(n)):
        raise ValueError('nonnegative local operators required')
    power=[[arb_poly([comb(W,k)*m[i,j] for k,m in enumerate(local)]) for j in range(n)] for i in range(n)]
    result=None;remaining=epochs
    while remaining:
        if remaining&1:result=power if result is None else polynomial_product(result,power)
        remaining>>=1
        if remaining:power=polynomial_product(power,power)
    groups=epochs*W
    return [arb_mat([[up(result[i][j][k]/comb(groups,k)) for j in range(n)] for i in range(n)])
            for k in range(groups+1)]


def propose_mgf(features,active,cell,interval,minimum_groups,groups,tilted_atom=False,tilts=TILTS,tilted_variance=False):
    """The LP supplies only dual coefficients; outward evaluation checks them."""
    fs=np.array(list(map(float,features)));vs=fs*(1-fs);a=np.array(active,dtype=float)
    mean=float(sum(cell)/2);lo,hi=map(float,interval);witnesses=[]
    for tilt in tilts:
        with np.errstate(divide='ignore'):
            g=np.logaddexp(np.log1p(-fs),np.log(fs)+float(tilt))
        fit=linprog(-g,A_ub=np.array([vs,-vs,-a]),b_ub=[hi,-lo,-minimum_groups/groups],
            A_eq=np.array([np.ones(len(fs)),fs]),b_eq=[1.,mean],bounds=(0,None),method='highs')
        if not fit.success:continue
        v=fit.ineqlin.marginals;e=fit.eqlin.marginals
        coefficients=(-e[1],-v[0]+v[1],min(0.,v[2]))
        dual=[Q(round(float(x)*10**9),10**9) for x in coefficients]
        row=dict(tilt=str(tilt),dual=list(map(str,dual)))
        if tilted_atom:row['tilted_atom']=True
        if tilted_variance:
            exponential=np.exp(float(tilt));values=vs*exponential/(1-fs+fs*exponential)**2
            lower=linprog(values,A_ub=np.array([vs,-vs,-a]),b_ub=[hi,-lo,-minimum_groups/groups],
                A_eq=np.array([np.ones(len(fs)),fs]),b_eq=[1.,mean],bounds=(0,None),method='highs')
            if lower.success:
                v=lower.ineqlin.marginals;e=lower.eqlin.marginals
                slopes=(e[1],v[0]-v[1],max(0.,-v[2]))
                row['tilted_variance_dual']=[str(Q(round(float(x)*10**9),10**9)) for x in slopes]
        witnesses.append(row)
    return witnesses


def mgf_upper(features,active,cell,interval,minimum_groups,groups,witness):
    """Check an affine majorant of every finite-feature log MGF."""
    if (len(features)!=len(active) or not features or any(not 0<=f<=1 for f in features)
            or any(a not in (0,1) for a in active) or type(groups) is not int or groups<1
            or not 0<=minimum_groups<=groups or not 0<cell[0]<=cell[1]<1
            or not 0<=interval[0]<=interval[1]<=Q(1,4) or len(witness['dual'])!=3):
        raise ValueError('valid feature family, count, and interior moment intervals required')
    tilt=Q(witness['tilt']);b,c,d=map(Q,witness['dual'])
    if d>0:raise ValueError('occupancy coefficient must be nonpositive')
    exponential=aq(tilt).exp()
    intercept=max(up((aq(1-f)+aq(f)*exponential).log()-aq(b*f+c*f*(1-f)+d*a))
        for f,a in zip(features,active))
    return up(groups*(intercept+aq(max(b*x for x in cell)+max(c*v for v in interval)))+aq(d*minimum_groups))


def binomial_masses(groups,p):
    p=Q(p)
    if type(groups) is not int or groups<0 or not 0<p<1:raise ValueError('interior binomial probability required')
    values=[aq(1-p)**groups];ratio=aq(p/(1-p))
    for j in range(groups):values.append(values[-1]*ratio*(groups-j)/(j+1))
    return values


def count_ratios(features,active,cell,interval,minimum_groups,groups,cap,witnesses):
    """Outward PB/binomial caps, uniform throughout the two intervals."""
    if not cap>0:raise ValueError('positive baseline ratio cap required')
    bounds=[up(cap) for _ in range(groups+1)]
    lower=[min(arb(a.lower()),arb(b.lower())) for a,b in zip(
        binomial_masses(groups,cell[0]),binomial_masses(groups,cell[1]))]
    if any(not x>0 for x in lower):raise ArithmeticError('binomial denominator is not strictly positive')
    # Multiplicative recurrences avoid an exponential per count per dual.
    for witness in witnesses:
        moment=mgf_upper(features,active,cell,interval,minimum_groups,groups,witness).exp()
        enabled=witness.get('tilted_atom',False)
        if type(enabled) is not bool:raise ValueError('boolean tilted-atom option required')
        if enabled:
            variance_lower=(tilted_variance_lower(features,active,cell,interval,minimum_groups,groups,witness)
                            if 'tilted_variance_dual' in witness else arb(0))
            moment*=tilted_atom_upper(features,interval[0],groups,Q(witness['tilt']),variance_lower)
        ratio=(-aq(Q(witness['tilt']))).exp()
        for j in range(groups+1):
            bounds[j]=min(bounds[j],up(moment/lower[j]))
            moment*=ratio
    return bounds


def binomial_interval_upper(groups,cell):
    """Maximum binomial atom over an interior interval of probabilities."""
    if type(groups) is not int or groups<1 or not 0<cell[0]<=cell[1]<1:
        raise ValueError('positive count and interior mean interval required')
    left=binomial_masses(groups,cell[0]);right=binomial_masses(groups,cell[1]);result=[]
    for j in range(groups+1):
        mode=Q(j,groups)
        if mode<=cell[0]:value=left[j]
        elif mode>=cell[1]:value=right[j]
        else:value=comb(groups,j)*aq(mode)**j*aq(1-mode)**(groups-j)
        result.append(up(value))
    return result


def count_mass_caps(features,active,cell,interval,minimum_groups,groups,cap,witnesses):
    """Bound count probabilities directly, without a reference denominator."""
    if not cap>0:raise ValueError('positive baseline ratio cap required')
    # Validate the domain even when the finite list of MGF witnesses is empty.
    mgf_upper(features,active,cell,interval,minimum_groups,groups,{'tilt':'0','dual':['0','0','0']})
    bounds=[min(arb(1),up(cap*mass)) for mass in binomial_interval_upper(groups,cell)]
    for witness in witnesses:
        moment=mgf_upper(features,active,cell,interval,minimum_groups,groups,witness).exp()
        enabled=witness.get('tilted_atom',False)
        if type(enabled) is not bool:raise ValueError('boolean tilted-atom option required')
        if enabled:
            lower=(tilted_variance_lower(features,active,cell,interval,minimum_groups,groups,witness)
                   if 'tilted_variance_dual' in witness else arb(0))
            moment*=tilted_atom_upper(features,interval[0],groups,Q(witness['tilt']),lower)
        ratio=(-aq(Q(witness['tilt']))).exp()
        for j in range(groups+1):
            bounds[j]=min(bounds[j],up(moment));moment*=ratio
    return bounds


def tilted_atom_upper(features,variance_lower,groups,tilt,checked_lower=arb(0)):
    """Bound any atom of the exponentially tilted Bernoulli count."""
    if (not features or any(not 0<=f<=1 for f in features)
            or not 0<=variance_lower<=Q(1,4) or type(groups) is not int or groups<1):
        raise ValueError('valid features and average variance lower bound required')
    exponential=aq(tilt).exp()
    ratios=[exponential/(aq(1-f)+aq(f)*exponential)**2 for f in features if 0<f<1]
    # A universal fallback is valid even for an empty set of interior types.
    factor=min(arb(r.lower()) for r in ratios) if ratios else (-abs(aq(tilt))).exp()
    lower=max(arb(0),arb((aq(groups*variance_lower)*factor).lower()),arb(checked_lower.lower()))
    if not lower:return arb(1)
    return min(arb(1),up((-lower).exp()*lower.bessel_i(0)))


def tilted_variance_lower(features,active,cell,interval,minimum_groups,groups,witness):
    """Check an affine minorant for variance after the count MGF tilt."""
    # Reuse the MGF domain checks; its value is not used in this minorant.
    mgf_upper(features,active,cell,interval,minimum_groups,groups,witness)
    if len(witness['tilted_variance_dual'])!=3:raise ValueError('three tilted variance slopes required')
    b,c,d=map(Q,witness['tilted_variance_dual'])
    if d<0:raise ValueError('nonnegative active-count coefficient required for a variance minorant')
    exponential=aq(Q(witness['tilt'])).exp()
    intercept=min(arb((aq(f*(1-f))*exponential/(aq(1-f)+aq(f)*exponential)**2
                      -aq(b*f+c*f*(1-f)+d*a)).lower()) for f,a in zip(features,active))
    value=groups*(intercept+aq(min(b*x for x in cell)+min(c*v for v in interval)))+aq(d*minimum_groups)
    return max(arb(0),arb(value.lower()))


def feedback_class_interval(witness,windows):
    """Validate an optional complete GF-feedback class allocation selector."""
    keys=('regional_feedback_classes_from','regional_feedback_classes_through')
    uniform=witness.get('regional_feedback_uniform_classes',False)
    replace=witness.get('regional_feedback_uniform_replace',False)
    if type(uniform) is not bool or type(replace) is not bool or replace and not uniform:
        raise ValueError('boolean uniform-source options required; replacement requires class allocation')
    if not any(key in witness for key in keys):
        if uniform:raise ValueError('uniform-source refinement requires a feedback-class interval')
        return None
    if not all(key in witness for key in keys):raise ValueError('both feedback-class endpoints required')
    start,through=(witness[key] for key in keys)
    if type(start) is not int or type(through) is not int or not 1<=start<=through<=windows:
        raise ValueError('valid positive feedback-class packet interval required')
    return start,through


def prepare_witness(model,cell,witness):
    """Propose missing MGF duals, or validate the saved partition exactly."""
    import variance_partition as variance
    if (Q(witness['tilt'])!=model.tilt or model.data['windows']!=32
            or not model.variance_shuffle or not model.variance_bins):
        raise ValueError('base-tilt variance partition and actual 32-packet steps required')
    if type(witness.get('regional_direct_counts',False)) is not bool:
        raise ValueError('boolean direct-count option required')
    if type(witness.get('regional_exact_zero',False)) is not bool:
        raise ValueError('boolean exact-zero option required')
    feedback_class_interval(witness,model.data['windows'])
    if 'regional_lazy_density_through' in witness:
        cutoff=witness['regional_lazy_density_through']
        if type(cutoff) is not int or not 0<=cutoff<=model.data['windows']:
            raise ValueError('valid integer lazy-density occupancy cutoff required')
    if 'regional_joint_return_through' in witness:
        cutoff=witness['regional_joint_return_through']
        if type(cutoff) is not int or not 0<=cutoff<=min(4,model.data['windows']):
            raise ValueError('valid integer joint-return occupancy cutoff required')
    parts=variance.validate(cell,witness['variance_partition'])
    saved=witness.get('regional_count_parts')
    if saved is None:
        atom=witness.get('regional_tilted_atom',False)
        fine=witness.get('regional_fine_tilts',False)
        tilted_variance=witness.get('regional_tilted_variance',False)
        if (any(type(v) is not bool for v in (atom,fine,tilted_variance)) or tilted_variance and not atom):
            raise ValueError('boolean regional refinement options required; tilted variance requires the atom bound')
        saved=[dict(interval=list(map(str,interval)),dual=list(map(str,dual)),
                    mgf_witnesses=propose_mgf(model.features,model.active,cell,interval,model.q_min,sc.G,
                        atom,FINE_TILTS if fine else TILTS,tilted_variance))
               for interval,dual in parts]
    if not isinstance(saved,list) or len(saved)!=len(parts):
        raise ValueError('one MGF family for every variance part required')
    for part,(interval,dual) in zip(saved,parts):
        if (tuple(map(Q,part['interval']))!=interval or tuple(map(Q,part['dual']))!=dual
                or not isinstance(part['mgf_witnesses'],list) or len(part['mgf_witnesses'])>128):
            raise ValueError('MGF witnesses must match the complete variance partition')
    return parts,dict(witness,regional_count_parts=saved)


def proposal_cache(model):
    """Bounded in-memory search scratch; never consulted by outward()."""
    context=(id(model.data),ctx.prec)
    saved=getattr(model,'regional_proposal_scratch',None)
    if saved is None or saved[0]!=context:
        saved=(context,{})
        model.regional_proposal_scratch=saved
    return saved[1]


def local_operators(model,witness,*,_proposal_scratch=None):
    import occupancy_birth_classes
    feedback=feedback_class_interval(witness,model.data['windows'])
    lam=Q(witness['parameters'][0])
    exact_zero=witness.get('regional_exact_zero',False)
    if type(exact_zero) is not bool:raise ValueError('boolean exact-zero option required')
    key=(lam,exact_zero)
    saved=None if _proposal_scratch is None else _proposal_scratch.get('local_base')
    if saved is not None and saved[0]==key:local=saved[1]
    else:
        local=occupancy_birth_classes.outward(model.data,lam)
        if exact_zero:local=occupancy_birth_classes.refine_zero(model.data,local,(-aq(lam)).exp())
        if _proposal_scratch is not None:_proposal_scratch['local_base']=(key,local)
    if 'regional_joint_return_through' in witness:
        import return_moment
        through=witness['regional_joint_return_through']
        # Integer counts do not depend on lambda. Each fresh model checks
        # them once; no census or acceptance label is read from a receipt.
        if not hasattr(model,'joint_return_cache'):model.joint_return_cache={}
        if through not in model.joint_return_cache:
            model.joint_return_cache[through]=return_moment.actual_census(model.data,through)
        local=return_moment.refine_class_returns(model.data,local,model.joint_return_cache[through],(-aq(lam)).exp())
    feedback_source=local
    if 'regional_lazy_density_through' in witness:
        import lazy_density
        through=witness['regional_lazy_density_through'];z=(-aq(lam)).exp()
        if _proposal_scratch is None:
            local=lazy_density.refine_actual(model.data,local,z,through)
        else:
            key=(lam,through)
            saved=_proposal_scratch.get('lazy_caps')
            if saved is None or saved[0]!=key:
                saved=(key,lazy_density.actual_caps(model.data,z,through))
                _proposal_scratch['lazy_caps']=saved
            local=lazy_density.candidate(model.data,local,z,saved[1],through)
    if feedback is not None:
        import fiber_density
        # Only integer profile counts are reused within this model.
        # The weighted fiber bounds are rebuilt at the current precision
        # and output tilt. No numerical cap is loaded from a receipt.
        if not hasattr(model,'fiber_density_data'):
            model.fiber_density_data=fiber_density.attach(model.data)
        uniform=witness.get('regional_feedback_uniform_classes',False)
        replace=witness.get('regional_feedback_uniform_replace',False)
        options={'include_uniform':True} if uniform and not replace else {}
        local=fiber_density.candidate(model.fiber_density_data,local,(-aq(lam)).exp(),
            *feedback,allocation='classes',**options)
        if replace:local=fiber_density.replace_uniform_classes(model.fiber_density_data,local,
            feedback_source,(-aq(lam)).exp(),*feedback)
    return local


def proposal_count_weights(model,cell,witness,parts,checked,scratch):
    """Reuse count caps across inner/tilt trials for exactly the same cell."""
    import variance_partition as variance
    direct=witness.get('regional_direct_counts',False)
    # These weights do not depend on the output tilt or the inner envelope.
    # Retain one family only; sibling cells, changed duals, geometry, and
    # precision rebuild it. No numerical scratch is written to a receipt.
    key=(ctx.prec,sc.G,model.q_min,model.tilt,tuple(model.features),tuple(model.active),tuple(cell),direct,
         tuple(witness['variance_dual']),json.dumps(checked['regional_count_parts'],sort_keys=True))
    saved=scratch.get('count_weights')
    if saved is not None and saved[0]==key:return saved[1]
    weights,_=model.weights(cell,model.tilt);scale=Q(1) if direct else sum(weights)
    log_masses=(-np.arange(sc.G+1)*sc.logq(model.tilt) if direct else
                np.array([float(m.log()) for m in binomial_masses(sc.G,weights[1]/scale)]))
    result=[]
    for (interval,dual),part in zip(parts,checked['regional_count_parts']):
        cap=variance.factor(model,cell,interval[0],witness['variance_dual'])
        ratios=(count_mass_caps if direct else count_ratios)(model.features,model.active,cell,interval,
                    model.q_min,sc.G,aq(cap),part['mgf_witnesses'])
        result.append(log_masses+np.array([float(r.log()) for r in ratios]))
    value=(scale,result);scratch['count_weights']=(key,value)
    return value


def propose(model,cell,witness):
    """Whole-cell floating proposal; the saved rational duals are replayable."""
    import regional_count_probe as floating
    parts,checked=prepare_witness(model,cell,witness)
    lam=Q(witness['parameters'][0])
    scratch=proposal_cache(model)
    local=local_operators(model,witness,_proposal_scratch=scratch);size=local[0].nrows()
    arrays=np.array([[[float(m[i,j]) for j in range(size)] for i in range(size)] for m in local])
    region,region_logs=floating.scaled_placement(arrays,64)
    scale,count_weights=proposal_count_weights(model,cell,witness,parts,checked,scratch)
    _,_,logs=model.family(model.tilt);features=np.array(list(map(float,model.features)))
    active=np.array(model.active);terms=[]
    for (interval,dual),log_weights in zip(parts,count_weights):
        matrix,shift=floating.weighted_region(region,log_weights+region_logs)
        moment=sc.log_power(matrix,sc.REGIONS)+sc.REGIONS*shift
        eta,mu,gamma=map(float,dual)
        count=sc.G*logsumexp(logs+eta*features+mu*active+gamma*features*(1-features))
        count-=sc.G*(min(eta*float(x) for x in cell)+min(gamma*float(v) for v in interval))+mu*model.q_min
        terms.append(count+moment)
    score=float((logsumexp(terms)+sc.PACKETS*sc.logq(scale)+float(lam)*model.threshold)/np.log(2))
    if not np.isfinite(score):raise ArithmeticError('nonfinite regional proposal')
    return score,checked


def outward(model,cell,witness):
    """Check one whole mean interval and every variance part, not all means."""
    import variance_partition as variance
    parts,checked=prepare_witness(model,cell,witness)
    lam=Q(witness['parameters'][0])
    if lam<=0:raise ValueError('positive output tilt required')
    precision=ctx.prec
    print('REGIONAL COUNT constructing ordered polynomial operators',str(lam),'precision',precision,flush=True)
    region=placement(local_operators(model,witness))
    if len(region)!=sc.G+1:raise ArithmeticError('regional geometry mismatch')
    direct=witness.get('regional_direct_counts',False)
    weights,_=model.weights(cell,model.tilt);scale=Q(1) if direct else sum(weights)
    masses=([aq(model.tilt)**(-j) for j in range(sc.G+1)] if direct else
            binomial_masses(sc.G,weights[1]/scale))
    cs,_,_=model.family(model.tilt)
    total=arb(0);size=region[0].nrows()
    for index,((interval,dual),part) in enumerate(zip(parts,checked['regional_count_parts'])):
        eta,mu,gamma=dual
        mgfs=part['mgf_witnesses']
        cap=variance.factor(model,cell,interval[0],witness['variance_dual'])
        ratios=(count_mass_caps if direct else count_ratios)(model.features,model.active,cell,interval,
                    model.q_min,sc.G,aq(cap),mgfs)
        matrix=sum((mass*ratio*value for mass,ratio,value in zip(masses,ratios,region)),arb_mat(size,size))
        power=matrix**sc.REGIONS;moment=sum((power[0,j] for j in range(size)),arb(0))
        count=sum((aq(c)*aq(eta*f+mu*a+gamma*f*(1-f)).exp()
            for c,f,a in zip(cs,model.features,model.active)),arb(0))
        exponent=sc.G*(min(eta*x for x in cell)+min(gamma*v for v in interval))+mu*model.q_min
        term=(sc.G*count.log()-aq(exponent)+moment.log()).exp()
        total=up(total+term)
        print('REGIONAL COUNT variance part',index+1,'/',len(parts),'checked',flush=True)
    upper=up((total.log()+sc.PACKETS*aq(scale).log()+aq(lam)*model.threshold).exp())
    if ctx.prec!=precision:raise ArithmeticError('working precision changed during regional replay')
    return upper,checked


def main():
    import birth_classes
    import variance_partition as variance
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distance',default='37/400')
    parser.add_argument('--updates',type=int,choices=(2,3,4),default=2,help='Actual independent inner updates; not a proof-only knob or a production default')
    parser.add_argument('--minimum-groups',type=int,default=97,help='First active-group occupancy included in the dense comparison')
    parser.add_argument('--activities',nargs='+',default=['.15','.175','.2'])
    parser.add_argument('--radius',default='1/1000000')
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--variance-bins',type=int,default=8)
    parser.add_argument('--tilted-atom',action='store_true',help='Retain atom anti-concentration after the count MGF tilt')
    parser.add_argument('--fine-tilts',action='store_true',help='Add count MGF tilts at spacing 1/8 between -2 and 2')
    parser.add_argument('--tilted-variance',action='store_true',help='Use a checked finite-feature variance minorant after tilting')
    parser.add_argument('--direct-counts',action='store_true',help='Bound count masses directly, avoiding mean-cell input inflation')
    parser.add_argument('--exact-zero',action='store_true',help='Use the exact weighted Fourier count for zero-to-zero transitions')
    parser.add_argument('--lazy-density-through',type=int,help='Retain density through this packet occupancy, using a freshly checked one-packet census')
    parser.add_argument('--joint-return-through',type=int,choices=(0,1,2,3,4),help='Use freshly cross-checked exact joint return counts through this occupancy')
    parser.add_argument('--feedback-classes-from',type=int,help='Retain both lazy mass and GF-feedback density in outgoing state classes from this packet count')
    parser.add_argument('--feedback-classes-through',type=int,default=32)
    parser.add_argument('--feedback-uniform-classes',action='store_true')
    parser.add_argument('--feedback-uniform-replace',action='store_true')
    parser.add_argument('--tilt-scales',nargs='+',default=['1'],help='Positive scales of the proposed output Chernoff tilt')
    parser.add_argument('--screen-only',action='store_true',help='Floating proposals only; no outward interval certificate')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.precision<128 or not 0<Q(args.radius)<Q(1,100):parser.error('valid precision and small positive radius required')
    if not 1<=args.minimum_groups<=sc.G:parser.error('minimum groups must be in 1..2048')
    if args.lazy_density_through is not None and not 0<=args.lazy_density_through<=32:parser.error('lazy-density cutoff must lie in 0..32')
    if args.feedback_classes_from is not None and not 1<=args.feedback_classes_from<=args.feedback_classes_through<=32:
        parser.error('feedback-class packet interval must lie in 1..32')
    if args.feedback_uniform_classes and args.feedback_classes_from is None:
        parser.error('uniform-source refinement requires --feedback-classes-from')
    if args.feedback_uniform_replace and not args.feedback_uniform_classes:
        parser.error('uniform-source replacement requires --feedback-uniform-classes')
    if any(Q(scale)<=0 for scale in args.tilt_scales):parser.error('positive output tilt scales required')
    if args.tilted_variance and not args.tilted_atom:parser.error('--tilted-variance requires --tilted-atom')
    ctx.prec=args.precision
    model=sc.Model(sc.actual_components(129,Q(21,50),Q(1001,1000),True),birth_classes.actual(args.updates),
        int(Q(args.distance)*sc.N),args.minimum_groups,tilt=Q(3,16),inner=birth_classes,
        variance_shuffle=True,variance_bins=args.variance_bins)
    ctx.prec=args.precision
    rows=[]
    for p in map(Q,args.activities):
        x=model.tilt*p/(1-p+model.tilt*p);cell=(x-Q(args.radius),x+Q(args.radius))
        _,witness=variance.propose(model,cell,model.propose_with(cell,model.tilt))
        if args.tilted_atom:witness['regional_tilted_atom']=True
        if args.fine_tilts:witness['regional_fine_tilts']=True
        if args.tilted_variance:witness['regional_tilted_variance']=True
        if args.direct_counts:witness['regional_direct_counts']=True
        if args.exact_zero:witness['regional_exact_zero']=True
        if args.lazy_density_through is not None:witness['regional_lazy_density_through']=args.lazy_density_through
        if args.joint_return_through is not None:witness['regional_joint_return_through']=args.joint_return_through
        if args.feedback_classes_from is not None:
            witness.update(regional_feedback_classes_from=args.feedback_classes_from,
                regional_feedback_classes_through=args.feedback_classes_through)
        if args.feedback_uniform_classes:witness['regional_feedback_uniform_classes']=True
        if args.feedback_uniform_replace:witness['regional_feedback_uniform_replace']=True
        _,witness=prepare_witness(model,cell,witness)
        lam=Q(witness['parameters'][0])
        for scale in map(Q,args.tilt_scales):
            trial=dict(witness,parameters=[str(lam*scale),*witness['parameters'][1:]])
            row=dict(activity=str(p),cell=list(map(str,cell)),tilt_scale=str(scale))
            if args.screen_only:
                score,checked=propose(model,cell,trial)
                row.update(proposal=score,witness=checked)
                print('REGIONAL COUNT proposal',str(p),'tilt scale',str(scale),'log2',score,flush=True)
            else:
                upper,checked=outward(model,cell,trial);m,e=upper.upper().man_exp()
                row.update(upper=[int(m),int(e)],witness=checked)
                print('REGIONAL COUNT checked interval',list(map(str,cell)),'tilt scale',str(scale),
                      'log2 upper',upper.log()/arb(2).log(),flush=True)
            rows.append(row)
            if args.output:args.output.write_text(json.dumps(dict(
                schema='gf16-regional-count-screen-1' if args.screen_only else 'gf16-regional-count-selected-1',
                threshold=model.threshold,minimum_groups=model.q_min,precision=args.precision,updates=args.updates,
                partial_only=True,diagnostic_only=args.screen_only,tilted_atom=args.tilted_atom,
                fine_tilts=args.fine_tilts,tilted_variance=args.tilted_variance,direct_counts=args.direct_counts,
                exact_zero=args.exact_zero,
                lazy_density_through=args.lazy_density_through,joint_return_through=args.joint_return_through,
                feedback_classes_from=args.feedback_classes_from,feedback_classes_through=args.feedback_classes_through,
                feedback_uniform_classes=args.feedback_uniform_classes,feedback_uniform_replace=args.feedback_uniform_replace,
                rows=rows),indent=2)+'\n')
    print('Floating diagnostic only.' if args.screen_only else
          'Selected mean intervals only; not a complete dense or whole-code certificate.',flush=True)


if __name__=='__main__':main()
