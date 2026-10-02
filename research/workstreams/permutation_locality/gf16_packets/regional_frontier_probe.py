"""Counterfactual regional diagnostics; never distance certificates.

Delete selected contributions from the bound at fixed saved witnesses.
These deliberately optimistic matrices need not bound any actual encoder.
The optional exact-return comparison instead refines valid local bounds;
its whole-length evaluation is still floating and is not a certificate.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from flint import ctx
import birth_classes
import frontier_probe
import occupancy_birth_classes as occupancy
import regional_count as regional
import regional_count_probe as floating
import scalar_cover as sc
import variance_partition as variance


scaled_placement=floating.scaled_placement


def evaluate(local,parts,offset,epochs=64,regions=256):
    """Floating ordered regional placement under fixed count envelopes."""
    placed,scales=scaled_placement(local,epochs)
    terms=[]
    for log_weights,count in parts:
        matrix,shift=floating.weighted_region(placed,log_weights+scales)
        terms.append(count+sc.log_power(matrix,regions)+regions*shift)
    return float((logsumexp(terms)+offset)/np.log(2))


def variants(local,data):
    rows=[frontier_probe.counterfactual_matrices(matrix,data) for matrix in local]
    result={'baseline':local}
    for name in rows[0]:result[name]=np.array([row[name] for row in rows])
    quiet=local.copy();quiet[:,0,1:]=0
    result['no_births_from_zero']=quiet
    combined=result['no_returns'].copy();combined[:,1:,1]=0
    result['no_returns_or_lazy_mass_to_arbitrary']=combined
    empty=np.zeros_like(local);empty[0,0,0]=1
    result['entirely_empty_input_path']=empty
    return result


def weighted_parts(model,cell,witness,parts,checked):
    """Reconstruct this witness's count measure, including direct masses."""
    direct=witness.get('regional_direct_counts',False)
    weights,_=model.weights(cell,model.tilt);scale=Q(1) if direct else sum(weights)
    log_masses=(-np.arange(sc.G+1)*sc.logq(model.tilt) if direct else
                np.array([float(m.log()) for m in regional.binomial_masses(sc.G,weights[1]/scale)]))
    _,_,logs=model.family(model.tilt);fs=np.array(list(map(float,model.features)));active=np.array(model.active)
    weighted=[];without_ratios=[]
    for (interval,dual),part in zip(parts,checked['regional_count_parts']):
        cap=variance.factor(model,cell,interval[0],witness['variance_dual'])
        masses=(regional.count_mass_caps if direct else regional.count_ratios)(
            model.features,model.active,cell,interval,model.q_min,sc.G,regional.aq(cap),part['mgf_witnesses'])
        eta,mu,gamma=map(float,dual)
        count=sc.G*logsumexp(logs+eta*fs+mu*active+gamma*fs*(1-fs))
        count-=sc.G*(min(eta*float(x) for x in cell)+min(gamma*float(v) for v in interval))+mu*model.q_min
        weighted.append((log_masses+np.array([float(r.log()) for r in masses]),count))
        if not direct:without_ratios.append((log_masses,count))
    return weighted,without_ratios,scale


def probe(source,row,data):
    precision=ctx.prec
    cell=tuple(map(Q,row['cell']));witness=row['witness']
    model=sc.Model(sc.actual_components(129,Q(21,50),Q(1001,1000),True),data,
        source['threshold'],source['minimum_groups'],tilt=Q(witness['tilt']),inner=birth_classes,
        variance_shuffle=True,variance_bins=len(witness['regional_count_parts']),regional_count=True)
    ctx.prec=precision
    parts,checked=regional.prepare_witness(model,cell,witness)
    lam=Q(witness['parameters'][0]);matrices=regional.local_operators(model,witness);size=matrices[0].nrows()
    local=np.array([[[float(m[i,j]) for j in range(size)] for i in range(size)] for m in matrices])
    weighted,without_ratios,scale=weighted_parts(model,cell,witness,parts,checked)
    offset=sc.PACKETS*sc.logq(scale)+float(lam)*model.threshold
    options={'baseline':local} if data.get('joint_only',False) else variants(local,data)
    if 'regional_lazy_density_through' in witness:
        # The U entry now includes lazy density, so U/L no longer
        # isolates only freshly refreshed returns. Omit that ablation.
        options.pop('no_lazy_returns',None)
    refined_source=any(key in witness for key in ('regional_joint_return_through','regional_lazy_density_through'))
    if refined_source or 'joint_return_census' in data or data.get('lazy_density_through'):
        # Alternative allocations start from the unrefined operator,
        # never from a matrix whose U entry already contains lazy mass.
        original=occupancy.outward(data,lam)
        if refined_source:
            options['without_local_refinements']=np.array([
                [[float(m[i,j]) for j in range(size)] for i in range(size)] for m in original])
    if 'joint_return_census' in data:
        import fixed_joint_return_probe as joint
        exact=joint.refine(data,original,data['joint_return_census'],(-regional.aq(lam)).exp())
        options['exact_small_occupancy_returns']=joint.arrays(exact)
    if data.get('lazy_density_through'):
        import lazy_density_probe as lazy
        z=(-regional.aq(lam)).exp();single=None
        if data.get('lazy_density_one_packet'):
            import single_packet
            single=single_packet.census(data['lazy_density_images'],data['columns'],z)
        caps=lazy.density_caps(data,z,single)
        for through in data['lazy_density_through']:
            options[f'lazy_density_through_{through}']=lazy.arrays(lazy.candidate(data,original,z,caps,through))
            if 'joint_return_census' in data:
                options[f'joint_and_lazy_through_{through}']=lazy.arrays(lazy.candidate(data,exact,z,caps,through))
    if data.get('feedback_density_starts'):
        import fiber_density
        import lazy_density_probe as lazy
        z=(-regional.aq(lam)).exp();caps=fiber_density.profile_caps(data,z)
        replace_uniform=data.get('feedback_uniform_replace',False)
        if replace_uniform:
            source_witness={key:value for key,value in witness.items() if key not in (
                'regional_lazy_density_through','regional_feedback_classes_from','regional_feedback_classes_through',
                'regional_feedback_uniform_classes','regional_feedback_uniform_replace')}
            source_operators=regional.local_operators(model,source_witness)
        for start in data['feedback_density_starts']:
            allocation=data.get('feedback_density_allocation','density')
            value=fiber_density.candidate(
                data,matrices,z,start,data['feedback_density_through'],caps,allocation,
                include_uniform=data.get('feedback_uniform_classes',False) and not replace_uniform)
            if replace_uniform:value=fiber_density.replace_uniform_classes(
                data,value,source_operators,z,start,data['feedback_density_through'],caps)
            suffix='_uniform_replace' if replace_uniform else ''
            options[f'feedback_{allocation}_from_{start}{suffix}']=lazy.arrays(value)
    scores={}
    for name,value in options.items():
        scores[name]=evaluate(value,weighted,offset)
        print('REGIONAL FRONTIER',row['activity'],name,scores[name],flush=True)
    if abs(scores['baseline']-row['proposal'])>1e-5:
        raise ArithmeticError('fixed-witness reconstruction does not match the saved diagnostic')
    if not data.get('joint_only',False) and without_ratios:
        scores['count_ratios_replaced_by_one']=evaluate(local,without_ratios,offset)
    return dict(activity=row['activity'],cell=row['cell'],tilt=witness['parameters'][0],scores=scores)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path,help='Saved regional-count floating screen; used only as diagnostic input')
    parser.add_argument('--activities',nargs='+',default=['.175','.2'])
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--joint-return-through',type=int,choices=(0,1,2,3,4),default=0)
    parser.add_argument('--joint-only',action='store_true',help='Compare only baseline and checked local refinements, omitting deletions')
    parser.add_argument('--lazy-density-through',nargs='+',type=int,default=[])
    parser.add_argument('--lazy-density-one-packet',action='store_true')
    parser.add_argument('--feedback-density-starts',nargs='+',type=int,default=[])
    parser.add_argument('--feedback-density-through',type=int,default=32)
    parser.add_argument('--feedback-density-allocation',choices=('density','classes'),default='density')
    parser.add_argument('--feedback-uniform-classes',action='store_true',help='Also retain class masses from U when its lazy branch still enters M')
    parser.add_argument('--feedback-uniform-replace',action='store_true',help='Rebuild complete U rows before prior lazy-density allocations')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();source=json.loads(args.input.read_text())
    if args.precision<128:parser.error('precision >=128 required')
    if args.joint_only and not (args.joint_return_through or args.lazy_density_through or args.feedback_density_starts):parser.error('--joint-only requires a local refinement')
    if any(not 0<=j<=32 for j in args.lazy_density_through):parser.error('lazy-density cutoff must lie in 0..32')
    if args.lazy_density_one_packet and not args.lazy_density_through:parser.error('one-packet density requires a lazy-density refinement')
    if any(not 1<=start<=args.feedback_density_through<=32 for start in args.feedback_density_starts):
        parser.error('feedback-density intervals must lie in 1..32')
    if args.feedback_uniform_classes and (args.feedback_density_allocation!='classes' or not args.feedback_density_starts):
        parser.error('--feedback-uniform-classes requires a class allocation interval')
    if args.feedback_uniform_replace and not args.feedback_uniform_classes:
        parser.error('--feedback-uniform-replace requires --feedback-uniform-classes')
    if source.get('schema')!='gf16-regional-count-screen-1' or not source.get('diagnostic_only'):
        parser.error('a regional floating screen is required; certificates are not imported')
    data=birth_classes.actual(source['updates']);ctx.prec=args.precision;rows=[]
    data['joint_only']=args.joint_only
    data['lazy_density_through']=args.lazy_density_through
    data['lazy_density_one_packet']=args.lazy_density_one_packet
    data['feedback_density_starts']=args.feedback_density_starts
    data['feedback_density_through']=args.feedback_density_through
    data['feedback_density_allocation']=args.feedback_density_allocation
    data['feedback_uniform_classes']=args.feedback_uniform_classes
    data['feedback_uniform_replace']=args.feedback_uniform_replace
    if args.feedback_density_starts:
        import fiber_density
        data=fiber_density.attach(data)
    if args.lazy_density_one_packet:
        images,columns,_=sc.kernel.maps()
        if columns!=data['columns']:raise ArithmeticError('feedback map mismatch')
        data['lazy_density_images']=images
    if args.joint_return_through:
        import return_moment
        images,columns,_=sc.kernel.maps()
        if columns!=data['columns']:raise ArithmeticError('feedback map mismatch')
        data['joint_return_census']=return_moment.census(images,columns,data['bits'],args.joint_return_through,verbose=True)
        ctx.prec=args.precision
        print('EXACT JOINT RETURN census checked through',args.joint_return_through,flush=True)
    for p in map(Q,args.activities):
        matching=[row for row in source['rows'] if Q(row['activity'])==p]
        if not matching:parser.error('requested activity is missing')
        row=min(matching,key=lambda row:row['proposal'])
        rows.append(probe(source,row,data))
        if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,
            fixed_witness_counterfactuals=not args.joint_only,joint_return_through=args.joint_return_through,
            lazy_density_through=args.lazy_density_through,lazy_density_one_packet=args.lazy_density_one_packet,
            feedback_density_starts=args.feedback_density_starts,feedback_density_through=args.feedback_density_through,
            feedback_density_allocation=args.feedback_density_allocation,
            feedback_uniform_classes=args.feedback_uniform_classes,
            feedback_uniform_replace=args.feedback_uniform_replace,
            updates=source['updates'],threshold=source['threshold'],precision=args.precision,
            source=str(args.input),rows=rows),indent=2,allow_nan=False)+'\n')
    print('Floating local-refinement comparison; not a whole-code certificate.' if args.joint_only else
          'Counterfactual contributions at fixed witnesses; not bounds on an altered code.',flush=True)


if __name__=='__main__':main()
