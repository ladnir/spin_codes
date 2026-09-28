"""Outward replay of the opt-in mass-based lazy concentration bounds.

By default this replays selected homogeneous support vectors only. The
--full-cover mode instead uses the existing complete support-cover driver.
It retains the baseline grid and adds the requested refined alternatives.
An optional occupancy list reuses the local and region operators; every
occupancy still needs its own complete cover and outward replay.
Experimental operators use a separate optional exact local memo, never
the production operator memo. No encoder or setup distribution is changed.
"""
import argparse
from fractions import Fraction as Q
from types import SimpleNamespace
from flint import arb

from mass_density_screen import as_array, baseline
from mass_density import blend, conditioned_coefficients
from occupancy_memory import Z, C
from shell_cover import cover
from local_family import build as build_families


def requested_occupancies(groups,occupancies,full_cover):
    if occupancies is not None and not full_cover:
        raise ValueError('an occupancy batch requires --full-cover')
    result=tuple(occupancies) if occupancies is not None else (groups,)
    if (not result or len(set(result))!=len(result)
            or any(not isinstance(q,int) or isinstance(q,bool) or not 1<=q<=2048 for q in result)):
        raise ValueError('distinct occupancies in 1..2048 required')
    return result


def run_covers(args,occupancies,operators,counts,shells):
    results={}
    for q in occupancies:
        args.groups=q
        print('COMPLETE COVER with baseline and refined alternatives; q',q,flush=True)
        results[q]=cover(args,operators,counts,baseline.TAIL_TERMINAL,shells,prefix_rank=True)
    verified={q:value for q,value in results.items() if value is not None}
    if verified:
        total=baseline.up(sum(verified.values(),arb(0)))
        print('BATCH VERIFIED OCCUPANCIES',list(verified),'upper',total,
              'margin',-total.log()/arb(2).log(),flush=True)
    missing=[q for q,value in results.items() if value is None]
    if missing:
        print('BATCH OCCUPANCIES WITHOUT CERTIFICATE',missing,flush=True)
    return results


def alternative(base, census, tilt, penalty, start=6, end=12, through=6):
    if not 1 <= start <= end <= 32 or end >= len(base):
        raise ValueError('replacement occupancies must fit the local operators')
    candidate = conditioned_coefficients(census, through, max(through,end), tilt, penalty)
    fractions = {j:Q(1) for j in range(start,end+1)}
    zero = blend(base,candidate,fractions,target=Z)
    return blend(zero,candidate,fractions,target=C)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group_options=parser.add_mutually_exclusive_group()
    group_options.add_argument('--groups',type=int,default=64)
    group_options.add_argument('--occupancies',type=int,nargs='+',help='Complete-cover batch; each occupancy is checked separately')
    parser.add_argument('--tilt',default='.056')
    parser.add_argument('--tilts',nargs='+',help='Optional output-tilt grid; overrides --tilt')
    parser.add_argument('--penalty',default='.9')
    parser.add_argument('--supports',type=int,nargs='+',default=[192,200,208])
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--start',type=int,default=6)
    parser.add_argument('--end',type=int,default=12)
    parser.add_argument('--target-bits',type=int,default=52)
    parser.add_argument('--full-cover',action='store_true')
    parser.add_argument('--max-splits',type=int,default=400)
    parser.add_argument('--operator-cache',help='Optional exact-dyadic memo for unchanged baseline operators only')
    parser.add_argument('--epoch-cache',help='Optional source-bound exact memo for the experimental local families')
    feedback_options=parser.add_mutually_exclusive_group()
    feedback_options.add_argument('--spectral-feedback',type=int,default=0,choices=[0,7,8,9,10],help='Optional exact class counts and bounded feedback peaks beyond six windows')
    feedback_options.add_argument('--exact-feedback',type=int,default=0,choices=[0,7,8,9,10],help='Optional exact feedback peaks using guarded limb inversion')
    parser.add_argument('--joint-four',action='store_true',help='Enumerate joint lazy-return/output weights at four windows')
    parser.add_argument('--density-through',type=int,default=0,choices=[0,7,8,9,10],help='Extend positive density convolutions with outward integer counts')
    args=parser.parse_args()
    try:
        occupancies=requested_occupancies(args.groups,args.occupancies,args.full_cover)
    except ValueError as error:
        parser.error(str(error))
    maximum_groups=max(64,max(occupancies))
    tilts=args.tilts or [args.tilt]
    if (not 1 <= args.groups <= 2048 or any(Q(tilt)<=0 for tilt in tilts) or not 0<Q(args.penalty)<=1
            or args.precision<128 or not 1<=args.start<=args.end<=32
            or args.target_bits<1 or args.max_splits<0
            or any(not 38<=u<=256 for u in args.supports)):
        parser.error('invalid occupancy, witness, precision, support, or replacement range')
    baseline.folding_test(); baseline.geometry_test(); baseline.retained_test()
    baseline.placement_prefix_test(); baseline.mixing_test()
    caps=baseline.authenticated_caps()
    count_sets={}; shell_sets={}
    def measure(label,rho):
        count_sets[label]=baseline.integer_cdf(baseline.weighted_cdf_upper(caps,1<<128,full_weight=1/Q(rho)))
        shell_sets[label]=baseline.weighted_union_shells(caps,full_weight=1/Q(rho))
        shell_sets[label][0]-=1
    families=build_families(tilts,args.penalty,args.precision,max(16,args.end),
        exact_feedback=args.exact_feedback,spectral_feedback=args.spectral_feedback,
        joint_four=args.joint_four,density_through=args.density_through,directory=args.epoch_cache)
    label=f'{args.penalty}:mass-ZC-{args.start}-{args.end}'
    if args.spectral_feedback:
        label+=f':spectral-{args.spectral_feedback}'
    if args.exact_feedback:
        label+=f':exact-{args.exact_feedback}'
    if args.joint_four:
        label+=':joint-4'
    if args.density_through:
        label+=f':density-{args.density_through}'
    measure(label,args.penalty)
    operators={}
    fractions={j:Q(1) for j in range(args.start,args.end+1)}
    for tilt,(base,coefficients) in families.items():
        changed=blend(base,coefficients,fractions,target=Z)
        changed=blend(changed,coefficients,fractions,target=C)
        exact=baseline.placement(changed,rounding=baseline.rounded,maximum_groups=maximum_groups)
        operators[tilt,label]=exact,as_array(exact)
    if args.full_cover:
        settings=SimpleNamespace(groups=maximum_groups,precision=args.precision,
            operator_cache=args.operator_cache,full_feedback=6,window_histogram=8,
            joint_cancellation=True,column_density=True,weight_tilt='1')
        for tilts,penalties,density in ((['.016','.024'],['.75','1'],0),
                (['.032','.04','.048'],['.75','.9','1'],6),
                (['.052','.054','.056','.058','.06'],['.9'],6)):
            settings.tilts=tilts; settings.penalties=penalties; settings.feedback_density=density
            operators.update(baseline.build_operators(settings))
            for rho in penalties:
                if rho not in count_sets:
                    measure(rho,rho)
        args.probe_supports=[]; args.probe_vector=[]
        args.retain_parents=True; args.joint_witness=True; args.joint_top=2; args.screen_only=False
        run_covers(args,occupancies,operators,count_sets,shell_sets)
    else:
        args.probe_supports=args.supports
        baseline.shell_points(args,operators,count_sets,shell_sets)
    print('Other support vectors or occupancies remain; this is not the full-code certificate.',flush=True)


if __name__=='__main__':
    main()
