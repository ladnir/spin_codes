"""Selected dense diagnostics for independently shuffled row pairs.

An exact positive majorant bounds each expected support shell. Probes cover
only their explicit mean cells, never the whole dense range by themselves.
"""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path
from math import comb
from flint import ctx

import pairwise_support
import shared_mixture
import scalar_cover
import birth_classes


def probe_cell(value,radius=Q(0)):
    """An exact declared interval, without silently clipping its scope."""
    value,radius=Q(value),Q(radius)
    if radius<0 or not 0<=value-radius<=value+radius<=1:
        raise ValueError('nonnegative probe radius and a complete cell within [0,1] required')
    return value-radius,value+radius


def pair_majorant(caps, centers, zero_bits=80, cost_tilt=Q(1,4), mass_bits=None,
                  central_bits=None):
    """Reuse a central binomial component before bounding residual shells."""
    if central_bits is not None and (type(central_bits) is not int or not 0 <= central_bits <= 512):
        raise ValueError('bounded nonnegative central exponent required')
    if central_bits is None:
        return shared_mixture.envelope(caps,centers,zero_bits,cost_tilt,mass_bits)
    n = len(caps)-1
    c,p = Q(2)**central_bits,Q(3,4)
    if zero_bits is not None and c*(1-p)**n > Q(2)**-zero_bits:
        raise ValueError('central component exceeds empty-active budget')
    if mass_bits is not None and c > Q(2)**mass_bits:
        raise ValueError('central component exceeds mass budget')
    residual = [max(Q(0),v-c*comb(n,u)*p**u*(1-p)**(n-u)) for u,v in enumerate(caps)]
    result = [(c,p)] + shared_mixture.envelope(residual,centers,zero_bits,cost_tilt,mass_bits)
    shared_mixture.verify(caps,result)
    return result


def combine_pair_mixtures(pair_mixture):
    """Positive push-forward under union, excluding only zero-pair squared."""
    if any(Q(c) <= 0 or not 0 < Q(p) <= 1 for c,p in pair_mixture):
        raise ValueError('positive exact pair mixture required')
    merged = {}
    for c,p in pair_mixture:
        merged[p] = merged.get(p,Q(0))+2*c
    for c,p in pair_mixture:
        for d,q in pair_mixture:
            activity = p+q-p*q
            merged[activity] = merged.get(activity,Q(0))+c*d
    return [(c,p) for p,c in sorted(merged.items())]


def components(step=8, zero_bits=80, cost_tilt=Q(1,4), mass_bits=None,
               refined_counts=False, coupled_counts=False, pair_shells=False, central_bits=None,
               maximum_activity=Q(1), transform_shells=False,prune_mixture=False,joint_mixture=False,
               count_witness=None,union_cost_tilt=None):
    if type(step) is not int or not 1 <= step <= 64:
        raise ValueError('integer support-grid step in 1..64 required')
    maximum_activity = Q(maximum_activity)
    if not Q(3,4) <= maximum_activity <= 1:
        raise ValueError('maximum component activity must lie in [3/4,1]')
    if transform_shells and not (pair_shells and refined_counts and coupled_counts):
        raise ValueError('transform shells require pair-shell and both count-refinement options')
    if count_witness is not None and not transform_shells:
        raise ValueError('shell witness replay requires the transformed quaternary model')
    if (prune_mixture or joint_mixture) and not pair_shells:
        raise ValueError('pair mixture pruning requires pair shells')
    if union_cost_tilt is not None and (not (prune_mixture or joint_mixture) or not 0<Q(union_cost_tilt)<=1):
        raise ValueError('union objective requires mixture pruning and a tilt in (0,1]')
    if transform_shells:
        from pairwise_macwilliams import actual
        system = actual()
        pair_cdf = system.cdfs[0]
        exact = pairwise_support.combine_pair_cdfs(pair_cdf,pair_cdf)
        counts = [-(-v.numerator//v.denominator) for v in exact]
    else:
        counts, info = pairwise_support.pairwise_counts(refined_counts, coupled_counts)
        pair_cdf = info['pair_cdf']
    centers = sorted({Q(u,256) for u in range(38,257,step)
                      if Q(u,256) <= maximum_activity} | {maximum_activity})
    if central_bits is not None and not pair_shells:
        raise ValueError('central pair component requires pair-shell bounds')
    if pair_shells:
        from pairwise_moments import shell_caps
        from dual_shortening import verify_bch_premise
        verify_bch_premise()
        shell_bounds = shell_caps(256,128,30)
        caps = [0]+[min(c-1,m) for c,m in zip(pair_cdf[1:],shell_bounds[1:])]
        if transform_shells:
            caps = [0]+[min(caps[u],system.transform_cap(u)) for u in range(1,257)]
            print('PAIRWISE inverse-transform shell caps regenerated exactly',flush=True)
        if count_witness is not None:
            from pairwise_macwilliams import replay_shell_witnesses
            caps,checked=replay_shell_witnesses(system,json.loads(Path(count_witness).read_text()),caps)
            print('PAIRWISE exact quaternary shell witnesses replayed',checked,flush=True)
        pair_mixture = pair_majorant(caps,centers,zero_bits,cost_tilt,mass_bits,central_bits)
        if joint_mixture:
            from positive_prune import pool
            candidates=[pair_mixture]
            for cost in (Q(1,8),Q(1,4),Q(1,2),Q(1)):
                candidates.append(pair_majorant(caps,sorted(set(centers)|{Q(3,4)}),
                                                zero_bits,cost,mass_bits))
            pair_mixture=pool(caps,candidates)
            print('PAIRWISE pooled',len(candidates),'valid majorants;',len(pair_mixture),'centers',flush=True)
        if prune_mixture or joint_mixture:
            from positive_prune import prune
            pair_mixture,info=prune(caps,pair_mixture,cost_tilt,union_tilt=union_cost_tilt)
            print('PAIRWISE exact-repaired mixture pruning',info,flush=True)
        # Expand (zero + nonzero-pair-majorant)^2, removing only zero*zero.
        # All other comparison components retain the active-group label,
        # even if their Bernoulli support happens to be empty.
        mixture = combine_pair_mixtures(pair_mixture)
    else:
        mixture = shared_mixture.envelope(counts, centers, zero_bits, cost_tilt, mass_bits)
    result = shared_mixture.as_components(mixture)
    result = [(name.replace('shared-', 'pairwise-'), c, law, active)
              for name,c,law,active in result]
    return result, counts, mixture


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--step', type=int, default=8)
    parser.add_argument('--zero-bits', type=int, default=80)
    parser.add_argument('--cost-tilt', default='1/4')
    parser.add_argument('--mass-bits', type=int)
    parser.add_argument('--refined-counts', action='store_true')
    parser.add_argument('--coupled-counts', action='store_true')
    parser.add_argument('--pair-shells', action='store_true', help='Use the exact two-bit moment shell caps and combine positive pair mixtures')
    parser.add_argument('--transform-shells', action='store_true', help='Also use the solver-free quaternary inverse-transform shell bounds')
    parser.add_argument('--prune-mixture',action='store_true',help='Reduce component masses with exact shell repair')
    parser.add_argument('--joint-mixture',action='store_true',help='Pool several valid sets of centers before exact-repaired LP reduction')
    parser.add_argument('--retune-regional',action='store_true',help='Retune the output tilt of the winning regional bound variant')
    parser.add_argument('--count-witness',type=Path,help='Recheck exact dual multipliers for improved pair shells')
    parser.add_argument('--union-cost-tilt',help='Use the gradient of the pair-union tilted mass as the pruning objective')
    parser.add_argument('--central-bits', type=int, help='First use a central pair component of mass 2^bits and activity 3/4')
    parser.add_argument('--maximum-activity', default='1', help='Largest Bernoulli center before combining pairs; all shell bounds are still checked exactly')
    parser.add_argument('--updates', type=int, choices=(2,3,4), default=4)
    parser.add_argument('--minimum-groups', type=int, default=49)
    parser.add_argument('--threshold', type=int, default=209715)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--base-tilt', default='1/32')
    parser.add_argument('--variance-bins', type=int, default=4)
    parser.add_argument('--probe', nargs='+', default=['.0005','.004','.032'])
    parser.add_argument('--probe-radius',default='0',help='Exact half-width of each selected mean interval; zero gives singleton diagnostics')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    radius=Q(args.probe_radius)
    if (args.precision < 128 or radius<0 or any(not 0 <= Q(v)-radius <= Q(v)+radius <= 1 for v in args.probe)
            or Q(args.base_tilt) <= 0 or not 1 <= args.variance_bins <= 64):
        parser.error('precision >=128, probe means in [0,1], positive base tilt and 1..64 variance bins required')
    rows, caps, mixture = components(args.step, args.zero_bits, Q(args.cost_tilt),
                                    args.mass_bits, args.refined_counts, args.coupled_counts,
                                    args.pair_shells, args.central_bits,Q(args.maximum_activity),
                                    args.transform_shells,args.prune_mixture,args.joint_mixture,args.count_witness,
                                    args.union_cost_tilt)
    print('PAIRWISE exact expected-shell majorant;',len(mixture),'components',flush=True)
    for c,p in mixture:
        print('  activity',p,'log2 mass',scalar_cover.logq(c)/scalar_cover.np.log(2),flush=True)
    data = birth_classes.actual(args.updates)
    ctx.prec = args.precision
    model_type=scalar_cover.Model
    if args.retune_regional:
        from pairwise_search import Model
        model_type=Model
    model = model_type(rows, data, args.threshold, args.minimum_groups,
                               Q(args.base_tilt), inner=birth_classes, variance_shuffle=True,
                               variance_bins=args.variance_bins, regional_count=True)
    results = []
    for value in map(Q,args.probe):
        cell=probe_cell(value,radius)
        if model.empty(cell):
            print('Empty mean',value,flush=True)
            continue
        proposal, witness = model.proposal(cell)
        upper = model.outward(cell,witness)
        print('PAIRWISE PROBE mean',value,'cell',tuple(map(str,cell)),'proposal',proposal,'outward log2',
              upper.log()/scalar_cover.arb(2).log(),flush=True)
        results.append(dict(mean=str(value),cell=list(map(str,cell)),proposal=proposal,witness=witness,
                            upper=[int(v) for v in upper.upper().man_exp()]))
    if args.output:
        record = dict(schema='pairwise-gf16-dense-probes-1',
                      parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
                      mixture=[dict(mass=str(c),activity=str(p)) for c,p in mixture],
                      shell_caps=caps,probes=results,
                      note='Expected support-shell majorant with outward bounds on the explicitly listed cells only. '
                           'No whole-domain cover or full distance certificate.')
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(record,indent=2)+'\n')


if __name__ == '__main__':
    main()
