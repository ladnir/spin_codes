"""Diagnose the active-group constraint in a shared-row dense witness.

The softmax composition describes the positive comparison measure after
exponential tilting, not the actual code ensemble. Outward point replays
are separate from these floating diagnostics and never cover all means.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from flint import arb, ctx
import birth_classes
import regional_count as regional
import regional_count_probe as floating
import scalar_cover as sc
import shared_mixture
import variance_partition as variance
from shared_relaxed_dense import reuse_exact_census, save


def composition(model, dual):
    """Floating dual composition; never a probability certificate."""
    eta, mu, gamma = list(map(float, map(Q, dual)))
    _, _, logs = model.family(model.tilt)
    features = np.array(list(map(float, model.features)))
    active = np.array(model.active)
    scores = logs + eta*features + mu*active + gamma*features*(1-features)
    masses = np.exp(scores-logsumexp(scores))
    return dict(expected_active_groups=float(sc.G*(masses@active)),
                mean=float(masses@features),
                average_variance=float(masses@(features*(1-features))),
                occupancy_dual=mu,
                components=[dict(activity=str(p), feature=str(f), active=a,
                                 probability=float(weight))
                            for (_, p, a), f, weight in zip(model.components, model.features, masses)])


def regional_parts(model, cell, witness):
    """Decompose the floating regional bound to locate its dominant parts."""
    parts, checked = regional.prepare_witness(model, cell, witness)
    scratch = regional.proposal_cache(model)
    local = regional.local_operators(model, witness, _proposal_scratch=scratch)
    size = local[0].nrows()
    arrays = np.array([[[float(m[i,j]) for j in range(size)] for i in range(size)] for m in local])
    matrices, region_logs = floating.scaled_placement(arrays, 64)
    scale, count_weights = regional.proposal_count_weights(model, cell, witness, parts, checked, scratch)
    _, _, logs = model.family(model.tilt)
    features = np.array(list(map(float, model.features)))
    active = np.array(model.active)
    lam = Q(witness['parameters'][0])
    rows = []
    for index, ((interval, dual), log_weights) in enumerate(zip(parts, count_weights)):
        matrix, shift = floating.weighted_region(matrices, log_weights+region_logs)
        moment = sc.log_power(matrix, sc.REGIONS)+sc.REGIONS*shift
        eta, mu, gamma = map(float, dual)
        outside = sc.G*logsumexp(logs+eta*features+mu*active+gamma*features*(1-features))
        outside -= sc.G*(min(eta*float(x) for x in cell)+min(gamma*float(v) for v in interval))+mu*model.q_min
        value = (outside+moment+sc.PACKETS*sc.logq(scale)+float(lam)*model.threshold)/np.log(2)
        rows.append(dict(index=index+1, interval=list(map(str, interval)),
                         log2_upper_proposal=float(value), **composition(model, dual)))
    return rows


def retune_outer(model, cell, old):
    """Retune only counting duals; retain valid regional MGF inequalities."""
    witness = deepcopy(old)
    _, _, logs = model.family(model.tilt)
    features = np.array(list(map(float, model.features)))
    active = np.array(model.active)
    best = None
    for sign in (False, True):
        eta, mu = sc.outer_dual(logs, features, active, float(cell[0 if sign else 1]), model.q_min/sc.G, sign)
        dual = [Q(round(float(x)*10**8),10**8) for x in (eta, mu)]
        e, u = map(float, dual)
        score = logsumexp(logs+e*features+u*active)-min(e*float(x) for x in cell)-u*model.q_min/sc.G
        if best is None or score < best[0]:
            best = score, dual
    witness['parameters'] = [old['parameters'][0], *map(str, best[1])]
    for index, part in enumerate(witness['variance_partition']):
        interval = tuple(map(Q, part['interval']))
        _, dual = variance.outer_witness(logs, features, active, cell, model.q_min/sc.G, interval)
        part['dual'] = list(map(str, dual))
        witness['regional_count_parts'][index]['dual'] = list(map(str, dual))
    return witness


def cases(path):
    record = json.loads(path.read_text())
    if 'results' in record:
        mixture = [(Q(x['mass']), Q(x['activity'])) for x in record['mixture']]
        for row in record['results']:
            if Q(row['distance']) not in (Q('.06'), Q('.07')):
                continue
            for probe in row['probes']:
                if Q(probe['mean']) == Q('.032') and 'witness' in probe:
                    yield record, mixture, Q(row['distance']), probe['witness']
    else:
        for row in record['rows']:
            if row['name'] != 'pruned-1/4':
                continue
            mixture = [(Q(x['mass']), Q(x['activity'])) for x in row['mixture']]
            for probe in row['probes']:
                if Q(probe['mean']) == Q('.032') and 'witness' in probe:
                    yield record, mixture, Q(record['parameters']['distance']), probe['witness']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipts', nargs='+', type=Path, required=True)
    parser.add_argument('--minimum-groups', nargs='+', type=int, default=[65,97,129,257])
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.precision < 128 or any(not 1 <= q <= sc.G for q in args.minimum_groups):
        parser.error('fresh output path, precision>=128, and valid group minima required')
    _, caps, _ = shared_mixture.actual_components(coupled=True, zero_bits=64, cost_tilt=Q(1,4))
    data = birth_classes.actual(2)
    ctx.prec = args.precision
    result = dict(schema='shared-gf16-relaxed-bridge-1', precision=args.precision,
                  mean='4/125', updates=2, base_tilt='3/16',
                  note='Softmax compositions are comparison diagnostics, not actual message occupancies. Outward results cover one mean point only.',
                  cap_sha256=hashlib.sha256(json.dumps([str(x) for x in caps]).encode()).hexdigest(), cases=[])
    previous = None
    for path in args.receipts:
        for record, mixture, distance, witness in cases(path):
            shared_mixture.verify(caps, mixture)
            model = sc.Model(shared_mixture.as_components(mixture), data, int(distance*sc.N),33,Q(3,16),
                             inner=birth_classes,variance_shuffle=True,variance_bins=16,regional_count=True)
            if previous is not None:
                reuse_exact_census(previous, model)
            previous = model
            cell = (Q(4,125), Q(4,125))
            row = dict(receipt=str(path), receipt_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                       distance=str(distance), checked_shell_domination=True,
                       basic=composition(model,[*witness['parameters'][1:],'0']),
                       variance_parts=regional_parts(model,cell,witness), replays=[])
            result['cases'].append(row)
            print('BRIDGE COMPOSITION', path, distance, 'basicq', row['basic']['expected_active_groups'],
                  'dominantparts',[(r['index'],r['log2_upper_proposal'],r['expected_active_groups'])
                                   for r in sorted(row['variance_parts'],key=lambda x:x['log2_upper_proposal'],reverse=True)[:3]],flush=True)
            save(args.output,result)
            if distance != Q('.06'):
                continue
            for q in args.minimum_groups:
                target = sc.Model(shared_mixture.as_components(mixture),data,int(distance*sc.N),q,Q(3,16),
                                  inner=birth_classes,variance_shuffle=True,variance_bins=16,regional_count=True)
                reuse_exact_census(previous,target)
                previous=target
                tuned=retune_outer(target,cell,witness)
                score, checked=regional.propose(target,cell,tuned)
                ctx.prec=args.precision
                upper=target.outward(cell,checked)
                if not upper>0:
                    raise ArithmeticError('positive point upper required')
                point=dict(minimum_groups=q,proposal=score,witness=checked,
                           basic=composition(target,[*checked['parameters'][1:],'0']),
                           upper=[int(v) for v in upper.upper().man_exp()],
                           log2_upper=str(upper.log()/arb(2).log()))
                row['replays'].append(point)
                print('BRIDGE OUTWARD', path, 'qmin',q,point['log2_upper'],flush=True)
                save(args.output,result)
    save(args.output,result)


if __name__=='__main__':
    main()
