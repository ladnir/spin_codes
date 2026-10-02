"""One-point diagnostic: decompose and refine two dominant variance bins.

This is not a full cover or a new certificate. All source objects are rebuilt;
the one refined witness is checked by the unmodified outward implementation.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import numpy as np
from flint import arb, arb_mat, ctx
import shared_relaxed_strategy as proof
import shared_relaxed_region_cache as cache_module
import regional_count as regional
import variance_partition as variance


def decompose(model, cell, witness, cache):
    parts, checked = regional.prepare_witness(model, cell, witness)
    polynomial = cache.polynomial(model, witness)
    lam = Q(witness['parameters'][0]); direct = witness.get('regional_direct_counts', False)
    weights, _ = model.weights(cell, model.tilt)
    scale = Q(1) if direct else sum(weights)
    masses = ([regional.aq(model.tilt)**(-j) for j in range(regional.sc.G+1)] if direct
              else regional.binomial_masses(regional.sc.G, weights[1]/scale))
    cs, _, _ = model.family(model.tilt)
    cutoff_log = regional.aq(lam)*model.threshold
    scale_log = regional.sc.PACKETS*regional.aq(scale).log()
    total = arb(0); rows = []; log2 = arb(2).log()
    for index, ((interval, dual), part) in enumerate(zip(parts, checked['regional_count_parts'])):
        eta, mu, gamma = dual
        cap = variance.factor(model, cell, interval[0], witness['variance_dual'])
        ratios = (regional.count_mass_caps if direct else regional.count_ratios)(
            model.features, model.active, cell, interval, model.q_min, regional.sc.G,
            regional.aq(cap), part['mgf_witnesses'])
        size = polynomial[0].nrows()
        matrix = sum((mass*ratio*value for mass, ratio, value in zip(masses, ratios, polynomial)),
                     arb_mat(size, size))
        power = matrix**regional.sc.REGIONS
        moment = sum((power[0,j] for j in range(size)), arb(0))
        summands = [regional.aq(c)*regional.aq(eta*f+mu*a+gamma*f*(1-f)).exp()
                    for c, f, a in zip(cs, model.features, model.active)]
        count = sum(summands, arb(0))
        exponent = (regional.sc.G*(min(eta*x for x in cell)+min(gamma*v for v in interval))
                    + mu*model.q_min)
        outside = regional.sc.G*count.log()-regional.aq(exponent)
        term = (outside+moment.log()).exp()
        total = regional.up(total+term)
        composition = [float(v/count) for v in summands]
        rows.append(dict(index=index, interval=list(map(str, interval)), dual=list(map(str, dual)),
            outer_bits=str(outside/log2), inner_bits=str(moment.log()/log2),
            cutoff_bits=str(cutoff_log/log2), scale_bits=str(scale_log/log2),
            total_bits=str((outside+moment.log()+cutoff_log+scale_log)/log2),
            score=float((outside+moment.log()+cutoff_log+scale_log)/log2),
            groups=[dict(activity=str(p), count=regional.sc.G*v)
                    for (_,p,_), v in zip(model.components, composition) if v > 1e-8]))
        print('OBSTRUCTION part',index,'bits',rows[-1]['total_bits'],flush=True)
    upper = regional.up((total.log()+cutoff_log+scale_log).exp())
    return upper, rows, checked


def refine(model, cell, witness, indices, splits=4):
    old_parts = variance.validate(cell, witness['variance_partition'])
    if (type(splits) is not int or not 1 <= splits <= 16
            or not isinstance(indices,(list,tuple)) or not indices
            or len(set(indices)) != len(indices)
            or any(type(i) is not int or not 0 <= i < len(old_parts) for i in indices)):
        raise ValueError('distinct existing bins and integer split count in 1..16 required')
    _, _, logs = model.family(model.tilt)
    features = np.array(list(map(float, model.features)))
    active = np.array(model.active)
    parts = []; families = []
    for i, (interval, dual) in enumerate(old_parts):
        if i not in indices:
            parts.append(copy.deepcopy(witness['variance_partition'][i]))
            families.append(copy.deepcopy(witness['regional_count_parts'][i]))
            continue
        lo, hi = interval
        for j in range(splits):
            child = (lo+(hi-lo)*j/splits, lo+(hi-lo)*(j+1)/splits)
            _, proposed = variance.outer_witness(logs, features, active, cell,
                model.q_min/regional.sc.G, child)
            part = dict(interval=list(map(str, child)), dual=list(map(str, proposed)))
            families.append(dict(part, mgf_witnesses=regional.propose_mgf(
                model.features, model.active, cell, child, model.q_min, regional.sc.G,
                witness.get('regional_tilted_atom',False),
                regional.FINE_TILTS if witness.get('regional_fine_tilts',False) else regional.TILTS,
                witness.get('regional_tilted_variance',False))))
            parts.append(part)
    variance.validate(cell, parts)
    return dict(witness, variance_partition=parts, regional_count_parts=families)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('new diagnostic output path required')
    raw = args.source.read_bytes(); record = json.loads(raw)
    row = proof.validate_record(record)
    point = next(p for p in row['probes'] if Q(p['mean']) == Q('4/125'))
    model = proof.build_model(record,256); cell = (Q('4/125'),)*2
    result = dict(schema='shared-gf16-one-point-obstruction-1',
        source_sha256=hashlib.sha256(raw).hexdigest(), threshold=model.threshold,
        mean='4/125',precision=256,note='One point only; no complete certificate.')
    with cache_module.install() as cache:
        ctx.prec = 256
        upper, parts, witness = decompose(model,cell,point['witness'],cache)
        result['baseline'] = dict(upper=proof.endpoint(proof.dyadic([int(v) for v in upper.upper().man_exp()])),
            log2_upper=str(upper.log()/arb(2).log()),parts=parts)
        indices = [r['index'] for r in sorted(parts,key=lambda v:v['score'],reverse=True)[:2]]
        candidate = refine(model,cell,witness,indices)
        result['refined_indices'] = indices
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        ctx.prec = 256
        checked, _ = regional.outward(model,cell,candidate)
        if ctx.prec != 256 or not checked > 0:
            raise ArithmeticError('positive fresh refined bound required')
        result['refined'] = dict(upper=proof.endpoint(proof.dyadic([int(v) for v in checked.upper().man_exp()])),
            log2_upper=str(checked.log()/arb(2).log()),witness=candidate)
        result['cache'] = dict(hits=cache.hits,misses=cache.misses)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
        print('OBSTRUCTION baseline',result['baseline']['log2_upper'],
              'refined',result['refined']['log2_upper'],'parts',indices,flush=True)


if __name__ == '__main__':main()
