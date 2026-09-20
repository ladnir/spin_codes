"""Continuous component-witness refinement of the BCH sparse dominance grid.

Preserve the authenticated v1 arrays; refine only components that can affect
the aggregate. Every numerical candidate remains a valid fixed witness.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

import occupation_refresh_v1 as transfer
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


def fixed_component(regions, counts, block, bands, indices, cutoff, lam, shift):
    roots, active, inactive = transfer.old.density_roots(counts, block, bands, shift)
    indices = np.array(indices, dtype=np.int32)[None, :]
    law = transfer.composition.probability_logs(indices, active, inactive)[0]
    matrix = np.logaddexp.reduce(regions[:len(law)]+law[:, None, None], axis=0)
    moment = float(transfer.terminal_logs(matrix[None], block)[0])
    cost = block*float(roots[indices].sum())
    trivial = sum(math.log(sum(counts[w] for w in bands[g])) for g in indices[0])
    return min(trivial, moment+cost+cutoff*lam)


def refine(payload, counts, config, maximum_components=24):
    summary = payload['summary']; block, t, s, exponent = (summary[k] for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
    length = (1 << exponent)//(block//2); cutoff = block*length//10
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    epochs = transfer.Epochs(t, s, ac, config['kernel_counts'])
    region_cache = {}
    def regions(a):
        if a not in region_cache:
            if len(region_cache) > 128:
                region_cache.clear()
            lam = math.exp(a)/length
            region_cache[a] = transfer.region_logs(epochs.at(lam, 4), t, length, 4)
        return region_cache[a]
    detail = payload['occupations']; checks = []
    for q in (2, 3, 4):
        d = detail[str(q)]
        best = np.array(d['component_log_upper'])
        indices = np.array(d['compositions'])
        multiplicities = transfer.composition.log_multiplicities(indices)
        initial_terms = best+multiplicities
        # Bound the total discarded relevance budget by the number of
        # components, rather than testing each against a fixed gap alone.
        cutoff_term = max(initial_terms)-(30*LN2+math.log(len(best)))
        selected = [int(j) for j in np.argsort(initial_terms)[::-1] if initial_terms[j] >= cutoff_term][:maximum_components]
        for j in selected:
            tilt, shift = d['witnesses'][j]; a = tilt+math.log(length)
            baseline = fixed_component(regions(a), counts, block, d['bands'], indices[j], cutoff, math.exp(a)/length, shift)
            if abs(baseline-best[j]) > 2e-7:
                raise ArithmeticError('sparse component failed direct replay')
            calls = 0
            original = float(best[j]); winner = [tilt, shift]
            def objective(coordinates):
                nonlocal calls, winner
                log_scaled, eta = map(float, coordinates)
                lam = math.exp(log_scaled)/length
                value = fixed_component(regions(log_scaled), counts, block, d['bands'], indices[j], cutoff, lam, eta)
                calls += 1
                if value < best[j]:
                    best[j] = value; winner = [log_scaled-math.log(length), eta]
                return value/(q*block)
            result = minimize(objective, [a, shift], method='L-BFGS-B',
                              bounds=[(a-.7, a+.7), (max(-5., shift-1.5), min(9., shift+1.5))],
                              options={'maxiter': 80, 'ftol': 1e-14, 'gtol': 1e-10, 'maxls': 30})
            objective(result.x)
            d['witnesses'][j] = winner
            checks.append(dict(occupation=q, component=j, improvement_bits=(original-float(best[j]))/LN2,
                               evaluations=calls, optimizer_success=bool(result.success)))
        d['component_log_upper'] = best.tolist()
        d['log_upper'] = math.log(math.comb(length, q))+float(np.logaddexp.reduce(best+multiplicities))
        d['dominant_composition'] = int(np.argmax(best+multiplicities))
        d['dominant_log_tilt'], d['dominant_shift'] = d['witnesses'][d['dominant_composition']]
    q1 = -summary['q1_margin_bits']*LN2
    tail = float(np.logaddexp.reduce([detail[str(q)]['log_upper'] for q in (2, 3, 4)]))
    result = dict(summary)
    for q in (2, 3, 4):
        result[f'q{q}_margin_bits'] = -detail[str(q)]['log_upper']/LN2
    result.update(q2_q4_margin_bits=-tail/LN2, q2_q4_to_q1_log2_ratio=(tail-q1)/LN2,
                  q1_q4_penalty_bits=(float(np.logaddexp(q1, tail))-q1)/LN2,
                  q1_q4_margin_bits=-float(np.logaddexp(q1, tail))/LN2)
    return dict(summary=result, occupations=detail, refinements=checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint', default='b64_t64_s20_e26.json')
    parser.add_argument('--maximum-components', type=int, default=24)
    args = parser.parse_args()
    _, counts, maps, dependencies = study.load_inputs()
    receipt = json.loads((HERE/'bch_dominance_sparse.json').read_text())
    path = HERE/'bch_dominance_v1'/args.checkpoint
    if study.sha(path) != receipt['checkpoints'][args.checkpoint]:
        raise ValueError('changed sparse checkpoint')
    payload = json.loads(path.read_text()); row = payload['summary']
    result = refine(payload, counts[row['block_bits']], maps[row['step_bits'], row['state_bits']], args.maximum_components)
    for source in (Path(__file__), Path(transfer.__file__), path, HERE/'bch_dominance_sparse.json'):
        dependencies[source.relative_to(study.ROOT).as_posix()] = study.sha(source)
    result['source_sha256'] = dependencies
    result['status'] = 'BINARY64_SPARSE_COMPONENT_REFINEMENT'
    directory = HERE/'bch_dominance_v2'; directory.mkdir(exist_ok=True)
    (directory/args.checkpoint).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result['summary'], indent=2), flush=True)
    print('components refined', len(result['refinements']), flush=True)


if __name__ == '__main__':
    main()
