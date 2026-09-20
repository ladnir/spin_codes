"""Replay every selected interval and aggregate a BCH full-reference bound."""
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

import occupation_refresh_v1 as transfer
import occupation_positive_poly_v1 as positive
import refine_bch_dense_v1 as dense
import trim_bch_dense_v3 as partitions
import study_bch_dominance_v1 as study
import verify_bch_dominance_v1 as mp_reference

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


def main():
    _, spectra, maps, _ = study.load_inputs()
    block, t, s, exponent = 128, 64, 20, 20
    length = (1 << exponent)//(block//2); cutoff = block*length//10
    counts, config = spectra[block], maps[t, s]
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    epoch = transfer.Epochs(t, s, ac, config['kernel_counts'])
    sparse_receipt = json.loads((HERE/'bch_dominance_sparse.json').read_text())
    checkpoint = HERE/'bch_dominance_v1/b128_t64_s20_e20.json'
    if study.sha(checkpoint) != sparse_receipt['checkpoints'][checkpoint.name]:
        raise ValueError('Q2..4 checkpoint changed')
    sparse = json.loads(checkpoint.read_text())
    paths = [HERE/'bch_sparse_tail_v1_b128_t64_s20_e20/cover.json',
             HERE/'bch_sparse_tail_v2_b128_t64_s20_e20_q65_256/cover.json',
             HERE/'bch_dense_v3_b128_t64_s20_e20_q257/cover.json']
    payloads = [json.loads(path.read_text()) for path in paths]
    dependencies = dict(sparse_receipt['source_sha256'])
    for data in payloads:
        for name, digest in data['source_sha256'].items():
            if name in dependencies and dependencies[name] != digest:
                raise ValueError('incompatible source versions')
            dependencies[name] = digest
    for name, digest in dependencies.items():
        if study.sha(study.ROOT/name) != digest:
            raise ValueError(f'changed dependency: {name}')
    intervals = [(data['occupation_min'], data['occupation_max']) for data in payloads]
    if intervals != [(5, 64), (65, 256), (257, length)]:
        raise ValueError('incomplete or overlapping occupation cover')
    coefficient_checks = []
    for data in payloads[:2]:
        if (data['arguments']['block'], data['arguments']['step'], data['arguments']['state'], data['arguments']['exponent']) != (block, t, s, exponent):
            raise ValueError('incompatible sparse geometry')
        maximum = data['occupation_max']
        models = [transfer.CompositionBoxes(counts, block, data['bands'], p, maximum) for p in data['probability_banks']]
        cache = {}
        def regions(log_tilt):
            if log_tilt not in cache:
                raw = epoch.at(math.exp(log_tilt), maximum)
                cache[log_tilt] = positive.region_logs(raw, t, length, maximum)
            return cache[log_tilt]
        rows = data['occupations']
        if [r['occupation'] for r in rows] != list(range(data['occupation_min'], maximum+1)):
            raise ValueError('missing integer occupation')
        for row in rows:
            q = row['occupation']; log_tilt = row['witness_log_tilt']
            model = models[row['witness_probability_bank']]
            current = regions(log_tilt)[:q+1]
            for _ in range(q):
                current = model.envelope.apply(current[:-1], current[1:])
            adaptive = (math.log(math.comb(length, q))+q*math.log(len(data['bands']))
                        +cutoff*math.exp(log_tilt)+float(transfer.terminal_logs(current, block)[0]))
            if abs(adaptive-row['adaptive_log_upper']) > 2e-7:
                raise ArithmeticError('adaptive occupation replay failed')
            cover = row['composition_cover']; selected = adaptive
            if cover is not None:
                boxes = cover['boxes']
                partitions.root_partition(dict(lower=[0]*len(data['bands']), upper=[q]*len(data['bands'])), boxes, q)
                values = []
                for box in boxes:
                    witness = box['witness']; z = witness['log_tilt']; m = models[witness['probability_bank']]
                    value = m.bound(regions(z), box['lower'], box['upper'], q, length, cutoff, math.exp(z))
                    if abs(value-box['log_bound']) > 2e-7:
                        raise ArithmeticError('composition box replay failed')
                    values.append(value)
                combined = float(np.logaddexp.reduce(values))
                if abs(combined-cover['log_union_upper']) > 2e-7:
                    raise ArithmeticError('composition cover sum failed')
                selected = min(selected, combined)
            if abs(selected-row['log_upper']) > 2e-7:
                raise ArithmeticError('selected occupation bound failed')
        combined = float(np.logaddexp.reduce([r['log_upper'] for r in rows]))
        if abs(combined-data['log_upper']) > 2e-7:
            raise ArithmeticError('sparse interval sum failed')
        for z in list(cache)[:2]:
            expected = transfer.region_logs(epoch.at(math.exp(z), maximum), t, length, maximum)
            actual = cache[z]
            if not np.array_equal(np.isfinite(expected), np.isfinite(actual)):
                raise ArithmeticError('positive polynomial support changed')
            mask = np.isfinite(expected)
            error = float(max(abs(expected[mask]-actual[mask])))
            if error > 2e-7:
                raise ArithmeticError('positive coefficient replay failed')
            coefficient_checks.append(dict(maximum_occupation=maximum, log_tilt=z, maximum_log_error=error))
        print(f'Q{data["occupation_min"]}..{maximum}: every bound replayed', flush=True)
    data = payloads[2]
    if tuple(data[k] for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent')) != (block, t, s, exponent):
        raise ValueError('incompatible dense geometry')
    leaves = [leaf for root in data['roots'] for leaf in root['leaves']]
    coverage = dense.check_cover(leaves, length, 257)
    model = dense.Refiner(counts, block, t, s, ac, config['kernel_counts'], length, data['bands'])
    for root in data['roots']:
        values = []
        for leaf in root['leaves']:
            corners = transfer.typed.vertices(leaf['lower'], leaf['upper'], length)
            witness = leaf['witness']; p = np.array(witness['probabilities']); costs = model.costs(p)
            if max(abs(costs-witness['log_density_costs'])) > 2e-10:
                raise ValueError('dense counting measure changed')
            vertex = np.array(witness['maximum_vertex'])
            if np.any(vertex < leaf['lower']) or np.any(vertex > leaf['upper']) or vertex.sum() != length:
                raise ValueError('reported dense vertex escapes its box')
            value = float(max(model.value(corners, witness['log_surprisal'], np.array(witness['proposal']), p, costs)))
            value += transfer.typed.lattice_log_count(leaf['lower'], leaf['upper'])
            if abs(value-leaf['own_log_bound']) > 2e-6:
                raise ArithmeticError('dense leaf replay failed')
            values.append(value)
        if abs(float(np.logaddexp.reduce(values))-root['log_upper']) > 2e-6:
            raise ArithmeticError('dense root sum failed')
    total = float(np.logaddexp.reduce([r['log_upper'] for r in data['roots']]))
    if abs(total-data['log_upper']) > 2e-6:
        raise ArithmeticError('dense interval sum failed')
    print('complete dense integer cover replayed', coverage, flush=True)
    mp_checks = []
    selected = [max(leaves, key=lambda r: r['own_log_bound']),
                min(leaves, key=lambda r: r['lower'][0]), leaves[len(leaves)//2]]
    with mp.workdps(90):
        for leaf in selected:
            witness = leaf['witness']; lam = math.exp(witness['log_surprisal'])
            p = [mp.mpf(x) for x in witness['probabilities']]
            proposal = [mp.mpf(x) for x in witness['proposal']]
            mass = sum(proposal); proposal = [x/mass for x in proposal]
            theta = sum(a*b for a, b in zip(proposal, p))
            raw = mp_reference.epoch_matrices(t, s, ac, config['kernel_counts'], lam, t)
            matrix = mp.matrix(4)
            for j, coefficient in enumerate(raw):
                matrix += coefficient*theta**j*(1-theta)**(t-j)
            result = matrix**(block*length//t)
            moment = mp.log(sum(result[0, j] for j in range(4)))
            costs = [mp.mpf(0)]
            for band, probability in zip(data['bands'], p[1:]):
                costs.append(max(mp.log(counts[w])-mp.log(math.comb(block, w))-w*mp.log(probability)
                                 -(block-w)*mp.log1p(-probability) for w in band))
            corners = transfer.typed.vertices(leaf['lower'], leaf['upper'], length)
            values = []
            for corner in corners:
                multinomial = mp.loggamma(length+1)-sum(mp.loggamma(int(c)+1) for c in corner)
                values.append(cutoff*mp.mpf(lam)+moment+sum(int(c)*(cost-block*mp.log(pi)) for c, cost, pi in zip(corner, costs, proposal))
                              -(block-1)*multinomial)
            widths = [b-a+1 for a, b in zip(leaf['lower'], leaf['upper'])]
            penalty = sum(mp.log(w) for w in widths)-mp.log(max(widths))
            error = abs(float(max(values)+penalty)-leaf['own_log_bound'])
            if error > 2e-6:
                raise ArithmeticError('90-digit dense witness replay failed')
            mp_checks.append(dict(lower=leaf['lower'], upper=leaf['upper'], absolute_log_error=error))
    q1 = -sparse['summary']['q1_margin_bits']*LN2
    components = [dict(occupation_min=1, occupation_max=1, log_upper=q1)]
    components.extend(dict(occupation_min=q, occupation_max=q, log_upper=sparse['occupations'][str(q)]['log_upper']) for q in (2, 3, 4))
    components.extend(dict(occupation_min=d['occupation_min'], occupation_max=d['occupation_max'], log_upper=d['log_upper']) for d in payloads)
    tail = float(np.logaddexp.reduce([r['log_upper'] for r in components[1:]]))
    full = float(np.logaddexp(q1, tail))
    result = dict(status='VERIFIED_BINARY64_FULL_REFERENCE', block_bits=block, step_bits=t, state_bits=s,
                  message_exponent=exponent, components=components, full_margin_bits=-full/LN2,
                  q1_margin_bits=-q1/LN2, higher_occupation_margin_bits=-tail/LN2,
                  higher_to_q1_ratio=math.exp(tail-q1), margin_penalty_bits=(full-q1)/LN2,
                  coefficient_checks=coefficient_checks, dense_cover=coverage, dense_90_digit_checks=mp_checks,
                  authenticated_dependencies=len(dependencies),
                  source_sha256={str(p): study.sha(p) for p in [Path(__file__), checkpoint, *paths]},
                  limitations=['One geometry; no interpolation to other parameters.', 'No outward arithmetic certificate.',
                               'Dominance compares bound contributions, not true event probabilities.'])
    (HERE/'bch_full_reference_b128_t64_s20_e20.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k.endswith('bits') or k.endswith('ratio')}, indent=2), flush=True)


if __name__ == '__main__':
    main()
