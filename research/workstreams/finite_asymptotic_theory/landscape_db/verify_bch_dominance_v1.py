"""Independent positive high-precision replay of selected composition witnesses."""
import json
import math
from pathlib import Path

import mpmath as mp

import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent


def epoch_matrices(t, s, ac, kernel, lam, maximum):
    m = mp.mpf((1 << s)-1); kappa = m/(m-1)
    z = mp.exp(-mp.mpf(lam)); d = min(ac)
    result = []
    for j in range(maximum+1):
        total = math.comb(t, j)
        moments = {w: sum(mp.mpf(math.comb(w, v)*math.comb(t-w, j-v))*z**(w+j-2*v)
                           for v in range(max(0, j-t+w), min(w, j)+1))/total for w in ac}
        arbitrary = max(moments.values())
        uniform = sum(ac[w]*moments[w] for w in ac)/m
        matrix = mp.matrix(4)
        if j == 0:
            matrix[0, 0] = 1
            matrix[1, 2] = z**d
            matrix[2, 2] = uniform
            matrix[3, 2] = min(z**d, kappa*uniform)
        else:
            beta = mp.mpf(kernel[j])/total
            matrix[0, 0] = beta*z**j
            matrix[0, 1] = (1-beta)*z**j
            for state, moment in ((1, arbitrary), (2, uniform), (3, min(arbitrary, kappa*uniform))):
                zero_part = min(moment, (1-beta)*z**max(0, d-j))/m
                live_part = (1-1/m)*moment+min(moment, beta*z**max(0, d-j))/m
                matrix[state, 0] = zero_part
                matrix[state, 3] = live_part
        result.append(matrix*total)
    return result


def convolve(a, b, maximum):
    out = [mp.matrix(4) for _ in range(min(maximum+1, len(a)+len(b)-1))]
    for i, left in enumerate(a):
        for j, right in enumerate(b[:len(out)-i]):
            out[i+j] += left*right
    return out


def component(t, s, ac, kernel, length, block, counts, bands, indices, tilt, shift):
    q = len(indices); lam = math.exp(tilt)
    power = epoch_matrices(t, s, ac, kernel, lam, q)
    current = [mp.eye(4)]; exponent = length//t
    while exponent:
        if exponent & 1:
            current = convolve(current, power, q)
        exponent >>= 1
        if exponent:
            power = convolve(power, power, q)
    regions = [v/math.comb(length, j) for j, v in enumerate(current)]
    law = [mp.mpf(1)]; cost = mp.mpf(0)
    for index in indices:
        band = bands[index]
        if band == [block]:
            p = mp.mpf(1); gamma = mp.log(counts[block])
        else:
            midpoint = mp.mpf(min(band)+max(band))/(2*block)
            eta = mp.log(midpoint/(1-midpoint))+mp.mpf(shift)
            p = 1/(1+mp.exp(-eta))
            gamma = max(mp.log(counts[w])-mp.log(math.comb(block, w))-w*mp.log(p)
                        -(block-w)*mp.log1p(-p) for w in band)
        updated = [mp.mpf(0)]*(len(law)+1)
        for j, value in enumerate(law):
            updated[j] += value*(1-p); updated[j+1] += value*p
        law = updated; cost += gamma
    matrix = mp.matrix(4)
    for probability, region in zip(law, regions):
        matrix += probability*region
    product = matrix**block
    bound = mp.log(sum(product[0, j] for j in range(4)))+cost+lam*(block*length//10)
    trivial = sum(mp.log(sum(counts[w] for w in bands[g])) for g in indices)
    return min(trivial, bound)


def main():
    _, counts, maps, _ = study.load_inputs()
    receipt = json.loads((HERE/'bch_dominance_sparse.json').read_text())
    for name, digest in receipt['source_sha256'].items():
        if study.sha(study.ROOT/name) != digest:
            raise ValueError(f'changed source: {name}')
    if study.sha(HERE/'bch_dominance_sparse.csv') != receipt['csv_sha256']:
        raise ValueError('changed summary')
    checked = []
    with mp.workdps(90):
        for block, t, s, exponent in ((64, 64, 20, 20), (128, 64, 20, 20),
                                       (64, 64, 12, 20), (128, 256, 14, 20)):
            name = f'b{block}_t{t}_s{s}_e{exponent}.json'
            path = HERE/'bch_dominance_v1'/name
            if not path.exists():
                raise ValueError(f'missing requested replay checkpoint: {name}')
            if study.sha(path) != receipt['checkpoints'][name]:
                raise ValueError('changed checkpoint')
            data = json.loads(path.read_text()); config = maps[t, s]
            ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
            for q in (2, 3, 4):
                detail = data['occupations'][str(q)]; j = detail['dominant_composition']
                tilt, shift = detail['witnesses'][j]
                value = component(t, s, ac, config['kernel_counts'], (1 << exponent)//(block//2),
                                  block, counts[block], detail['bands'], detail['compositions'][j], tilt, shift)
                error = abs(float(value)-detail['component_log_upper'][j])
                if error > 2e-8:
                    raise ArithmeticError(f'composition replay mismatch: {name} Q{q} {error}')
                check = dict(checkpoint=name, occupation=q, absolute_log_error=error)
                checked.append(check); print(check, flush=True)
    result = dict(decimal_digits=90, checks=checked, status='SELECTED_WITNESS_REPLAY',
                  source_sha256={str(Path(__file__)): study.sha(Path(__file__)),
                                 str(HERE/'bch_dominance_sparse.json'): study.sha(HERE/'bch_dominance_sparse.json')},
                  note='Independent positive arithmetic; no outward certificate.')
    (HERE/'bch_dominance_verification.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
