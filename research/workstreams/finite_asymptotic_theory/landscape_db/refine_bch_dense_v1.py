"""Replay a complete BCH type cover and refine its continuous witnesses.

Finite witness searches only propose parameters. Every selected bound is
re-evaluated directly; the integer cover is checked independently.
"""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, logsumexp

import occupation_refresh_v1 as transfer
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


def lattice_count(lower, upper, total):
    """Exact bounded stars-and-bars count by inclusion-exclusion."""
    remaining = total-sum(lower); widths = [b-a+1 for a, b in zip(lower, upper)]
    groups = len(widths); result = 0
    for mask in range(1 << groups):
        left = remaining-sum(widths[g] for g in range(groups) if mask & (1 << g))
        if left >= 0:
            result += (-1 if mask.bit_count() & 1 else 1)*math.comb(left+groups-1, groups-1)
    return result


def check_cover(boxes, length, minimum):
    lo = np.array([b['lower'] for b in boxes], dtype=np.int64)
    hi = np.array([b['upper'] for b in boxes], dtype=np.int64)
    if np.any(lo < 0) or np.any(hi > length) or np.any(lo > hi) or np.any(hi[:, 0] > length-minimum):
        raise ValueError('box escapes the claimed occupation domain')
    groups = lo.shape[1]
    volumes = [lattice_count(list(map(int, a)), list(map(int, b)), length) for a, b in zip(lo, hi)]
    if any(v <= 0 for v in volumes):
        raise ValueError('empty selected type box')
    expected = math.comb(length+groups-1, groups-1)-math.comb(minimum+groups-2, groups-1)
    if sum(volumes) != expected:
        raise ValueError('selected boxes have incorrect total integer volume')
    # Batch the pair intersections to keep memory bounded. Distinct
    # integer hyperrectangles intersect the simplex iff their coordinate
    # intervals are nonempty and bracket the required total.
    for start in range(0, len(boxes), 64):
        lower = np.maximum(lo[start:start+64, None, :], lo[None, :, :])
        upper = np.minimum(hi[start:start+64, None, :], hi[None, :, :])
        overlap = np.all(lower <= upper, axis=2)&(lower.sum(axis=2) <= length)&(upper.sum(axis=2) >= length)
        for j in range(len(overlap)):
            overlap[j, start+j] = False
        if np.any(overlap):
            raise ValueError('selected type boxes overlap')
    return dict(boxes=len(boxes), exact_integer_type_count=expected, disjoint=True, complete=True)


class Refiner:
    def __init__(self, counts, block, t, s, ac, kernel, length, bands):
        self.counts, self.block, self.t, self.length = counts, block, t, length
        self.bands = bands; self.cutoff = block*length//10
        self.epochs = transfer.Epochs(t, s, ac, kernel)
        self.epoch_cache = {}

    def epoch(self, log_tilt):
        if log_tilt not in self.epoch_cache:
            if len(self.epoch_cache) >= 128:
                self.epoch_cache.clear()
            self.epoch_cache[log_tilt] = self.epochs.at(math.exp(log_tilt))
        return self.epoch_cache[log_tilt]

    def costs(self, probabilities):
        if len(probabilities) != len(self.bands)+1 or probabilities[0] != 0:
            raise ValueError('invalid band probability vector')
        costs = [0.]
        for band, p in zip(self.bands, probabilities[1:]):
            if band == [self.block] and p == 1:
                costs.append(math.log(self.counts[self.block])); continue
            if not 0 < p < 1:
                raise ValueError('invalid interior band probability')
            costs.append(max(math.log(self.counts[w])-math.log(math.comb(self.block, w))
                             -w*math.log(p)-(self.block-w)*math.log1p(-p) for w in band))
        return np.array(costs)

    def value(self, corners, log_tilt, proposal, probabilities, costs):
        theta = float(proposal@probabilities)
        moment = float(transfer.terminal_logs(transfer.epoch_mixture(self.epoch(log_tilt), self.t, [theta]),
                                               self.block*self.length//self.t)[0])
        return transfer.typed.point_logs(corners, self.length, self.block, costs, proposal,
                                         moment, self.cutoff, math.exp(log_tilt))

    def refine(self, box):
        corners = transfer.typed.vertices(box['lower'], box['upper'], self.length)
        penalty = transfer.typed.lattice_log_count(box['lower'], box['upper'])
        witness = box['witness']; probabilities = np.array(witness['probabilities'])
        costs = self.costs(probabilities)
        if np.max(np.abs(costs-witness['log_density_costs'])) > 2e-10:
            raise ValueError('stored counting costs fail replay')
        initial = np.array(witness['proposal'])
        if np.any(initial <= 0) or abs(float(initial.sum())-1) > 1e-12:
            raise ValueError('invalid stored proposal')
        tilt = float(witness['log_surprisal'])
        baseline = float(max(self.value(corners, tilt, initial, probabilities, costs)))+penalty
        if abs(baseline-box['own_log_bound']) > 2e-6:
            raise ArithmeticError('original type-box bound failed direct replay')
        anchor = int(np.argmax(initial)); free = np.array([j for j in range(len(initial)) if j != anchor])
        start = np.r_[tilt, np.log(initial/initial[anchor])[free]]
        best, selected = baseline, dict(witness)
        calls = 0
        def objective(coordinates):
            nonlocal best, selected, calls
            logits = np.zeros(len(initial)); logits[free] = coordinates[1:]
            proposal = np.exp(logits-logsumexp(logits)); z = float(coordinates[0])
            values = self.value(corners, z, proposal, probabilities, costs)
            bound = float(max(values))+penalty
            calls += 1
            if bound < best:
                best = bound
                selected = dict(log_surprisal=z, proposal=proposal.tolist(), probabilities=probabilities.tolist(),
                                log_density_costs=costs.tolist(), continuous_joint_search=True,
                                maximum_vertex=corners[int(np.argmax(values))].tolist())
            return (float(logsumexp(values))+penalty)/(self.block*self.length)
        # The positive transfer remains valid at every proposed witness;
        # optimizer convergence is not needed for validity of the bound.
        solution = minimize(objective, start, method='L-BFGS-B',
                            bounds=[(max(-18., tilt-.6), min(2., tilt+.6))]+[(-28., 28.)]*len(free),
                            options={'maxiter': 100, 'ftol': 1e-14, 'gtol': 1e-10, 'maxls': 30})
        objective(solution.x)
        return dict(lower=box['lower'], upper=box['upper'], own_log_bound=best, witness=selected,
                    replayed_log_bound=baseline, improvement_bits=(baseline-best)/LN2,
                    objective_evaluations=calls, optimizer_success=bool(solution.success))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=HERE/'bch_occupation_probe_b128_t64_s20_e20.json')
    parser.add_argument('--limit', type=int, default=8)
    args = parser.parse_args()
    _, counts, maps, dependencies = study.load_inputs()
    payload = json.loads(args.input.read_text()); settings = payload['arguments']
    block, t, s, exponent = (settings[k] for k in ('block', 'step', 'state', 'exponent'))
    length = (1 << exponent)//(block//2); config = maps[t, s]
    dense = payload['dense']; boxes = dense['selected_boxes']
    cover = check_cover(boxes, length, dense['occupation_min'])
    print('cover verified', cover, flush=True)
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    model = Refiner(counts[block], block, t, s, ac, config['kernel_counts'], length, dense['bands'])
    selected = sorted(enumerate(boxes), key=lambda pair: pair[1]['own_log_bound'], reverse=True)[:args.limit]
    results = []
    for index, box in selected:
        result = model.refine(box); result['source_box_index'] = index
        results.append(result)
        print(f'box {index}, Q={length-box["upper"][0]}..{length-box["lower"][0]}: '
              f'{-result["replayed_log_bound"]/LN2:.4f} -> {-result["own_log_bound"]/LN2:.4f} bits; '
              f'{result["objective_evaluations"]} evaluations', flush=True)
    for path in (Path(__file__), Path(transfer.__file__), args.input):
        dependencies[path.resolve().relative_to(study.ROOT).as_posix()] = study.sha(path)
    output = dict(status='SELECTED_TYPE_BOX_REFINEMENT', source_sha256=dependencies,
                  complete_input_cover=cover, boxes=results,
                  limitations=['Only the listed boxes have been replayed and refined.',
                               'No full-tail dominance follows from this subset.', 'Nearest binary64 diagnostic.'])
    (HERE/'bch_dense_continuous_probe.json').write_text(json.dumps(output, indent=2)+'\n')


if __name__ == '__main__':
    main()
