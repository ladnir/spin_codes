"""Resumable complete-tail refinement with continuous witnesses and type splits."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import refine_bch_dense_v1 as base
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


class CoverRefiner(base.Refiner):
    def __init__(self, *args):
        super().__init__(*args)
        self.warm = []

    def direct(self, lower, upper, witness):
        corners = base.transfer.typed.vertices(lower, upper, self.length)
        p = np.array(witness['probabilities']); costs = self.costs(p)
        proposal = np.array(witness['proposal'])
        if np.any(proposal <= 0) or abs(float(proposal.sum())-1) > 1e-12:
            raise ValueError('invalid positive type proposal')
        if max(abs(costs-witness['log_density_costs'])) > 2e-10:
            raise ValueError('density costs fail replay')
        values = self.value(corners, witness['log_surprisal'], proposal, p, costs)
        return float(max(values))+base.transfer.typed.lattice_log_count(lower, upper)

    def refine_box(self, box, target):
        original = self.direct(box['lower'], box['upper'], box['witness'])
        if abs(original-box['own_log_bound']) > 2e-6:
            raise ArithmeticError('input box bound failed replay')
        current = dict(box, own_log_bound=original)
        if original <= target:
            return current
        corners = base.transfer.typed.vertices(box['lower'], box['upper'], self.length)
        center = np.maximum(corners.mean(axis=0), .25)
        q = self.length-center[0]
        candidates = []
        for prior, old_center in self.warm[-12:]:
            old_q = self.length-old_center[0]
            if q <= 0 or old_q <= 0 or not .125 < q/old_q < 8:
                continue
            proposal = np.array(prior['proposal'])*center/old_center
            proposal /= proposal.sum()
            candidates.append(dict(prior, proposal=proposal.tolist(),
                                   log_surprisal=prior['log_surprisal']+math.log(q/old_q)))
        # Cheap phase-shifted witnesses bridge the coarse original tilt
        # grid before invoking a continuous optimizer.
        for shift in (.25, .5, .75):
            prior = box['witness']; proposal = np.array(prior['proposal'])
            proposal[1:] *= math.exp(shift); proposal /= proposal.sum()
            candidates.append(dict(prior, proposal=proposal.tolist(), log_surprisal=prior['log_surprisal']+shift))
        for witness in candidates:
            value = self.direct(box['lower'], box['upper'], witness)
            if value < current['own_log_bound']:
                current = dict(box, own_log_bound=value, witness=witness)
            if current['own_log_bound'] <= target:
                break
        if current['own_log_bound'] > target:
            result = super().refine(current)
            current = dict(box, own_log_bound=result['own_log_bound'], witness=result['witness'])
        self.warm.append((current['witness'], center))
        if len(self.warm) > 64:
            self.warm.pop(0)
        return current

    def root(self, box, target, maximum_nodes):
        remaining = maximum_nodes
        evaluated = 0
        def recurse(current):
            nonlocal remaining, evaluated
            remaining -= 1; evaluated += 1
            current = self.refine_box(current, target)
            if current['own_log_bound'] <= target or remaining < 2:
                return current['own_log_bound'], [current]
            corners = base.transfer.typed.vertices(current['lower'], current['upper'], self.length)
            children = base.transfer.typed.split_box(current['lower'], current['upper'], corners)
            if not children:
                return current['own_log_bound'], [current]
            leaves = []; values = []
            # Splitting is an exact disjoint integer partition. Replay
            # its volumes rather than relying on the search tree's intent.
            parent_volume = base.lattice_count(current['lower'], current['upper'], self.length)
            prepared = []
            for lo, hi in children:
                lo, hi = list(map(int, lo)), list(map(int, hi))
                if base.lattice_count(lo, hi, self.length):
                    prepared.append(dict(lower=lo, upper=hi, witness=current['witness'],
                                         own_log_bound=self.direct(lo, hi, current['witness'])))
            if sum(base.lattice_count(c['lower'], c['upper'], self.length) for c in prepared) != parent_volume:
                raise ArithmeticError('type split changed the integer volume')
            for child in prepared:
                if remaining:
                    value, chosen = recurse(child)
                else:
                    value, chosen = child['own_log_bound'], [child]
                values.append(value); leaves.extend(chosen)
            combined = float(np.logaddexp.reduce(values))
            if combined < current['own_log_bound']:
                return combined, leaves
            return current['own_log_bound'], [current]
        value, leaves = recurse(box)
        return dict(log_upper=value, leaves=leaves, evaluated_nodes=evaluated)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=HERE/'bch_occupation_probe_b128_t64_s20_e20.json')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--target-bits', type=float, default=60.)
    parser.add_argument('--nodes-per-root', type=int, default=31)
    args = parser.parse_args()
    _, counts, maps, dependencies = study.load_inputs()
    payload = json.loads(args.input.read_text()); settings = payload['arguments']
    block, t, s, exponent = (settings[k] for k in ('block', 'step', 'state', 'exponent'))
    length = (1 << exponent)//(block//2); config = maps[t, s]
    dense = payload['dense']; boxes = dense['selected_boxes']
    cover = base.check_cover(boxes, length, dense['occupation_min'])
    print('input cover verified', cover, flush=True)
    for path in (Path(__file__), Path(base.__file__), Path(base.transfer.__file__), args.input):
        dependencies[path.resolve().relative_to(study.ROOT).as_posix()] = study.sha(path)
    bound_arguments = dict(target_bits=args.target_bits, nodes_per_root=args.nodes_per_root)
    fingerprint = hashlib.sha256(json.dumps([dependencies, bound_arguments], sort_keys=True).encode()).hexdigest()
    directory = HERE/f'bch_dense_v2_b{block}_t{t}_s{s}_e{exponent}'
    directory.mkdir(exist_ok=True)
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    model = CoverRefiner(counts[block], block, t, s, ac, config['kernel_counts'], length, dense['bands'])
    selected = sorted(enumerate(boxes), key=lambda pair: pair[1]['own_log_bound'], reverse=True)
    output = []; processed = 0
    for ordinal, (index, box) in enumerate(selected, 1):
        path = directory/f'root_{index}.json'
        if path.exists():
            result = json.loads(path.read_text())
            if result['input_fingerprint'] != fingerprint:
                raise ValueError('stale dense root checkpoint')
            if result['source_box_index'] != index:
                raise ValueError('wrong dense root identity')
        elif args.limit is None or processed < args.limit:
            result = model.root(box, -args.target_bits*LN2, args.nodes_per_root)
            result.update(source_box_index=index, input_fingerprint=fingerprint)
            path.write_text(json.dumps(result, indent=2)+'\n')
            processed += 1
        else:
            # Retain a directly replayed original bound for every untouched
            # root, so the output is always a complete occupation cover.
            value = model.direct(box['lower'], box['upper'], box['witness'])
            if abs(value-box['own_log_bound']) > 2e-6:
                raise ArithmeticError('unrefined root failed replay')
            result = dict(source_box_index=index, log_upper=value, leaves=[box], refined=False)
        output.append(result)
        if ordinal <= 8 or ordinal % 100 == 0:
            print(f'{ordinal}/{len(boxes)} root {index}: {-result["log_upper"]/LN2:.3f} bits', flush=True)
    total = float(np.logaddexp.reduce([r['log_upper'] for r in output]))
    result = dict(status='BINARY64_COMPLETE_DENSE_COVER', input_fingerprint=fingerprint,
                  source_sha256=dependencies, arguments=bound_arguments, input_cover=cover,
                  occupation_min=dense['occupation_min'], occupation_max=length, block_bits=block,
                  step_bits=t, state_bits=s, message_exponent=exponent, log_upper=total,
                  margin_bits=-total/LN2, roots=output,
                  remaining_weak_roots=sum(r['log_upper'] > -args.target_bits*LN2 for r in output),
                  limitations=['Only the displayed occupation interval is covered.',
                               'No outward arithmetic or full Q1 dominance inference without the missing sparse occupations.'])
    (directory/'cover.json').write_text(json.dumps(result, indent=2)+'\n')
    print('complete dense interval margin', -total/LN2, 'weak roots', result['remaining_weak_roots'], flush=True)


if __name__ == '__main__':
    main()
