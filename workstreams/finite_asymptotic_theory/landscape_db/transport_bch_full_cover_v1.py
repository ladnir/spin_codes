"""Replay fixed BCH occupation witnesses at other state sizes, then verify.

B, T and K remain fixed. Each target has its own transfer coefficients,
component refinement and full replay; no monotonicity is assumed. Existing
transport receipts must match their authenticated inputs before reuse.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np

import occupation_refresh_v1 as transfer
import occupation_positive_poly_v1 as positive
import occupation_composition_positive_v1 as fast
import refine_bch_dense_v1 as dense
import syndrome_density_v1 as density
import study_bch_dominance_v1 as study
import trim_bch_dense_v3 as partitions

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


def authenticate(data):
    for name, digest in data['source_sha256'].items():
        if study.sha(study.ROOT/name) != digest:
            raise ValueError(f'Changed dependency: {name}')


def select_covers(reference):
    """Select the exact intervals aggregated in the verified reference."""
    geometry = tuple(reference[k] for k in
                     ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
    candidates = []
    for name in reference['source_sha256']:
        path = study.ROOT/name
        if path.name != 'cover.json':
            continue
        data = json.loads(path.read_text())
        if data.get('status') == 'BINARY64_COMPLETE_SPARSE_INTERVAL':
            current = tuple(data['arguments'][k] for k in ('block', 'step', 'state', 'exponent'))
        elif data.get('status') == 'BINARY64_COMPLETE_DENSE_INTERVAL':
            current = tuple(data[k] for k in
                            ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
        else:
            continue
        if current == geometry:
            candidates.append((path.resolve(), data))
    result = []
    for interval in reference['components']:
        if interval['occupation_min'] < 5:
            continue
        matches = [(p, d) for p, d in candidates
                   if all(d[k] == interval[k] for k in ('occupation_min', 'occupation_max', 'log_upper'))]
        if len(matches) != 1:
            raise ValueError('Cannot uniquely identify a selected source interval')
        result.append(matches[0])
    return result


def sparse_cover(source, counts, block, t, state, length, config):
    data = copy.deepcopy(source)
    maximum = data['occupation_max']; cutoff = block*length//10
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    epoch = transfer.Epochs(t, state, ac, config['kernel_counts'])
    models = [fast.CompositionBoxes(counts, block, data['bands'], p, maximum)
              for p in data['probability_banks']]
    regions = {}; adaptive = {}

    def region(z):
        if z not in regions:
            regions[z] = positive.region_logs(epoch.at(math.exp(z), maximum), t, length, maximum)
        return regions[z]

    rows = data['occupations']
    if [r['occupation'] for r in rows] != list(range(data['occupation_min'], maximum+1)):
        raise ValueError('Source sparse interval has missing occupations')
    for row in rows:
        q = row['occupation']; z = row['witness_log_tilt']; bank = row['witness_probability_bank']
        key = (z, bank)
        if key not in adaptive:
            current = region(z); matrices = []
            for _ in range(maximum):
                current = models[bank].envelope.apply(current[:-1], current[1:])
                matrices.append(current[0])
            adaptive[key] = transfer.terminal_logs(np.array(matrices), block)
        value = math.log(math.comb(length, q))+q*math.log(len(data['bands']))
        value += cutoff*math.exp(z)+float(adaptive[key][q-1])
        row['adaptive_log_upper'] = value
        cover = row['composition_cover']
        if cover is not None:
            boxes = cover['boxes']
            partitions.root_partition(dict(lower=[0]*len(data['bands']), upper=[q]*len(data['bands'])), boxes, q)
            for box in boxes:
                witness = box['witness']; tilt = witness['log_tilt']
                model = models[witness['probability_bank']]
                box['log_bound'] = model.bound(region(tilt), box['lower'], box['upper'], q,
                                               length, cutoff, math.exp(tilt))
            cover['log_union_upper'] = float(np.logaddexp.reduce([b['log_bound'] for b in boxes]))
            # Search statistics describe the source search, not a new target search.
            for key in list(cover):
                if key not in ('boxes', 'log_union_upper'):
                    del cover[key]
            value = min(value, cover['log_union_upper'])
        row['log_upper'] = value
        if q == data['occupation_min'] or q % 64 == 0:
            print(f'  transported Q{q}: {-value/LN2:.6f} bits', flush=True)
    data['arguments']['state'] = state
    return data


def dense_cover(source, counts, block, t, state, length, config):
    data = copy.deepcopy(source)
    if data.get('transfer_model') != 'syndrome_density_v1':
        raise ValueError('Transport requires the activation-density dense transfer')
    if data['occupation_max'] != length:
        raise ValueError('Dense source has incompatible K')
    leaves = [leaf for root in data['roots'] for leaf in root['leaves']]
    dense.check_cover(leaves, length, data['occupation_min'])
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    model = dense.Refiner(counts, block, t, state, ac, config['kernel_counts'], length, data['bands'])
    model.epochs = density.Epochs(t, state, ac, config['kernel_counts'])
    for root in data['roots']:
        for leaf in root['leaves']:
            witness = leaf['witness']; probabilities = np.array(witness['probabilities'])
            costs = model.costs(probabilities)
            if max(abs(costs-witness['log_density_costs'])) > 2e-10:
                raise ValueError('Source dense counting costs changed')
            corners = transfer.typed.vertices(leaf['lower'], leaf['upper'], length)
            values = model.value(corners, witness['log_surprisal'], np.array(witness['proposal']), probabilities, costs)
            witness['maximum_vertex'] = corners[int(np.argmax(values))].tolist()
            leaf['own_log_bound'] = float(max(values))+transfer.typed.lattice_log_count(leaf['lower'], leaf['upper'])
        root['log_upper'] = float(np.logaddexp.reduce([leaf['own_log_bound'] for leaf in root['leaves']]))
    data['state_bits'] = state; data['arguments']['state'] = state
    data['activation_density'] = model.epochs.activation
    data.pop('refinements', None)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--state', type=int, action='append', required=True)
    args = parser.parse_args()
    reference_path = args.reference.resolve()
    reference = json.loads(reference_path.read_text())
    if reference['status'] != 'VERIFIED_BINARY64_FULL_REFERENCE':
        raise ValueError('A verified full source reference is required')
    authenticate(reference)
    sources = select_covers(reference)
    block, t, source_state, exponent = (reference[k] for k in
                  ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
    length = (1 << exponent)//(block//2)
    _, spectra, maps, dependencies = study.load_inputs()
    dependencies.update(reference['source_sha256'])
    for path in (Path(__file__), reference_path, Path(fast.__file__), Path(dense.__file__),
                 Path(density.__file__), Path(positive.__file__), Path(transfer.__file__),
                 Path(partitions.__file__)):
        dependencies[str(path.resolve())] = study.sha(path)
    # Validate the entire request before any producer runs.
    for state in args.state:
        if state == source_state or (t, state) not in maps:
            raise ValueError('Target must be an available map distinct from the source')
        if not (HERE/f'bch_dominance_v1/b{block}_t{t}_s{state}_e{exponent}.json').exists():
            raise ValueError('Target needs the original Q2..4 checkpoint')
    for state in dict.fromkeys(args.state):
        tag = f'b{block}_t{t}_s{state}_e{exponent}'
        print(tag, flush=True)
        paths = []; dense_path = None; dense_minimum = length+1
        for path, source in sources:
            is_dense = source['status'] == 'BINARY64_COMPLETE_DENSE_INTERVAL'
            kind = 'dense' if is_dense else 'sparse'
            directory = HERE/f'bch_{kind}_transport_v1_{tag}_q{source["occupation_min"]}_{source["occupation_max"]}'
            directory.mkdir(exist_ok=True)
            target = directory/'cover.json'
            origin = dict(reference=str(reference_path), reference_sha256=study.sha(reference_path),
                          interval=str(path), interval_sha256=study.sha(path), source_state=source_state,
                          target_state=state)
            if target.exists():
                result = json.loads(target.read_text())
                if result.get('transport') != origin or result['source_sha256'] != dependencies:
                    raise ValueError('Existing transported cover has different inputs')
                authenticate(result)
            else:
                fn = dense_cover if is_dense else sparse_cover
                result = fn(source, spectra[block], block, t, state, length, maps[t, state])
                terms = ([r['log_upper'] for r in result['roots']] if is_dense else
                         [r['log_upper'] for r in result['occupations']])
                total = float(np.logaddexp.reduce(terms))
                result.update(log_upper=total, margin_bits=-total/LN2, source_sha256=dependencies,
                              transport=origin, limitations=['Fixed witnesses replayed at the target map.',
                                  'Only the displayed interval is covered.', 'No outward arithmetic.'])
                result.pop('input_fingerprint', None)
                temporary = target.with_suffix('.json.tmp')
                temporary.write_text(json.dumps(result, indent=2)+'\n'); temporary.replace(target)
            print(f'  {kind} interval: {result["margin_bits"]:.6f} bits', flush=True)
            if is_dense:
                if dense_path is not None:
                    raise ValueError('Multiple dense intervals are unsupported')
                dense_path = target; dense_minimum = source['occupation_min']
            else:
                paths.append(target)
        refined = HERE/f'bch_dominance_v2/{tag}.json'
        if not refined.exists():
            with (HERE/f'bch_transport_v1_{tag}_refine.log').open('w') as output:
                subprocess.run([sys.executable, '-u', str(HERE/'refine_bch_sparse_v2.py'),
                                '--checkpoint', f'{tag}.json'], cwd=HERE,
                               stdout=output, stderr=subprocess.STDOUT, check=True)
        command = [sys.executable, '-u', str(HERE/'verify_bch_full_reference_v5.py'),
                   '--block', str(block), '--step', str(t), '--state', str(state), '--exponent', str(exponent),
                   '--dense-minimum', str(dense_minimum), '--sparse-checkpoint', str(refined)]
        for path in paths:
            command.extend(['--sparse-cover', str(path)])
        if dense_path:
            command.extend(['--dense-cover', str(dense_path)])
        reference_target = HERE/f'bch_full_reference_{tag}.json'
        if reference_target.exists():
            history = HERE/'bch_full_reference_history'; history.mkdir(exist_ok=True)
            archive = history/f'{reference_target.stem}_{study.sha(reference_target)}.json'
            if not archive.exists():
                archive.write_bytes(reference_target.read_bytes())
        with (HERE/f'bch_transport_v1_{tag}_verify.log').open('w') as output:
            subprocess.run(command, cwd=HERE, stdout=output, stderr=subprocess.STDOUT, check=True)
        final = json.loads(reference_target.read_text())
        print(f'  VERIFIED full {final["full_margin_bits"]:.9f}; loss {final["margin_penalty_bits"]:.9g} bits', flush=True)


if __name__ == '__main__':
    main()
