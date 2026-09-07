"""Reconstruct and refine a complete dense tail after extending exact regions."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import close_bch_dense_v2 as previous
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


def root_partition(root, leaves, length):
    lo = np.array([x['lower'] for x in leaves]); hi = np.array([x['upper'] for x in leaves])
    if np.any(lo < root['lower']) or np.any(hi > root['upper']):
        raise ValueError('cached children escape their original root')
    volumes = [previous.base.lattice_count(list(map(int, a)), list(map(int, b)), length) for a, b in zip(lo, hi)]
    if any(v <= 0 for v in volumes) or sum(volumes) != previous.base.lattice_count(root['lower'], root['upper'], length):
        raise ValueError('cached root changed its integer volume')
    for i in range(len(leaves)):
        for j in range(i):
            a = np.maximum(lo[i], lo[j]); b = np.minimum(hi[i], hi[j])
            if np.all(a <= b) and a.sum() <= length <= b.sum():
                raise ValueError('cached root children overlap')


def canonical_witness(model, box):
    old = box['witness']
    witness = {key: old[key] for key in ('log_surprisal', 'proposal', 'probabilities', 'log_density_costs')}
    corners = previous.base.transfer.typed.vertices(box['lower'], box['upper'], model.length)
    values = model.value(corners, witness['log_surprisal'], np.array(witness['proposal']),
                         np.array(witness['probabilities']), np.array(witness['log_density_costs']))
    witness['maximum_vertex'] = corners[int(np.argmax(values))].tolist()
    bound = float(max(values))+previous.base.transfer.typed.lattice_log_count(box['lower'], box['upper'])
    if abs(bound-box['own_log_bound']) > 2e-6:
        raise ArithmeticError('selected trimmed witness failed replay')
    return dict(lower=box['lower'], upper=box['upper'], own_log_bound=bound, witness=witness)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=HERE/'bch_occupation_probe_b128_t64_s20_e20.json')
    parser.add_argument('--minimum', type=int, default=257)
    parser.add_argument('--target-bits', type=float, default=60.)
    parser.add_argument('--nodes-per-root', type=int, default=31)
    args = parser.parse_args()
    _, spectra, maps, dependencies = study.load_inputs()
    payload = json.loads(args.input.read_text()); settings = payload['arguments']
    block, t, s, exponent = (settings[k] for k in ('block', 'step', 'state', 'exponent'))
    length = (1 << exponent)//(block//2); config = maps[t, s]; dense = payload['dense']
    if not dense['occupation_min'] <= args.minimum <= length:
        raise ValueError('invalid subinterval')
    previous.base.check_cover(dense['selected_boxes'], length, dense['occupation_min'])
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    model = previous.CoverRefiner(spectra[block], block, t, s, ac, config['kernel_counts'], length, dense['bands'])
    cache = HERE/f'bch_dense_v2_b{block}_t{t}_s{s}_e{exponent}'
    seeds = []
    for index, root in enumerate(dense['selected_boxes']):
        path = cache/f'root_{index}.json'
        if path.exists():
            cached = json.loads(path.read_text()); leaves = cached['leaves']
            if cached['source_box_index'] != index:
                raise ValueError('cached root identity changed')
            root_partition(root, leaves, length)
            dependencies[path.relative_to(study.ROOT).as_posix()] = study.sha(path)
        else:
            leaves = [root]
        for leaf in leaves:
            lo, hi = list(leaf['lower']), list(leaf['upper'])
            hi[0] = min(hi[0], length-args.minimum)
            if any(a > b for a, b in zip(lo, hi)) or not previous.base.lattice_count(lo, hi, length):
                continue
            # The same fixed witness bounds every clipped box. Recompute
            # both its vertex maximum and its counting factor directly.
            value = model.direct(lo, hi, leaf['witness'])
            seeds.append(dict(lower=lo, upper=hi, own_log_bound=value, witness=leaf['witness']))
        if index % 500 == 0:
            print(f'rebuilt {index+1}/{len(dense["selected_boxes"])} roots', flush=True)
    cover = previous.base.check_cover(seeds, length, args.minimum)
    print('trimmed cover verified', cover, flush=True)
    for path in (Path(__file__), Path(previous.__file__), Path(previous.base.__file__),
                 Path(previous.base.transfer.__file__), args.input):
        dependencies[path.resolve().relative_to(study.ROOT).as_posix()] = study.sha(path)
    parameters = dict(minimum=args.minimum, target_bits=args.target_bits, nodes_per_root=args.nodes_per_root)
    fingerprint = hashlib.sha256(json.dumps([dependencies, parameters], sort_keys=True).encode()).hexdigest()
    directory = HERE/f'bch_dense_v3_b{block}_t{t}_s{s}_e{exponent}_q{args.minimum}'
    directory.mkdir(exist_ok=True)
    results = []
    for ordinal, (index, seed) in enumerate(sorted(enumerate(seeds), key=lambda x: x[1]['own_log_bound'], reverse=True), 1):
        path = directory/f'root_{index}.json'
        if path.exists():
            result = json.loads(path.read_text())
            if result['input_fingerprint'] != fingerprint:
                raise ValueError('stale trimmed checkpoint')
        else:
            result = model.root(seed, -args.target_bits*LN2, args.nodes_per_root)
            result['leaves'] = [canonical_witness(model, leaf) for leaf in result['leaves']]
            root_partition(seed, result['leaves'], length)
            result.update(input_fingerprint=fingerprint, source_seed_index=index)
            path.write_text(json.dumps(result, indent=2)+'\n')
        results.append(result)
        if ordinal <= 8 or ordinal % 250 == 0:
            print(f'trimmed {ordinal}/{len(seeds)}: {-result["log_upper"]/LN2:.3f} bits', flush=True)
    total = float(np.logaddexp.reduce([r['log_upper'] for r in results]))
    output = dict(status='BINARY64_COMPLETE_DENSE_INTERVAL', source_sha256=dependencies,
                  input_fingerprint=fingerprint, parameters=parameters, input_cover=cover,
                  block_bits=block, step_bits=t, state_bits=s, message_exponent=exponent,
                  occupation_min=args.minimum, occupation_max=length, bands=dense['bands'],
                  log_upper=total, margin_bits=-total/LN2, roots=results,
                  limitations=['Only the displayed occupation interval is covered.', 'No outward arithmetic.'])
    (directory/'cover.json').write_text(json.dumps(output, indent=2)+'\n')
    print('trimmed dense margin', -total/LN2, flush=True)


if __name__ == '__main__':
    main()
