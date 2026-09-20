"""Fine BCH occupation intervals with scaled positive region products."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

import occupation_refresh_v1 as transfer
import occupation_positive_poly_v1 as positive
import occupation_composition_lazy_v1 as lazy
import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--block', type=int, choices=(64, 128), default=128)
    parser.add_argument('--step', type=int, default=64)
    parser.add_argument('--state', type=int, default=20)
    parser.add_argument('--exponent', type=int, default=20)
    parser.add_argument('--maximum', type=int, default=64)
    parser.add_argument('--minimum', type=int, default=5)
    parser.add_argument('--nodes', type=int, default=127)
    parser.add_argument('--target-bits', type=float, default=55.)
    args = parser.parse_args()
    _, spectra, maps, dependencies = study.load_inputs()
    block, t, s, exponent = args.block, args.step, args.state, args.exponent
    counts, config = spectra[block], maps[t, s]
    length = (1 << exponent)//(block//2); maximum = min(length, args.maximum)
    if args.minimum < 2 or maximum < args.minimum:
        raise ValueError('empty requested sparse tail')
    cutoff = block*length//10
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    epoch = transfer.Epochs(t, s, ac, config['kernel_counts'])
    groups = [[w for w in counts if w <= min(counts)+1],
              [w for w in counts if w > min(counts)+1 and 16*w <= 13*block],
              [w for w in counts if 16*w > 13*block]]
    bands = [g for g in groups if g]
    models = [lazy.CompositionBoxes(counts, block, bands,
              [p for g, p in zip(groups, [low, middle, .88]) if g], maximum)
              for low in (.25, .35, .5, .65, .8, .9) for middle in (.5, .65)]
    tilts = np.unique(np.r_[np.arange(25, 226)/25-math.log(length), np.arange(-30, 6)/5])
    choose = np.array([math.log(math.comb(length, q))+q*math.log(len(bands)) for q in range(1, maximum+1)])
    best = np.full(maximum, np.inf); selected = np.zeros((maximum, 2), dtype=int)
    directory = HERE/f'bch_sparse_tail_v3_b{block}_t{t}_s{s}_e{exponent}_q{args.minimum}_{maximum}'
    directory.mkdir(exist_ok=True)
    for path in (Path(__file__), Path(transfer.__file__), Path(positive.__file__), Path(lazy.__file__)):
        dependencies[path.relative_to(study.ROOT).as_posix()] = study.sha(path)
    fingerprint = hashlib.sha256(json.dumps([dependencies, vars(args)], sort_keys=True).encode()).hexdigest()
    regions = []
    for ordinal, tilt in enumerate(tilts):
        path = directory/f'tilt_{ordinal}.json'
        if path.exists():
            cached = json.loads(path.read_text())
            if cached['input_fingerprint'] != fingerprint:
                raise ValueError('stale sparse-tail checkpoint')
            region = np.array(cached['regions']); values = np.array(cached['values'])
        else:
            lam = math.exp(tilt)
            region = positive.region_logs(epoch.at(lam, maximum), t, length, maximum)
            values = []
            for model in models:
                current = region; matrices = []
                for q in range(1, maximum+1):
                    current = model.envelope.apply(current[:-1], current[1:]); matrices.append(current[0])
                values.append(choose+cutoff*lam+transfer.terminal_logs(np.array(matrices), block))
            values = np.array(values)
            path.write_text(json.dumps(dict(input_fingerprint=fingerprint, log_tilt=float(tilt),
                                           regions=region.tolist(), values=values.tolist()))+'\n')
        regions.append(region)
        banks = np.argmin(values, axis=0); candidate = values[banks, np.arange(maximum)]
        improve = candidate < best; best[improve] = candidate[improve]
        selected[improve, 0] = ordinal; selected[improve, 1] = banks[improve]
        if ordinal % 25 == 0:
            print(f'fine sparse tilt {ordinal+1}/{len(tilts)}', flush=True)
    results = []
    for q in range(args.minimum, maximum+1):
        tilt_index, bank = selected[q-1]
        value = float(best[q-1]); cover = None
        if value > -args.target_bits*LN2:
            candidate_indices = [i for i, z in enumerate(tilts) if abs(z-tilts[tilt_index]) <= .081]
            # Keep multiple probability choices fixed within each box.
            banks = sorted(set((int(bank), 0, 2, 4, 6, 8, 10)))
            witnesses = [(models[k], regions[i], math.exp(tilts[i]),
                          dict(log_tilt=float(tilts[i]), probability_bank=k))
                         for i in candidate_indices for k in banks]
            cover = transfer.boxes.search(witnesses, q, length, cutoff,
                                           maximum_nodes=args.nodes, target_bits=args.target_bits)
            value = min(value, cover['log_union_upper'])
        results.append(dict(occupation=q, log_upper=value, adaptive_log_upper=float(best[q-1]),
                            witness_log_tilt=float(tilts[tilt_index]), witness_probability_bank=int(bank),
                            composition_cover=cover))
        print(f'Q{q}: {-value/LN2:.4f} bits', flush=True)
    total = float(np.logaddexp.reduce([r['log_upper'] for r in results]))
    payload = dict(status='BINARY64_COMPLETE_SPARSE_INTERVAL', source_sha256=dependencies,
                   input_fingerprint=fingerprint, arguments=vars(args), bands=bands,
                   probability_banks=[m.probabilities.tolist() for m in models],
                   occupation_min=args.minimum, occupation_max=maximum, log_upper=total, margin_bits=-total/LN2,
                   occupations=results, limitations=['Only the displayed occupation interval is covered.',
                                                   'No outward arithmetic.'])
    (directory/'cover.json').write_text(json.dumps(payload, indent=2)+'\n')
    print('sparse interval margin', -total/LN2, flush=True)


if __name__ == '__main__':
    main()
