"""Best-first continuous refinement across an entire saved BCH dense cover."""
import argparse
import heapq
import json
import math
import hashlib
from pathlib import Path

import numpy as np

import close_bch_dense_v2 as previous
import trim_bch_dense_v3 as checks
import study_bch_dominance_v1 as study
import density_dense_witness_v1 as fast
import entropy_type_split_v1 as splitting

HERE = Path(__file__).resolve().parent
LN2 = math.log(2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--block', type=int, default=64)
    parser.add_argument('--step', type=int, default=64)
    parser.add_argument('--state', type=int, default=20)
    parser.add_argument('--exponent', type=int, default=20)
    parser.add_argument('--minimum', type=int, default=257)
    parser.add_argument('--target-bits', type=float, default=25.)
    parser.add_argument('--maximum-refinements', type=int, default=8192)
    args = parser.parse_args()
    _, spectra, maps, dependencies = study.load_inputs()
    block, t, s, e = args.block, args.step, args.state, args.exponent
    length = (1 << e)//(block//2); config = maps[t, s]; counts = spectra[block]
    seed_path = HERE/f'bch_dense_seed_split_one_b{block}_t{t}_s{s}_e{e}_q{args.minimum}.json'
    seed = json.loads(seed_path.read_text()); raw = seed['dense']
    previous.base.check_cover(raw['selected_boxes'], length, args.minimum)
    ac = {w: n for w, n in enumerate(config['a_counts']) if w and n}
    model = fast.DensityCoverRefiner(counts, block, t, s, ac, config['kernel_counts'], length, raw['bands'])
    directory = HERE/f'bch_dense_v10_b{block}_t{t}_s{s}_e{e}_q{args.minimum}'
    directory.mkdir(exist_ok=True)
    old_directory = None
    resume = directory/'checkpoint.json'
    if not resume.exists():
        resume = HERE/f'bch_dense_v9_b{block}_t{t}_s{s}_e{e}_q{args.minimum}/checkpoint.json'
    if resume.exists():
        data = resume.read_bytes(); prior = json.loads(data)
        for name, value in (('block', block), ('step', t), ('state', s), ('exponent', e), ('minimum', args.minimum)):
            if prior['arguments'][name] != value:
                raise ValueError('resume checkpoint geometry changed')
        raw = dict(raw, selected_boxes=prior['leaves'])
        previous.base.check_cover(raw['selected_boxes'], length, args.minimum)
        saved = directory/f'resume_{hashlib.sha256(data).hexdigest()}.json'
        if not saved.exists():
            saved.write_bytes(data)
        dependencies[saved.relative_to(study.ROOT).as_posix()] = study.sha(saved)
        old_directory = None
    leaves = {}; heap = []; serial = 0
    for index, root in enumerate(raw['selected_boxes']):
        path = old_directory/f'root_{index}.json' if old_directory is not None else None
        selected = [root]
        if path is not None and path.exists():
            cached = json.loads(path.read_text())
            if cached['source_seed_index'] != index:
                raise ValueError('old root identity changed')
            checks.root_partition(root, cached['leaves'], length)
            selected = cached['leaves']
            dependencies[path.relative_to(study.ROOT).as_posix()] = study.sha(path)
        for leaf in selected:
            rebound = dict(leaf, own_log_bound=model.direct(leaf['lower'], leaf['upper'], leaf['witness']))
            current = checks.canonical_witness(model, rebound)
            leaves[serial] = current; heapq.heappush(heap, (-current['own_log_bound'], serial)); serial += 1
    previous.base.check_cover(list(leaves.values()), length, args.minimum)
    for path in (Path(__file__), Path(fast.__file__), Path(splitting.__file__), Path(fast.density.__file__), Path(fast.joint.__file__), Path(fast.joint.fast.__file__), Path(previous.__file__), Path(checks.__file__), Path(previous.base.__file__),
                 Path(previous.base.transfer.__file__), seed_path):
        dependencies[path.relative_to(study.ROOT).as_posix()] = study.sha(path)
    def total():
        return float(np.logaddexp.reduce([leaf['own_log_bound'] for leaf in leaves.values()]))
    target = -args.target_bits*LN2; used = 0
    while heap and used < args.maximum_refinements and total() > target:
        _, index = heapq.heappop(heap)
        if index not in leaves:
            continue
        current = leaves[index]
        # A stronger per-leaf target makes the total insensitive to the
        # final number of boxes, but the stopping rule uses the actual sum.
        improved = model.refine_box(current, -(args.target_bits+20)*LN2)
        improved = checks.canonical_witness(model, improved)
        leaves[index] = improved; used += 1
        if total() <= target:
            break
        if improved['own_log_bound'] <= -(args.target_bits+20)*LN2:
            continue
        corners = previous.base.transfer.typed.vertices(improved['lower'], improved['upper'], length)
        children = splitting.split_box(improved['lower'], improved['upper'], corners)
        if children:
            prepared = []
            for lo, hi in children:
                lo, hi = list(map(int, lo)), list(map(int, hi))
                if not previous.base.lattice_count(lo, hi, length):
                    continue
                value = model.direct(lo, hi, improved['witness'])
                prepared.append(dict(lower=lo, upper=hi, own_log_bound=value, witness=improved['witness']))
            checks.root_partition(improved, prepared, length)
            del leaves[index]
            for child in prepared:
                leaves[serial] = child; heapq.heappush(heap, (-child['own_log_bound'], serial)); serial += 1
        if used % 50 == 0:
            print(f'{used} refinements, {len(leaves)} boxes: {-total()/LN2:.4f} bits', flush=True)
            snapshot = dict(transfer_model='syndrome_density_v1', status='PARTIAL_CONTINUOUS_DENSE_REFINEMENT', source_sha256=dependencies,
                            arguments=vars(args), refinements=used, leaves=list(leaves.values()), log_upper=total())
            (directory/'checkpoint.json').write_text(json.dumps(snapshot)+'\n')
    canonical = [checks.canonical_witness(model, leaf) for leaf in leaves.values()]
    coverage = previous.base.check_cover(canonical, length, args.minimum)
    bound = float(np.logaddexp.reduce([x['own_log_bound'] for x in canonical]))
    # Keep the root/leaves interface consumed by the full-reference replay.
    roots = [dict(log_upper=leaf['own_log_bound'], leaves=[leaf]) for leaf in canonical]
    payload = dict(transfer_model='syndrome_density_v1', activation_density=model.epochs.activation, status='BINARY64_COMPLETE_DENSE_INTERVAL', source_sha256=dependencies, arguments=vars(args),
                   block_bits=block, step_bits=t, state_bits=s, message_exponent=e,
                   occupation_min=args.minimum, occupation_max=length, bands=raw['bands'],
                   log_upper=bound, margin_bits=-bound/LN2, roots=roots, input_cover=coverage,
                   refinements=used, limitations=['Only the displayed occupation interval is covered.', 'No outward arithmetic.'])
    (directory/'cover.json').write_text(json.dumps(payload, indent=2)+'\n')
    print('best-first complete dense margin', -bound/LN2, flush=True)


if __name__ == '__main__':
    main()
