"""Apply the same Q2..4 refinement to every currently useful full BCH reference.

Archive the preceding reference, preserve its occupation intervals, and
replay the complete union after refinement. All numerical jobs are serial.
"""
import json
from pathlib import Path
import subprocess
import sys

import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent


def main():
    originals = []
    for path in sorted(HERE.glob('bch_full_reference_b*_t*_s*_e*.json')):
        data = json.loads(path.read_text())
        if data['status'] != 'VERIFIED_BINARY64_FULL_REFERENCE' or data['full_margin_bits'] <= 0:
            continue
        for source, digest in data['source_sha256'].items():
            if study.sha(Path(source)) != digest:
                raise ValueError(f'stale full reference: {path.name}')
        originals.append((path, data))
    history = HERE/'bch_full_reference_history'
    history.mkdir(exist_ok=True)
    rows = []; inputs = [Path(__file__), HERE/'refine_bch_sparse_v2.py', HERE/'verify_bch_full_reference_v4.py']
    for path, original in originals:
        b, t, s, e = (original[k] for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
        if b not in (64, 128):
            raise ValueError('unexpected block size')
        length = (1 << e)//(b//2)
        tag = f'b{b}_t{t}_s{s}_e{e}'
        archive = history/f'{path.stem}_{study.sha(path)}.json'
        if not archive.exists():
            archive.write_bytes(path.read_bytes())
        inputs.append(archive)
        has_dense = original['dense_cover'] is not None
        tail = [c for c in original['components'] if c['occupation_min'] >= 5]
        minimum = tail[-1]['occupation_min'] if has_dense else max(257, length+1)
        sparse_parts = tail[:-1] if has_dense else tail
        candidates = []
        for source in original['source_sha256']:
            candidate = Path(source)
            if candidate.name != 'cover.json':
                continue
            cover = json.loads(candidate.read_text())
            if cover.get('status') == 'BINARY64_COMPLETE_SPARSE_INTERVAL':
                geometry = cover['arguments']
                if tuple(geometry[k] for k in ('block', 'step', 'state', 'exponent')) == (b, t, s, e):
                    candidates.append((candidate, cover))
        selected = []
        for part in sparse_parts:
            matches = [candidate for candidate, cover in candidates if
                       all(cover[k] == part[k] for k in ('occupation_min', 'occupation_max', 'log_upper'))]
            if not matches:
                raise ValueError('could not reconstruct original sparse interval')
            selected.append(sorted(matches)[0])

        def run(script, arguments, stage):
            print(f'{tag}: {stage}', flush=True)
            with (HERE/f'bch_anchor_refinement_v1_{tag}_{stage}.log').open('w') as log:
                subprocess.run([sys.executable, '-u', str(HERE/script), *arguments], cwd=HERE,
                               stdout=log, stderr=subprocess.STDOUT, check=True)

        checkpoint = HERE/f'bch_dominance_v2/{tag}.json'
        if not checkpoint.exists():
            run('refine_bch_sparse_v2.py', ['--checkpoint', f'{tag}.json'], 'refine')
        arguments = ['--block', str(b), '--step', str(t), '--state', str(s), '--exponent', str(e),
                     '--dense-minimum', str(minimum), '--sparse-checkpoint', str(checkpoint)]
        for interval in selected:
            arguments += ['--sparse-cover', str(interval)]
        run('verify_bch_full_reference_v4.py', arguments, 'verify')
        current = json.loads(path.read_text())
        if current['full_margin_bits']+1e-8 < original['full_margin_bits']:
            raise ArithmeticError('sparse refinement weakened the full bound')
        rows.append(dict(block_bits=b, step_bits=t, state_bits=s, message_exponent=e,
                         previous_margin_bits=original['full_margin_bits'],
                         full_margin_bits=current['full_margin_bits'],
                         margin_penalty_bits=current['margin_penalty_bits']))
        inputs += [path, checkpoint]
        print(f"{tag}: full margin {current['full_margin_bits']:.9f}; "
              f"Q1 loss {current['margin_penalty_bits']:.9g} bits", flush=True)
        result = dict(status='PARTIAL_ANCHOR_REFINEMENT', expected_references=len(originals), rows=rows,
                      source_sha256={str(p): study.sha(p) for p in inputs})
        (HERE/'bch_full_anchor_refinement_v1.json').write_text(json.dumps(result, indent=2)+'\n')
    if not rows:
        raise ValueError('no useful full references to refine')
    result['status'] = 'VERIFIED_FULL_ANCHOR_REFINEMENT'
    (HERE/'bch_full_anchor_refinement_v1.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
