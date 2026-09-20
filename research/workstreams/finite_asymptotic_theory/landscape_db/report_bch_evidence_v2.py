"""Join BCH contribution surfaces with complete bounds and independent obstructions."""
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import study_bch_dominance_v1 as study

HERE = Path(__file__).resolve().parent


def main():
    lower_path = HERE/'bch_zero_state_lower_grid_v2.json'
    lower = json.loads(lower_path.read_text())
    verification_path = HERE/'bch_zero_state_verification_v2.json'
    verification = json.loads(verification_path.read_text())
    if verification['source_sha256'][str(lower_path)] != study.sha(lower_path):
        raise ValueError('zero-state verification is stale')
    receipt = json.loads((HERE/'bch_dominance_sparse.json').read_text())
    if study.sha(HERE/'bch_dominance_sparse.csv') != receipt['csv_sha256']:
        raise ValueError('sparse CSV changed')
    with (HERE/'bch_dominance_sparse.csv').open(newline='') as handle:
        sparse = list(csv.DictReader(handle))
    def key(r):
        return tuple(int(r[k]) for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent'))
    lower_index = {key(r): r for r in lower['rows']}
    extra_inputs = []
    exact_path = HERE/'bch_zero_state_exact_grid_v3.json'
    exact_verification = HERE/'bch_zero_state_verification_v3.json'
    if exact_path.exists() and exact_verification.exists():
        replay = json.loads(exact_verification.read_text())
        if replay['source_sha256'][str(exact_path)] != study.sha(exact_path):
            raise ValueError('exact-coefficient lower-bound verification is stale')
        for row in json.loads(exact_path.read_text())['rows']:
            old = lower_index.get(key(row), {}).get('first_moment_lower_bits')
            new = row.get('first_moment_lower_bits')
            if new is not None and (old is None or new > old):
                lower_index[key(row)] = row
        extra_inputs = [exact_path, exact_verification]
    full_index = {}
    full_paths = sorted(HERE.glob('bch_full_reference_b*_t*_s*_e*.json'))
    for path in full_paths:
        data = json.loads(path.read_text())
        if data['status'] != 'VERIFIED_BINARY64_FULL_REFERENCE':
            raise ValueError('unverified full reference')
        for name, digest in data['source_sha256'].items():
            if study.sha(Path(name)) != digest:
                raise ValueError('full-reference verification is stale')
        full_index[key(data)] = data
    rows = []
    for row in sparse:
        geometry = key(row); low = lower_index[geometry]; full = full_index.get(geometry)
        lower_bits = low.get('first_moment_lower_bits')
        obstructed = lower_bits is not None and lower_bits > 0
        if obstructed and full and full['full_margin_bits'] > 0:
            raise ArithmeticError('positive upper and lower evidence contradict each other')
        result = {k: int(row[k]) for k in ('block_bits', 'step_bits', 'state_bits', 'message_exponent')}
        result.update(q1_margin_bits=float(row['q1_margin_bits']), q2_q4_margin_bits=float(row['q2_q4_margin_bits']),
                      sparse_aggregation_penalty_bits=float(row['q1_q4_penalty_bits']),
                      first_moment_lower_bits=lower_bits, full_margin_bits=full['full_margin_bits'] if full else None,
                      full_higher_to_q1_ratio=full['higher_to_q1_ratio'] if full else None,
                      full_aggregation_penalty_bits=full['margin_penalty_bits'] if full else None,
                      evidence=('FULL_BOUND_Q1_DOMINANT' if full and full['higher_to_q1_ratio'] <= .1
                                else 'FULL_BOUND_TAIL_MATERIAL' if full else 'FIRST_MOMENT_OBSTRUCTION' if obstructed else 'SPARSE_ONLY'))
        rows.append(result)
    path = HERE/'bch_engineering_evidence_v2.csv'
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout='constrained')
    colors = {64: '#2563a6', 128: '#ce782e', 256: '#a33f55'}
    for column, block in enumerate((64, 128)):
        selected = [r for r in rows if r['block_bits'] == block and r['message_exponent'] == 20]
        top, bottom = axes[:, column]
        for t in (64, 128, 256):
            data = sorted([r for r in selected if r['step_bits'] == t], key=lambda r: r['state_bits'])
            top.plot([r['state_bits'] for r in data], [r['q1_margin_bits'] for r in data], '.-', color=colors[t], label=f't={t}')
            available = [r for r in data if r['first_moment_lower_bits'] is not None]
            if available:
                bottom.plot([r['state_bits'] for r in available],
                            [r['first_moment_lower_bits']/(1 << 21) for r in available], '.-', color=colors[t], label=f't={t}')
        for r in selected:
            if r['full_margin_bits'] is not None:
                top.scatter([r['state_bits']], [r['full_margin_bits']], marker='*', s=160, color='#25835c', zorder=5,
                            label=f'Full bound: t={r["step_bits"]}, s={r["state_bits"]}')
        top.set(title=f'BCH-{block}: Q1 barely changes with t', xlabel='State bits s', ylabel='Q1 margin (bits)')
        top.grid(alpha=.2); top.legend(fontsize=9)
        bottom.axhline(0., color='#555', linewidth=1)
        bottom.axhspan(0., .1, color='#a33f55', alpha=.08)
        bottom.set(title='Zero-state trajectories expose a different tradeoff', xlabel='State bits s',
                   ylabel='Zero-state lower exponent / output bit')
        bottom.grid(alpha=.2); bottom.legend(fontsize=9)
        bottom.text(.98, .04, 'Below zero: inconclusive', transform=bottom.transAxes, ha='right', fontsize=9, bbox=dict(facecolor='white', edgecolor='none', alpha=.9))
        bottom.text(.03, .96, 'Above zero: first-moment obstruction', transform=bottom.transAxes, va='top', fontsize=10)
    fig.suptitle('BCH engineering audit at K=2²⁰ and relative distance 10%\nQ1 dominance needs full-occupation evidence; an unavailable lower bound does not imply closure.', fontsize=14)
    fig.savefig(HERE/'bch_q1_vs_dense_tradeoff_v2.png', dpi=180)
    inputs = [Path(__file__), lower_path, verification_path, HERE/'bch_dominance_sparse.csv',
              HERE/'bch_dominance_sparse.json', *full_paths, *extra_inputs]
    summary = dict(rows=len(rows), evidence_counts={label: sum(r['evidence'] == label for r in rows) for label in sorted({r['evidence'] for r in rows})},
                   csv_sha256=study.sha(path), source_sha256={str(p): study.sha(p) for p in inputs},
                   limitations=['Full-bound conclusions apply only to geometries with full-reference evidence.',
                                'Lower bounds concern expected bad-word counts, not failure probabilities.'])
    (HERE/'bch_engineering_evidence_v2.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary['evidence_counts'], indent=2))


if __name__ == '__main__':
    main()

