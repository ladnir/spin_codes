"""Plot coverage-qualified BCH occupation margins and aggregation costs."""
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
    receipt = json.loads((HERE/'bch_dominance_sparse.json').read_text())
    if study.sha(HERE/'bch_dominance_sparse.csv') != receipt['csv_sha256']:
        raise ValueError('changed sparse CSV')
    with (HERE/'bch_dominance_sparse.csv').open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    integer = ('block_bits', 'step_bits', 'state_bits', 'message_exponent')
    for row in rows:
        for key in integer:
            row[key] = int(row[key])
        for key in row.keys()-set(integer)-{'full_dominance_status'}:
            row[key] = float(row[key])
    def select(**values):
        return [r for r in rows if all(r[k] == v for k, v in values.items())]
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), layout='constrained')
    for index, block in enumerate((64, 128)):
        ax = axes[index, 0]
        data = sorted(select(block_bits=block, step_bits=64, state_bits=20), key=lambda r: r['message_exponent'])
        x = [r['message_exponent'] for r in data]
        for name, label, color in (('q1_margin_bits', 'Q1', '#2563a6'),
                                    ('q2_q4_margin_bits', 'Sum of Q2–Q4 bounds', '#cc742c')):
            ax.plot(x, [r[name] for r in data], '.-', label=label, color=color)
        ax.set(title=f'BCH-{block}: message-length slice', xlabel='log₂ K   (t=64, s=20)', ylabel='Contribution margin (bits)')
        ax.legend(fontsize=9); ax.grid(alpha=.2)
        ax = axes[index, 1]
        for t, color in ((64, '#2563a6'), (128, '#cc742c'), (256, '#3c9071')):
            data = sorted(select(block_bits=block, step_bits=t, message_exponent=20), key=lambda r: r['state_bits'])
            ax.plot([r['state_bits'] for r in data], [r['q2_q4_to_q1_log2_ratio'] for r in data], '.-', label=f't={t}', color=color)
        ax.axhline(math.log2(.1), color='#555', linestyle='--', linewidth=1, label='10% of Q1 bound')
        ax.set(title=f'BCH-{block}: sparse dominance', xlabel='State bits s   (K=2²⁰)', ylabel='log₂[(Q2–Q4 bound) / (Q1 bound)]')
        ax.legend(fontsize=8); ax.grid(alpha=.2)
        ax = axes[index, 2]
        exponents = (16, 18, 20, 22, 24); states = (10, 12, 16, 20)
        values = np.full((len(states), len(exponents)), np.nan)
        for i, s in enumerate(states):
            for j, e in enumerate(exponents):
                match = select(block_bits=block, step_bits=64, state_bits=s, message_exponent=e)
                if match:
                    values[i, j] = match[0]['q1_q4_penalty_bits']
        mesh = ax.pcolormesh(np.arange(15, 26, 2), [9, 11, 14, 18, 22], values, cmap='YlOrRd', vmin=0, shading='flat')
        for i, s in enumerate(states):
            center = ([9, 11, 14, 18][i]+[11, 14, 18, 22][i])/2
            for j, e in enumerate(exponents):
                if np.isfinite(values[i, j]):
                    ax.text(e, center, f'{values[i,j]:.2g}', ha='center', va='center', fontsize=9)
        ax.set(title=f'BCH-{block}: cost of adding Q2–Q4', xlabel='log₂ K   (t=64)', ylabel='State bits s', xticks=exponents)
        ax.set_yticks([10, 12.5, 16, 20], labels=states)
        fig.colorbar(mesh, ax=ax, label='Added margin loss (bits)', shrink=.8)
    fig.suptitle('BCH occupation audit: Q2–Q4 across the engineering grid\nQ5 and higher remain a separate obligation; these panels do not establish full-tail dominance.', fontsize=15)
    fig.savefig(HERE/'bch_sparse_dominance.png', dpi=180)
    summary = dict(rows=len(rows), sparse_ratio_at_most_tenth=sum(r['q2_q4_to_q1_log2_ratio'] <= math.log2(.1) for r in rows),
                   maximum_sparse_penalty_bits=max(r['q1_q4_penalty_bits'] for r in rows),
                   reference_points=select(step_bits=64, state_bits=20, message_exponent=20),
                   full_tail_status='UNRESOLVED', source_receipt_sha256=study.sha(HERE/'bch_dominance_sparse.json'))
    (HERE/'bch_dominance_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
