"""Join BCH contribution surfaces with complete bounds and independent obstructions."""
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
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
        q1 = float(row['q1_margin_bits'])
        sparse_margin = float(row['q2_q4_margin_bits'])
        sparse_penalty = float(row['q1_q4_penalty_bits'])
        if full:
            logs = [component['log_upper'] for component in full['components']
                    if component['occupation_min'] == component['occupation_max']
                    and component['occupation_min'] in (2, 3, 4)]
            if len(logs) != 3 or abs(q1-full['q1_margin_bits']) > 1e-8:
                raise ValueError('incompatible full reference')
            tail = float(np.logaddexp.reduce(logs))
            sparse_margin = -tail/math.log(2)
            sparse_penalty = float(np.logaddexp(-q1*math.log(2), tail))/math.log(2)+q1
        result.update(q1_margin_bits=q1, q2_q4_margin_bits=sparse_margin,
                      sparse_aggregation_penalty_bits=sparse_penalty,
                      coarse_sparse_aggregation_penalty_bits=float(row['q1_q4_penalty_bits']),
                      sparse_source='full_reference' if full else 'original_grid',
                      first_moment_lower_bits=lower_bits, full_margin_bits=full['full_margin_bits'] if full else None,
                      full_higher_to_q1_ratio=full['higher_to_q1_ratio'] if full else None,
                      full_higher_to_q1_log2_ratio=full['higher_to_q1_log2_ratio'] if full else None,
                      full_aggregation_penalty_bits=full['margin_penalty_bits'] if full else None,
                      evidence=('FIRST_MOMENT_OBSTRUCTION' if obstructed else
                                'FULL_BOUND_UNINFORMATIVE' if full and full['full_margin_bits'] <= 0 else
                                'FULL_BOUND_Q1_DOMINANT' if full and full['margin_penalty_bits'] <= math.log2(1.1) else
                                'FULL_BOUND_LARGER_PENALTY' if full else 'SPARSE_ONLY'))
        rows.append(result)
    path = HERE/'bch_engineering_evidence_v7.csv'
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
        marked = set()
        for r in selected:
            if r['full_margin_bits'] is not None:
                if r['full_margin_bits'] <= 0:
                    continue
                t = r['step_bits']
                top.scatter([r['state_bits']], [r['full_margin_bits']], marker='*', s=140, color=colors[t],
                            edgecolor='white', linewidth=.5, zorder=5,
                            label=f'Full bound: t={t}' if t not in marked else None)
                marked.add(t)
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
    fig.savefig(HERE/'bch_q1_vs_dense_tradeoff_v7.png', dpi=180)
    plt.close(fig)
    # A categorical coverage map makes untested states visible without
    # interpolating a full-margin curve through them.
    states = list(range(7, 21))
    geometries = [(b, t) for b in (64, 128) for t in (64, 128, 256)]
    entries = {key(r): r for r in rows}
    matrix = np.full((len(geometries), len(states)), np.nan)
    codes = {'SPARSE_ONLY': 0, 'FIRST_MOMENT_OBSTRUCTION': 1,
             'FULL_BOUND_UNINFORMATIVE': 2, 'FULL_BOUND_LARGER_PENALTY': 3,
             'FULL_BOUND_Q1_DOMINANT': 3}
    palette = ['#e1e5e9', '#db8791', '#edc477', '#89c9a7']
    fig, ax = plt.subplots(figsize=(13, 5), layout='constrained')
    for i, (b, t) in enumerate(geometries):
        for j, s in enumerate(states):
            r = entries.get((b, t, s, 20))
            if r is None:
                continue
            matrix[i, j] = codes[r['evidence']]
            if codes[r['evidence']] == 3:
                penalty = r['full_aggregation_penalty_bits']
                label = f'{penalty:.3f}' if penalty >= .001 else f'{penalty:.1e}'
                ax.text(j, i, label, ha='center', va='center', fontsize=9)
            elif codes[r['evidence']] == 2:
                ax.text(j, i, '?', ha='center', va='center', fontsize=12)
    ax.imshow(matrix, cmap=ListedColormap(palette), vmin=-.5, vmax=3.5, aspect='auto', interpolation='none')
    ax.set(xticks=range(len(states)), xticklabels=states,
           yticks=range(len(geometries)), yticklabels=[f'BCH-{b}, T={t}' for b, t in geometries],
           xlabel='State bits S', title='Where the complete BCH bound is understood at K=2²⁰\nNumbers give the full-bound margin loss relative to Q1, in bits')
    ax.set_xticks(np.arange(-.5, len(states)), minor=True)
    ax.set_yticks(np.arange(-.5, len(geometries)), minor=True)
    ax.grid(which='minor', color='white', linewidth=2)
    ax.tick_params(which='minor', bottom=False, left=False)
    labels = ['Only sparse evidence', 'First-moment obstruction', 'Full upper bound remains weak', 'Useful full bound']
    ax.legend(handles=[Patch(facecolor=c, label=l) for c, l in zip(palette, labels)],
              loc='upper center', bbox_to_anchor=(.5, -.18), ncol=2, frameon=False)
    fig.savefig(HERE/'bch_full_bound_coverage_v7.png', dpi=180)
    plt.close(fig)
    model_checks = []
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), layout='constrained')
    for column, block in enumerate((64, 128)):
        curve = sorted((r for r in rows if r['block_bits'] == block and r['step_bits'] == 64
                        and r['state_bits'] == 20), key=lambda r: r['message_exponent'])
        verified = [r for r in curve if r['full_margin_bits'] is not None and r['full_margin_bits'] > 0]
        top, bottom = axes[:, column]
        top.plot([r['message_exponent'] for r in curve], [r['q1_margin_bits'] for r in curve],
                 '.-', color='#9a9fa6', label='Q1 grid')
        top.scatter([r['message_exponent'] for r in verified], [r['full_margin_bits'] for r in verified],
                    marker='o', s=55, facecolor='none', edgecolor='#176c52', linewidth=1.5,
                    label='Verified useful full bound', zorder=3)
        weak = [r for r in curve if r['full_margin_bits'] is not None and r['full_margin_bits'] <= 0]
        if weak:
            top.scatter([r['message_exponent'] for r in weak], [r['q1_margin_bits'] for r in weak],
                        marker='x', color='#b57924', label='Full bound weak (marked at Q1)')
        top.set(title=f'BCH-{block}: margin versus message size', ylabel='Margin (bits)')
        bottom.scatter([r['message_exponent'] for r in verified],
                       [r['full_aggregation_penalty_bits'] for r in verified], s=50, color='#176c52')
        bottom.set(yscale='log', ylabel='Full-bound loss relative to Q1 (bits)', xlabel='log₂ K')
        for r in verified:
            bottom.annotate(f"{r['full_aggregation_penalty_bits']:.2g}",
                            (r['message_exponent'], r['full_aggregation_penalty_bits']),
                            xytext=(0, 7), textcoords='offset points', ha='center', fontsize=9)
        anchor = full_index.get((block, 64, 20, 20))
        endpoint = full_index.get((block, 64, 20, 26))
        if anchor is not None and endpoint is not None:
            def q2_log(reference):
                return next(c['log_upper'] for c in reference['components']
                            if c['occupation_min'] == c['occupation_max'] == 2)
            ln2 = math.log(2)
            length0 = (1 << 20)//(block//2)
            h1 = -anchor['q1_margin_bits']*ln2-math.log(length0)
            ratio21 = q2_log(anchor)+anchor['q1_margin_bits']*ln2-math.log((length0-1)/2)
            exponents = np.arange(20, 27)
            lengths = np.exp2(exponents)/(block//2)
            q1_predictions = h1+np.log(lengths)
            ratio_predictions = ratio21+np.log((lengths-1)/2)
            margins = -np.logaddexp(q1_predictions, q1_predictions+ratio_predictions)/ln2
            losses = np.logaddexp(0., ratio_predictions)/ln2
            top.plot(exponents, margins, '--', color='#2875aa',
                     label='Q1+Q2 model from K=2²⁰')
            bottom.plot(exponents, losses, '--', color='#2875aa', label='Q2 aggregation model')
            for comparison in verified:
                exponent = comparison['message_exponent']
                if exponent <= 20:
                    continue
                observed = full_index[(block, 64, 20, exponent)]
                current_length = (1 << exponent)//(block//2)
                observed_q1 = -observed['q1_margin_bits']*ln2
                observed_q2 = q2_log(observed)
                observed_ratio = observed_q2-observed_q1-math.log((current_length-1)/2)
                predicted_q1 = h1+math.log(current_length)
                predicted_ratio = ratio21+math.log((current_length-1)/2)
                predicted_margin = -float(np.logaddexp(predicted_q1, predicted_q1+predicted_ratio))/ln2
                two_term_log = float(np.logaddexp(observed_q1, observed_q2))
                two_term_margin = -two_term_log/ln2
                remaining_log = float(np.logaddexp.reduce([c['log_upper'] for c in observed['components']
                                                          if c['occupation_min'] >= 3]))
                # Subtracting margins loses corrections far below one ulp.
                remaining_penalty = float(np.logaddexp(0., remaining_log-two_term_log))/ln2
                model_checks.append(dict(block_bits=block, step_bits=64, state_bits=20,
                                         calibration_exponent=20, comparison_exponent=exponent,
                                         log2_h1=h1/ln2, log2_h2_over_h1=ratio21/ln2,
                                         comparison_log2_h2_over_h1=observed_ratio/ln2,
                                         normalized_ratio_change_bits=(observed_ratio-ratio21)/ln2,
                                         predicted_two_term_margin_bits=predicted_margin,
                                         observed_two_term_margin_bits=two_term_margin,
                                         verified_full_margin_bits=observed['full_margin_bits'],
                                         count_model_error_bits=predicted_margin-two_term_margin,
                                         q3_and_higher_penalty_bits=remaining_penalty,
                                         prediction_minus_full_bits=predicted_margin-observed['full_margin_bits']))
            bottom.legend(fontsize=9, loc='best')
        top.legend(fontsize=9)
        for ax in (top, bottom):
            ax.set_xticks(range(12, 27, 2)); ax.set_xlim(11.5, 26.5); ax.grid(alpha=.2)
    fig.suptitle('BCH message-size scaling: T=64, S=20, relative distance 10%\nCircles: replayed full bounds. Dashed curves: a two-term estimate calibrated at K=2²⁰.', fontsize=13)
    fig.savefig(HERE/'bch_full_bound_k_scaling_v7.png', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout='constrained')
    for column, block in enumerate((64, 128)):
        curve = sorted((r for r in rows if r['block_bits'] == block and r['step_bits'] == 64
                        and r['message_exponent'] == 20 and r['full_margin_bits'] is not None
                        and r['full_margin_bits'] > 0), key=lambda r: r['state_bits'])
        states = [r['state_bits'] for r in curve]
        top, bottom = axes[:, column]
        top.plot(states, [r['q1_margin_bits'] for r in curve], '.-', color='#9a9fa6', label='Q1')
        top.plot(states, [r['full_margin_bits'] for r in curve], 'o-', color='#176c52',
                 markersize=4, label='Replayed full bound')
        top.set(title=f'BCH-{block}', ylabel='Margin (bits)')
        top.legend(fontsize=9)
        bottom.plot(states, [r['full_aggregation_penalty_bits'] for r in curve], 'o-',
                    markersize=4, color='#176c52')
        bottom.set(yscale='log', ylabel='Full-bound loss relative to Q1 (bits)', xlabel='State bits S')
        for r in curve[::max(1, len(curve)-1)]:
            bottom.annotate(f"{r['full_aggregation_penalty_bits']:.3g}",
                            (r['state_bits'], r['full_aggregation_penalty_bits']),
                            xytext=(0, 8), textcoords='offset points', ha='center', fontsize=9)
        for ax in (top, bottom):
            ax.set_xticks(states); ax.grid(alpha=.2); ax.margins(x=.08, y=.15)
    fig.suptitle('BCH state-size tradeoff at K=2²⁰, T=64, relative distance 10%\nEach marker replays the complete occupation sum; only useful full bounds are shown.', fontsize=13)
    fig.savefig(HERE/'bch_full_bound_state_scaling_v7.png', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout='constrained')
    for column, block in enumerate((64, 128)):
        top, bottom = axes[:, column]
        for state, color in ((16, '#ce782e'), (20, '#176c52')):
            selected = [next((r for r in rows if r['block_bits'] == block and r['step_bits'] == 64
                             and r['state_bits'] == state and r['message_exponent'] == e), None)
                        for e in (20, 22, 24)]
            margins = [r['full_margin_bits'] if r is not None and r['full_margin_bits'] is not None
                       and r['full_margin_bits'] > 0 else np.nan for r in selected]
            losses = [r['full_aggregation_penalty_bits'] if math.isfinite(margin) else np.nan
                      for r, margin in zip(selected, margins)]
            top.plot([20, 22, 24], margins, 'o-', color=color, label=f'S={state}')
            bottom.plot([20, 22, 24], losses, 'o-', color=color, label=f'S={state}')
        top.set(title=f'BCH-{block}', ylabel='Full margin (bits)')
        bottom.set(yscale='log', ylabel='Full-bound loss relative to Q1 (bits)', xlabel='log₂ K')
        top.legend()
        for ax in (top, bottom):
            ax.set_xticks([20, 22, 24]); ax.grid(alpha=.2)
    fig.suptitle('BCH message/state interaction: T=64, relative distance 10%\nEvery marker has a replayed complete bound; weak or missing bounds break the curves.', fontsize=13)
    fig.savefig(HERE/'bch_full_bound_k_state_scaling_v7.png', dpi=180)
    plt.close(fig)
    inputs = [Path(__file__), lower_path, verification_path, HERE/'bch_dominance_sparse.csv',
              HERE/'bch_dominance_sparse.json', *full_paths, *extra_inputs]
    summary = dict(rows=len(rows), evidence_counts={label: sum(r['evidence'] == label for r in rows) for label in sorted({r['evidence'] for r in rows})},
                   csv_sha256=study.sha(path), source_sha256={str(p): study.sha(p) for p in inputs},
                   limitations=['Full-bound conclusions apply only to geometries with full-reference evidence.',
                                'Lower bounds concern expected bad-word counts, not failure probabilities.'])
    (HERE/'bch_engineering_evidence_v7.json').write_text(json.dumps(summary, indent=2)+'\n')
    model = dict(status='ENGINEERING_MODEL_WITH_FIXED_CALIBRATION_CHECKS', checks=model_checks,
                 source_sha256={str(p): study.sha(p) for p in [Path(__file__), *full_paths]},
                 limitations=['The model is an estimate, not a full-tail upper bound.',
                              'Calibration uses K=2^20 only; comparison points do not refit the coefficients.',
                              'No uniform claim between anchors or extrapolation beyond log2 K=26.',
                              'Coefficients describe selected bounds, not true failure probabilities.'])
    (HERE/'bch_k_counting_model_v3.json').write_text(json.dumps(model, indent=2)+'\n')
    print(json.dumps(summary['evidence_counts'], indent=2))


if __name__ == '__main__':
    main()
