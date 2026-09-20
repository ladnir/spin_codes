"""Audit the new Q1 observations and export parameter comparisons and a plot."""
import csv
import hashlib
import json
import sqlite3
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent.parent


def main():
    db = sqlite3.connect(f'{(HERE / "spin_landscape.sqlite3").as_uri()}?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    sources = db.execute('SELECT * FROM sources').fetchall()
    for source in sources:
        path = HERE.parent / source['path']
        if hashlib.sha256(path.read_bytes()).hexdigest() != source['sha256']:
            raise ValueError(f'changed source: {path}')
    rows = [dict(r) for r in db.execute(
        "SELECT * FROM landscape WHERE transfer_review_status='activation_aware' AND comparison_eligible=1 AND occupation_min=1 AND occupation_max=1 ORDER BY outer_id,message_exponent,step_bits,state_bits")]
    if not rows:
        raise ValueError('no activation-aware observations')
    for row in rows:
        if row['occupation_min'] != 1 or row['occupation_max'] != 1 or row['coverage_kind'] != 'single_occupation':
            raise ValueError('this report is specifically Q1')
        if row['outer_model_kind'] == 'fixed_certified_constraints':
            raise ValueError('bounded-spectrum BCH must remain outside the primary pilot')
        if row['outer_rows'] % row['step_bits'] or row['bad_weight'] != row['output_bits'] // 10:
            raise ValueError('noncomparable pilot geometry or cutoff')
    identities = [(r['outer_id'], r['message_exponent'], r['step_bits'], r['state_bits'], r['map_tag']) for r in rows]
    if len(set(identities)) != len(identities):
        raise ValueError('duplicate pilot observation')
    groups = {}
    for r in rows:
        groups.setdefault((r['outer_id'], r['message_exponent'], r['step_bits']), []).append(r)
    frontier = []
    for group in groups.values():
        for target in (0, 20, 40, 60):
            passing = [r for r in group if r['margin_bits'] >= target]
            best = min(passing, key=lambda r: r['state_bits']) if passing else None
            r = group[0]
            frontier.append(dict(outer_label=r['outer_label'], message_exponent=r['message_exponent'],
                                 step_bits=r['step_bits'], target_q1_margin=target,
                                 smallest_tested_passing_state=best['state_bits'] if best else '',
                                 attained_margin=best['margin_bits'] if best else '',
                                 witness_at_grid_edge=best['dominant_witness_at_grid_edge'] if best else '',
                                 status='Q1 screen only; no full-occupation guarantee'))
    with (HERE / 'activation_q1_frontier.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(frontier[0]))
        writer.writeheader()
        writer.writerows(frontier)
    maps = {}
    for r in rows:
        key = r['inner_id']
        if key in maps:
            continue
        # Locate the snapshotted generators through the source CSV, whose hash
        # was checked above. These counts use independent XOR reductions only.
        with (HERE.parent / r['source_path']).open(newline='', encoding='utf-8') as handle:
            source_row = next(x for x in csv.DictReader(handle) if x['map_tag'] == r['map_tag'])
        payload = json.loads((ROOT / source_row['map_source']).read_text())
        generators = [int(g, 16) for g in payload['generator_words_hex']]
        t, s = r['step_bits'], r['state_bits']
        a_cost = sum(max(0, sum((g >> j) & 1 for g in generators) - 1) for j in range(t))
        b_cost = sum(max(0, g.bit_count() - 1) for g in generators)
        maps[key] = dict(t=t, s=s, map_tag=r['map_tag'],
                         a_naive_xors=a_cost, b_naive_xors=b_cost,
                         a_b_and_output_xors_per_bit=(a_cost + b_cost + t) / t,
                         state_updates_per_output_bit=1 / t,
                         scope='Direct independent XOR reductions; excludes field multiply, update addition, outer, routing, SIMD and circuit sharing.')
    audit = dict(status='AUDITED_Q1_PILOT_INDEX', observation_count=len(rows),
                 dominant_grid_edges=sum(bool(r['dominant_witness_at_grid_edge']) for r in rows),
                 source_hashes_checked=len(sources),
                 message_exponents=sorted({r['message_exponent'] for r in rows}),
                 configurations=len(maps), exact_bch_observations=sum(r['outer_family'] == 'bch' for r in rows),
                 exact_rm_observations=sum(r['outer_family'] == 'rm' for r in rows),
                 random_reference_observations=sum(r['outer_family'] == 'random' for r in rows),
                 historical_rows_under_review=db.execute('SELECT COUNT(*) FROM results_under_review').fetchone()[0],
                 current_certificate_rows=db.execute('SELECT COUNT(*) FROM certified_results').fetchone()[0],
                 superseded_grid_rows=db.execute("SELECT COUNT(*) FROM landscape WHERE transfer_review_status='activation_aware' AND comparison_eligible=0").fetchone()[0],
                 cost_inputs=list(maps.values()),
                 limitations=['No BCH256 bounded-spectrum observations in primary pilot.',
                              'The completed BCH256 certificate is still in its owning worktree.',
                              'These Q1 frontiers are not full-certificate frontiers or runtime rankings.'])
    (HERE / 'activation_pilot_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=False)
    colors = {64: '#176b87', 128: '#b65325', 256: '#7451a8'}
    for ax, family, size, title in zip(axes, ['bch', 'rm'], [128, 512], ['Exact BCH [128,64,22]', 'Exact RM(4,9) [512,256,32]']):
        for e, style in [(16, '--'), (20, '-')]:
            for t in (64, 128, 256):
                selected = sorted([r for r in rows if r['outer_family'] == family and r['block_bits'] == size and
                                   r['message_exponent'] == e and r['step_bits'] == t], key=lambda r: r['state_bits'])
                if selected:
                    ax.plot([r['state_bits'] for r in selected], [r['margin_bits'] for r in selected],
                            color=colors[t], linestyle=style, marker='o', markersize=3,
                            label=f't={t}, log2(k)={e}')
        ax.axhline(40, color='#555555', linewidth=1, alpha=.6)
        ax.set_title(title)
        ax.set_xlabel('State bits s')
        ax.set_ylabel('Q1 diagnostic margin (bits)')
        ax.set_xticks([7, 9, 12, 14, 16, 18, 20])
        ax.grid(alpha=.18)
        ax.legend(fontsize=8)
    fig.suptitle('Activation-aware RM2Sub: epoch/state pilot', fontsize=14)
    fig.text(.5, .015, 'One unselected nested chain per epoch size. Binary64 Q1 bounds; higher occupations are not covered.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .96))
    fig.savefig(HERE / 'activation_ts_tradeoff.png', dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, family, size, title in zip(axes, ['bch', 'rm'], [128, 512], ['Exact BCH [128,64,22]', 'Exact RM(4,9) [512,256,32]']):
        for s, style in [(14, '--'), (20, '-')]:
            for t in (64, 128, 256):
                selected = sorted([r for r in rows if r['outer_family'] == family and r['block_bits'] == size and
                                   r['state_bits'] == s and r['step_bits'] == t], key=lambda r: r['message_exponent'])
                ax.plot([r['message_exponent'] for r in selected], [r['margin_bits'] for r in selected],
                        color=colors[t], linestyle=style, marker='o', markersize=3, label=f't={t}, s={s}')
        ax.axhline(40, color='#555555', linewidth=1, alpha=.6)
        ax.set_title(title)
        ax.set_xlabel('Message exponent log2(k)')
        ax.set_ylabel('Q1 diagnostic margin (bits)')
        ax.set_xticks(audit['message_exponents'])
        ax.grid(alpha=.18)
        ax.legend(fontsize=8)
    fig.suptitle('Fixed-map scaling: epoch size matters more when state is small', fontsize=13)
    fig.text(.5, .015, 'One unselected chain per epoch size. Binary64 Q1 bounds; no complete distance or runtime claim.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .96))
    fig.savefig(HERE / 'activation_k_scaling.png', dpi=170)
    plt.close(fig)
    db.close()
    print(json.dumps({k: v for k, v in audit.items() if k not in ('cost_inputs', 'limitations')}, indent=2))


if __name__ == '__main__':
    main()
