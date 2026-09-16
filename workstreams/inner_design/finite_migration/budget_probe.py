"""Diagnose centers of the largest remaining cover leaves; never a full proof."""
import argparse
from pathlib import Path

from flint import arb, ctx
import budget_dense as budget


def run(seed_path, output, count):
    assert count > 0 and not output.exists()
    ctx.prec = 256
    saved = budget.model.base.read(seed_path)
    budget.model.authenticate(saved)
    checker = budget.short_dense.Checker(saved['exponent'])
    assert checker.engine.identity() == saved['instance']
    selected = sorted(saved['leaves'].items(), key=lambda item: item[1]['power'], reverse=True)[:count]
    points = []
    for key, node in selected:
        lo, hi, a, b = budget.geometry.geometry(node)
        q, coordinate = (lo + hi) // 2, (a + b) / 2
        old = checker.bound(q, q, coordinate, coordinate, node['witness'])
        witness = checker.witness(q, q, coordinate, coordinate)
        current = checker.bound(q, q, coordinate, coordinate, witness)
        if old < current:
            current, witness = old, node['witness']
        margin = float(-current / arb(2).log())
        points.append(dict(key=key, leaf_power=node['power'], q=q, coordinate=str(coordinate),
                           witness=witness, log_upper=budget.model.base.encode(budget.model.exact(current)),
                           margin_bits=margin))
        print('center', key, q, str(coordinate), 'leaf', node['power'], 'margin', margin, flush=True)
    sources = budget.short_dense.mixed_dense.ladder.candidate.sources()
    sources[seed_path.relative_to(budget.model.ROOT).as_posix()] = budget.model.base.sha(seed_path)
    budget.model.base.write_new(output, dict(status='IMT_COVER_CENTER_DIAGNOSTICS',
        instance=checker.engine.identity(), points=points, full_distance_proved=False,
        source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--count', type=int, default=5)
    a = p.parse_args()
    run(a.seed.resolve(), a.output.resolve(), a.count)
