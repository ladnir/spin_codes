"""Bounded dense refinement: fractional ceilings and one density bisection.

Only selected leaves of an authenticated complete cover are replaced. Each
replacement bounds the same parent rectangle; untouched leaves retain their
accepted ceilings. This is an overlay, never a rewrite of the frozen cover.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path

from flint import ctx
import coupled_input_dense
import budget_dense
import q1_slack as audit

model = audit.model
HERE = Path(__file__).resolve().parent


def run(output, limit=8, verify=False):
    ctx.prec = 512 if verify else 256
    if not verify and output.exists():
        raise FileExistsError('Use a fresh output path')
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
    dense_path = HERE / 'DENSE_M16_coupled_v1.json'
    full_path = HERE / 'FULL_M16_VERIFIED.json'
    q1_path = HERE / 'Q1_INJECTION_SHELLS_v1.json'
    q1_replay_path = q1_path.with_name(q1_path.stem + '_replay.json')
    dense, full, q1, replay = map(model.base.read, (dense_path, full_path, q1_path, q1_replay_path))
    assert replay['producer_sha256'] == model.base.sha(q1_path)
    assert replay['status'] == 'IMT_K16_Q1_INJECTION_SHELL_512_BIT_LINEAR_REPLAY_PASSED'
    model.authenticate(q1)
    old_dense = budget_dense.checked_union(dense, 64, 512)
    parts = [model.base.decode(v) for v in full['component_upper']]
    assert old_dense == parts[2]
    assert sum(parts, F()) == model.base.decode(full['union_upper'])
    checker = coupled_input_dense.Checker(16)
    assert checker.engine.identity() == dense['instance'] == full['instance'] == q1['instance']
    keys = (list(saved['replacements']) if saved else
            sorted(dense['leaves'], key=lambda k: dense['leaves'][k]['power'], reverse=True)[:limit])
    replacements = {}
    total = old_dense
    for key in keys:
        node = dense['leaves'][key]
        lo, hi, a, b = budget_dense.geometry.geometry(node)
        mid = (a+b)/2
        def bound(left, right):
            log_value = checker.bound(lo, hi, left, right, node['witness'])
            return model.exact(model.unpack(model.pack(model.up(log_value.exp()))))
        parent = bound(a, b)
        children = bound(a, mid) + bound(mid, b)
        old = F(2)**node['power']
        value = min(old, parent, children)
        row = dict(upper=model.base.encode(value), parent_upper=model.base.encode(parent),
                   children_upper=model.base.encode(children), original_power=node['power'])
        if saved:
            accepted = saved['replacements'][key]
            assert accepted['original_power'] == node['power']
            assert parent <= model.base.decode(accepted['parent_upper'])
            assert children <= model.base.decode(accepted['children_upper'])
            value = min(old, model.base.decode(accepted['parent_upper']),
                        model.base.decode(accepted['children_upper']))
            assert value == model.base.decode(accepted['upper'])
            row = accepted
        replacements[key] = row
        total += value-old
        print('leaf', key, 'old_bits', -node['power'], 'parent_bits', audit.bits(parent),
              'split_bits', audit.bits(children), 'dense_bits', audit.bits(total), flush=True)
    q1_upper = model.base.decode(q1['q1_upper'])
    union = q1_upper + parts[1] + total
    diagnostic = dict(original_dense_bits=audit.bits(old_dense), dense_bits=audit.bits(total),
                      q1_bits=audit.bits(q1_upper), original_full_bits=audit.bits(sum(parts, F())),
                      full_bits=audit.bits(union), replaced_leaves=len(replacements),
                      original_leaf_count=len(dense['leaves']))
    result = dict(status='IMT_K16_DENSE_SLACK_OVERLAY', instance=checker.engine.identity(),
                  replacements=replacements, dense_upper=model.base.encode(total),
                  q1_upper=model.base.encode(q1_upper), sparse_upper=model.base.encode(parts[1]),
                  union_upper=model.base.encode(union), diagnostic=diagnostic)
    if saved:
        for key in ('instance', 'dense_upper', 'q1_upper', 'sparse_upper', 'union_upper', 'diagnostic'):
            assert result[key] == saved[key], key
        model.base.write_new(output.with_name(output.stem + '_replay.json'), dict(
            status='IMT_K16_DENSE_SLACK_512_BIT_REPLAY_PASSED',
            producer_sha256=model.base.sha(output), diagnostic=diagnostic))
    else:
        result['source_sha256'] = {**audit.ladder.candidate.sources(),
            **{p.relative_to(model.ROOT).as_posix(): model.base.sha(p)
               for p in (dense_path, full_path, q1_path, q1_replay_path)}}
        model.base.write_new(output, result)
    print(diagnostic, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--limit', type=int, default=8)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.limit < 1:
        parser.error('--limit must be positive')
    run(args.output.resolve(), args.limit, args.verify)
