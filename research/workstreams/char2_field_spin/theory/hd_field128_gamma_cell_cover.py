"""Bounded adaptive coverage of one selected weight using directed cell receipts.

The public CLI builds a fresh all-pair bank from exact witness literals only.
The core accepts an evaluator interface for tiny synthetic tests. Every split
is an exact integer partition. A complete union requires accepted leaves over
all canonical roots; partial results expose their whole pending frontier and
have no selected-weight bound. This module never claims a whole-code result.
"""

import argparse
import ast
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from time import monotonic

from flint import arb, ctx

import hd_field128_gamma_cell_certificate as numeric
from hd_field128_gamma_cell_geometry import Cell


base, tail, provider = numeric.base, numeric.tail, numeric.provider


def source_pins():
    directory, pending, pins = Path(__file__).resolve().parent, [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in pins:
            continue
        raw = (directory/name).read_bytes()
        pins[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw)):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in pins:
                    pending.append(target)
    return dict(sorted(pins.items()))


def check_source_pins(pins):
    directory = Path(__file__).resolve().parent
    if any(sha256((directory/name).read_bytes()).hexdigest() != digest for name, digest in pins.items()):
        raise RuntimeError('A pinned adaptive-cover source changed')


def write_new(path, value):
    raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')
    with Path(path).open('xb') as stream:
        stream.write(raw)
    return dict(path=str(Path(path).resolve()), sha256=sha256(raw).hexdigest())


def selected_weight_roots(weight, pairs=4096):
    """Disjoint integer roots: interior and, at even weight, both pure faces."""
    if (type(pairs) is not int or pairs < 1 or type(weight) is not int
            or not 1 <= weight < 2*pairs):
        raise ValueError('Expected a nonempty, non-full selected weight')
    roots = [Cell('interior', ((weight, weight), (0, weight//2), (0, weight//2)), pairs)]
    if not weight % 2:
        u = weight//2
        ranges = ((u, u), (1, min(u, pairs-u)))
        roots.extend(Cell(kind, ranges, pairs) for kind in ('even_face', 'odd_face'))
    return tuple(roots)


def choose_split(cell):
    """Largest exact histogram-impact span; no numerical bound is inferred."""
    vertices = cell.vertices()
    if not vertices:
        return None
    impacts = (2, 3, 3) if cell.kind == 'interior' else (2, 3)
    choices = []
    for axis, impact in enumerate(impacts):
        low, high = min(v[axis] for v in vertices), max(v[axis] for v in vertices)
        low, high = -((-low.numerator)//low.denominator), high.numerator//high.denominator
        if low < high:
            choices.append(((high-low)*impact, -axis, axis, (low+high)//2))
    if not choices:
        return None
    unused, unused_axis, axis, cut = max(choices)
    return axis, cut


def _receipt_pass(receipt, cell, failure_bits):
    if not isinstance(receipt, dict) or receipt.get('cell') != cell.descriptor():
        raise ValueError('Cell receipt identity does not match the requested cell')
    if receipt.get('complete_cell') is not True:
        return False
    if receipt.get('directed') is not True or receipt.get('status') != 'complete':
        raise ValueError('A complete cell receipt must be directed and complete')
    if receipt.get('zero') is True:
        if receipt.get('log2_upper_rational') is not None:
            raise ValueError('A zero cell receipt must not supply a finite logarithm')
        passing = True
    elif receipt.get('zero') is False:
        value = receipt.get('log2_upper_rational')
        if type(value) not in (str, int):
            raise ValueError('Expected an exact rational cell ceiling')
        passing = F(value) <= -failure_bits
    else:
        raise ValueError('A complete receipt must explicitly classify its zero bound')
    if receipt.get('failure_bits') != failure_bits or receipt.get('meets_cell_budget') is not passing:
        raise ValueError('Cell budget metadata disagrees with its exact ceiling')
    return passing


def audit_tree(nodes, weight, pairs, failure_bits):
    """Check root coverage and every stored lattice split; return leaf id lists."""
    by_id = {node['id']: node for node in nodes}
    if len(by_id) != len(nodes):
        raise ValueError('Repeated tree node identity')
    roots = selected_weight_roots(weight, pairs)
    expected_ids = tuple(f'n{i:06d}' for i in range(len(roots)))
    seen, accepted, pending = set(), [], []

    def visit(node_id, expected, parent):
        if node_id in seen or node_id not in by_id:
            raise ValueError('Cyclic, shared, or absent child in the cover tree')
        seen.add(node_id)
        node = by_id[node_id]
        if node['cell'] != expected.descriptor() or node['parent'] != parent:
            raise ValueError('Tree geometry or parent link does not match its checked split')
        status = node['status']
        if status == 'split':
            split, children = node['split'], node['children']
            if len(children) != 2:
                raise ValueError('A split must retain both integer children')
            left, right = expected.split(split['axis'], split['integer_cut'])
            visit(children[0], left, node_id)
            visit(children[1], right, node_id)
        elif status == 'accepted':
            if not _receipt_pass(node['receipt'], expected, failure_bits):
                raise ValueError('An accepted leaf lacks a passing complete receipt')
            accepted.append(node_id)
        elif status in ('pending', 'unresolved'):
            pending.append(node_id)
        else:
            raise ValueError('Unknown tree node state')

    for node_id, root in zip(expected_ids, roots):
        visit(node_id, root, None)
    if seen != set(by_id):
        raise ValueError('Unreachable extra node in the cover tree')
    return accepted, pending


def union_ceiling(receipts, precision):
    """Directed sum of accepted leaf ceilings, called only after full coverage."""
    with ctx.workprec(precision):
        total = sum((arb(2)**base.rational(F(row['log2_upper_rational']))
                     for row in receipts if not row['zero']), arb(0))
        if total.is_zero():
            return dict(zero=True, log2_upper_rational=None, log2_upper_enclosure=None)
        value = tail.log2arb(total)
        return dict(zero=False, log2_upper_rational=str(base.rational_ceiling(value)),
                    log2_upper_enclosure=str(value))


def checked_limits(max_nodes, max_seconds, cell_seconds, failure_bits, aggregate_bits, precision):
    if (type(max_nodes) is not int or not 1 <= max_nodes <= 100000
            or any(type(value) is not int or not 0 <= value <= 4096 for value in (failure_bits, aggregate_bits))
            or type(precision) is not int or not 128 <= precision <= 4096
            or any(type(value) not in (int, float) or not isfinite(value) or value <= 0
                   for value in (max_seconds, cell_seconds))):
        raise ValueError('Invalid explicit adaptive-cover resource or probability budget')


def cover(evaluator, weight, *, pairs=4096, max_nodes=256, max_seconds=180.,
          cell_seconds=60., failure_bits=64, aggregate_bits=40, precision=256,
          marker_pairs=None, checkpoint=None, input_check=None):
    """Evaluate at most max_nodes cells; checkpoints are append-only tree deltas.

    A delta contains replaced/new node records plus the complete pending id
    list. Replaying the initial root records and later replacements reconstructs
    the tree. Final output includes that entire tree. No resume is attempted.
    """
    roots = selected_weight_roots(weight, pairs)
    checked_limits(max_nodes, max_seconds, cell_seconds, failure_bits, aggregate_bits, precision)
    if marker_pairs is not None:
        marker_pairs = tuple(tuple(row) for row in marker_pairs)
    started, pins = monotonic(), source_pins()
    production = isinstance(evaluator, numeric.Evaluator)
    evaluator.check_pins()
    if input_check:
        input_check()
    nodes = [dict(id=f'n{i:06d}', parent=None, cell=cell.descriptor(), status='pending')
             for i, cell in enumerate(roots)]
    cells = list(roots)
    queue, blocked, accepted = deque(range(len(nodes))), [], []
    evaluations, sequence = 0, 0

    def publish(changed, action):
        nonlocal sequence
        check_source_pins(pins)
        if input_check:
            input_check()
        value = dict(sequence=sequence, action=action, changed_nodes=changed,
                     pending_frontier=[nodes[i]['id'] for i in tuple(queue)+tuple(blocked)],
                     accepted_leaves=[nodes[i]['id'] for i in accepted], evaluated_nodes=evaluations,
                     active_groups=weight, pairs=pairs, complete_selected_weight=False,
                     directed=False, whole_code_certificate=False, log2_upper_rational=None)
        if checkpoint:
            checkpoint(value)
        sequence += 1

    publish(list(nodes), 'initial-roots')
    stop_reason = None
    while queue:
        remaining = max_seconds-(monotonic()-started)
        if evaluations >= max_nodes or remaining <= 0:
            stop_reason = 'max_nodes' if evaluations >= max_nodes else 'max_seconds'
            break
        index = queue.popleft()
        cell, node = cells[index], nodes[index]
        try:
            receipt = evaluator.receipt(cell, marker_pairs=marker_pairs,
                                        max_seconds=min(cell_seconds, remaining), failure_bits=failure_bits)
        except (TimeoutError, ArithmeticError) as error:
            receipt = dict(cell=cell.descriptor(), status='pending', complete_cell=False,
                           directed=False, log2_upper_rational=None, meets_cell_budget=False,
                           failure=dict(error_type=type(error).__name__, message=str(error)))
        evaluations += 1
        passing = _receipt_pass(receipt, cell, failure_bits)
        if passing:
            nodes[index] = dict(node, status='accepted', receipt=receipt)
            accepted.append(index)
            changed, action = [nodes[index]], 'accepted-leaf'
        elif receipt.get('complete_cell') is not True:
            nodes[index] = dict(node, status='unresolved', receipt=receipt, pending_reason='incomplete-cell-receipt')
            blocked.append(index)
            changed, action = [nodes[index]], 'unresolved-leaf'
        else:
            split = choose_split(cell)
            if split is None:
                nodes[index] = dict(node, status='unresolved', receipt=receipt, pending_reason='unsplittable-over-budget')
                blocked.append(index)
                changed, action = [nodes[index]], 'unresolved-leaf'
            else:
                axis, cut = split
                children = cell.split(axis, cut)
                indices = (len(nodes), len(nodes)+1)
                records = [dict(id=f'n{i:06d}', parent=node['id'], cell=child.descriptor(), status='pending')
                           for i, child in zip(indices, children)]
                nodes[index] = dict(node, status='split', receipt=receipt,
                                    split=dict(axis=axis, integer_cut=cut), children=[row['id'] for row in records])
                nodes.extend(records)
                cells.extend(children)
                # Breadth first; higher a/b (or lower face k) goes first at equal depth.
                order = tuple(reversed(indices)) if cell.kind == 'interior' and axis else indices
                queue.extend(order)
                changed, action = [nodes[index]]+records, 'integer-split'
        publish(changed, action)
    evaluator.check_pins()
    check_source_pins(pins)
    if input_check:
        input_check()
    checked_accepted, pending = audit_tree(nodes, weight, pairs, failure_bits)
    complete = not pending
    result = dict(status='complete' if complete else 'partial', active_groups=weight, pairs=pairs,
                  root_cells=[cell.descriptor() for cell in roots], nodes=nodes,
                  accepted_leaves=checked_accepted, pending_frontier=pending,
                  evaluated_nodes=evaluations, tree_nodes=len(nodes), checkpoints=sequence,
                  stop_reason=None if complete else (stop_reason or 'unresolved-leaves'),
                  complete_selected_weight=complete, covers_all_supports_at_weight=complete,
                  complete_requested_roots=complete, whole_code_certificate=False,
                  directed=complete and production, numerical_backend='directed-cell-evaluator' if production else 'injected-test-evaluator',
                  max_nodes=max_nodes, max_seconds=max_seconds, cell_seconds=cell_seconds,
                  cell_failure_bits=failure_bits, aggregate_failure_bits=aggregate_bits, precision=precision,
                  source_sha256=pins, old_numerical_bounds_reused=False,
                  zero=None, log2_upper_rational=None, log2_upper_enclosure=None,
                  meets_requested_selected_weight_target=False, elapsed_seconds=monotonic()-started,
                  scope='one complete selected weight only if all roots pass; injected evaluators are test interfaces, never production certificates')
    if complete:
        by_id = {node['id']: node for node in nodes}
        result.update(union_ceiling([by_id[node_id]['receipt'] for node_id in checked_accepted], precision))
        result['meets_requested_selected_weight_target'] = result['zero'] or F(result['log2_upper_rational']) <= -aggregate_bits
    check_source_pins(pins)
    if input_check:
        input_check()
    return result


def witness_literals(artifact):
    """Import positive exact literals only, never a prior bound or marker log."""
    if not isinstance(artifact, dict) or not isinstance(artifact.get('witness'), dict):
        raise ValueError('Expected a saved directed summary with a witness object')
    witness = artifact['witness']
    def exact(value):
        if type(value) not in (str, int):
            raise ValueError('Witness literals must be rational strings or integers')
        return F(value)
    z, alpha = exact(witness['z']), exact(witness['alpha'])
    markers = tuple(tuple(exact(value) for value in row) for row in witness['markers'])
    if (not 0 < z < 1 or not 0 < alpha <= 1 or not 1 <= len(markers) <= 128
            or any(len(row) != 2 or min(row) <= 0 for row in markers)):
        raise ValueError('Invalid exact positive witness literals')
    return z, alpha, markers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--witness-source', type=Path, required=True)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--state', type=int, choices=provider.STATES, required=True)
    parser.add_argument('--distance', required=True)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--grid', type=int, default=128)
    parser.add_argument('--bits', type=int, default=80)
    parser.add_argument('--max-nodes', type=int, default=256)
    parser.add_argument('--max-seconds', type=float, default=180.)
    parser.add_argument('--cell-seconds', type=float, default=60.)
    parser.add_argument('--build-seconds', type=float, default=600.)
    parser.add_argument('--failure-bits', type=int, default=64)
    parser.add_argument('--aggregate-bits', type=int, default=40)
    parser.add_argument('--max-table-bytes', type=int, default=384 << 20)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('Output directory must be unused')
    checked_limits(args.max_nodes, args.max_seconds, args.cell_seconds,
                   args.failure_bits, args.aggregate_bits, args.precision)
    base.Budget(args.build_seconds)
    roots = selected_weight_roots(args.r)
    raw, pins = args.witness_source.read_bytes(), source_pins()
    z, alpha, markers = witness_literals(json.loads(raw))
    provenance = dict(path=str(args.witness_source.resolve()), sha256=sha256(raw).hexdigest(),
                      role='exact witness literals only; no old ceiling or coverage imported')
    def input_check():
        if args.witness_source.read_bytes() != raw:
            raise RuntimeError('Pinned exact witness source changed')
    args.output_dir.mkdir(exist_ok=False)
    setup = dict(witness_source=provenance, source_sha256=pins, active_groups=args.r,
                 construction=provider.local.construction(args.state), state_dimension=args.state,
                 target_distance=str(F(args.distance)), field_order=str(provider.FIELD_ORDER),
                 exact_witness=dict(z=str(z), alpha=str(alpha), markers=[[str(x) for x in row] for row in markers]),
                 root_cells=[cell.descriptor() for cell in roots], directed=False,
                 complete_selected_weight=False, whole_code_certificate=False,
                 log2_upper_rational=None, old_numerical_bounds_reused=False)
    write_new(args.output_dir/'setup.json', setup)
    try:
        bank = provider.AliasBank(z, alpha, markers, state=args.state, distance=F(args.distance),
                                  grid=args.grid, precision=args.precision, bits=args.bits,
                                  max_seconds=args.build_seconds, max_table_bytes=args.max_table_bytes,
                                  progress=lambda row: print(json.dumps(row), flush=True))
    except (TimeoutError, ArithmeticError) as error:
        input_check()
        check_source_pins(pins)
        write_new(args.output_dir/'summary.json', dict(setup, status='bank-build-incomplete',
                  pending_root_cells=[cell.descriptor() for cell in roots],
                  failure=dict(error_type=type(error).__name__, message=str(error))))
        raise SystemExit(2)
    check_source_pins(pins)
    input_check()
    evaluator = numeric.Evaluator(bank)
    def checkpoint(row):
        check_source_pins(pins)
        input_check()
        write_new(args.output_dir/f"checkpoint_{row['sequence']:06d}.json", row)
        print(json.dumps({key: row[key] for key in ('sequence', 'action', 'evaluated_nodes')}
                         | dict(accepted=len(row['accepted_leaves']), pending=len(row['pending_frontier']))), flush=True)
    result = cover(evaluator, args.r, max_nodes=args.max_nodes, max_seconds=args.max_seconds,
                   cell_seconds=args.cell_seconds, failure_bits=args.failure_bits,
                   aggregate_bits=args.aggregate_bits, precision=args.precision,
                   checkpoint=checkpoint, input_check=input_check)
    result.update(witness_source=provenance, bank_metadata=bank.metadata(),
                  bank_source_sha256=bank.pins, setup_source_sha256=pins)
    evaluator.check_pins()
    input_check()
    check_source_pins(pins)
    write_new(args.output_dir/'summary.json', result)
    print(json.dumps({key: result[key] for key in ('status', 'evaluated_nodes', 'complete_selected_weight',
                                                  'log2_upper_rational', 'stop_reason')}), flush=True)
    if not result['complete_selected_weight']:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
