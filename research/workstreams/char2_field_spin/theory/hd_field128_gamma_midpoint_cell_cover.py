"""Refine the frozen regional-table cell cover with midpoint box counts.

The regional transfer is unchanged. Each cell is freshly bounded by the
frozen evaluator, then rescaled by the exact ratio of two valid integer
support-count caps. The selected weight remains fixed throughout.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

from flint import arb, ctx

import hd_field128_gamma_cell_pinned_batch as frozen
import hd_field128_gamma_midpoint_boxes as boxes


base, tail = frozen.regional.base, frozen.regional.tail


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
    return dict(sorted({**pins, **boxes.source_pins()}.items()))


def refine_receipt(receipt, midpoint_count, precision):
    """Adjust only a completed fresh count-times-transfer cell bound."""
    if (receipt.get('complete_cell') is not True or receipt.get('directed') is not True
            or receipt.get('status') != 'complete' or type(midpoint_count) is not int
            or midpoint_count < 0):
        raise ValueError('Expected a complete directed cell receipt and nonnegative integer cap')
    old_count = receipt['count']['support_count_upper']
    if type(old_count) is not int or old_count <= 0:
        raise ValueError('Expected a positive unrefined cell support count')
    new_count = min(old_count, midpoint_count)
    result = dict(receipt, unrefined_count=receipt['count'],
                  unrefined_log2_upper_rational=receipt['log2_upper_rational'],
                  midpoint_rectangle_count_upper=midpoint_count,
                  count=dict(support_count_upper=new_count,
                             exact_support_count=False,
                             method='minimum of frozen projected count and raw midpoint rectangle sum'),
                  bound_method='frozen normalized regional transfer with midpoint rectangle count refinement',
                  exact_count_rescaling=dict(numerator=new_count, denominator=old_count))
    if not new_count:
        result.update(zero=True, log2_upper_rational=None, log2_upper_enclosure=None,
                      meets_cell_budget=True)
    elif not receipt['zero'] and new_count < old_count:
        with ctx.workprec(precision):
            value = (base.rational(F(receipt['log2_upper_rational']))
                     + tail.log2arb(arb(new_count))-tail.log2arb(arb(old_count)))
            ceiling = base.rational_ceiling(value)
        result.update(log2_upper_rational=str(ceiling), log2_upper_enclosure=str(value),
                      meets_cell_budget=ceiling <= -receipt['failure_bits'])
    return result


class Evaluator(frozen.PinnedBatchEvaluator):
    """Fresh directed regional tables; O(r) count preparation per weight."""

    def __init__(self, bank, weight, *, max_seconds=120.):
        expected_pins = source_pins()
        super().__init__(bank, weight, max_seconds=max_seconds)
        budget = base.Budget(max_seconds)
        self.boxes = boxes.MidpointBoxes(base.PAIRS, weight, budget=budget)
        self.pins = expected_pins
        self._force_check_pins()

    def receipt(self, cell, *, marker_pairs=None, max_seconds=60., failure_bits=64):
        result = super().receipt(cell, marker_pairs=marker_pairs, max_seconds=max_seconds,
                                 failure_bits=failure_bits)
        if result.get('complete_cell') is not True or result.get('integer_empty'):
            return result
        even, odd = result['full_count_intervals']
        count = self.boxes.rectangle(*even, *odd)
        result = refine_receipt(result, count, self.bank.precision)
        self.check_pins()
        return result
