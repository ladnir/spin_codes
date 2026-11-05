"""Directed bounds on native profile cells, not a whole-code certificate.

Use one fixed positive marker per parity, an exact support-count envelope,
and every rational vertex of the cell. The inherited scalar is called with
count_matching=None: the support count is supplied separately, exactly once.
The coefficient provider and its frozen receipts are not modified.
"""

import ast
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path

from flint import arb, ctx

import hd_field128_gamma_all_pair_certificate as provider
import hd_field128_gamma_profile_certificate as profiles
import hd_field128_gamma_cell_counts as counts
from hd_field128_gamma_cell_geometry import Cell

base, tail = provider.base, provider.tail


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


def signed_upper(value):
    if not value.is_finite():
        raise ArithmeticError('Nonfinite scalar or logarithmic enclosure')
    return value.upper()


def representative(cell):
    """Feasible integer point near the box center, for witness selection only."""
    if cell.integer_empty():
        return None
    p = cell.pairs
    if cell.kind == 'interior':
        (rl, rh), (al, ah), (bl, bh) = cell.ranges
        al, ah, bl, bh = max(0, al), min(p, ah), max(0, bl), min(p, bh)
        rl, rh = max(1, rl, 2*al+1, 2*bl+1), min(2*p-1, rh, p+ah, p+bh)
        r = (rl+rh)//2
        a = (max(al, r-p)+min(ah, (r-1)//2))//2
        b = (max(bl, r-p)+min(bh, (r-1)//2))//2
        point = r, a, b
    else:
        (ul, uh), (kl, kh) = cell.ranges
        kl = max(1, kl)
        u = (max(1, ul, kl)+min(p-1, uh, p-kl))//2
        point = u, (kl+min(kh, u, p-u))//2
    if not cell.contains(point):
        raise ArithmeticError('Invalid integer witness for marker selection')
    return point


def fixed_bound(bank, cell, marker_indices, support_count, budget):
    """Arithmetic core; the public Evaluator checks the actual provider scope.

    No numerical claim is returned unless both numerator rectangles, the
    complete vertex set, and the chronological transfer finish successfully.
    """
    budget.check()
    if (cell.pairs != base.PAIRS or cell.integer_empty()
            or type(support_count) is not int or support_count <= 0
            or len(marker_indices) != 2
            or any(type(i) is not int or not 0 <= i < len(bank.markers) for i in marker_indices)):
        raise ValueError('Expected a nonempty 4096-pair cell, positive count, and two fixed markers')
    vertices = cell.vertices()
    rectangles = tuple(cell.histogram_bounds(parity) for parity in (0, 1))
    if (not vertices or any(rect is None or any(lo < 0 or lo > hi for lo, hi in rect)
                            for rect in rectangles)):
        raise ArithmeticError('Incomplete cell geometry or invalid histogram envelope')
    numerators = tuple(tuple(base.upper(value) for value in bank.numerator(index, rect))
                       for index, rect in zip(marker_indices, rectangles))
    values = []
    for vertex in vertices:
        budget.check()
        values.append(signed_upper(bank.scalar(cell.profile(vertex), marker_indices, count_matching=None)))
    envelope = max(values)
    transfer = profiles.pair_row_sum128(*numerators)
    if transfer.is_zero():
        # A zero positive upper bound is a valid zero-event certificate.
        result = dict(zero=True, log2_upper_rational=None, log2_upper_enclosure=None)
    else:
        value = tail.log2arb(arb(support_count))+tail.log2arb(transfer)+envelope
        result = dict(zero=False, log2_upper_rational=str(base.rational_ceiling(value)),
                      log2_upper_enclosure=str(value))
    budget.check()
    return dict(result, marker_indices=list(marker_indices), histogram_rectangles=rectangles,
                vertices=[[str(x) for x in point] for point in vertices],
                scalar_log2_upper_rational=str(base.rational_ceiling(envelope)),
                transfer=profiles.support_log(transfer), fixed_marker_per_parity=True,
                count_matching=None, chronological_order='even-odd-to-power128')


class Evaluator:
    """Batch evaluator with a bounded count cache and immutable source checks.

    Call check_pins before publishing a batch. A cell receipt alone makes no
    coverage or union claim about any other cell, weight, or construction.
    """

    def __init__(self, bank):
        if not isinstance(bank, provider.AliasBank):
            raise ValueError('Expected a freshly built directed gamma all-pair AliasBank')
        bank.check_pins()
        self.bank, self.pins = bank, source_pins()
        self.cache = counts.PrefixCache()

    def check_pins(self):
        profiles.check_source_pins(self.pins)
        self.bank.check_pins()

    def candidate_indices(self, cell, budget):
        """Heuristic shortlist only; every chosen pair receives a full cell bound.

        Keep each entry's best marked numerator at one feasible integer point.
        Its matching normalizer is marker-independent and may be omitted here.
        This selection is not used as a bound on the rest of the cell.
        """
        point = representative(cell)
        if point is None:
            return ((), ())
        r, a, b = cell.profile(point)
        choices = []
        for full in (a, b):
            single = r-2*full
            best = [None]*4
            for index, logarithms in enumerate(self.bank.log_markers):
                budget.check()
                row = self.bank.numerator(index, ((single, single), (full, full)))
                denominator_log = single*logarithms[0]+full*logarithms[1]
                for entry, cap in enumerate(row):
                    if cap.is_zero():
                        score = None
                    else:
                        score = signed_upper(tail.log2arb(cap)-denominator_log)
                    old = best[entry]
                    if old is None or (score is None and old[0] is not None) or (
                            score is not None and old[0] is not None and score < old[0]):
                        best[entry] = score, index
            choices.append(tuple(sorted({value[1] for value in best})))
        return tuple(choices)

    def receipt(self, cell, *, marker_pairs=None, max_seconds=60., failure_bits=64):
        if not isinstance(cell, Cell) or cell.pairs != base.PAIRS:
            raise ValueError('Expected the production 4096-pair native cell domain')
        if type(failure_bits) is not int or failure_bits < 0:
            raise ValueError('Expected a nonnegative integer probability budget')
        budget, bank = base.Budget(max_seconds), self.bank
        identity = dict(status='pending', cell=cell.descriptor(), construction=bank.construction,
                        state=bank.state, q=str(bank.q), target_distance=str(bank.distance),
                        forbidden_output_weight=bank.forbidden, precision=bank.precision,
                        witness=dict(z=str(bank.z), alpha=str(bank.alpha), grid=bank.grid, bits=bank.bits),
                        directed=False, complete_cell=False, whole_code_certificate=False,
                        log2_upper_rational=None, failure_bits=failure_bits,
                        cell_source_sha256=self.pins, bank_source_sha256=bank.pins)
        try:
            with ctx.workprec(bank.precision):
                count = counts.count_diagnostics(cell, budget, cache=self.cache)
                identity['count'] = count
                if not count['support_count_upper']:
                    best = dict(zero=True, log2_upper_rational=None, log2_upper_enclosure=None,
                                integer_empty=cell.integer_empty(), marker_indices=None)
                    tried = 0
                else:
                    pairs = (tuple(product(*self.candidate_indices(cell, budget)))
                             if marker_pairs is None else tuple(dict.fromkeys(tuple(x) for x in marker_pairs)))
                    if not pairs:
                        raise ValueError('At least one fixed marker pair is required')
                    best, tried = None, 0
                    for indices in pairs:
                        value = fixed_bound(bank, cell, indices, count['support_count_upper'], budget)
                        tried += 1
                        if best is None or value['zero'] or (not best['zero'] and
                                F(value['log2_upper_rational']) < F(best['log2_upper_rational'])):
                            best = value
                        if best['zero']:
                            break
                budget.check()
            identity.update(best, status='complete', complete_cell=True, directed=True,
                            marker_pairs_evaluated=tried,
                            meets_cell_budget=best['zero'] or F(best['log2_upper_rational']) <= -failure_bits)
            if best['marker_indices'] is not None:
                identity['witness']['selected_markers'] = [[str(v) for v in bank.markers[i]]
                                                            for i in best['marker_indices']]
        except (TimeoutError, ArithmeticError) as error:
            identity.update(failure=dict(error_type=type(error).__name__, message=str(error)),
                            meets_cell_budget=False)
        self.check_pins()
        return identity
