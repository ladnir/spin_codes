"""Directed entrywise-marker fallback for one gamma profile cell.

Every marker gives a uniform upper bound on each regional matrix entry.
Taking their entrywise minimum preserves that upper bound. The optional
balanced mode distributes the density slope among256 regional factors;
the control mode maximizes the outer prefix separately. Both retain mixed
markers at a point, and may lose correlations on broad cells.
No frozen provider or cell checker is modified.
"""

import ast
from hashlib import sha256
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


def _signed_upper(value):
    if not value.is_finite():
        raise ArithmeticError('Nonfinite cell normalization or outer prefix')
    return value.upper()


def _entrywise_bound(bank, cell, support_count, budget, *, balanced_density=False):
    """Numerical core for synthetic tests; public entry verifies the bank/count.

    The caller supplies an upper bound on the represented support count.
    A deadline or arithmetic failure raises without returning a cell claim.
    """
    budget.check()
    if (not isinstance(cell, Cell) or cell.pairs != base.PAIRS
            or type(support_count) is not int or support_count < 0
            or type(balanced_density) is not bool):
        raise ValueError('Expected a production native cell and nonnegative integer support count')
    common = dict(cell=cell.descriptor(), support_count_upper=support_count,
                  complete_cell=True, whole_code_certificate=False,
                  bound_method='entrywise-regional-all-marker-minimum',
                  chronological_order='even-odd-to-power128',
                  balanced_density=balanced_density,
                  regional_matrix='2^(slope*r/256) R' if balanced_density else 'R',
                  outer_prefix_maximized_separately=not balanced_density,
                  common_state_potential=False)
    if cell.integer_empty():
        budget.check()
        return dict(common, zero=True, integer_empty=True, marker_entries_evaluated=0,
                    log2_upper_rational=None, log2_upper_enclosure=None)
    if (support_count == 0 or not bank.markers or len(bank.log_markers) != len(bank.markers)
            or any(len(row) != 2 for row in bank.log_markers)):
        raise ValueError('A nonempty cell needs a positive count and a complete marker bank')
    vertices = cell.vertices()
    rectangles = tuple(cell.histogram_bounds(parity) for parity in (0, 1))
    if (not vertices or any(rect is None or any(lo < 0 or lo > hi for lo, hi in rect)
                            for rect in rectangles)):
        raise ArithmeticError('Incomplete vertices or invalid histogram envelope')
    logtwo, logfactorial = arb(2).log(), arb(base.PAIRS+1).lgamma()
    histograms, prefix_values = ([], []), []
    # Marker-independent Gamma terms are evaluated only once per vertex/parity.
    for vertex in vertices:
        budget.check()
        weight, even, odd = cell.profile(vertex)
        density = bank.slope*base.rational(weight)
        prefix_values.append(_signed_upper(bank.constant if balanced_density else bank.constant+density))
        for parity, full in enumerate((even, odd)):
            single, empty = weight-2*full, base.PAIRS-weight+full
            minus_logh = (sum(((base.rational(value)+1).lgamma()
                              for value in (empty, single, full)), arb(0))-logfactorial)/logtwo
            # There are256 regional matrices in (S_even S_odd)^128.
            # The added affine term preserves convexity on the entire cell.
            scale_log = minus_logh+density/256 if balanced_density else minus_logh
            histograms[parity].append((base.rational(single), base.rational(full), scale_log))
    prefix = max(prefix_values)
    regions, choices, normalization = [], [], []
    for parity in (0, 1):
        best, selected, normalizers = [None]*4, [None]*4, [None]*4
        for index, (logx, logy) in enumerate(bank.log_markers):
            budget.check()
            numerator = tuple(base.upper(value) for value in bank.numerator(index, rectangles[parity]))
            if len(numerator) != 4:
                raise ArithmeticError('A regional numerator must have four entries')
            exponent = max(_signed_upper(minus_logh-single*logx-full*logy)
                           for single, full, minus_logh in histograms[parity])
            scale = base.upper(arb(2)**exponent)
            for entry in range(4):
                cap = base.upper(numerator[entry]*scale)
                if best[entry] is None or cap < best[entry]:
                    best[entry], selected[entry], normalizers[entry] = cap, index, exponent
        regions.append(tuple(best))
        choices.append(selected)
        normalization.append([str(base.rational_ceiling(value)) for value in normalizers])
    budget.check()
    transfer = profiles.pair_row_sum128(*regions)
    if transfer.is_zero():
        result = dict(zero=True, log2_upper_rational=None, log2_upper_enclosure=None)
    else:
        value = tail.log2arb(arb(support_count))+prefix+tail.log2arb(transfer)
        result = dict(zero=False, log2_upper_rational=str(base.rational_ceiling(value)),
                      log2_upper_enclosure=str(value))
    budget.check()
    return dict(common, **result, integer_empty=False,
                marker_entries_evaluated=8*len(bank.markers),
                marker_indices_by_parity_entry=choices,
                regional_log2_upper_rational=[[
                    None if value.is_zero() else str(base.rational_ceiling(tail.log2arb(value)))
                    for value in row] for row in regions],
                regional_scale_log2_upper_rational=normalization,
                prefix_log2_upper_rational=str(base.rational_ceiling(prefix)),
                transfer=profiles.support_log(transfer),
                histogram_rectangles=rectangles,
                vertices=[[str(value) for value in vertex] for vertex in vertices])


def entrywise_bound(bank, cell, support_count, budget, *, balanced_density=False):
    """Standalone directed cell receipt, checking both provider and count.

    The supplied count must dominate the frozen exact projection count.
    This routine evaluates every marker, not a heuristic shortlist. With
    balanced_density=True it bounds S(r,f)=2^(slope*r/256) R(r,f) entrywise
    and leaves only bank.constant outside the chronological transfer.
    """
    if (not isinstance(bank, provider.AliasBank) or not isinstance(cell, Cell)
            or cell.pairs != base.PAIRS or type(support_count) is not int or support_count < 0
            or type(balanced_density) is not bool):
        raise ValueError('Expected the frozen gamma all-pair bank, production cell, and integer count')
    budget.check()
    bank.check_pins()
    pins = source_pins()
    checked_count = counts.support_count(cell, budget)
    if support_count < checked_count:
        raise ValueError('Supplied count is below the verified projection support count')
    with ctx.workprec(bank.precision):
        result = _entrywise_bound(bank, cell, support_count, budget, balanced_density=balanced_density)
    budget.check()
    profiles.check_source_pins(pins)
    bank.check_pins()
    receipt = dict(bank.metadata())
    receipt.update(result, source_sha256=pins, cell_source_sha256=pins,
                   bank_source_sha256=bank.pins, checked_projection_support_count=checked_count,
                   support_count_verified=True, directed=True,
                   scope='this complete integer cell only; not a coverage or whole-code certificate')
    return receipt
