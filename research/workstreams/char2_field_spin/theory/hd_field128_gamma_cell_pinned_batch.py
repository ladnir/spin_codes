"""Amortized evaluator source checks around the frozen selected-weight cover.

Standalone cell receipts check sources normally. During run_cover, evaluator
checks are deferred until the guarded batch exits; the frozen driver's own
checks still run. Checkpoints are partial, and nested cell receipts therein
remain conditional on the final batch source check. No final result returns
unless the complete evaluator closure and bank closure still match.
"""

import ast
from hashlib import sha256
from pathlib import Path

import hd_field128_gamma_cell_cover as frozen_cover
import hd_field128_gamma_cell_regional_table as regional


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


class PinnedBatchEvaluator(regional.Evaluator):
    """One normalized regional table with guarded batch pin amortization.

    Only run_cover opens the deferred-check scope. Its returned final record
    validates the batch as a whole, including receipts whose local metadata
    records that their check was deferred. Published checkpoints never acquire
    a completed batch claim retroactively.
    """

    def __init__(self, bank, weight, *, max_seconds=120.):
        self._batch_active = False
        super().__init__(bank, weight, max_seconds=max_seconds)
        self.pins = source_pins()
        self._force_check_pins()

    def _force_check_pins(self):
        # Explicit base dispatch bypasses only this subclass's deferral guard.
        # It validates self.pins, which includes this wrapper and the driver.
        regional.Evaluator.check_pins(self)

    def check_pins(self):
        if not self._batch_active:
            self._force_check_pins()

    def receipt(self, cell, *, marker_pairs=None, max_seconds=60., failure_bits=64):
        result = super().receipt(cell, marker_pairs=marker_pairs, max_seconds=max_seconds,
                                 failure_bits=failure_bits)
        result['source_verification'] = ('pending-batch-exit' if self._batch_active
                                         else 'verified-standalone')
        return result

    def run_cover(self, weight=None, **options):
        """Return only after before/after checks; exceptions also close/check.

        All ordinary cover options retain the frozen driver's interpretation.
        The callback receives partial checkpoint envelopes, never the final
        result. A source failure at exit prevents any completed return.
        """
        if self._batch_active:
            raise RuntimeError('A pinned cover batch is already active; nested or reentrant batches are forbidden')
        weight = self.weight if weight is None else weight
        if (type(weight) is not int or weight != self.weight
                or options.get('pairs', regional.base.PAIRS) != regional.base.PAIRS):
            raise ValueError('The batch must use the locked production selected weight')
        callback = options.pop('checkpoint', None)
        if callback is not None and not callable(callback):
            raise ValueError('Checkpoint must be callable or None')
        self._force_check_pins()

        def checkpoint(row):
            if (row.get('directed') is not False or row.get('complete_selected_weight') is not False
                    or row.get('whole_code_certificate') is not False
                    or row.get('log2_upper_rational') is not None):
                raise RuntimeError('A batch checkpoint must remain an incomplete, non-directed envelope')
            if callback is not None:
                callback(dict(row, batch_source_sha256=dict(self.pins),
                              batch_source_verified_at_exit=False,
                              cell_receipts_conditional_on_final_source_check=True))

        self._batch_active = True
        try:
            result = frozen_cover.cover(self, weight, checkpoint=checkpoint, **options)
        finally:
            self._batch_active = False
            self._force_check_pins()
        # Do not mutate nested checkpoint/node objects previously published by
        # the frozen driver. This final envelope validates the whole batch.
        return dict(result, batch_source_sha256=dict(self.pins),
                    batch_source_verified_at_entry=True, batch_source_verified_at_exit=True,
                    batch_pin_policy='full evaluator and bank closure before and after frozen cover',
                    checkpoint_source_claim='partial checkpoints require this successful final source check',
                    old_numerical_bounds_reused=False)
