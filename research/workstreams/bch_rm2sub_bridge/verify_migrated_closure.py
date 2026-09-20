"""Replay the frozen closure audit without replacing its original receipt."""
import argparse
import json
from pathlib import Path

import audit_larger_state_closure as original


def verify(output=None):
    base = original.base
    receipt = base.HERE / 'generated/larger_t64_s20_full_closure_audit.json'
    expected = base.read(receipt)
    completed = []
    writer = base.write_new
    directory = base.HERE

    class ReplayPath(type(directory)):
        def exists(self):
            # The producer requires an unused output filename. Virtualize only
            # that precondition; reads still use the untouched original receipt,
            # and compare() replaces the sole output operation below.
            return False if self == receipt else super().exists()

    def compare(path, actual):
        assert Path(path) == receipt
        # Subprocess test timings vary; their success is checked by the original
        # audit with check=True. All remaining receipt content must match exactly.
        assert [r['file'] for r in actual['tests']] == [r['file'] for r in expected['tests']]
        assert {k: v for k, v in actual.items() if k != 'tests'} == {
            k: v for k, v in expected.items() if k != 'tests'}
        completed.append(actual)

    base.write_new = compare
    base.HERE = ReplayPath(directory)
    try:
        original.run()
    finally:
        base.write_new = writer
        base.HERE = directory
    assert len(completed) == 1
    result = dict(status='MIGRATED_FROZEN_CLOSURE_AUDIT_REPLAY_PASSED',
                  original_receipt_sha256=base.sha(receipt),
                  migration_manifest_sha256=base.sha(base.HERE / 'MIGRATION_MANIFEST.json'),
                  verifier_sha256=base.sha(Path(__file__)),
                  frozen_checker_sha256=base.sha(Path(original.__file__)),
                  coverage=[1, 8192], exact_shell_caps_replayed=expected['exact_shell_caps_replayed'],
                  margin_bits_diagnostic=expected['margin_bits_diagnostic'],
                  failure_below_2_to_minus_50=expected['full_failure_below_2_to_minus_50'],
                  scope='Frozen final audit replay; not a fresh 512-bit recomputation of every inner coefficient.',
                  tests_passed=[r['file'] for r in completed[0]['tests']])
    if output:
        writer(output, result)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    verify(parser.parse_args().output)
