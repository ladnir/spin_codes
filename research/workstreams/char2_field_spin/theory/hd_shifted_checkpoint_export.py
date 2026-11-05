"""Recompute established sparse or dense checkpoints with complete receipts.

This is an export of fresh arithmetic, not an importer of saved bounds.
It preserves the established witnesses and records every accepted leaf,
so later disjoint composition need not rely on summary-only console logs.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from hd_shifted_cover_ledger import Box
from hd_shifted_swap_ledger import SwapLedger, replay_swap_cover
from hd_shifted_swap_sparse import BATCH_SCHEDULE, BATCH_SCHEDULE_256, batch_receipt, family_receipt


PLANS = (
    ('unsigned', 'hd_shifted_bulk_plan4.json', 5120, 6143,
     'c8ce97119f9c1422a95f724bba1ca7000e9e4282caa5b09bf1d2ce52dbbec4a2'),
    ('unsigned', 'hd_shifted_bulk_plan2.json', 6144, 7167,
     '563aacd1f3d90fea2368edce7a6f5d3b2cee2dac8d03401450cf11352c9c997b'),
    ('unsigned', 'hd_shifted_bulk_plan.json', 7168, 7648,
     '384071145a74b64e99e0cd81d474b096e441d5937a4a24dd9db9ed3a9a338314'),
)


def sparse_checkpoint(precision=256, progress=None):
    ledger = SwapLedger(Box(((1, 256), (0, 4096), (0, 4096))))
    first = family_receipt(64, precision=precision, progress=progress)
    ledger.accept(first)
    for maximum, schedule in ((128, BATCH_SCHEDULE), (256, BATCH_SCHEDULE_256)):
        for low, high, z, alpha in schedule:
            receipt = batch_receipt(low, high, z, alpha, maximum=maximum,
                                    precision=precision, progress=progress)
            ledger.accept(receipt)
    summary = ledger.summary(precision)
    if not summary['complete_local_cover'] or F(summary['accepted_union_log2_upper']) > -40:
        raise ArithmeticError('The established sparse checkpoint did not replay')
    return dict(kind='sparse', precision=precision, construction=ledger.construction,
                accepted=ledger.accepted, summary=summary, fresh_arithmetic=True,
                whole_code_certificate=False)


def tail_checkpoint(directory, precision=256, progress=None):
    inputs, hashes = [], {}
    for kind, name, low, high, expected in PLANS:
        raw = (Path(directory)/name).read_bytes()
        actual = sha256(raw).hexdigest()
        if actual != expected:
            raise ValueError(f'Established proposal changed: {name}')
        inputs.append((name, json.loads(raw), low, high))
        hashes[name] = actual
    ledger, result = replay_swap_cover(inputs, precision=precision, threshold=-40,
                                       progress=progress)
    # replay_swap_cover also adds a separately proved r=4096 endpoint.
    # Export only the requested completed tail, with whole leaf bounds.
    tail = SwapLedger(Box(((5120, 8192), (0, 4096), (0, 4096))))
    tail.import_band(ledger, 5120, 8192, label='established-unsigned-tail')
    summary = tail.summary(precision)
    if not summary['complete_local_cover'] or F(summary['accepted_union_log2_upper']) > -40:
        raise ArithmeticError('The established dense tail did not replay')
    return dict(kind='tail', precision=precision, construction=tail.construction,
                accepted=tail.accepted, summary=summary, proposal_sha256=hashes,
                sources=result['sources'], source_accepted=ledger.accepted,
                fresh_arithmetic=True, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kind', choices=('sparse', 'tail'), required=True)
    parser.add_argument('--precision', type=int, choices=(256, 384), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Frozen artifacts are never overwritten')
    def progress(row):
        if 'active_groups' not in row or row['active_groups'] % 16 == 0:
            print(json.dumps(row, allow_nan=False), flush=True)
    result = (sparse_checkpoint(args.precision, progress) if args.kind == 'sparse'
              else tail_checkpoint(Path(__file__).parent, args.precision, progress))
    result['exporter_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, allow_nan=False, separators=(',', ':'))
        target.write('\n')
    print(json.dumps(dict(path=str(args.output.resolve()), kind=args.kind,
                          precision=args.precision, summary=result['summary'])), flush=True)
