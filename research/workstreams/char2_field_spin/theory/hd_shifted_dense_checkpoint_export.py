"""Fresh complete exports of the two established swapped-route dense bands."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

from hd_shifted_cover_ledger import Box
from hd_shifted_swap_ledger import SwapLedger, replay_swap_plan


PLANS = (
    ('hd_shifted_swap_bulk_plan3.json', 4096, 4351,
     '401bce5f667df40cc8431ec50519cfa6af02e450c10df53dda11839360aa7339'),
    ('hd_shifted_swap_bulk_plan4.json', 4352, 5119,
     'f023bfd9c50f0ddb9a1bf14d86630d667cacc5ec751233a22e9f511b53c3c972'),
)


def dense_checkpoint(directory, precision=256, progress=None):
    ledger = SwapLedger(Box(((4096, 5119), (0, 4096), (0, 4096))))
    hashes = {}
    for name, low, high, expected in PLANS:
        raw = (Path(directory)/name).read_bytes()
        actual = sha256(raw).hexdigest()
        if actual != expected:
            raise ValueError(f'Established proposal changed: {name}')
        source, result = replay_swap_plan(json.loads(raw), precision=precision,
                                           threshold=-40, progress=progress)
        if result['failed']:
            raise ArithmeticError(f'An established dense leaf failed replay: {name}')
        ledger.import_band(source, low, high, label=name)
        hashes[name] = actual
    summary = ledger.summary(precision)
    if not summary['complete_local_cover'] or F(summary['accepted_union_log2_upper']) > -40:
        raise ArithmeticError('The established dense checkpoint did not replay')
    return dict(kind='dense', precision=precision, construction=ledger.construction,
                accepted=ledger.accepted, summary=summary, proposal_sha256=hashes,
                fresh_arithmetic=True, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, choices=(256, 384), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Frozen artifacts are never overwritten')
    result = dense_checkpoint(Path(__file__).parent, args.precision,
                               lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['exporter_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, allow_nan=False, separators=(',', ':'))
        target.write('\n')
    print(json.dumps(dict(path=str(args.output.resolve()), precision=args.precision,
                          summary=result['summary'])), flush=True)
