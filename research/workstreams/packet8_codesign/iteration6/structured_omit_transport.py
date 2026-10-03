"""Transport the wider certificate for one identity symbol and15 MDS sandwiches.

The premise is aggregate expected group-measure domination, not a uniform
image for each fixed message at the identity symbol. See the adjacent note.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

import structured_randomizers as s


def run(args):
    if args.output.exists():
        raise ValueError('fresh output path required')
    result = s.correct_receipt(args.receipt, F((1 << 32)-1, 255**4), symbols=15)
    result['schema'] = 'packet8-wider24-one-identity-mds-sandwich-transport-1'
    result['proof_premise'] = ('complete expected group measure dominated by gamma^15 times '
        'the old measure; one identity symbol; other15 independent MDS sandwiches; '
        'group setup and routing independent')
    result['fixed_message_transitivity_at_identity_symbol_claimed'] = False
    result['identity_symbol_position'] = 0
    result['other_symbol_family'] = 'normalized seven-diagonal GF256 MDS sandwich from STRUCTURED_RANDOMIZERS.md'
    for path in (Path(__file__).resolve(), Path(__file__).with_name('STRUCTURED_RANDOMIZERS.md')):
        result['source_sha256'][str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    for name, digest in result['source_sha256'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ArithmeticError('source changed during identity-symbol transport')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('combined_margin_bits_display',
        'q1_margin_bits_display', 'weakest_q', 'exact_union_at_most_2_minus40')}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, default=Path(__file__).resolve().parents[1]/'iteration5'/'whole_wider24_v1.json')
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
