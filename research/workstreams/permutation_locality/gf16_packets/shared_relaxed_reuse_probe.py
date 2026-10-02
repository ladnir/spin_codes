"""Check a retained parent witness and its two restrictions on the real model.

This is a search diagnostic, not a full cover or a whole-code certificate.
No numerical endpoint from the source is used in an inequality.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

from flint import arb, ctx
import shared_relaxed_parallel as parallel
import shared_relaxed_strategy as proof


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--path', required=True)
    parser.add_argument('--extra-cutoff', type=int, default=0)
    parser.add_argument('--target-bits', type=int, default=40)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if (args.output.exists() or args.output.resolve() == args.source.resolve()
            or args.extra_cutoff < 0 or args.target_bits < 1 or args.precision < 256):
        parser.error('new output, nonnegative cutoff increment, and valid precision/margin required')
    raw = args.source.read_bytes()
    record = json.loads(raw)
    row = proof.validate_record(record)
    original_cutoff = row['threshold']
    if args.path not in row['cover']['leaves']:
        parser.error('select a saved leaf with a rational witness')
    witness = row['cover']['leaves'][args.path]['witness']
    source_sha = hashlib.sha256(raw).hexdigest()
    if args.extra_cutoff:
        record = parallel.retarget_record(record,
            Q(original_cutoff+args.extra_cutoff, proof.sc.N), 0, source_sha)
        row = record['results'][0]
    model = proof.build_model(record, args.precision)
    cells = proof.sc.partition(model, row['cover']['leaves'], row['cover']['unresolved'])
    parent = cells[args.path]
    midpoint = sum(parent)/2
    scopes = [('parent', parent, witness)]
    for name, cell in (('left', (parent[0], midpoint)), ('right', (midpoint, parent[1]))):
        scopes.append((name, cell, parallel.restrict_witness(parent, cell, witness)))
    result = dict(schema='shared-gf16-witness-reuse-probe-1', source_sha256=source_sha,
        source_threshold=original_cutoff, threshold=row['threshold'], path=args.path,
        precision=args.precision, target_bits=args.target_bits, rows=[],
        note='Only this parent and its children were freshly checked; no complete certificate.')
    for name, cell, candidate in scopes:
        ctx.prec = args.precision
        value = model.outward(cell, candidate)
        if ctx.prec != args.precision or not value > 0:
            raise ArithmeticError('positive bound at the requested precision required')
        result['rows'].append(dict(kind=name, cell=list(map(str, cell)),
            upper=[int(v) for v in value.upper().man_exp()],
            margin_passed=bool(value < arb(2)**-args.target_bits), witness=candidate))
        print('REUSE', name, 'log2 upper', value.log()/arb(2).log(),
              'passed', result['rows'][-1]['margin_passed'], flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
