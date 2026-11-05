"""Replay one preserved full marker bank for the new full-field ensemble."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic

from hd_field128_lift_alias_certificate import AliasBank
from hd_field128_lift_certificate import Budget


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--center', type=int, required=True)
    parser.add_argument('--ratio', required=True, help='Exact positive common state potential')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--max-seconds', type=float, default=180.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists')
    directory = Path(__file__).resolve().parent
    manifest = directory/'hd_field128_lift_tail_seed_manifest.json'
    seed = next(row for row in json.loads(manifest.read_text())['seeds']
                if row['central_weight']==args.center)
    source = directory/seed['source_path']
    digest = sha256(source.read_bytes()).hexdigest()
    if digest != seed['source_sha256']:
        raise RuntimeError('Preserved source pin mismatch')
    original = json.loads(source.read_text())['fixed_witness']['source_exact_witness']
    markers = [tuple(map(F, row.split())) for row in original['markers']]
    # Decimal rational proposals avoid needlessly huge denominators from
    # raising a binary64 rational to the sixteenth power. These exact values,
    # not the old witness, are passed to the directed verifier and recorded.
    z = F(str(float(F(seed['source_exact_z'])**16)))
    alpha = F(str(float(F(seed['source_exact_alpha'])/16)))
    budget = Budget(args.max_seconds)
    bank = AliasBank(z,alpha,markers,precision=args.precision,max_seconds=args.max_seconds,
                     progress=lambda row:print(json.dumps(row),flush=True))
    budget.check()
    receipt = bank.weight_receipt(args.center,F(args.ratio),method='runs',
                                  max_seconds=args.max_seconds-(monotonic()-budget.started))
    receipt.update(source_witness_path=source.name,source_witness_sha256=digest,
                   driver_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                   seed_manifest_sha256=sha256(manifest.read_bytes()).hexdigest(),
                   elapsed_seconds=monotonic()-budget.started)
    with args.output.open('x',encoding='utf-8') as out:
        json.dump(receipt,out,indent=2,allow_nan=False)
        out.write('\n')
    print(json.dumps({key:receipt[key] for key in ('active_groups','log2_upper_rational',
                     'log2_upper_enclosure','elapsed_seconds','whole_code_certificate')},indent=2))


if __name__=='__main__':
    main()
