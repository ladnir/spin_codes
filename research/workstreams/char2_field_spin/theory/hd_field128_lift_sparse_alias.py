"""Propose nearby sparse alias witnesses and replay every weight in each band.

Scaling a saved witness is only a proposal, never an inference of coverage.
Every reported bound is recomputed for the full-field ensemble by the
frozen directed verifier. Each completed band is saved independently.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import hd_field128_lift_alias_cover as cover
from hd_field128_lift_alias_certificate import AliasBank
from hd_field128_lift_certificate import FIELD_ORDER, CONSTRUCTION
from hd_shifted_swap_potential import SAVED_RUN_WITNESSES

BANDS = {80:(65,96), 96:(81,112), 160:(140,184), 192:(160,224),
         224:(192,256), 256:(224,279)}

def proposal(center):
    if center not in BANDS:
        raise ValueError('Use an explicit bounded center')
    old = SAVED_RUN_WITNESSES[128]
    z = F(str(float(F(old['z'])**(center//8))))
    alpha = F(str(float(F(old['alpha'])/16)))
    scale = F(center,128)
    markers = [tuple(F(x)*scale for x in pair.split()) for pair in old['markers']]
    ratio = F(old['potential'])/scale
    return z,alpha,markers,tuple(ratio*f for f in (F(1,2),F(1),F(3,2),F(2)))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--centers',type=int,nargs='+',default=list(BANDS))
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--build-seconds',type=float,default=180)
    parser.add_argument('--band-seconds',type=float,default=180)
    args = parser.parse_args()
    if args.output_dir.exists() or any(c not in BANDS for c in args.centers):
        parser.error('Use a fresh directory and listed centers')
    args.output_dir.mkdir(parents=False)
    directory = Path(__file__).resolve().parent
    pins = {name:sha256((directory/name).read_bytes()).hexdigest() for name in
            ('hd_field128_lift_sparse_alias.py','hd_shifted_swap_potential.py')}
    for center in args.centers:
        z,alpha,markers,ratios = proposal(center)
        witness = dict(construction=CONSTRUCTION,field_order=FIELD_ORDER,
                       full_field_random_coefficients=True,pair_swaps=True,
                       status='proposal only; no numerical bound',directed=False,
                       whole_code_certificate=False,central_weight=center,
                       state_potential_ratio=str(ratios[1]),
                       witness=dict(z=str(z),alpha=str(alpha),markers=[list(map(str,p)) for p in markers],grid=128),
                       source_sha256=pins)
        with (args.output_dir/f'center{center}_proposal.json').open('x') as out:
            json.dump(witness,out,indent=2,allow_nan=False)
        print(json.dumps(dict(stage='build',center=center)),flush=True)
        bank = AliasBank(z,alpha,markers,max_seconds=args.build_seconds)
        first,last = BANDS[center]
        receipt = cover.band_receipt(bank,first,last,ratios,max_weights=last-first+1,
                                     max_seconds=args.band_seconds)
        receipt.update(proposal_source_sha256=pins,proposal_center=center)
        with (args.output_dir/f'center{center}_receipt.json').open('x') as out:
            json.dump(receipt,out,indent=2,allow_nan=False)
        good = [row['active_groups'] for row in receipt['weight_receipts']
                if F(row['log2_upper_rational'])<=-80]
        print(json.dumps(dict(stage='completed',center=center,
                              good_at_80_bits=good,complete=receipt['complete_requested_range'],
                              failure=receipt['failure'])),flush=True)

if __name__=='__main__':
    main()
