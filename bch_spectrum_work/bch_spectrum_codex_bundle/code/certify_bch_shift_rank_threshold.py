"""Audit another retained Fourier-case cover using the frozen independent checker."""
import argparse
import json
from pathlib import Path
from certify_bch_shift_rank import audit_folder
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]


def build(distance,folder_name):
    value=audit_folder(distance,folder_name)
    return dict(classification='Independent exact Fourier-case certificate; only the minimum checked rank is claimed',
        code=value,extended_dual_distance_lower=value['extended_dual_distance_lower'],
        literature_table_used=False,original_M22_target_closed=False,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/certify_bch_shift_rank.py',ROOT/'code/bch_quotient.py',ROOT/'code/affine_wambach.py')})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('distance',type=int,choices=[37,39])
    parser.add_argument('folder_name')
    parser.add_argument('receipt_name')
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    assert args.receipt_name.isidentifier()
    value=build(args.distance,args.folder_name)
    path=ROOT/'generated'/(args.receipt_name+'.json')
    if args.verify:
        assert value==json.loads(path.read_text())
    else:
        write_new(path,value)
    print(json.dumps(dict(classification=value['classification'],
        extended_dual_distance_lower=value['extended_dual_distance_lower'],
        cases=len(value['code']['cases']),literature_table_used=False),indent=2))
