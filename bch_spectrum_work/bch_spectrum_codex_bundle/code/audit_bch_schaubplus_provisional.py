"""Verify exact LP arithmetic without promoting unaudited literature inputs."""
import json
import sys
from pathlib import Path
from audit_bch_extended_hull_cap import build as exact_audit
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]


def build():
    folder,value=exact_audit('prepare_bch_schaubplus_provisional','oa21_hull_dual18_probe')
    value['classification']='Exact LP sensitivity result CONDITIONAL on unaudited SchaubPlus literature inputs; NOT a certified BCH bound'
    value['primary_literature_input_audit_pending']=True
    value['conditional_aggregate_below_target']=value.pop('original_M22_target_closed')
    value['original_M22_target_closed']=False
    value['conditional_A38_cap']=value.pop('A38_cap_without_lattice_rounding')
    value['conditional_full_M22_upper']=value.pop('full_M22_first_moment_upper')
    value['conditional_margin_bits']=value.pop('full_M22_margin_bits_diagnostic')
    value['source_sha256'][str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    return folder,value


if __name__=='__main__':
    folder,value=build()
    path=folder/'provisional_audit.json'
    if '--verify' in sys.argv:
        assert value==json.loads(path.read_text())
    else:
        write_new(path,value)
    print(json.dumps({k:v for k,v in value.items() if k not in
        ('source_sha256','h38_optimum','conditional_full_M22_upper','feasible_primal_true_tail_lower','nonzero_dual_row_names')},indent=2))
