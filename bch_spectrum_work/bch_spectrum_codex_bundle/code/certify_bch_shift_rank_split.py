"""Independently audit a binary refinement of the last rank-28 Q-dual case."""
import json
import sys
from pathlib import Path
from certify_bch_shift_rank import audit_folder,check_path,orbit
from bch_quotient import generator_polynomial,binary_poly_divmod
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]


def build():
    base=audit_folder(39,'shift_rank_q30_complete')
    folder=ROOT/'generated/shift_rank_q15_split23'
    parent_path=ROOT/'generated/shift_rank_q30_complete/case_015.json'
    parent=json.loads(parent_path.read_text())
    summary=json.loads((folder/'summary.json').read_text())
    assert summary['source_folder']=='shift_rank_q30_complete'
    assert summary['source_case_pivot']==15 and summary['split_coset_representative']==23
    zeros=set(parent['zero_indices'])
    assert 15 not in zeros and orbit(23).isdisjoint(zeros|orbit(15))
    leaves=[]
    for filename,zero_branch in (('nonzero.json',False),('zero.json',True)):
        leaf=json.loads((folder/filename).read_text())
        expected=zeros|orbit(23) if zero_branch else zeros
        pivot=15 if zero_branch else 23
        assert leaf['split_coset_is_zero']==zero_branch and leaf['pivot']==pivot
        assert leaf['zero_indices']==sorted(expected)
        state=check_path(expected,pivot,leaf['steps'])
        assert len(state)==leaf['rank']
        assert sorted((e-pivot)%255 for e in state)==leaf['relative_independent_set']
        leaves.append(dict(split_coset_is_zero=zero_branch,pivot=pivot,rank=len(state),steps_checked=len(leaf['steps'])))
    others=[c for c in base['cases'] if c['pivot']!=15]
    assert len(others)==16
    rank=min([c['rank'] for c in others]+[c['rank'] for c in leaves])
    assert rank>=29
    # Q is contained in P, so Pdual is contained in Qdual.
    _,rem=binary_poly_divmod(generator_polynomial(39),generator_polynomial(37))
    assert rem==0
    return dict(classification='Independent full Fourier-case proof with a checked binary refinement',
        base_cover=base,refined_case_pivot=15,split_coset_representative=23,refinement_leaves=leaves,
        complete_refined_leaf_count=18,minimum_certified_rank=rank,
        Qdual_extended_distance_lower=rank+rank%2,Pdual_extended_distance_lower_from_containment=rank+rank%2,
        Q_and_cosets_OA_strength_at_least=rank+rank%2-1,literature_table_used=False,
        original_M22_target_closed=False,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/certify_bch_shift_rank.py',ROOT/'code/bch_quotient.py',
             parent_path,folder/'summary.json',folder/'nonzero.json',folder/'zero.json')})


if __name__=='__main__':
    value=build()
    path=ROOT/'generated/bch256_shift_rank_q30_refined.json'
    if '--verify' in sys.argv:
        assert value==json.loads(path.read_text())
    else:
        write_new(path,value)
    print(json.dumps({k:v for k,v in value.items() if k not in ('base_cover','source_sha256')},indent=2))
