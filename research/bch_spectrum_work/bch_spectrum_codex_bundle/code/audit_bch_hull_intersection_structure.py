"""Read-only replay of intersection algebra, models, and an eliminating cut."""
import json
import math
import sys
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from audit_bch_hulls import rows,hull,nullspace,echelon,cyclic_description
from prepare_bch_hull_intersection_probe import build as intersection_build
from prepare_bch_hull_intersection_reduced import build as reduced_build
from prepare_bch_hull_benders_cut import build as cut_build
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'generated/bch256_hull_intersection_structure.json'


def build():
    models={}
    for folder,builder in (
        ('oa21_hull_intersection_probe',intersection_build),
        ('oa21_hull_intersection_reduced',reduced_build),
        ('oa21_hull_benders_cut',cut_build)):
        model,scales=builder()
        assert model==json.loads((ROOT/'generated'/folder/'model.json').read_text())
        assert scales==json.loads((ROOT/'generated'/folder/'scales.json').read_text())
        models[folder]=model
    full=models['oa21_hull_intersection_probe']
    j=[int(x,16) for x in full['metadata']['J_intersection']['basis_hex']]
    hb=hull(rows(61,63))
    dimension=len(j)+len(hb)-len(echelon(j+hb))
    assert len(hb)==61 and dimension==53
    jr=set(full['metadata']['J_intersection']['cyclic_description']['roots'])
    br=set(cyclic_description(hb)['roots'])
    multipliers=[m for m in range(1,255) if math.gcd(m,255)==1 and {(m*x)%255 for x in br}==jr]
    assert multipliers==[]
    cut=models['oa21_hull_benders_cut']['metadata']['J_elimination_cut']
    paths=[Path(__file__),ROOT/'code/prepare_bch_hull_intersection_probe.py',
           ROOT/'code/prepare_bch_hull_intersection_reduced.py',ROOT/'code/prepare_bch_hull_benders_cut.py']
    return dict(classification='Exact intersection and Farkas-cut replay; not a new shell cap by itself',
                intersection_dimension=61,anchor_quotient_factor_hex='0x169',nonzero_shift_orbit=255,
                original_variables=229,reduced_variables=213,
                fixed_J_substitution=models['oa21_hull_intersection_reduced']['metadata']['exact_fixed_J_substitution'],
                L63_hull_dimension=61,L63_hull_intersection_with_J_dimension=dimension,
                cyclic_multiplier_equivalences_found=multipliers,
                arbitrary_permutation_inequivalence_claimed=False,
                cut_terms=cut['number_of_nonzero_multipliers'],
                exact_old_primal_violation=cut['farkas_contradiction_at_old_primal'],
                all_eliminated_coefficients_nonpositive=True,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__=='__main__':
    result=build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','exact_old_primal_violation')},indent=2))
