"""Verify a Farkas witness and eliminate the unknown J spectrum exactly."""
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from prepare_bch_hull_intersection_probe import build as full_build
from prepare_bch_hull_fixed_extension import build as extension_build
from prepare_bch_hull_reciprocal_probe import build as previous_build
from verify_scaled_rational_solution import parse_assignments
from prepare_bch_hull_probe import export
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/oa21_hull_benders_cut'
FARKAS=ROOT/'generated/oa21_hull_extension_farkas'


def build():
    extension=extension_build()
    assert extension==json.loads((ROOT/'generated/oa21_hull_fixed_extension/model.json').read_text())
    source=json.loads((FARKAS/'model.json').read_text())
    text=(FARKAS/'h_38.sol').read_text()
    assert 'status = OPTIMAL' in text
    raw=parse_assignments(text,'VARS:','REDUCED COST:')
    assert set(raw)<=set(source['variables'])
    y={n:raw.get(n,Fraction(0)) for n in source['variables']}
    assert all(v>=0 for v in y.values()) and sum(y.values())<=1
    expected_terms=[]
    for i,row in enumerate(extension['constraints']):
        for sign in ((1,-1) if row['sense']=='eq' else ((1,) if row['sense']=='ge' else (-1,))):
            expected_terms.append(dict(variable=f'f_{len(expected_terms)}',row_index=i,sign=sign))
    assert source['terms']==expected_terms
    combined={n:Fraction(0) for n in extension['variables']}
    contradiction=Fraction(0)
    for term in expected_terms:
        row=extension['constraints'][term['row_index']]
        mult=term['sign']*y[term['variable']]
        contradiction+=mult*Fraction(row['rhs'])
        for n,c in row['coeffs'].items():
            combined[n]+=mult*Fraction(c)
    assert all(v<=0 for v in combined.values()) and contradiction>0

    full,fullscales=full_build()
    fullrows={r['name']:r for r in full['constraints']}
    fixed={n:Fraction(v) for n,v in extension['fixed_physical'].items()}
    old,scales=previous_build()
    kept=set(old['metadata']['variables'])
    cut={n:Fraction(0) for n in kept}
    residual={n:Fraction(0) for n in extension['variables']}
    cut_rhs=Fraction(0)
    used=[]
    for term in expected_terms:
        signed=term['sign']*y[term['variable']]
        if not signed:
            continue
        row=extension['constraints'][term['row_index']]
        original=fullrows[row['name']]
        rhs_at_point=Fraction(int(original['rhs']))-sum((int(c)*fixed[n] for n,c in original['coeffs'].items() if n in fixed),Fraction(0))
        remaining={n:int(c)*fullscales[n] for n,c in original['coeffs'].items() if n not in fixed and int(c)}
        norm=max(Fraction(1),abs(rhs_at_point),*(abs(v) for v in remaining.values()))
        assert row['rhs']==str(rhs_at_point/norm)
        assert row['coeffs']=={n:str(Fraction(c)/norm) for n,c in remaining.items()}
        mult=signed/norm
        # Unfix all old spectrum coordinates, but retain t_0=1 and zero J shells.
        rhs=Fraction(int(original['rhs']))-sum((int(c)*fixed[n] for n,c in original['coeffs'].items() if n in fixed and n not in kept),Fraction(0))
        cut_rhs+=mult*rhs
        for n,c in original['coeffs'].items():
            if n in kept:
                cut[n]+=mult*int(c)
            elif n in residual:
                residual[n]+=mult*int(c)
        used.append(dict(row=row['name'],physical_multiplier=str(mult)))
    assert all(v<=0 for v in residual.values())
    violation=cut_rhs-sum((c*fixed[n] for n,c in cut.items()),Fraction(0))
    assert violation==contradiction
    denominator=math.lcm(cut_rhs.denominator,*(v.denominator for v in cut.values()))
    coefficients={n:int(v*denominator) for n,v in cut.items() if v}
    rhs=int(cut_rhs*denominator)
    divisor=math.gcd(abs(rhs),*(abs(v) for v in coefficients.values()))
    coefficients={n:str(v//divisor) for n,v in sorted(coefficients.items())}
    rhs//=divisor
    old['constraints'].append(dict(name='J_Farkas_elimination_cut_1',coeffs=coefficients,sense='ge',rhs=str(rhs)))
    old['metadata']['J_elimination_cut']={
        'farkas_contradiction_at_old_primal':str(contradiction),
        'all_eliminated_coefficients_nonpositive':True,
        'eliminated_physical_coefficients':{n:str(v) for n,v in residual.items()},
        'used_original_rows':used,'number_of_nonzero_multipliers':len(used),
        'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),ROOT/'code/prepare_bch_hull_intersection_probe.py',
             ROOT/'code/prepare_bch_hull_fixed_extension.py',ROOT/'code/prepare_bch_hull_reciprocal_probe.py',
             FARKAS/'model.json',FARKAS/'h_38.sol')}}
    return old,scales


if __name__=='__main__':
    model,scales=build()
    lp,count=export(model,scales)
    FOLDER.mkdir(exist_ok=False)
    write_new(FOLDER/'model.json',model)
    write_new(FOLDER/'scales.json',scales)
    with (FOLDER/'h_38.lp').open('x') as f:
        f.write(lp)
    print(json.dumps(dict(variables=len(scales),rows=count,
          farkas_terms=model['metadata']['J_elimination_cut']['number_of_nonzero_multipliers'])))
