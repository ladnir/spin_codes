"""Read-only exact audit of the modulo-four refinement and M22 aggregation."""
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from audit_bch_quadratic_sums import build as quadratic_audit
from audit_bch_q1_full_arb import encode,decode
from verify_scaled_rational_solution import normalized_rows,parse_assignments
from export_scaled_rational_lp import variable_scales
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
FOLDER=GEN/'oa21_mod4_probe'
OUTPUT=FOLDER/'audit.json'


def build():
    sums=quadratic_audit()
    assert sums==json.loads((GEN/'bch256_quadratic_weight_sums.json').read_text())
    old=json.loads((GEN/'oa21_closure_probe/model.json').read_text())
    model=json.loads((FOLDER/'model.json').read_text())
    expected=old['constraints'].copy()
    for prefix,label in (('q','Q'),('h','H')):
        expected.append(dict(name=f'exact_mod4_{prefix}',
                            coeffs={f'{prefix}_{w}':str((1 if w==128 else 2)*(-1)**(w//2)) for w in range(0,129,2)},
                            sense='eq',rhs=sums['results'][label]['signed_sum']))
    assert model['constraints']==expected
    assert model['metadata']['variables']==old['metadata']['variables']
    scales=variable_scales(model)
    rows=normalized_rows(model,scales)
    variables=model['metadata']['variables']
    text=(FOLDER/'h_38.sol').read_text()
    assert 'status = OPTIMAL' in text
    raw=parse_assignments(text,'VARS:','REDUCED COST:')
    prices=parse_assignments(text,'PI:','SLACK:')
    assert set(raw)<=set(variables) and set(prices)<={f'c{i}' for i in range(1,len(rows)+1)}
    x={n:raw.get(n,Fraction(0)) for n in variables}
    assert all(v>=0 for v in x.values())
    combined={n:Fraction(0) for n in variables}
    objective=Fraction(0)
    for i,row in enumerate(rows,1):
        actual=sum((x[n]*v for n,v in row['coeffs'].items()),Fraction(0))
        rhs=row['rhs']
        assert {'eq':actual==rhs,'le':actual<=rhs,'ge':actual>=rhs}[row['sense']],row['name']
        price=prices.get(f'c{i}',Fraction(0))
        assert row['sense']=='eq' or (price>=0 if row['sense']=='le' else price<=0)
        objective+=price*rhs
        for n,v in row['coeffs'].items():
            combined[n]+=price*v
    assert all(combined[n]>=int(n=='h_38') for n in variables)
    assert objective==x['h_38']
    physical=objective*scales['h_38']
    cap=31*(physical.numerator//physical.denominator)
    oa=json.loads((GEN/'oa21_closure_probe/audit.json').read_text())
    cells=json.loads((GEN/'bch256_q1_activity_cells_outward.json').read_text())
    transfer=json.loads((GEN/'bch256_q1_full_arb_transfer.json').read_text())
    previous=oa['shells'][0]['cap']
    assert cap<previous
    pair=decode(cells['shells'][0]['coefficient_upper'])+decode(transfer['coefficient_rows']['218']['coefficient_upper'])
    full=decode(cells['full_first_moment_upper'])-(previous-cap)*pair
    lower=31*physical*decode(cells['shells'][0]['coefficient_lower'])
    assert lower>Fraction(1,1<<40)
    paths=[Path(__file__),ROOT/'code/audit_bch_quadratic_sums.py',ROOT/'code/verify_scaled_rational_solution.py',
           ROOT/'code/export_scaled_rational_lp.py',ROOT/'code/export_lp.py',ROOT/'code/prepare_bch_mod4_probe.py',
           GEN/'bch256_quadratic_weight_sums.json',GEN/'oa21_closure_probe/model.json',GEN/'oa21_closure_probe/audit.json',
           GEN/'bch256_q1_activity_cells_outward.json',GEN/'bch256_q1_full_arb_transfer.json',
           FOLDER/'model.json',FOLDER/'h_38.lp',FOLDER/'h_38.sol']
    return dict(classification='Exact OA21 plus signed-weight identities LP certificate',
                rows_checked=len(rows),variables_checked=len(variables),all_primal_dual_checks_passed=True,
                h38_optimum=encode(physical),A38_cap_without_lattice_rounding=cap,
                previous_A38_cap=previous,improvement_bits_diagnostic=math.log2(previous)-math.log2(cap),
                full_M22_first_moment_upper=encode(full),full_M22_margin_bits_diagnostic=-encode(full)['log2_diagnostic'],
                feasible_primal_true_tail_lower=encode(lower),
                OA21_plus_mod4_still_insufficient_even_with_exact_inner_tails=True,
                A38_cap_at_most_10_to_13=cap<=10**13,original_M22_target_closed=False,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__=='__main__':
    result=build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps({k:result[k] for k in ('rows_checked','all_primal_dual_checks_passed',
                     'A38_cap_without_lattice_rounding','improvement_bits_diagnostic','full_M22_margin_bits_diagnostic',
                     'OA21_plus_mod4_still_insufficient_even_with_exact_inner_tails','original_M22_target_closed')},indent=2))
