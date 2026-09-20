"""Read-only exact rational audit of the coupled hull spectrum LP."""
import json
import sys
from fractions import Fraction
from pathlib import Path
from flint import fmpq
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from prepare_bch_hull_probe import build as model_build, export
from check_bch_hull_fourier import build as toy_build
from verify_scaled_rational_solution import normalized_rows,parse_assignments
from audit_bch_q1_full_arb import encode,decode
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
FOLDER=GEN/'oa21_hull_probe'
OUTPUT=FOLDER/'audit.json'


def fq(value):
    return fmpq(value.numerator,value.denominator)


def build():
    toy=toy_build()
    assert toy==json.loads((GEN/'bch256_hull_fourier_toy.json').read_text())
    model,scales=model_build()
    assert model==json.loads((FOLDER/'model.json').read_text())
    assert scales==json.loads((FOLDER/'scales.json').read_text())
    assert export(model,scales)[0]==(FOLDER/'h_38.lp').read_text()
    rows=normalized_rows(model,scales)
    variables=model['metadata']['variables']
    text=(FOLDER/'h_38.sol').read_text()
    assert 'status = OPTIMAL' in text
    raw=parse_assignments(text,'VARS:','REDUCED COST:')
    prices=parse_assignments(text,'PI:','SLACK:')
    assert set(raw)<=set(variables)
    assert set(prices)<={f'c{i}' for i in range(1,len(rows)+1)}
    x={n:fq(raw.get(n,Fraction(0))) for n in variables}
    assert all(v>=0 for v in x.values())
    combined={n:fmpq(0) for n in variables}
    objective=fmpq(0)
    for i,row in enumerate(rows,1):
        coeffs={n:fq(v) for n,v in row['coeffs'].items()}
        actual=sum((x[n]*v for n,v in coeffs.items()),fmpq(0))
        rhs=fq(row['rhs'])
        assert {'eq':actual==rhs,'le':actual<=rhs,'ge':actual>=rhs}[row['sense']],row['name']
        price=fq(prices.get(f'c{i}',Fraction(0)))
        assert row['sense']=='eq' or (price>=0 if row['sense']=='le' else price<=0)
        objective+=price*rhs
        for n,v in coeffs.items():
            combined[n]+=price*v
    assert all(combined[n]>=int(n=='h_38') for n in variables)
    assert objective==x['h_38']
    physical=Fraction(int(objective.numerator),int(objective.denominator))*scales['h_38']
    cap=31*(physical.numerator//physical.denominator)
    previous=json.loads((GEN/'oa21_mod4_probe/audit.json').read_text())
    cells=json.loads((GEN/'bch256_q1_activity_cells_outward.json').read_text())
    transfer=json.loads((GEN/'bch256_q1_full_arb_transfer.json').read_text())
    oldcap=previous['A38_cap_without_lattice_rounding']
    assert cap<=oldcap
    pair=decode(cells['shells'][0]['coefficient_upper'])+decode(transfer['coefficient_rows']['218']['coefficient_upper'])
    full=decode(previous['full_M22_first_moment_upper'])-(oldcap-cap)*pair
    lower=31*physical*decode(cells['shells'][0]['coefficient_lower'])
    paths=[Path(__file__),ROOT/'code/prepare_bch_hull_probe.py',ROOT/'code/check_bch_hull_fourier.py',
           ROOT/'code/audit_bch_hulls.py',GEN/'bch256_hull_structure.json',GEN/'bch256_hull_fourier_toy.json',
           GEN/'oa21_mod4_probe/audit.json',GEN/'bch256_q1_activity_cells_outward.json',
           GEN/'bch256_q1_full_arb_transfer.json',FOLDER/'model.json',FOLDER/'scales.json',
           FOLDER/'h_38.lp',FOLDER/'h_38.sol']
    return dict(classification='Exact coupled hull LP certificate; anchor and OA21 inputs inherited',
                rows_checked=len(rows),variables_checked=len(variables),all_primal_dual_checks_passed=True,
                h38_optimum=encode(physical),A38_cap_without_lattice_rounding=cap,previous_A38_cap=oldcap,
                full_M22_first_moment_upper=encode(full),full_M22_margin_bits_diagnostic=-encode(full)['log2_diagnostic'],
                feasible_primal_true_tail_lower=encode(lower),
                relaxation_insufficient_even_with_exact_inner_tails=lower>Fraction(1,1<<40),
                original_M22_target_closed=full<Fraction(1,1<<40),
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__=='__main__':
    result=build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','h38_optimum',
            'full_M22_first_moment_upper','feasible_primal_true_tail_lower')},indent=2))
