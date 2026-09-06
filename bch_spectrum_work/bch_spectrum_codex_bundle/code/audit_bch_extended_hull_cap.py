"""Reconstruct an extended hull LP and independently check its rational optimum."""
import argparse
import importlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from flint import fmpq
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from prepare_bch_hull_probe import export
from verify_scaled_rational_solution import normalized_rows,parse_assignments
from audit_bch_q1_full_arb import encode,decode
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'


def fq(x):
    return fmpq(x.numerator,x.denominator)


def build(module_name,previous_folder):
    assert module_name.startswith('prepare_bch_') and module_name.isidentifier()
    assert previous_folder.isidentifier()
    module=importlib.import_module(module_name)
    folder=module.FOLDER
    assert folder.parent==GEN
    model,scales=module.build()
    assert model==json.loads((folder/'model.json').read_text())
    assert scales==json.loads((folder/'scales.json').read_text())
    assert export(model,scales)[0]==(folder/'h_38.lp').read_text()
    rows=normalized_rows(model,scales)
    variables=model['metadata']['variables']
    text=(folder/'h_38.sol').read_text()
    assert 'status = OPTIMAL' in text
    primal=parse_assignments(text,'VARS:','REDUCED COST:')
    prices=parse_assignments(text,'PI:','SLACK:')
    assert set(primal)<=set(variables) and set(prices)<={f'c{i}' for i in range(1,len(rows)+1)}
    x={n:fq(primal.get(n,Fraction(0))) for n in variables}
    assert all(v>=0 for v in x.values())
    combined={n:fmpq(0) for n in variables}
    objective=fmpq(0)
    used_rows=[]
    for i,row in enumerate(rows,1):
        coefficients={n:fq(v) for n,v in row['coeffs'].items()}
        lhs=sum((x[n]*v for n,v in coefficients.items()),fmpq(0))
        rhs=fq(row['rhs'])
        assert {'eq':lhs==rhs,'ge':lhs>=rhs,'le':lhs<=rhs}[row['sense']],row['name']
        price=fq(prices.get(f'c{i}',Fraction(0)))
        assert row['sense']=='eq' or (price>=0 if row['sense']=='le' else price<=0)
        if price:
            used_rows.append(row['name'])
        objective+=price*rhs
        for n,v in coefficients.items():
            combined[n]+=price*v
    assert all(combined[n]>=int(n=='h_38') for n in variables)
    assert objective==x['h_38']
    physical=Fraction(int(objective.numerator),int(objective.denominator))*scales['h_38']
    cap=31*(physical.numerator//physical.denominator)
    previous=json.loads((GEN/previous_folder/'audit.json').read_text())
    cells=json.loads((GEN/'bch256_q1_activity_cells_outward.json').read_text())
    transfer=json.loads((GEN/'bch256_q1_full_arb_transfer.json').read_text())
    oldcap=previous['A38_cap_without_lattice_rounding']
    assert cap<=oldcap
    pair=decode(cells['shells'][0]['coefficient_upper'])+decode(transfer['coefficient_rows']['218']['coefficient_upper'])
    full=decode(previous['full_M22_first_moment_upper'])-(oldcap-cap)*pair
    lower=31*physical*decode(cells['shells'][0]['coefficient_lower'])
    paths=[Path(__file__),Path(module.__file__),ROOT/'code/prepare_bch_hull_probe.py',
           ROOT/'code/verify_scaled_rational_solution.py',ROOT/'code/export_lp.py',
           GEN/previous_folder/'audit.json',GEN/'bch256_q1_activity_cells_outward.json',
           GEN/'bch256_q1_full_arb_transfer.json',folder/'model.json',folder/'scales.json',
           folder/'h_38.lp',folder/'h_38.sol']
    return folder,dict(classification='Exact reconstructed extended hull LP certificate',
                model_module=module_name,previous_folder=previous_folder,
                rows_checked=len(rows),variables_checked=len(variables),all_primal_dual_checks_passed=True,
                h38_optimum=encode(physical),A38_cap_without_lattice_rounding=cap,previous_A38_cap=oldcap,
                full_M22_first_moment_upper=encode(full),full_M22_margin_bits_diagnostic=-encode(full)['log2_diagnostic'],
                feasible_primal_true_tail_lower=encode(lower),
                relaxation_insufficient_even_with_exact_inner_tails=lower>Fraction(1,1<<40),
                original_M22_target_closed=full<Fraction(1,1<<40),nonzero_dual_row_names=used_rows,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('module')
    parser.add_argument('previous_folder')
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    folder,result=build(args.module,args.previous_folder)
    if args.verify:
        assert json.loads((folder/'audit.json').read_text())==result
    else:
        write_new(folder/'audit.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','h38_optimum',
            'full_M22_first_moment_upper','feasible_primal_true_tail_lower','nonzero_dual_row_names')},indent=2))
