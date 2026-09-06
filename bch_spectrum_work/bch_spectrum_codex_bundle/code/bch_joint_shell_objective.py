"""Jointly bound the six paired dominant shells with exact objective coefficients."""
import argparse
import importlib
import json
import sys
from fractions import Fraction as Q
from pathlib import Path
from flint import fmpq
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from prepare_bch_hull_probe import export
from verify_scaled_rational_solution import normalized_rows,parse_assignments
from export_scaled_rational_lp import expression
from bch_hull_cut_iteration import solve
from certify_bch_m25_closure import deterministic_caps,old_higher_occupations
from audit_bch_q1_full_arb import encode,decode
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'


def read(p):
    return json.loads(p.read_text())


def build(module_name):
    assert module_name.isidentifier() and module_name.startswith('prepare_bch_')
    module=importlib.import_module(module_name)
    model,scales=module.build()
    assert model==read(module.FOLDER/'model.json') and scales==read(module.FOLDER/'scales.json')
    assert read(module.FOLDER/'audit.json')['all_primal_dual_checks_passed']
    cells=read(GEN/'bch256_q1_activity_cells_outward.json')
    transfer=read(GEN/'bch256_q1_full_arb_transfer.json')
    oa=read(GEN/'oa21_closure_probe/audit.json')
    co={int(w):decode(row['coefficient_upper']) for w,row in transfer['coefficient_rows'].items()}
    co.update({r['weight']:decode(r['coefficient_upper']) for r in cells['shells']})
    caps=deterministic_caps()
    for r in oa['shells']:
        caps[r['weight']]=caps[256-r['weight']]=r['cap']
    higher=Q(128,127)*sum(old_higher_occupations(),Q(0))
    assert sum((caps[w]*c for w,c in co.items()),Q(0))+higher==decode(cells['full_first_moment_upper'])
    paired={w:co[w]+co[256-w] for w in (38,40,42)}
    rest=sum((caps[w]*c for w,c in co.items() if w not in (38,40,42,214,216,218)),Q(0))+higher
    objective={}
    rounded={}
    scale_bits=192
    for w,c in paired.items():
        rounded[w]=Q(-(-(c.numerator<<scale_bits)//c.denominator),1<<scale_bits)
        assert c<=rounded[w]<c+Q(1,1<<scale_bits)
        for a,mult in (('q',1),('h',31)):
            objective[f'{a}_{w}']=rounded[w]*mult*scales[f'{a}_{w}']
    normalization=max(objective.values())
    normalized={n:v/normalization for n,v in objective.items()}
    lp=export(model,scales)[0]
    assert lp.count(' obj: h_38\n')==1
    lp=lp.replace(' obj: h_38\n',' obj: '+expression([(v,n) for n,v in normalized.items()])+'\n')
    data=dict(model_module=module_name,classification='Joint dominant-shell objective, NOT an h38-cap LP',
        physical_objective_normalization=encode(normalization),scaled_objective={n:encode(v) for n,v in normalized.items()},
        exact_pair_coefficients={str(w):encode(v) for w,v in paired.items()},
        dyadic_upper_pair_coefficients={str(w):encode(v) for w,v in rounded.items()},
        constant_rest_upper=encode(rest),higher_occupation_upper=encode(higher),rounding_fractional_bits=scale_bits)
    return module,model,scales,lp,data


def fq(q):
    return fmpq(q.numerator,q.denominator)


def audit(module_name,folder):
    module,model,scales,lp,data=build(module_name)
    assert read(folder/'model.json')==model and read(folder/'scales.json')==scales
    assert read(folder/'objective.json')==data and (folder/'h_38.lp').read_text()==lp
    text=(folder/'h_38.sol').read_text()
    assert 'status = OPTIMAL' in text
    variables=model['metadata']['variables']; rows=normalized_rows(model,scales)
    primal=parse_assignments(text,'VARS:','REDUCED COST:')
    prices=parse_assignments(text,'PI:','SLACK:')
    assert set(primal)<=set(variables) and set(prices)<={f'c{i}' for i in range(1,len(rows)+1)}
    x={n:fq(primal.get(n,Q(0))) for n in variables}
    assert all(v>=0 for v in x.values())
    objective={n:fq(decode(v)) for n,v in data['scaled_objective'].items()}
    combined={n:fmpq(0) for n in variables}; bound=fmpq(0)
    for i,row in enumerate(rows,1):
        coeffs={n:fq(v) for n,v in row['coeffs'].items()}
        lhs=sum((x[n]*v for n,v in coeffs.items()),fmpq(0)); rhs=fq(row['rhs'])
        assert {'eq':lhs==rhs,'ge':lhs>=rhs,'le':lhs<=rhs}[row['sense']],row['name']
        price=fq(prices.get(f'c{i}',Q(0)))
        assert row['sense']=='eq' or (price>=0 if row['sense']=='le' else price<=0)
        bound+=price*rhs
        for n,v in coeffs.items():
            combined[n]+=price*v
    assert all(combined[n]>=objective.get(n,fmpq(0)) for n in variables)
    assert bound==sum((x[n]*v for n,v in objective.items()),fmpq(0))
    physical=Q(int(bound.numerator),int(bound.denominator))*decode(data['physical_objective_normalization'])
    full=physical+decode(data['constant_rest_upper'])
    used={str(w):encode((primal.get(f'q_{w}',Q(0))*scales[f'q_{w}']+31*primal.get(f'h_{w}',Q(0))*scales[f'h_{w}'])) for w in (38,40,42)}
    return dict(classification='Exact joint dominant-shell and full M22 first-moment certificate',
        model_module=module_name,rows_checked=len(rows),variables_checked=len(variables),all_primal_dual_checks_passed=True,
        paired_shells_upper=encode(physical),constant_rest_upper=data['constant_rest_upper'],
        full_M22_first_moment_upper=encode(full),full_M22_margin_bits_diagnostic=-encode(full)['log2_diagnostic'],
        original_M22_target_closed=full<Q(1,1<<40),statistical_shell_assumptions_used=False,
        parameters=dict(memory_bits=22,outer_rows=8192,outer_length=256,message_bits=1<<20,output_bits=1<<21,distance_cutoff=209716),
        diagnostic_primal_shell_values=used,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
            (Path(__file__),Path(module.__file__),module.FOLDER/'audit.json',
             ROOT/'code/verify_scaled_rational_solution.py',ROOT/'code/certify_bch_m25_closure.py',
             GEN/'bch256_q1_activity_cells_outward.json',GEN/'bch256_q1_full_arb_transfer.json',
             GEN/'oa21_closure_probe/audit.json',folder/'model.json',folder/'scales.json',folder/'objective.json',
             folder/'h_38.lp',folder/'h_38.sol')})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['run','audit','verify'])
    parser.add_argument('module')
    parser.add_argument('folder_name')
    args=parser.parse_args()
    assert args.folder_name.isidentifier()
    folder=GEN/args.folder_name
    if args.mode=='run':
        module,model,scales,lp,data=build(args.module)
        folder.mkdir(exist_ok=False)
        write_new(folder/'model.json',model); write_new(folder/'scales.json',scales); write_new(folder/'objective.json',data)
        for name,text in (('h_38.lp',lp),('warm.bas',(module.FOLDER/'h_38.bas').read_text())):
            with (folder/name).open('x') as f:
                f.write(text)
        print(solve(folder,120,True),flush=True)
    else:
        value=audit(args.module,folder)
        if args.mode=='audit':
            write_new(folder/'audit.json',value)
        else:
            assert value==read(folder/'audit.json')
        print(json.dumps({k:v for k,v in value.items() if k not in ('source_sha256','paired_shells_upper',
            'constant_rest_upper','full_M22_first_moment_upper','diagnostic_primal_shell_values')},indent=2))
