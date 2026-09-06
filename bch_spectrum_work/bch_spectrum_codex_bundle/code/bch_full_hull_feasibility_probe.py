"""Configurable exact feasibility test, never a BCH upper-bound model."""
import argparse
import copy
import importlib
import json
import sys
from fractions import Fraction as Q
from pathlib import Path
from flint import fmpq
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from prepare_bch_hull_intersection_probe import build as full_build
from prepare_bch_hull_intersection_reduced import build as reduced_build
from prepare_bch_hull_probe import export
from bch_hull_cut_iteration import warm_basis,solve,witness,check_master
from audit_bch_q1_full_arb import encode,decode
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'


def read(p):
    return json.loads(p.read_text())


def build(module_name):
    assert module_name.isidentifier() and module_name.startswith('prepare_bch_')
    module=importlib.import_module(module_name)
    previous=module.FOLDER
    old,scales=module.build()
    assert old==read(previous/'model.json') and scales==read(previous/'scales.json')
    reduced,rscales=reduced_build()
    model=copy.deepcopy(old)
    oldnames={r['name']:r for r in old['constraints']}
    newvars=[n for n in reduced['metadata']['variables'] if n not in old['metadata']['variables']]
    model['metadata']['variables']+=newvars
    newscales=scales|{n:rscales[n] for n in newvars}
    for row in reduced['constraints']:
        if row['name'] not in oldnames:
            model['constraints'].append(row)
        else:
            other=oldnames[row['name']]
            assert row['sense']==other['sense'] and int(row['rhs'])==int(other['rhs'])
            assert {n:int(c) for n,c in row['coeffs'].items() if int(c)}=={n:int(c) for n,c in other['coeffs'].items() if int(c)}
    model['constraints'].append(dict(name='AUXILIARY_h38_lower_for_relaxation_only',coeffs={'h_38':'1'},sense='ge',rhs=str(5*10**12)))
    model['metadata']['classification']='Auxiliary feasibility test, NOT a valid BCH upper-bound model'
    model['metadata']['fixed_J_substitution']=reduced['metadata']['exact_fixed_J_substitution']
    model['metadata']['source_master_module']=module_name
    lp=export(model,newscales)[0]
    assert lp.count(' obj: h_38\n')==1
    lp=lp.replace(' obj: h_38\n',' obj: 0 h_38\n')
    basis,mapping=warm_basis(old,scales,model,newscales,previous/'h_38.bas')
    return previous,model,newscales,lp,basis,mapping


def audit(module_name,folder):
    previous,model,scales,lp,basis,mapping=build(module_name)
    assert model==read(folder/'model.json') and scales==read(folder/'scales.json')
    assert lp==(folder/'h_38.lp').read_text() and basis==(folder/'warm.bas').read_text()
    assert mapping==read(folder/'warm_mapping.json')
    raw=witness(folder,model['metadata']['variables'])
    physical={n:x*scales[n] for n,x in raw.items()}
    assert all(x>=0 for x in physical.values())
    def check(rows,x):
        for row in rows:
            lhs=sum((fmpq(int(c))*fmpq(x[n].numerator,x[n].denominator) for n,c in row['coeffs'].items()),fmpq(0))
            rhs=fmpq(int(row['rhs']))
            assert {'eq':lhs==rhs,'ge':lhs>=rhs,'le':lhs<=rhs}[row['sense']],row['name']
    check(model['constraints'],physical)
    full,_=full_build()
    restored=physical|{n:Q(v) for n,v in model['metadata']['fixed_J_substitution'].items()}
    assert set(restored)==set(full['metadata']['variables'])
    check(full['constraints'],restored)
    cells=read(GEN/'bch256_q1_activity_cells_outward.json')
    lower=31*physical['h_38']*decode(cells['shells'][0]['coefficient_lower'])
    assert lower>Q(1,1<<40)
    paths=[Path(__file__),ROOT/'code'/f'{module_name}.py',
        ROOT/'code/prepare_bch_hull_intersection_probe.py',ROOT/'code/prepare_bch_hull_intersection_reduced.py',
        previous/'model.json',previous/'scales.json',previous/'audit.json',
        GEN/'bch256_q1_activity_cells_outward.json',folder/'model.json',folder/'scales.json',folder/'h_38.sol']
    return dict(classification='Exact hypothetical spectrum proves insufficiency of the specified relaxation',
        source_master_module=module_name,physical_A38=encode(31*physical['h_38']),
        weight38_true_tail_first_moment_lower=encode(lower),
        all_nonnegative_and_all_rows_passed=True,original_intersection_raw_rows_checked=len(full['constraints']),
        augmented_raw_rows_checked=len(model['constraints']),relaxation_insufficient_even_with_exact_inner_tails=True,
        actual_BCH_spectrum_claimed=False,actual_failure_probability_lower_bound_claimed=False,
        original_M22_target_closed=False,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['run','audit','verify'])
    parser.add_argument('module')
    parser.add_argument('folder')
    args=parser.parse_args()
    assert args.folder.isidentifier()
    folder=GEN/args.folder
    if args.mode=='run':
        previous,model,scales,lp,basis,mapping=build(args.module)
        folder.mkdir(exist_ok=False)
        write_new(folder/'model.json',model)
        write_new(folder/'scales.json',scales)
        write_new(folder/'warm_mapping.json',mapping)
        for name,text in (('h_38.lp',lp),('warm.bas',basis)):
            with (folder/name).open('x') as f:
                f.write(text)
        print(solve(folder,120,True),flush=True)
    else:
        value=audit(args.module,folder)
        if args.mode=='audit':
            write_new(folder/'audit.json',value)
        else:
            assert value==read(folder/'audit.json')
        print(json.dumps({k:v for k,v in value.items() if k not in ('source_sha256','physical_A38','weight38_true_tail_first_moment_lower')},indent=2))
