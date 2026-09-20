"""Seek a counterexample to sufficiency of the relaxation, not an actual code.

The auxiliary h38 lower bound is NOT a BCH constraint and must never produce
a claimed spectrum upper bound. A feasible rational point only limits this LP.
"""
import argparse
import copy
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
from bch_hull_cut_iteration import warm_basis,solve,witness
from audit_bch_q1_full_arb import encode,decode
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
FOLDER=GEN/'hull_intersection_feasibility'
PREVIOUS=GEN/'hull_cut_iterations/round_002/master'


def read(p):
    return json.loads(p.read_text())


def text_new(p,text):
    with p.open('x',encoding='ascii') as f:
        f.write(text)


def build():
    old,scales=read(PREVIOUS/'model.json'),read(PREVIOUS/'scales.json')
    reduced,rscales=reduced_build()
    oldnames={r['name']:r for r in old['constraints']}
    model=copy.deepcopy(old)
    newvars=[n for n in reduced['metadata']['variables'] if n not in model['metadata']['variables']]
    model['metadata']['variables']+=newvars
    newscales=scales|{n:rscales[n] for n in newvars}
    for row in reduced['constraints']:
        if row['name'] not in oldnames:
            model['constraints'].append(row)
        else:
            other=oldnames[row['name']]
            assert row['sense']==other['sense'] and int(row['rhs'])==int(other['rhs'])
            assert {n:int(c) for n,c in row['coeffs'].items() if int(c)}=={n:int(c) for n,c in other['coeffs'].items() if int(c)}
    model['constraints'].append(dict(name='AUXILIARY_h38_lower_for_relaxation_test',coeffs={'h_38':'1'},sense='ge',rhs=str(5*10**12)))
    model['metadata']['classification']='Auxiliary feasibility test of a relaxation, NOT a valid BCH upper-bound model'
    model['metadata']['fixed_J_substitution']=reduced['metadata']['exact_fixed_J_substitution']
    model['metadata']['purpose']='Find a feasible hypothetical spectrum whose true-tail lower first moment exceeds 2^-40'
    lp=export(model,newscales)[0]
    assert lp.count(' obj: h_38\n')==1
    lp=lp.replace(' obj: h_38\n',' obj: 0 h_38\n')
    basis,mapping=warm_basis(old,scales,model,newscales,PREVIOUS/'h_38.bas')
    return model,newscales,lp,basis,mapping


def audit():
    model,scales,lp,basis,mapping=build()
    assert model==read(FOLDER/'model.json') and scales==read(FOLDER/'scales.json')
    assert lp==(FOLDER/'h_38.lp').read_text() and basis==(FOLDER/'warm.bas').read_text()
    assert mapping==read(FOLDER/'warm_mapping.json')
    raw=witness(FOLDER,model['metadata']['variables'])
    physical={n:x*scales[n] for n,x in raw.items()}
    assert all(x>=0 for x in physical.values())
    def verify_rows(rows,x):
        for row in rows:
            lhs=sum((fmpq(int(c))*fmpq(x[n].numerator,x[n].denominator) for n,c in row['coeffs'].items()),fmpq(0))
            rhs=fmpq(int(row['rhs']))
            assert {'eq':lhs==rhs,'ge':lhs>=rhs,'le':lhs<=rhs}[row['sense']],row['name']
    verify_rows(model['constraints'],physical)
    full,_=full_build()
    fullphysical=physical|{n:Q(v) for n,v in model['metadata']['fixed_J_substitution'].items()}
    assert set(fullphysical)==set(full['metadata']['variables'])
    verify_rows(full['constraints'],fullphysical)
    cells=read(GEN/'bch256_q1_activity_cells_outward.json')
    weight38_lower=31*physical['h_38']*decode(cells['shells'][0]['coefficient_lower'])
    assert weight38_lower>Q(1,1<<40)
    paths=[Path(__file__),ROOT/'code/bch_hull_cut_iteration.py',
           ROOT/'code/prepare_bch_hull_intersection_probe.py',ROOT/'code/prepare_bch_hull_intersection_reduced.py',
           PREVIOUS/'model.json',PREVIOUS/'scales.json',PREVIOUS/'h_38.bas',PREVIOUS/'audit.json',
           GEN/'bch256_q1_activity_cells_outward.json',FOLDER/'model.json',FOLDER/'scales.json',FOLDER/'h_38.sol']
    return dict(classification='Exact feasible hypothetical spectrum disproves sufficiency of this relaxation',
                physical_A38=encode(31*physical['h_38']),weight38_true_tail_first_moment_lower=encode(weight38_lower),
                full_intersection_rows_verified=len(full['constraints']),augmented_rows_verified=len(model['constraints']),
                all_nonnegative_and_all_rows_passed=True,relaxation_insufficient_even_with_exact_inner_tails=True,
                actual_BCH_spectrum_claimed=False,actual_failure_probability_lower_bound_claimed=False,
                original_M22_target_closed=False,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['run','audit','verify'])
    args=parser.parse_args()
    if args.mode=='run':
        model,scales,lp,basis,mapping=build()
        FOLDER.mkdir(exist_ok=False)
        write_new(FOLDER/'model.json',model)
        write_new(FOLDER/'scales.json',scales)
        fresh=[('h_38.lp',lp),('warm.bas',basis)]
        for n,t in fresh:
            text_new(FOLDER/n,t)
        write_new(FOLDER/'warm_mapping.json',mapping)
        print(solve(FOLDER,120,True),flush=True)
    else:
        result=audit()
        if args.mode=='verify':
            assert read(FOLDER/'audit.json')==result
        else:
            write_new(FOLDER/'audit.json',result)
        print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','physical_A38','weight38_true_tail_first_moment_lower')},indent=2))
