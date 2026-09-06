"""Sequential exact Benders iterations with verified cuts and mapped warm bases.

All outputs are new files. Replay uses saved rational witnesses, not an optimizer.
"""
import argparse
import copy
import json
import math
import re
import subprocess
import sys
import time
from fractions import Fraction as Q
from pathlib import Path
from flint import fmpq
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from prepare_bch_hull_intersection_probe import build as full_build
from prepare_bch_hull_probe import export
from export_scaled_rational_lp import expression
from verify_scaled_rational_solution import normalized_rows,parse_assignments
from audit_bch_q1_full_arb import encode,decode
from run_higher_endpoint_preflight import write_new,sha
ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'


def read(path):
    return json.loads(path.read_text())


def fresh_text(path,text):
    with path.open('x',encoding='ascii') as f:
        f.write(text)


def fq(v):
    v=Q(v)
    return fmpq(v.numerator,v.denominator)


def witness(folder,variables):
    text=(folder/'h_38.sol').read_text()
    assert 'status = OPTIMAL' in text
    values=parse_assignments(text,'VARS:','REDUCED COST:')
    assert set(values)<=set(variables)
    return {n:values.get(n,Q(0)) for n in variables}


def check_master(model,scales,folder):
    variables=model['metadata']['variables']
    raw=witness(folder,variables)
    x={n:fq(v) for n,v in raw.items()}
    assert all(v>=0 for v in x.values())
    rows=normalized_rows(model,scales)
    prices=parse_assignments((folder/'h_38.sol').read_text(),'PI:','SLACK:')
    assert set(prices)<={f'c{i}' for i in range(1,len(rows)+1)}
    combined={n:fmpq(0) for n in variables}
    objective=fmpq(0)
    for i,row in enumerate(rows,1):
        coefficients={n:fq(v) for n,v in row['coeffs'].items()}
        lhs=sum((v*x[n] for n,v in coefficients.items()),fmpq(0))
        rhs=fq(row['rhs'])
        assert {'eq':lhs==rhs,'ge':lhs>=rhs,'le':lhs<=rhs}[row['sense']],row['name']
        price=fq(prices.get(f'c{i}',Q(0)))
        assert row['sense']=='eq' or (price>=0 if row['sense']=='le' else price<=0)
        objective+=price*rhs
        for n,v in coefficients.items():
            combined[n]+=price*v
    assert all(combined[n]>=int(n=='h_38') for n in variables)
    assert objective==x['h_38']
    return raw,dict(rows_checked=len(rows),variables_checked=len(variables),
                    all_primal_dual_checks_passed=True)


def extension_model(full,scales,oldmodel,oldscales,oldraw):
    fixed={n:oldraw[n]*oldscales[n] for n in oldmodel['metadata']['variables']}
    for row in full['constraints']:
        if row['name'].startswith('J_support_'):
            n,c=next(iter(row['coeffs'].items()))
            assert row['sense']=='eq' and int(c)==1
            fixed[n]=Q(int(row['rhs']))
    variables=[n for n in full['metadata']['variables'] if n not in fixed]
    constraints=[]
    for row in full['constraints']:
        rhs=Q(int(row['rhs']))-sum((int(c)*fixed[n] for n,c in row['coeffs'].items() if n in fixed),Q(0))
        coefficients={n:int(c)*scales[n] for n,c in row['coeffs'].items() if n not in fixed and int(c)}
        if not coefficients:
            assert {'eq':rhs==0,'ge':rhs<=0,'le':rhs>=0}[row['sense']],row['name']
            continue
        norm=max(Q(1),abs(rhs),*(abs(c) for c in coefficients.values()))
        constraints.append(dict(name=row['name'],coeffs={n:str(Q(c)/norm) for n,c in coefficients.items()},
                                sense=row['sense'],rhs=str(rhs/norm),normalization=str(norm)))
    return dict(variables=variables,scales={n:scales[n] for n in variables},
                fixed_physical={n:str(v) for n,v in fixed.items()},constraints=constraints,
                objective={'t_64':'1'})


def auxiliary_lp(model):
    def wrapped(coeffs):
        pieces=[expression([(Q(v),n)]) for n,v in coeffs.items() if Q(v)]
        return [' '.join(pieces[j:j+4]) for j in range(0,len(pieces),4)] or ['0 '+model['variables'][0]]
    chunks=wrapped(model['objective'])
    chunks[0]=' obj: '+chunks[0]
    lines=['Maximize']+chunks+['Subject To']
    for i,row in enumerate(model['constraints'],1):
        chunks=wrapped(row['coeffs'])
        chunks[0]=f' c{i}: '+chunks[0]
        chunks[-1]+=' '+{'eq':'=','ge':'>=','le':'<='}[row['sense']]+' '+row['rhs']
        lines+=chunks
    return '\n'.join(lines+['Bounds']+[f' 0 <= {n}' for n in model['variables']]+['End'])+'\n'


def farkas_model(extension):
    terms=[]
    for i,row in enumerate(extension['constraints']):
        for sign in ((1,-1) if row['sense']=='eq' else ((1,) if row['sense']=='ge' else (-1,))):
            terms.append(dict(variable=f'f_{len(terms)}',row_index=i,sign=sign))
    variables=[t['variable'] for t in terms]
    objective={t['variable']:str(t['sign']*Q(extension['constraints'][t['row_index']]['rhs'])) for t in terms}
    constraints=[]
    for n in extension['variables']:
        constraints.append(dict(name=n,coeffs={t['variable']:str(t['sign']*Q(extension['constraints'][t['row_index']]['coeffs'].get(n,'0'))) for t in terms},sense='le',rhs='0'))
    constraints.append(dict(name='mass',coeffs={n:'1' for n in variables},sense='le',rhs='1'))
    return dict(variables=variables,objective=objective,constraints=constraints,terms=terms)


def derive_cut(full,oldmodel,extension,farkas,y,number):
    assert all(v>=0 for v in y.values()) and sum(y.values())<=1
    combined={n:Q(0) for n in extension['variables']}
    contradiction=Q(0)
    original={r['name']:r for r in full['constraints']}
    fixed={n:Q(v) for n,v in extension['fixed_physical'].items()}
    kept=set(oldmodel['metadata']['variables'])
    cut={n:Q(0) for n in kept}
    physical_residual={n:Q(0) for n in extension['variables']}
    cut_rhs=Q(0)
    used=[]
    for term in farkas['terms']:
        signed=term['sign']*y[term['variable']]
        if not signed:
            continue
        row=extension['constraints'][term['row_index']]
        contradiction+=signed*Q(row['rhs'])
        for n,v in row['coeffs'].items():
            combined[n]+=signed*Q(v)
        source=original[row['name']]
        mult=signed/Q(row['normalization'])
        rhs=Q(int(source['rhs']))-sum((int(c)*fixed[n] for n,c in source['coeffs'].items() if n in fixed and n not in kept),Q(0))
        cut_rhs+=mult*rhs
        for n,c in source['coeffs'].items():
            if n in kept:
                cut[n]+=mult*int(c)
            elif n in physical_residual:
                physical_residual[n]+=mult*int(c)
        used.append(dict(row=row['name'],physical_multiplier=str(mult)))
    assert contradiction>0 and all(v<=0 for v in combined.values())
    assert all(v<=0 for v in physical_residual.values())
    assert cut_rhs-sum((c*fixed[n] for n,c in cut.items()),Q(0))==contradiction
    denominator=math.lcm(cut_rhs.denominator,*(v.denominator for v in cut.values()))
    integer={n:int(v*denominator) for n,v in cut.items() if v}
    rhs=int(cut_rhs*denominator)
    divisor=math.gcd(abs(rhs),*(abs(v) for v in integer.values()))
    row=dict(name=f'J_Farkas_iteration_{number}',coeffs={n:str(v//divisor) for n,v in sorted(integer.items())},sense='ge',rhs=str(rhs//divisor))
    proof=dict(contradiction=str(contradiction),used_original_rows=used,
               eliminated_physical_coefficients={n:str(v) for n,v in physical_residual.items()})
    return row,proof


def check_extension(full,extension,x):
    assert all(v>=0 for v in x.values())
    physical={n:Q(v) for n,v in extension['fixed_physical'].items()}
    physical.update({n:v*extension['scales'][n] for n,v in x.items()})
    assert set(physical)==set(full['metadata']['variables'])
    for row in full['constraints']:
        lhs=sum((fq(c)*fq(physical[n]) for n,c in row['coeffs'].items()),fmpq(0))
        rhs=fmpq(row['rhs'])
        assert {'eq':lhs==rhs,'ge':lhs>=rhs,'le':lhs<=rhs}[row['sense']],row['name']
    return {n:str(v) for n,v in physical.items()}


def warm_basis(oldmodel,oldscales,newmodel,newscales,source):
    oldrows=normalized_rows(oldmodel,oldscales)
    newrows=normalized_rows(newmodel,newscales)
    indices={row['name']:i for i,row in enumerate(newrows,1)}
    assert len(indices)==len(newrows)
    mapping={}
    for i,row in enumerate(oldrows,1):
        j=indices[row['name']]
        assert newrows[j-1]==row
        mapping[f'c{i}']=f'c{j}'
    text=source.read_text()
    def replace(match):
        return mapping[match.group(0)]
    mapped=re.sub(r'\bc\d+\b',replace,text)
    return mapped,dict(old_rows=len(oldrows),new_rows=len(newrows),
                       unchanged_row_coefficients_verified=True,row_name_mapping=mapping)


def solve(folder,seconds,warm=False):
    assert not (folder/'attempt.json').exists()
    cwd='/mnt/c/'+str(folder.resolve())[3:].replace('\\','/')
    runtime=GEN/'oa21_closure_probe/tools/runtime'
    wruntime='/mnt/c/'+str(runtime.resolve())[3:].replace('\\','/')
    command=['wsl.exe','-d','Ubuntu-24.04','--cd',cwd,'--','env',
        'LD_LIBRARY_PATH='+wruntime+'/usr/lib/x86_64-linux-gnu',wruntime+'/usr/bin/esolver',
        '-S','-d','7','-P','256','-L','-R',str(seconds),'-b','h_38.bas','-O','h_38.sol']
    if warm:
        command+=['-B','warm.bas']
    command+=['h_38.lp']
    start=time.monotonic()
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    with (folder/'esolver.log').open('xb') as f:
        f.write(result.stdout)
    value=dict(command=command,elapsed_seconds=time.monotonic()-start,exit_code=result.returncode,
                solution_exists=(folder/'h_38.sol').exists(),warm_start=warm)
    write_new(folder/'attempt.json',value)
    print(json.dumps(dict(stage=folder.name,**{k:v for k,v in value.items() if k!='command'})),flush=True)
    return (folder/'h_38.sol').read_text().splitlines()[0] if (folder/'h_38.sol').exists() else 'NO_SOLUTION'


def prepare_aux(folder,model):
    folder.mkdir()
    write_new(folder/'model.json',model)
    fresh_text(folder/'h_38.lp',auxiliary_lp(model))


def aggregate(previous,physical,checks):
    old=read(previous/'audit.json')
    cap=31*(physical.numerator//physical.denominator)
    assert cap<=old['A38_cap_without_lattice_rounding']
    cells=read(GEN/'bch256_q1_activity_cells_outward.json')
    transfer=read(GEN/'bch256_q1_full_arb_transfer.json')
    pair=decode(cells['shells'][0]['coefficient_upper'])+decode(transfer['coefficient_rows']['218']['coefficient_upper'])
    full=decode(old['full_M22_first_moment_upper'])-(old['A38_cap_without_lattice_rounding']-cap)*pair
    lower=31*physical*decode(cells['shells'][0]['coefficient_lower'])
    return dict(**checks,A38_cap_without_lattice_rounding=cap,previous_A38_cap=old['A38_cap_without_lattice_rounding'],
        h38_optimum=encode(physical),full_M22_first_moment_upper=encode(full),
        full_M22_margin_bits_diagnostic=-encode(full)['log2_diagnostic'],
        feasible_primal_true_tail_lower=encode(lower),
        relaxation_insufficient_even_with_exact_inner_tails=lower>Q(1,1<<40),
        original_M22_target_closed=full<Q(1,1<<40))


def execute(previous,folder,number):
    folder.mkdir(parents=True,exist_ok=False)
    write_new(folder/'inputs.json',dict(previous=str(previous.relative_to(GEN)),number=number,source_sha256=sha(Path(__file__))))
    full,fullscales=full_build()
    assert full==read(GEN/'oa21_hull_intersection_probe/model.json')
    oldmodel,oldscales=read(previous/'model.json'),read(previous/'scales.json')
    oldraw,oldchecks=check_master(oldmodel,oldscales,previous)
    extension=extension_model(full,fullscales,oldmodel,oldscales,oldraw)
    ef=folder/'extension'
    prepare_aux(ef,extension)
    status=solve(ef,45)
    if status=='status = OPTIMAL':
        x=witness(ef,extension['variables'])
        physical=check_extension(full,extension,x)
        write_new(folder/'exhaustion.json',dict(full_intersection_feasible=True,
            physical_spectrum=physical,old_objective_unchanged=True,original_M22_target_closed=False))
        print('EXHAUSTED: exact feasible intersection extension',flush=True)
        return
    if status!='status = INFEASIBLE':
        print('STOP: extension has no decisive saved status',flush=True)
        return
    farkas=farkas_model(extension)
    ff=folder/'farkas'
    prepare_aux(ff,farkas)
    if solve(ff,45)!='status = OPTIMAL':
        print('STOP: no Farkas witness',flush=True)
        return
    y=witness(ff,farkas['variables'])
    cut,proof=derive_cut(full,oldmodel,extension,farkas,y,number)
    write_new(folder/'cut.json',dict(row=cut,proof=proof))
    model=copy.deepcopy(oldmodel)
    model['constraints'].append(cut)
    model['metadata'].setdefault('iterated_J_cuts',[]).append(dict(number=number,proof_source=str((folder/'cut.json').relative_to(GEN)),sha256=sha(folder/'cut.json')))
    mf=folder/'master'
    mf.mkdir()
    write_new(mf/'model.json',model)
    write_new(mf/'scales.json',oldscales)
    fresh_text(mf/'h_38.lp',export(model,oldscales)[0])
    basis,mapping=warm_basis(oldmodel,oldscales,model,oldscales,previous/'h_38.bas')
    fresh_text(mf/'warm.bas',basis)
    write_new(mf/'warm_mapping.json',mapping)
    if solve(mf,90,True)!='status = OPTIMAL':
        print('STOP: no exact master optimum',flush=True)
        return
    raw,checks=check_master(model,oldscales,mf)
    audit=aggregate(previous,raw['h_38']*oldscales['h_38'],checks)
    audit['source_sha256']={str(p.relative_to(ROOT)):sha(p) for p in
        (Path(__file__),folder/'cut.json',mf/'model.json',mf/'scales.json',mf/'h_38.sol',previous/'audit.json')}
    write_new(mf/'audit.json',audit)
    print(json.dumps({k:audit[k] for k in ('A38_cap_without_lattice_rounding','full_M22_margin_bits_diagnostic','original_M22_target_closed')}),flush=True)


def replay(folder):
    inputs=read(folder/'inputs.json')
    assert inputs['source_sha256']==sha(Path(__file__))
    previous=GEN/inputs['previous']
    full,fullscales=full_build()
    assert full==read(GEN/'oa21_hull_intersection_probe/model.json')
    oldmodel,oldscales=read(previous/'model.json'),read(previous/'scales.json')
    raw,checks=check_master(oldmodel,oldscales,previous)
    extension=extension_model(full,fullscales,oldmodel,oldscales,raw)
    ef=folder/'extension'
    assert extension==read(ef/'model.json') and auxiliary_lp(extension)==(ef/'h_38.lp').read_text()
    if (folder/'exhaustion.json').exists():
        x=witness(ef,extension['variables'])
        assert read(folder/'exhaustion.json')==dict(full_intersection_feasible=True,
            physical_spectrum=check_extension(full,extension,x),old_objective_unchanged=True,original_M22_target_closed=False)
        print('PASS: exact feasible extension; intersection relaxation exhausted')
        return
    farkas=farkas_model(extension)
    ff=folder/'farkas'
    assert farkas==read(ff/'model.json') and auxiliary_lp(farkas)==(ff/'h_38.lp').read_text()
    y=witness(ff,farkas['variables'])
    cut,proof=derive_cut(full,oldmodel,extension,farkas,y,inputs['number'])
    assert read(folder/'cut.json')==dict(row=cut,proof=proof)
    model=copy.deepcopy(oldmodel)
    model['constraints'].append(cut)
    model['metadata'].setdefault('iterated_J_cuts',[]).append(dict(number=inputs['number'],proof_source=str((folder/'cut.json').relative_to(GEN)),sha256=sha(folder/'cut.json')))
    mf=folder/'master'
    assert model==read(mf/'model.json') and oldscales==read(mf/'scales.json')
    assert export(model,oldscales)[0]==(mf/'h_38.lp').read_text()
    basis,mapping=warm_basis(oldmodel,oldscales,model,oldscales,previous/'h_38.bas')
    assert basis==(mf/'warm.bas').read_text() and mapping==read(mf/'warm_mapping.json')
    raw,checks=check_master(model,oldscales,mf)
    expected=aggregate(previous,raw['h_38']*oldscales['h_38'],checks)
    actual=read(mf/'audit.json')
    assert all(actual[k]==v for k,v in expected.items())
    assert all(sha(ROOT/p)==v for p,v in actual['source_sha256'].items())
    print('PASS: exact cut, mapped basis, primal/dual and full M22 aggregation')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['run','verify'])
    parser.add_argument('number',type=int)
    parser.add_argument('--previous',default='oa21_hull_benders_cut')
    args=parser.parse_args()
    assert 2<=args.number<=1000
    folder=GEN/'hull_cut_iterations'/f'round_{args.number:03d}'
    previous=(GEN/args.previous).resolve()
    assert previous.is_relative_to(GEN.resolve())
    if args.mode=='run':
        execute(previous,folder,args.number)
    else:
        replay(folder)
