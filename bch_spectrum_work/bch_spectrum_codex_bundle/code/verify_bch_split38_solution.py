"""Exact primal/dual check for the code-specific split LP, without solver calls."""
import json
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from flint import fmpq
from audit_bch_q1_full_arb import encode
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'generated/split38_local_probe'


def assignments(text,start,end):
    result={}
    for line in text.split(start,1)[1].split(end,1)[0].splitlines():
        if ' = ' in line:
            name,value=line.split(' = ',1)
            result[name.strip()]=fmpq(value.strip())
    return result


def verify():
    model=json.loads((FOLDER/'model.json').read_text())
    solution=FOLDER/'h38_512.sol'
    text=solution.read_text()
    assert 'status = OPTIMAL' in text
    variables=model['variables']
    scales=model['scales']
    rows=model['constraints']
    raw=assignments(text,'VARS:','REDUCED COST:')
    prices=assignments(text,'PI:','SLACK:')
    assert set(raw)<=set(variables)
    assert set(prices)<={f'c{i}' for i in range(1,len(rows)+1)}
    x={name:raw.get(name,fmpq(0))*scales[name] for name in variables}
    assert all(value>=0 for value in x.values())
    objective_scale=model['objective']['physical_multiplier']
    combined={name:fmpq(0) for name in variables}
    dual=fmpq(0)
    for i,row in enumerate(rows,1):
        co={n:int(v) for n,v in row['coeffs'].items()}
        rhs=int(row['rhs'])
        actual=sum((x[n]*v for n,v in co.items()),fmpq(0))
        assert {'eq':actual==rhs,'le':actual<=rhs,'ge':actual>=rhs}[row['sense']],row['name']
        norm=max(1,abs(rhs),*(abs(v*scales[n]) for n,v in co.items()))
        price=prices.get(f'c{i}',fmpq(0))
        assert row['sense']=='eq' or (price>=0 if row['sense']=='le' else price<=0)
        factor=price*objective_scale/norm
        dual+=factor*rhs
        for name,value in co.items():
            combined[name]+=factor*value
    assert all(combined[n]>=int(n=='h_38') for n in variables)
    assert dual==x['h_38']
    omitted_failures=[]
    for row in json.loads((FOLDER/'all_rows_after_substitution.json').read_text()):
        actual=sum((x[n]*int(v) for n,v in row['coeffs'].items()),fmpq(0))
        rhs=int(row['rhs'])
        if not {'eq':actual==rhs,'le':actual<=rhs,'ge':actual>=rhs}[row['sense']]:
            omitted_failures.append(row['name'])
    integer=int(dual.numerator)//int(dual.denominator)
    paths=[Path(__file__),FOLDER/'model.json',FOLDER/'all_rows_after_substitution.json',FOLDER/'h38.lp',solution,
           ROOT/'generated/bch256_wambach_shortening.json']
    return dict(classification='Exact split-LP primal and dual checks; model applicability is stated separately',
                rows_checked=len(rows),variables_checked=len(variables),all_primal_and_dual_checks_passed=True,
                h38_optimum=encode(Fraction(int(dual.numerator),int(dual.denominator))),
                A38_cap_without_lattice_rounding=31*integer,omitted_row_primal_failures=omitted_failures,
                A38_cap_at_most_10_to_13=31*integer<=10**13,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__=='__main__':
    result=verify()
    output=FOLDER/'certificate.json'
    if '--verify' in sys.argv:
        assert json.loads(output.read_text())==result
    else:
        write_new(output,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','h38_optimum')},indent=2))
