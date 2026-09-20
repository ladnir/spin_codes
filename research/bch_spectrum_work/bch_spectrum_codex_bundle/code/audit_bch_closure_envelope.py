"""Re-derive BCH LP rows and remove all empirical orbit-search lower bounds.

The dual witnesses remain feasible when those lower bounds are weakened to
zero. The resulting slightly larger caps require no retained search orbit.
Published spectra remain explicit external mathematical inputs.
"""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from pypdf import PdfReader
sys.dont_write_bytecode = True
from anchors import L71_HALF, U187_DUAL_HALF
from bch_quotient import (generator_polynomial, binary_poly_degree,
                          binary_poly_divmod, verify_quotient_algebra)
from run_higher_endpoint_preflight import sha, write_new

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / 'generated'
WEIGHTS = list(range(0, 129, 2))


def kraw_table():
    columns = []
    for w in range(257):
        col = [1, 256-2*w]
        for j in range(1, 256):
            a, rem = divmod((256-2*w)*col[-1]-(257-j)*col[-2], j+1)
            assert rem == 0
            col.append(a)
        columns.append(col)
    return [[columns[w][j] for w in range(257)] for j in range(257)]


def audit():
    paper = ROOT / 'paper/e104-a_9_1321.pdf'
    pdf = PdfReader(paper)
    table7 = pdf.pages[5].extract_text().split('Table 7', 1)[1].split('Table 8', 1)[0]
    table10 = pdf.pages[6].extract_text().split('Table 10', 1)[1]
    def extract(text, column):
        values = {}
        for line in text.splitlines():
            fields = line.split()
            if len(fields)>1 and fields[0].isdigit() and fields[column].isdigit():
                values[int(fields[0])] = int(fields[column])
        return values
    assert extract(table7, -1) == L71_HALF
    assert extract(table10, 1) == U187_DUAL_HALF
    kt = kraw_table()
    def expand(half):
        return [half.get(min(w, 256-w), 0) for w in range(257)]
    def dual(a, dimension):
        assert sum(a) == 1 << dimension
        result = []
        for row in kt:
            value, rem = divmod(sum(x*y for x,y in zip(a,row)), 1<<dimension)
            assert rem == 0 and value >= 0
            result.append(value)
        assert sum(result) == 1 << (256-dimension)
        return result
    low, upper_dual = expand(L71_HALF), expand(U187_DUAL_HALF)
    low_dual, upper = dual(low, 71), dual(upper_dual, 69)
    assert next(w for w in range(1,257) if low_dual[w]) == 16
    quotient = verify_quotient_algebra()
    generators = {d: generator_polynomial(d) for d in (19,37,39,59)}
    assert [255-binary_poly_degree(generators[d]) for d in generators] == [187,131,123,71]
    for small, large in ((19,37),(37,39),(39,59)):
        assert binary_poly_divmod(generators[large], generators[small])[1] == 0
    assert all(binary_poly_divmod((1<<255)-1,g)[1] == 0 for g in generators.values())

    model = json.loads((GEN/'coupled_lp_exact.json').read_text())
    variables = [f'{prefix}_{w}' for prefix in ('q','h') for w in WEIGHTS]
    assert model['metadata']['variables'] == variables
    scales = {name: max(1, math.comb(256,int(name.split('_')[1]))//(1<<132)) for name in variables}
    rows = [r for r in model['constraints'] if not r['name'].startswith(('Q_OA_t','H_OA_t','q_nonneg_','h_nonneg_'))]
    def coeffs(prefix, degree):
        return {f'{prefix}_{w}': kt[degree][w]+(kt[degree][256-w] if w!=128 else 0) for w in WEIGHTS}
    for label, prefix in (('Q','q'),('H','h')):
        for degree in range(0,16,2):
            rows.append(dict(name=f'{label}_OA_Kraw_{degree}',coeffs=coeffs(prefix,degree),
                             sense='eq',rhs=(1<<123) if degree==0 else 0))
    checked_rows = []
    for row in rows:
        name = row['name']
        actual = {k:int(v) for k,v in row['coeffs'].items() if int(v)}
        rhs = int(row['rhs'])
        weakened = False
        w = int(name.split('_')[-1]) if name.split('_')[-1].isdigit() else None
        if name in ('q_0','h_0'):
            expected, sense, value = {name:1}, 'eq', int(name=='q_0')
        elif name.startswith(('q_support_','h_support_')):
            prefix = name[0]
            assert w%2==0 and (0<w<40 if prefix=='q' else 0<=w<38)
            expected, sense, value = {f'{prefix}_{w}':1}, 'eq', 0
        elif name.startswith(('Wambach_','orbit_search_')):
            part = name.removeprefix('Wambach_').removeprefix('orbit_search_').removesuffix('_lower')
            assert part in variables and rhs>=0
            expected, sense, value = {part:1}, 'ge', rhs
            weakened = True
        elif name.startswith('Q_ge_L_'):
            expected, sense, value = {f'q_{w}':1}, 'ge', low[w]
        elif name.startswith('P_le_U_'):
            expected, sense, value = {f'q_{w}':1,f'h_{w}':255}, 'le', upper[w]
        elif name.startswith('Bp_ge_Udual_'):
            expected = coeffs('q',w) | {n:255*v for n,v in coeffs('h',w).items()}
            sense, value = 'ge', upper_dual[w]*(1<<131)
        elif name.startswith('Bq_le_Ldual_'):
            expected, sense, value = coeffs('q',w), 'le', low_dual[w]*(1<<123)
        elif name.startswith('Bq_ge_Bp_'):
            expected = coeffs('q',w) | {n:-v for n,v in coeffs('h',w).items()}
            sense, value = 'ge', 0
        elif name.startswith(('Q_OA_Kraw_','H_OA_Kraw_')):
            expected, sense, value = coeffs('q' if name[0]=='Q' else 'h',w), 'eq', (1<<123) if w==0 else 0
        else:
            raise AssertionError(name)
        expected = {k:v for k,v in expected.items() if v}
        assert actual == expected and row['sense']==sense and rhs==value, name
        normalization = max(1,abs(rhs),*(abs(v*scales[n]) for n,v in actual.items()))
        checked_rows.append(dict(name=name,coeffs=actual,sense=sense,
                                 original_rhs=rhs,rhs=0 if weakened else rhs,
                                 normalization=normalization,weakened=weakened))
    results = []
    paths = [Path(__file__),paper,GEN/'coupled_lp_exact.json',GEN/'exact_bch_lp_caps.json',
             ROOT/'code/anchors.py',ROOT/'code/bch_quotient.py',ROOT/'code/affine_wambach.py']
    old_caps = {r['weight']:r['exact_Aw_C_cap'] for r in json.loads((GEN/'exact_bch_lp_caps.json').read_text())['shells']}
    for w in range(38,52,2):
        objective = 'h_38' if w==38 else f'c_{w}'
        path = GEN/f'max_{objective}_scaled_rational.sol'
        prices = {}
        for line in path.read_text().split('PI:',1)[1].split('SLACK:',1)[0].splitlines():
            if ' = ' in line:
                name, val = line.split(' = ',1)
                prices[name.strip()] = Fraction(val.strip())
        assert set(prices) <= {f'c{i}' for i in range(1,len(rows)+1)}
        combined = {name:Fraction(0) for name in variables}
        upper_bound = Fraction(0)
        for i,row in enumerate(checked_rows,1):
            price = prices.get(f'c{i}',Fraction(0))
            assert row['sense']=='eq' or (price>=0 if row['sense']=='le' else price<=0)
            physical = price*scales[f'q_{w}']/row['normalization']
            for name, value in row['coeffs'].items():
                combined[name] += physical*value
            upper_bound += physical*row['rhs']
        target = {f'h_{w}':1} if w==38 else {f'q_{w}':1,f'h_{w}':31}
        assert all(combined[n]>=target.get(n,0) for n in variables)
        cap = upper_bound.numerator//upper_bound.denominator
        if w==38:
            cap *= 31  # q_38=0. No orbit/divisibility rounding is used.
        assert cap <= old_caps[w]*Fraction((1<<20)+1,1<<20)
        results.append(dict(weight=w,cap=cap,old_cap=old_caps[w],
                            ratio_diagnostic=cap/old_caps[w]))
        paths.append(path)
    return dict(classification='Exact BCH model and dual-certificate audit with all orbit-search lower bounds replaced by zero',
        anchor_tables_checked=[7,10],containment_dimensions=[71,123,131,187],
        quotient=quotient,OA_strength=15,linear_rows_rederived=len(rows),
        orbit_lower_bounds_used=False,statistical_caps_used=False,
        lattice_rounding_used=False,shells=results,
        each_new_cap_at_most_old_times_one_plus_2_to_minus_20=True,
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


def main():
    value = audit()
    output = GEN/'bch256_closure_deterministic_envelope.json'
    if '--verify' in sys.argv:
        assert json.loads(output.read_text())==value
    else:
        write_new(output,value)
    print(json.dumps({k:v for k,v in value.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':
    main()
