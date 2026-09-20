"""Read-only exact LP checks and deterministic M22 aggregation for OA21.

The published dual-distance bound is an external mathematical input. The
extension argument is in BCH_OA21_REFINEMENT.md, not inferred from solver output.
No orbit-search lower bounds, sampling caps, or lattice rounding are used.
"""
from __future__ import annotations
import copy
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
from audit_bch_closure_envelope import audit, kraw_table
from audit_bch_q1_full_arb import decode, encode
from certify_bch_m25_closure import deterministic_caps, old_higher_occupations
from export_scaled_rational_lp import variable_scales
from verify_scaled_rational_solution import normalized_rows, parse_assignments
from run_higher_endpoint_preflight import sha, write_new

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT/'generated'
FOLDER = GEN/'oa21_closure_probe'
OUTPUT = FOLDER/'audit.json'


def build():
    baseline = audit()
    assert baseline == json.loads((GEN/'bch256_closure_deterministic_envelope.json').read_text())
    old = json.loads((GEN/'coupled_lp_exact.json').read_text())
    model = json.loads((FOLDER/'model.json').read_text())
    expected = copy.deepcopy(old['constraints'])
    for row in expected:
        if row['name'].startswith(('Wambach_', 'orbit_search_')):
            row['rhs'] = '0'
    kt = kraw_table()
    for prefix in ('q', 'h'):
        for degree in (16, 18, 20):
            coeffs = {f'{prefix}_{w}': str(kt[degree][w]+(kt[degree][256-w] if w != 128 else 0))
                      for w in range(0,129,2)}
            expected.append(dict(name=f'provisional_OA21_{prefix}_{degree}',
                                 coeffs={n:v for n,v in coeffs.items() if int(v)}, sense='eq', rhs='0'))
    assert model['constraints'] == expected
    assert model['metadata']['variables'] == old['metadata']['variables']
    scales = variable_scales(model)
    rows = normalized_rows(model, scales)
    variables = model['metadata']['variables']
    caps = deterministic_caps()
    transfer = json.loads((GEN/'bch256_q1_full_arb_transfer.json').read_text())
    coefficients = {int(w):decode(r['coefficient_upper']) for w,r in transfer['coefficient_rows'].items()}
    results = []
    paths = [Path(__file__), ROOT/'code/prepare_bch_oa21_probe.py',
             ROOT/'code/verify_scaled_rational_solution.py', ROOT/'code/export_scaled_rational_lp.py',
             ROOT/'code/export_lp.py', ROOT/'code/certify_bch_m25_closure.py',
             ROOT/'code/audit_bch_closure_envelope.py', ROOT/'paper/augot_levy_dual_bch_1996.pdf',
             FOLDER/'model.json', GEN/'coupled_lp_exact.json', GEN/'bch256_q1_full_arb_transfer.json',
             GEN/'bch256_closure_deterministic_envelope.json', GEN/'johnson_n256_w52_d38.json',
             GEN/'johnson_n256_w54_d38.json', GEN/'bch256_q2_positive_outward.json',
             GEN/'bch256_q3_capped_outward.json', GEN/'bch256_occupation_tail_outward.json']
    for w in (38,40,42):
        label = 'h_38' if w == 38 else f'c_{w}'
        solution = FOLDER/('h_38_dual.sol' if w == 38 else label+'.sol')
        metadata = json.loads((FOLDER/(label+'.json')).read_text())
        assert metadata['objective_label'] == label
        objective = {f'h_{w}':1} if w == 38 else {f'q_{w}':1,f'h_{w}':31}
        assert metadata['objective_scaled_coefficients'] == objective
        assert int(metadata['objective_physical_multiplier']) == scales[f'q_{w}']
        text = solution.read_text()
        assert 'status = OPTIMAL' in text
        raw = parse_assignments(text,'VARS:','REDUCED COST:')
        x = {n:raw.get(n,Fraction(0)) for n in variables}
        prices = parse_assignments(text,'PI:','SLACK:')
        assert set(raw) <= set(variables)
        assert set(prices) <= {f'c{i}' for i in range(1,len(rows)+1)}
        assert all(v >= 0 for v in x.values())
        combined = {n:Fraction(0) for n in variables}
        dual = Fraction(0)
        for i,row in enumerate(rows,1):
            value = sum((v*x[n] for n,v in row['coeffs'].items()),Fraction(0))
            assert {'eq':value==row['rhs'],'le':value<=row['rhs'],'ge':value>=row['rhs']}[row['sense']]
            price = prices.get(f'c{i}',Fraction(0))
            assert row['sense']=='eq' or (price >= 0 if row['sense']=='le' else price <= 0)
            dual += price*row['rhs']
            for n,v in row['coeffs'].items():
                combined[n] += price*v
        assert all(combined[n] >= objective.get(n,0) for n in variables)
        primal = sum((v*x[n] for n,v in objective.items()),Fraction(0))
        assert primal == dual
        physical = primal*scales[f'q_{w}']
        integer = physical.numerator//physical.denominator
        cap = integer*(31 if w == 38 else 1)
        prior = caps[w]
        assert cap < prior
        caps[w] = caps[256-w] = cap
        # This is a feasible LP pseudoenumerator, not an actual code spectrum.
        pseudo_q1 = sum((coefficients[v]*(x[f'q_{min(v,256-v)}']*scales[f'q_{min(v,256-v)}']
                        +31*x[f'h_{min(v,256-v)}']*scales[f'h_{min(v,256-v)}'])
                        for v in coefficients),Fraction(0))
        results.append(dict(weight=w,cap=cap,prior_cap=prior,cap_log2_diagnostic=math.log2(cap),
                            improvement_bits_diagnostic=math.log2(prior)-math.log2(cap),
                            physical_optimum=encode(physical),all_primal_dual_checks_passed=True,
                            feasible_pseudoenumerator_Q1=encode(pseudo_q1),
                            fixed_transfer_LP_relaxation_cannot_certify_target=pseudo_q1>Fraction(1,1<<40)))
        paths.extend([solution,FOLDER/(label+'.lp'),FOLDER/(label+'.json')])
    q1 = sum((coefficients[w]*caps[w] for w in coefficients),Fraction(0))
    q2,q3,tail = old_higher_occupations()
    full = q1+Fraction(128,127)*(q2+q3+tail)
    return dict(classification='Exact OA21 LP refinement using a published dual-distance bound; M22 remains open',
                rows_verified=len(rows),variables_verified=len(variables),OA_strength=21,
                external_mathematical_input='Augot and Levy-dit-Vehel, author manuscript p.10: binary length 255, designed distance 39, dual minimum distance at least 22',
                published_rank_certificate_independently_reproduced=False,
                extension_argument='BCH_OA21_REFINEMENT.md',
                orbit_lower_bounds_used=False,statistical_shell_caps_used=False,lattice_rounding_used=False,
                shells=results,Q1_upper=encode(q1),full_first_moment_upper=encode(full),
                full_margin_bits_diagnostic=-encode(full)['log2_diagnostic'],
                exact_full_bound_below_2_to_minus_40=full<Fraction(1,1<<40),
                original_M22_target_closed=False,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__ == '__main__':
    result = build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text()) == result
    else:
        write_new(OUTPUT,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','Q1_upper','full_first_moment_upper','shells')},indent=2))
    for row in result['shells']:
        print(json.dumps({k:v for k,v in row.items() if k not in ('physical_optimum','feasible_pseudoenumerator_Q1')},indent=2))
