"""Verify saved local primal/dual witnesses without invoking an LP optimizer."""
import json
from fractions import Fraction
from pathlib import Path
from certify_bch_local_slack import residuals
from prepare_bch_split38_probe import kraw
from audit_bch_q1_full_arb import decode
from run_higher_endpoint_preflight import sha
import math
import sys
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1]
receipt=json.loads((ROOT/'generated/bch256_local_slack_certificate.json').read_text())
assert all(sha(ROOT/path)==value for path,value in receipt['source_sha256'].items())
r,upper=residuals()
kt=kraw(218)
total=Fraction(0)
for block in receipt['blocks']:
    i=block['inside_weight']
    nodes=[j for j in range(219) if (i+j)%2==0 and (i+j in (0,256) or 40<=i+j<=216) and 38<=38-i+j<=218]
    assert nodes==block['outside_nodes']
    basis=[[Fraction(kt[j][t],math.comb(218,t)) for j in nodes] for t in range(7)]
    costs=[]
    for j in nodes:
        a,b=i+j,38-i+j
        costs.append(r[f'q_{min(a,256-a)}']/(1 if a==128 else 2)+r[f'h_{min(b,256-b)}']/(1 if b==128 else 2))
    polynomial=[decode(v) for v in block['normalized_kraw_polynomial_coefficients']]
    assert all(sum((polynomial[t]*basis[t][j] for t in range(7)),Fraction(0))<=costs[j] for j in range(len(nodes)))
    witness=block['exact_zero_cost_primal']
    assert witness is not None
    masses=[(v['node_index'],decode(v['mass'])) for v in witness]
    assert all(x>=0 and costs[j]==0 for j,x in masses)
    assert all(sum((x*basis[t][j] for j,x in masses),Fraction(0))==int(t==0) for t in range(7))
    assert polynomial[0]==0
    assert decode(block['certified_slack_lower'])==0
    total+=decode(block['certified_slack_lower'])
assert [b['inside_weight'] for b in receipt['blocks']]==list(range(20))
assert total==decode(receipt['certified_total_slack'])==0
assert upper==decode(receipt['old_h38_upper'])==decode(receipt['new_h38_upper'])
print('PASS: all 20 exact local zero-cost primal witnesses and polynomial inequalities; no optimizer invoked')
