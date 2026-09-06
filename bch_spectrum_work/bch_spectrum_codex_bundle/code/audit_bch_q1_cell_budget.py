"""Independent integer lower-tail check and exact one-shell M22 budget.

This does not certify the remaining A38 cap. The integer check establishes
that OA21 alone remains insufficient even with exact one-row transfer tails.
"""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from flint import arb,ctx
from certify_bch_q1_activity_cells import region_polynomials
from audit_bch_q1_full_arb import decode,encode,rational
from run_higher_endpoint_preflight import sha,write_new

ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
OUTPUT=GEN/'bch256_q1_cell_budget.json'


def integer_lower():
    ctx.prec=256
    z,a=region_polynomials(8192,4,22,'ceil',arb)
    bits=64
    scale=1<<bits
    def integer(poly):
        result=[]
        for matrix in poly:
            row=[]
            for entries in matrix:
                vals=[]
                for x in entries:
                    lower=max(Fraction(0),rational(x.lower()))
                    vals.append((lower.numerator*scale)//lower.denominator)
                row.append(vals)
            result.append(row)
        return result
    zi,ai=integer(z),integer(a)
    # Fixed-point masses use 192 fractional bits; every product is truncated
    # downward using an integer shift, and every accumulation is exact.
    fraction_bits=192
    maximum=204
    current=[[[0]*(maximum+1) for _ in range(2)] for _ in range(39)]
    current[0][0][0]=1<<fraction_bits
    for completed in range(256):
        updated=[[[0]*(maximum+1) for _ in range(2)] for _ in range(39)]
        for w in range(min(completed,38)+1):
            for start in range(2):
                source=current[w][start]
                for end in range(2):
                    for charge in range(5):
                        width=min(completed*4+1,maximum+1-charge)
                        coeff=zi[start][end][charge]
                        if coeff:
                            dest=updated[w][end]
                            for k in range(width):
                                dest[k+charge]+=(source[k]*coeff)>>bits
                        coeff=ai[start][end][charge]
                        if coeff and w<38:
                            dest=updated[w+1][end]
                            for k in range(width):
                                dest[k+charge]+=(source[k]*coeff)>>bits
        current=updated
    # At most 204 charged cells means active time <=417792. The output
    # is stochastically bounded above by Bin(417792,1/2).
    h=204*2048
    t=arb(209717)-arb(h)/2
    binomial_lower=max(Fraction(0),rational((1-(-2*t*t/h).exp()).lower()))
    mass=Fraction(sum(sum(row) for row in current[38]),(1<<fraction_bits)*math.comb(256,38))
    coefficient=8192*mass*binomial_lower
    return dict(classification='Independent integer-truncated lower bound; no binary64 recurrence',
                cells_per_region=4,cell_length=2048,maximum_charge=maximum,
                coefficient_fraction_bits=bits,accumulator_fraction_bits=fraction_bits,
                charged_time_event_probability_lower=encode(mass),
                binomial_factor_lower=encode(binomial_lower),coefficient_lower=encode(coefficient))


def build():
    cells=json.loads((GEN/'bch256_q1_activity_cells_outward.json').read_text())
    oa=json.loads((GEN/'oa21_closure_probe/audit.json').read_text())
    transfer=json.loads((GEN/'bch256_q1_full_arb_transfer.json').read_text())
    assert [r['weight'] for r in cells['shells']]==[38,40,42]
    pair=decode(cells['shells'][0]['coefficient_upper'])+decode(transfer['coefficient_rows']['218']['coefficient_upper'])
    cap38=oa['shells'][0]['cap']
    rest=decode(cells['full_first_moment_upper'])-cap38*pair
    target=Fraction(1,1<<40)
    assert 0<rest<target
    threshold=(target-rest)/pair
    largest=(threshold.numerator-1)//threshold.denominator
    assert rest+largest*pair<target<=rest+(largest+1)*pair
    readable_cap=10_000_000_000_000
    total=rest+readable_cap*pair
    assert total<target
    check=integer_lower()
    feasible=31*decode(oa['shells'][0]['physical_optimum'])*decode(check['coefficient_lower'])
    assert feasible>target
    paths=[Path(__file__),ROOT/'code/certify_bch_q1_activity_cells.py',ROOT/'code/audit_bch_q1_full_arb.py',
           GEN/'bch256_q1_activity_cells_outward.json',GEN/'bch256_q1_activity_cells_outward.npz',
           GEN/'oa21_closure_probe/audit.json',GEN/'bch256_q1_full_arb_transfer.json']
    return dict(classification='Exact one-shell conditional M22 budget plus independent integer obstruction',
                remaining_assumption=dict(weight=38,A38_cap=readable_cap,deterministically_proved=False),
                statistical_weight40_or_42_assumptions_used=False,
                weight40_deterministic_cap=oa['shells'][1]['cap'],weight42_deterministic_cap=oa['shells'][2]['cap'],
                largest_integer_A38_cap_sufficient_for_saved_bound=largest,
                rest_without_weights38_and218=encode(rest),rest_fraction_of_target_diagnostic=float(rest/target),
                paired_weight38_coefficient_upper=encode(pair),conditional_full_upper=encode(total),
                conditional_margin_bits_diagnostic=-encode(total)['log2_diagnostic'],
                exact_conditional_full_bound_below_2_to_minus_40=True,
                integer_check=check,feasible_OA21_weight38_contribution_lower=encode(feasible),
                exact_true_tail_OA21_obstruction_verified=True,original_M22_target_closed=False,
                source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})


if __name__=='__main__':
    result=build()
    if '--verify' in sys.argv:
        assert json.loads(OUTPUT.read_text())==result
    else:
        write_new(OUTPUT,result)
    print(json.dumps({k:result[k] for k in ('remaining_assumption','largest_integer_A38_cap_sufficient_for_saved_bound',
                     'rest_fraction_of_target_diagnostic','conditional_margin_bits_diagnostic',
                     'exact_true_tail_OA21_obstruction_verified')},indent=2))
    print('INTEGER LOWER COEFFICIENT LOG2',result['integer_check']['coefficient_lower']['log2_diagnostic'])
    print('INTEGER PRIMAL OBSTRUCTION MARGIN',-result['feasible_OA21_weight38_contribution_lower']['log2_diagnostic'])
