"""Independent Arb count recurrence and direct-power checks of tail bounds."""
from __future__ import annotations
import json
import math
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from flint import arb,ctx
from audit_bch_q1_full_arb import decode,rational
from run_higher_endpoint_preflight import sha,write_new
ROOT=Path(__file__).resolve().parents[1]
GEN=ROOT/'generated'
ctx.prec=256


def aa(x):
    return arb(x.numerator)/x.denominator if isinstance(x,Fraction) else arb(x)


def multiply(x,y):
    return tuple(sum((x[2*i+k]*y[2*k+j] for k in (0,1)),arb(0)) for i in (0,1) for j in (0,1))


def power(matrix,n):
    result=(arb(1),arb(0),arb(0),arb(1))
    for bit in bin(n)[2:]:
        result=multiply(result,result)
        if bit=='1':
            result=multiply(result,matrix)
    return result


def count_regions(zero,candidate,length,maximum):
    coefficients=[(arb(0),)*4 for _ in range(maximum+1)]
    coefficients[0]=(arb(1),arb(0),arb(0),arb(1))
    for n in range(1,length+1):
        for j in range(min(n,maximum),0,-1):
            first=multiply(coefficients[j],zero)
            second=multiply(coefficients[j-1],candidate)
            coefficients[j]=tuple(x+y for x,y in zip(first,second))
        coefficients[0]=multiply(coefficients[0],zero)
    return [tuple(x/math.comb(length,j) for x in row) for j,row in enumerate(coefficients)]


def main():
    output=GEN/'bch256_tail_transfer_crosschecks.json'
    assert not output.exists()
    source=GEN/'bch256_occupation_tail_outward.json'
    cert=json.loads(source.read_text())
    for name,expected in cert['source_sha256'].items():
        assert sha(ROOT/name)==expected
    envelope_path=GEN/'bch256_tail_moment_envelopes.json'
    envelopes=json.loads(envelope_path.read_text())
    factors={r['rho']:int(r['factor']) for r in envelopes['references']}
    groups=defaultdict(list)
    dense=[]
    for row in cert['rows']:
        if row['method']=='outward_positive_region':
            groups[(row['rho'],decode(row['s_dyadic']))].append(row)
        else:
            assert row['method']=='Arb_Perron_conditioned'
            dense.append(row)
    sparse_checks=0
    for (rho_text,s_fraction),rows in groups.items():
        rho,s=aa(Fraction(rho_text)),aa(s_fraction)
        b=(1+(-s).exp())/2
        end=b/(1<<22)
        stay=b-end
        zero=(arb(1),arb(0),end,stay)
        candidate=(1-rho+rho*end,rho*stay,end,stay)
        maximum=max(r['q'] for r in rows)
        regions=count_regions(zero,candidate,8192,maximum)
        for row in rows:
            q=row['q']
            full=power(regions[q],256)
            bound=(arb(math.comb(8192,q)).log()+q*arb(factors[rho_text]).log()+
                   (full[0]+full[1]).log()+209716*s)/arb(2).log()
            assert rational(bound.upper())<=decode(row['log2_upper'])
            sparse_checks+=1
        print('INDEPENDENT sparse',rho_text,maximum,flush=True)
    for row in dense:
        q=row['q']
        rho=aa(Fraction(row['rho']))
        p,s=aa(decode(row['p_dyadic'])),aa(decode(row['s_dyadic']))
        b=(1+(-s).exp())/2
        end=b/(1<<22)
        stay=b-end
        eta=p*rho
        # Direct full-length matrix power: no eigenvalues or Perron prefactor.
        matrix=(1-eta+eta*end,eta*stay,end,stay)
        complete=power(matrix,1<<21)
        choose=arb(math.comb(8192,q)).log()
        point=choose+q*p.log()
        if q<8192:
            point+=(8192-q)*(1-p).log()
        bound=(choose+q*arb(factors[row['rho']]).log()+(complete[0]+complete[1]).log()+
               209716*s-256*point)/arb(2).log()
        assert rational(bound.upper())<=decode(row['log2_upper']),q
        if q%1024==0:
            print('INDEPENDENT dense',q,flush=True)
    assert sparse_checks+len(dense)==8189
    result=dict(classification='Independent full-Arb count and direct-power comparison against all tail certificates',
        Arb_precision_bits=ctx.prec,unnormalized_region_count_checks=sparse_checks,
        direct_full_length_power_checks=len(dense),all_8189_occupations_checked=True,
        limitations='Checks transfer arithmetic; the mathematical prefix-domination and coefficient-bound arguments remain explicit proof steps.',
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__),source,envelope_path,
            ROOT/'code/audit_bch_q1_full_arb.py')})
    write_new(output,result)
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':
    main()
