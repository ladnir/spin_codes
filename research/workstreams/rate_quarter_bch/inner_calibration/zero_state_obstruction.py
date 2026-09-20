"""Exact-integer lower certificates for the bad-message first moment.

Restrict to an even number Q of outer words of one fixed even weight w,
with Q*w <= the bad-weight cutoff. Uniform row-coordinate permutations
give even region counts with probability at least 2^(1-B). Exact rational
Chernoff comparisons leave probability at least 2^-B inside a count interval.
An explicit two-adjacent-weight kernel coefficient proves beta_j >= 2^-R
for every even region count in that interval. Thus E[Z] >= C(L,Q)*A_w^Q /
2^(B+B*R). This is NOT a lower bound on setup-failure probability.
"""
import argparse
import json
import math
from pathlib import Path

import check_candidates as check
import bch_zero_state_lower_v1 as concentration
import bch_zero_state_lower_v2 as lower


def exact_tail_check(q,a,b,boundary,exponent,lower_tail):
    """Prove the optimized binomial Chernoff bound <= 2^-exponent."""
    r = boundary
    assert 0<r<q and ((r*b<q*a) if lower_tail else (r*b>q*a))
    numerator = pow(a,r)*pow(b-a,q-r)*pow(q,q)
    denominator = pow(b,q)*pow(r,r)*pow(q-r,q-r)
    assert (numerator<<exponent)<=denominator,'tail charge too small'


def coefficient_term(kernel,epochs,j):
    a = 2*(j//(2*epochs))
    h = (j-a*epochs)//2
    assert a*epochs+2*h==j and 0<=h<epochs
    return math.comb(epochs,h)*pow(kernel[a],epochs-h)*(pow(kernel[a+2],h) if h else 1)


def certify(record,counts,distance):
    t = record['step_bits']; epochs = check.L//t
    assert check.L%t==0
    w = 52
    cutoff = check.N*distance.numerator//distance.denominator
    q = 2*(cutoff//w//2)
    assert q>0 and q<=check.L and w%2==0 and q*w<=cutoff
    lo,hi,_ = concentration.concentrated_range(q,w/check.B,check.B)
    # B union-bound terms, each with two tails <=2^-(B+8) for B=128,
    # consume at most 2^-B from the even-parity lower bound 2^(1-B).
    tail_bits = check.B+(2*check.B).bit_length()-1
    assert 1<<(tail_bits-check.B)==2*check.B
    exact_tail_check(q,w,check.B,lo-1,tail_bits,True)
    exact_tail_check(q,w,check.B,hi+1,tail_bits,False)
    first,last = lo+lo%2,hi-hi%2
    kernel = record['kernel_counts']
    logs = lower.two_weight_logs(kernel,epochs,last)
    proposed = math.ceil(-float(min(logs[first//2:last//2+1]))/math.log(2))+1
    for j in range(first,last+1,2):
        numerator = coefficient_term(kernel,epochs,j)
        denominator = math.comb(check.L,j)
        assert numerator>0 and (numerator<<proposed)>=denominator,'regional lower bound failed'
    family = math.comb(check.L,q)*pow(counts[w],q)
    family_bits = family.bit_length()-1
    assert family >= 1<<family_bits
    exponent = family_bits-check.B-check.B*proposed
    result = dict(distance_target=str(distance),bad_weight=cutoff,occupation=q,outer_weight=w,
                  shell_multiplicity=counts[w],actual_output_weight_on_zero_path=q*w,
                  region_count_interval=[lo,hi],even_region_interval=[first,last],
                  individual_tail_upper_bits=tail_bits,good_parity_lower_bits=check.B,
                  region_kernel_probability_lower_bits=proposed,epochs_per_region=epochs,
                  regional_coefficients_checked=(last-first)//2+1,
                  family_size_log2_floor=family_bits,expected_bad_messages_lower_power_of_two=exponent,
                  exact_integer_checks_passed=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tags',nargs='+',default=['t256_s18_nested','t256_s19_nested','t256_s20_nested'])
    args = parser.parse_args()
    folder = check.calibration.HERE
    screen_path = folder/'Q1_SCREEN.json'
    screen = json.loads(screen_path.read_text())
    for name,digest in screen['source_sha256'].items():
        assert check.fixed.sha(check.fixed.ROOT/name)==digest,name
    construction = check.calibration.smaller_outer.construction()
    counts = check.calibration.smaller_outer.spectrum()
    paths = [Path(__file__),Path(check.__file__),Path(concentration.__file__),Path(lower.__file__),screen_path,
             Path(check.calibration.smaller_outer.__file__),check.fixed.HERE/'sources/EBCH128_29.wd',
             check.fixed.HERE/'sources/EBCH128_36.wd',Path(check.fixed.outer.__file__)]
    rows = []
    for tag in args.tags:
        path = folder/'maps'/f'{tag}.json'
        entry = next(r for r in screen['results'] if r['tag']==tag)
        assert check.fixed.sha(path)==entry['map_sha256']
        record = json.loads(path.read_text())
        # Reconstruct the generator and spectra before using its kernel.
        audited = check.calibration.audit(record['step_bits'],record['state_bits'],
                                         [int(x,16) for x in record['quadratic_masks_hex']])
        assert all(audited[k]==record[k] for k in ('columns','generator_rows_hex','a_counts','kernel_counts'))
        paths.append(path)
        for distance in check.DISTANCES:
            row = dict(tag=tag,**certify(record,counts,distance))
            rows.append(row); print(row,flush=True)
    payload = dict(status='EXACT_INTEGER_BAD_MESSAGE_FIRST_MOMENT_LOWER_CERTIFICATES',outer=construction,
                   message_exponent=20,output_bits=check.N,results=rows,
                   limitations=['A large expected bad-message count does not imply a large setup-failure probability.',
                                'Nonpositive lower exponents do not establish an upper bound.',
                                'Float arithmetic only proposes intervals and exponents; every bound used is checked with integers.'],
                   source_sha256={p.relative_to(check.fixed.ROOT).as_posix():check.fixed.sha(p) for p in paths})
    (folder/'ZERO_STATE_OBSTRUCTION.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__': main()
