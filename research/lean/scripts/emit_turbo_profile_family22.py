"""One exact-rational physical-family candidate; Lean checks remain separate."""
from __future__ import annotations
import argparse
import json
import re
from fractions import Fraction as Q
from math import comb, isqrt
from pathlib import Path

STATE_COUNT = 1 << 22
LIVE_COUNT = STATE_COUNT - 1
LEVELS = [0, 0, 0, 16, 24, 28, 32, 36, 40, 48, 64]


def coefficients(profile, factors):
    values = [Q(1)]
    for count, (inactive, active) in zip(profile, factors):
        for _ in range(count):
            new = [Q(0) for _ in range(len(values)+1)]
            for j, value in enumerate(values):
                new[j] += inactive*value
                new[j+1] += active*value
            values = new
    return values


def upper(value, denominator, exact_boundary=False):
    assert value >= 0
    if exact_boundary and value in (0, 1):
        return int(value)*denominator
    return (value.numerator*denominator + value.denominator-1)//value.denominator + 2


def sqrt_quotient_ceiling(product, divisor):
    result = (isqrt(product)+divisor-1)//divisor
    if (result*divisor)**2 < product:
        result += 1
    return result


def lean_list(values):
    if values and isinstance(values[0], list):
        return '[\n' + ',\n'.join('  '+lean_list(row) for row in values) + '\n]'
    return '['+','.join(map(str,values))+']'


DEFAULT_NAMESPACE = 'Spin.Turbo.ProfileFamilySmoke22'


def emit_lean(data, path, namespace=None, hi=False):
    namespace = namespace or data.get('lean_namespace', DEFAULT_NAMESPACE)
    if not re.fullmatch(r'[A-Za-z_][A-Za-z_0-9]*(\.[A-Za-z_][A-Za-z_0-9]*)*', namespace):
        raise ValueError('The Lean namespace must consist of dotted ASCII identifiers.')
    lines = [('import SpinCodes.Turbo.Hi.PhysicalProfileNumericData22' if hi else
              'import SpinCodes.Turbo.PhysicalProfileNumericData22'),
             f'namespace {namespace}', ('open Spin.NumericHi' if hi else 'open Spin.Numeric'),
             '/- Untrusted candidates until the integer evaluator accepts them. -/',
             f'def zNum : Int := {data["z_num"]}', f'def zDen : Int := {data["z_den"]}',
             'def zbox : Fix := Fix.ofFrac zNum zDen',
             f'def endpointDen : Int := {data["endpoint_den"]}']
    for name, key, ty in [
        ('firstCaps','first_caps','List (List Int)'),
        ('secondCaps','second_caps','List (List Int)'),
        ('returnCaps','return_caps','List (List Int)'),
        ('feedbackActiveCaps','feedback_active_caps','List Int'),
        ('zeroHoldCaps','zero_hold_caps','List Int'),
        ('birthOneCaps','birth_one_caps','List Int'),
        ('birthTwoCaps','birth_two_caps','List Int')]:
        lines.append(f'def {name} : {ty} := {lean_list(data[key])}')
    lines.append(f'end {namespace}')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')


def produce_family(candidate, z_num, z_den, namespace=DEFAULT_NAMESPACE, endpoint_bits=80):
    if not 0 < z_num <= z_den:
        raise ValueError('Expected 0 < z numerator <= positive denominator.')
    profiles, counts = candidate['profiles'], candidate['multiplicities']
    assert len(profiles) == 25 and sum(counts) == STATE_COUNT
    assert all(len(p)==5 and sum(p)==16 for p in profiles)
    weights = [sum(r*n for r,n in enumerate(p)) for p in profiles]
    z = Q(z_num,z_den)
    endpoint_den = 1 << endpoint_bits
    divisors = [comb(16,j)*15**j for j in range(17)]
    moments = []
    for point in (z,z*z):
        bank = [(point**r,(1+point)**4-point**r) for r in range(5)]
        moments.append([[v/divisors[j] for j,v in enumerate(coefficients(p,bank))] for p in profiles])
    feedback, zero_hold = [],[]
    for point,target in ((Q(1),feedback),(z,zero_hold)):
        bank = [(Q(1),(1+point)**(4-r)*(1-point)**r-1) for r in range(5)]
        totals = [Q(0) for _ in range(17)]
        for p,count in zip(profiles,counts):
            for j,v in enumerate(coefficients(p,bank)):
                totals[j] += count*v
        target.extend(totals[j]/(STATE_COUNT*divisors[j]) for j in range(17))
    assert all(0<=v<=1 for group in moments for row in group for v in row)
    assert all(0<=v<=1 for v in feedback+zero_hold)
    indices = [[k for k,w in enumerate(weights) if w==0],
               [k for k,w in enumerate(weights) if w!=0],
               [k for k,w in enumerate(weights) if w!=0]]
    indices += [[k for k,w in enumerate(weights) if w==level] for level in LEVELS[3:]]
    class_counts = [sum(counts[k] for k in ks) for ks in indices]
    assert class_counts == [1,4194303,4194303,2604,291648,888832,1828134,888832,291648,2604,1]
    caps = []
    for values in moments:
        rows = []
        for j in range(17):
            row = []
            for i,ks in enumerate(indices):
                value = (sum((counts[k]*values[k][j] for k in ks),Q(0))/class_counts[i]
                         if i==2 else max(values[k][j] for k in ks))
                row.append(upper(value,endpoint_den,exact_boundary=(i==0 and j==0)))
            rows.append(row)
        caps.append(rows)
    feedback_caps = [upper(1-v,endpoint_den,exact_boundary=True) for v in feedback]
    assert all(0<=n<=endpoint_den for n in feedback_caps)
    returns = []
    for j in range(17):
        row = [0]
        for i in range(1,11):
            first,second,active = caps[0][j][i],caps[1][j][i],feedback_caps[j]
            r = min((first+LIVE_COUNT-1)//LIVE_COUNT,
                    sqrt_quotient_ceiling(second*active,LIVE_COUNT))
            assert first<=r*LIVE_COUNT or second*active<=(r*LIVE_COUNT)**2
            row.append(r)
        returns.append(row)
    birth_vectors = {}
    positive_birth_values = []
    for birth in candidate['small_births']:
        values = [sum((Q(n)*z**w for w,n in enumerate(poly)),Q(0))/birth['denominator']
                  for poly in birth['coefficients']]
        j = birth['occupancy']
        assert sum(values)==moments[0][0][j]
        assert values[0]==zero_hold[j]
        positive_birth_values.extend(value for value in values if value > 0)
        vector = [0]*11
        for k,value in enumerate(values):
            vector[0 if k==0 else k+2] = upper(value,endpoint_den,exact_boundary=True)
        birth_vectors[j] = vector
    data = dict(schema='turbo-rational-profile-family-candidate-1',
                scope='Exact rational producer candidate; no row certificate until Lean checks accept it.',
                lean_namespace=namespace,
                z_num=z_num,z_den=z_den,endpoint_den=endpoint_den,
                first_caps=caps[0],second_caps=caps[1],return_caps=returns,
                feedback_active_caps=feedback_caps,
                zero_hold_caps=[upper(v,endpoint_den,exact_boundary=(j==0)) for j,v in enumerate(zero_hold)],
                birth_one_caps=birth_vectors[1],birth_two_caps=birth_vectors[2],
                feedback_zero_exact=[dict(num=v.numerator,den=v.denominator) for v in feedback])
    # These are exact producer diagnostics, not additional proof assumptions.
    audit = {
        'moment_uniform_lower': z**128,
        'minimum_positive_feedback_active': min(1-v for v in feedback if v < 1),
        'minimum_positive_zero_hold': min(v for v in zero_hold if v > 0),
        'minimum_positive_small_birth': min(positive_birth_values),
    }
    data['precision_audit'] = {key: dict(num=value.numerator,den=value.denominator)
                               for key,value in audit.items()}
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',type=Path,default=Path('.lake/turbo/coarse_packet_profiles.json'))
    parser.add_argument('--output',type=Path,default=Path('.lake/turbo/profile_family_3_20_rational.json'))
    parser.add_argument('--lean-output',type=Path,default=Path('.lake/turbo/ProfileFamilySmoke22Data.lean'))
    parser.add_argument('--z-num',type=int)
    parser.add_argument('--z-den',type=int)
    parser.add_argument('--namespace',help='Dotted Lean namespace; retained JSON supplies it in emit-only mode.')
    parser.add_argument('--emit-only',action='store_true')
    args = parser.parse_args()
    if args.emit_only:
        if args.z_num is not None or args.z_den is not None:
            parser.error('--emit-only uses retained rational data; do not supply --z-num/--z-den.')
        emit_lean(json.loads(args.output.read_text(encoding='utf-8')),args.lean_output,args.namespace)
        print(f'Copied retained candidate definitions to {args.lean_output}.')
        return
    if (args.z_num is None) != (args.z_den is None):
        parser.error('--z-num and --z-den must be supplied together.')
    z_num = args.z_num if args.z_num is not None else 8607079764250578
    z_den = args.z_den if args.z_den is not None else 10000000000000000
    candidate = json.loads(args.candidate.read_text(encoding='utf-8'))
    data = produce_family(candidate,z_num,z_den,args.namespace or DEFAULT_NAMESPACE)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    emit_lean(data,args.lean_output)
    print(f'Emitted one exact-rational 25-profile/17-occupation family to {args.output}.')
    print(f'Birth j=1,2 vectors and shared 2^80 upper endpoints: {args.lean_output}.')
    print('No DP, FWHT, raw-state enumeration, floating powers, or logarithms were used.')


if __name__=='__main__':
    main()
