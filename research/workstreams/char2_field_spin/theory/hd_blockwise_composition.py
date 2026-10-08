"""Capped support-wise composition and directed one-active-row receipt.

No whole-code certificate. Fractional powers act on MDS prescribed-support
counts and local transfer entries, never on route probabilities or the
number of support patterns. The code also records a route-only type bound.
"""

from collections import Counter
from fractions import Fraction as F
from math import comb,factorial

from flint import arb,ctx

from hd_blockwise_maps import local_transfer


def prescribed_count(length,dimension,weight,q):
    if weight==0:
        return 1
    parity=length-dimension
    if not parity<weight<=length:
        return 0
    return (q-1)*sum((-1)**j*comb(weight-1,j)*q**(weight-parity-1-j)
                     for j in range(weight-parity))


def rational(x):
    x=F(x)
    return arb(x.numerator)/x.denominator


def zero_matrix(size,zero=F(0)):
    return [[zero for _ in range(size)] for _ in range(size)]


def identity_matrix(size,zero=F(0),one=F(1)):
    return [[one if i==j else zero for j in range(size)] for i in range(size)]


def matmul(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),a[0][0]*0)
             for j in range(len(b[0]))] for i in range(len(a))]


def add(a,b):
    return [[x+y for x,y in zip(left,right)] for left,right in zip(a,b)]


def marked_power(base,mark,exponent):
    """Return base^B and sum_j base^j mark base^(B-1-j), exactly."""
    zero=base[0][0]*0
    one=zero+1
    result=(identity_matrix(len(base),zero,one),zero_matrix(len(base),zero))
    current=(base,mark)
    while exponent:
        if exponent&1:
            result=(matmul(result[0],current[0]),
                    add(matmul(result[0],current[1]),matmul(result[1],current[0])))
        current=(matmul(current[0],current[0]),
                 add(matmul(current[0],current[1]),matmul(current[1],current[0])))
        exponent>>=1
    return result


def apply_row(row,matrix):
    return [sum((row[i]*matrix[i][j] for i in range(len(row))),row[0]*0)
            for j in range(len(row))]


def region_coefficients(absent,present,length):
    """Positive coefficients of e0*(absent+z present)^length*1."""
    zero=absent[0][0]*0
    coefficients=[[zero+1,zero,zero]]
    for _ in range(length):
        following=[[zero,zero,zero] for _ in range(len(coefficients)+1)]
        for weight,row in enumerate(coefficients):
            no=apply_row(row,absent)
            yes=apply_row(row,present)
            following[weight]=[x+y for x,y in zip(following[weight],no)]
            following[weight+1]=[x+y for x,y in zip(following[weight+1],yes)]
        coefficients=following
    return [sum(row,zero) for row in coefficients]


def entry_power(matrix,rho):
    exponent=rational(rho)
    return [[arb(0) if value==0 else rational(value)**exponent for value in row] for row in matrix]


def singleton_receipt(state=16,width=64,length=128,dimension=64,rows=16384,
                      q=65536,theta=F(9993,10000),rho=F(3,4),distance=F(2,5),
                      uniform_lazy='density',precision=192):
    """All nonzero messages supported in exactly one outer row.

    The bound unions every row and every exact MDS coordinate support, and
    averages its one marked position per occupied region exactly. It does
    not cover even the two-active-row class.
    """
    if (rows%width or not 0<rho<=1 or not 0<theta<1 or state*2>width
            or not 0<dimension<length or q<length):
        raise ValueError('invalid construction/witness')
    old=ctx.prec
    try:
        ctx.prec=precision
        n=length*rows
        threshold=(n*distance.numerator-1)//distance.denominator
        quiet=entry_power(local_transfer(0,theta,width,state,q,uniform_lazy=uniform_lazy),rho)
        active=entry_power(local_transfer(1,theta,width,state,q,uniform_lazy=uniform_lazy),rho)
        packets=rows//width
        absent,insertion=marked_power(quiet,active,packets)
        present=[[x/packets for x in row] for row in insertion]
        coefficients=region_coefficients(absent,present,length)
        total=arb(0)
        terms=[]
        for weight in range(length-dimension+1,length+1):
            count=prescribed_count(length,dimension,weight,q)
            term=arb(count)**rational(rho)*coefficients[weight]
            total+=term
            terms.append((weight,term))
        prefactor=rows*arb(q-1)**(-rational(rho))*rational(theta)**(-rational(rho)*threshold)
        failure=total*prefactor
        bits=-failure.log()/arb(2).log()
        worst=max(terms,key=lambda pair:float(pair[1].log()))[0]
        return dict(construction='blockwise-residue-rankone-no-input-diagonal',
                    outer=[length,dimension],rows=rows,n=n,width=width,state=state,q=q,
                    distance=str(distance),output_threshold=threshold,
                    source_class='exactly one nonzero outer row',theta=str(theta),rho=str(rho),
                    uniform_lazy=uniform_lazy,all_coordinate_supports=True,
                    all_source_rows=True,all_mixer_setups_averaged=True,
                    route_probabilities_unpowered=True,support_pattern_multiplicities_unpowered=True,
                    source_row_choices_unpowered=True,prescribed_support_count_power=str(rho),
                    largest_support_weight=worst,margin_bits=str(bits),
                    margin_bits_floor=int(bits.lower().fmpq().p//bits.lower().fmpq().q),
                    rigorous=True,whole_code=False,precision=precision)
    finally:
        ctx.prec=old


def type_probability(values):
    """Probability one uniform permutation equals one feasible target vector."""
    counts=Counter(values)
    numerator=1
    for count in counts.values():
        numerator*=factorial(count)
    return F(numerator,factorial(len(values)))


def information_set_type_bound(columns,dimension):
    """Fixed inverse word; choose parity columns with greatest type entropy."""
    probabilities=sorted(type_probability(column) for column in columns)
    result=F(1)
    for probability in probabilities[:len(columns)-dimension]:
        result*=probability
    return result


def type_gap_receipt(rows=16384,length=128,dimension=64,q=65536,
                     distance=F(2,5),precision=160):
    """Diagnostic: even all-distinct columns cannot make this bound 40%."""
    old=ctx.prec
    try:
        ctx.prec=precision
        n=rows*length
        e=(n*distance.numerator-1)//distance.denominator
        entropy=(length-dimension)*arb(factorial(rows)).log()/arb(2).log()
        # One ball shell is a lower bound. The binomial type is modal at p=e/n,
        # so its probability is >=1/(n+1), giving this directed entropy bound.
        p=arb(e)/n
        binomial_lower=-e*p.log()-(n-e)*(1-p).log()-arb(n+1).log()
        ball_lower=(binomial_lower+e*arb(q-1).log())/arb(2).log()
        projective_lower=ball_lower-arb(q-1).log()/arb(2).log()
        return dict(type_entropy_upper_bits=str(entropy),ball_shell_lower_bits=str(ball_lower),
                    projective_shell_lower_bits=str(projective_lower),
                    unavoidable_log2_gap=str(projective_lower-entropy),
                    relation='limitation of this route-only upper-bound strategy',
                    actual_failure_lower_bound=False,rigorous=True)
    finally:
        ctx.prec=old


def fractional_obstruction(active=64,state=16,rho=F(3,4),precision=160):
    """Lower bound on the uncapped support-wise fractional moment itself.

    On the previously proved route event, every mixer setup has a large
    low-output subspace. Concavity of the power gives sum_P Z_P^rho >=
    (sum_P Z_P)^rho, irrespective of how its exact supports are distributed.
    This is not a lower bound on actual distance-failure probability.
    """
    from hd_blockwise_global import moment_obstruction
    receipt=moment_obstruction(active,state,precision)
    old=ctx.prec
    try:
        ctx.prec=precision
        lower=arb(receipt['route_event_log2'])+16*(receipt['message_dimension_lower']-2)*rational(rho)
        return dict(active_rows=active,state=state,rho=str(rho),
                    source_dimension_lower=receipt['message_dimension_lower'],
                    fractional_sum_log2_lower=str(lower),
                    rules_out_common_rho_bound=bool(lower>0),
                    actual_failure_lower_bound=False,rigorous=True)
    finally:
        ctx.prec=old


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser()
    parser.add_argument('--state',type=int,default=16)
    parser.add_argument('--rho',type=F,default=F(3,4))
    parser.add_argument('--theta',type=F,default=F(9993,10000))
    args=parser.parse_args()
    print(json.dumps(singleton_receipt(state=args.state,rho=args.rho,theta=args.theta),indent=2))
