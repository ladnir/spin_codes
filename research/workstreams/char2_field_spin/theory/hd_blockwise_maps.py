"""Fixed residue maps and a composable three-component local envelope.

No whole-code certificate or native timing. The native profile uses embedded
AES GF256 maps on GF65536 symbols, with full-field rank-one refreshes.
"""

from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import comb

from hd_blockwise_inner import law_weights


def mul(a,b,q=256):
    modulus={2:3,4:7,8:11,16:19,256:0x11b}[q]
    answer=0
    while b:
        if b&1:
            answer^=a
        b>>=1
        a<<=1
        if a&q:
            a^=modulus
    return answer


def power(a,n,q=256):
    answer=1
    while n:
        if n&1:
            answer=mul(answer,a,q)
        a=mul(a,a,q)
        n>>=1
    return answer


def xor_sum(values):
    answer=0
    for value in values:
        answer^=value
    return answer


def matvec(matrix,vector,q=256):
    return tuple(xor_sum(mul(a,b,q) for a,b in zip(row,vector)) for row in matrix)


def matmul(left,right,q=256):
    columns=tuple(zip(*right))
    return tuple(tuple(xor_sum(mul(a,b,q) for a,b in zip(row,column))
                       for column in columns) for row in left)


@lru_cache(None)
def maps(width=64,state=16,base=128,q=256):
    """Lagrange interpolation V=base+[0,s), U=[0,t); C is A transpose."""
    if (state<1 or width<2*state or width&(width-1)
            or state&(state-1) or base%state or base<width or base+state>q):
        raise ValueError('require disjoint additive U,V and width>=2*state')
    nodes=tuple(base^j for j in range(state))
    denominators=tuple(_product((v^other for other in nodes if other!=v),q=q)
                       for v in nodes)
    assert len(set(denominators))==1
    rows=[]
    for u in range(width):
        numerator=_product((u^v for v in nodes),q=q)
        rows.append(tuple(mul(numerator,power(mul(denominators[j],u^v,q),q-2,q),q)
                          for j,v in enumerate(nodes)))
    a=tuple(rows)
    return a,tuple(zip(*a))


def _product(values,q):
    result=1
    for value in values:
        result=mul(result,value,q)
    return result


def native_receipt():
    a,c=maps()
    ca=matmul(c,a)
    assert not any(x for row in ca for x in row)
    return dict(width=64,state=16,coefficient_subfield=256,field_size=65536,
                input_nodes=list(range(64)),state_nodes=list(range(128,144)),
                A_first_row=list(a[0]),CA_zero=True,C_is_A_transpose=True,
                A_sha256=sha256(bytes(x for row in a for x in row)).hexdigest(),
                emission_minimum=49,feedback_kernel_minimum=17,
                whole_certificate=False)


def emission_minimum(weight,width,state,syndrome_nonzero=None):
    """Guaranteed wt(x+Az) for z!=0, optionally using Cx and CA=0."""
    if not 0<=weight<=width:
        raise ValueError('invalid input weight')
    floor=max(0,width-state+1-weight)
    if syndrome_nonzero is True:
        return max(1,floor)
    if syndrome_nonzero is False and weight<width-state+1:
        return max(state+1,floor)
    return floor


@lru_cache(None)
def emission_upper(weight,theta,width=64,state=16,q=65536):
    """Fixed x, uniform nonzero state; MDS forced-zero expansion, all exact.

    Outside supp(x), s zero equations force z=0 and contribute nothing.
    With at least one selected nonzero x coordinate, a consistent system
    has q**max(s-number_of_equations,0) states; consistency is relaxed.
    """
    theta=F(theta)
    if not (0<=weight<=width and 0<theta<=1 and state<=width):
        raise ValueError('invalid local parameters')
    total=F(0)
    for inside in range(weight+1):
        for outside in range(min(width-weight,state-1)+1):
            count=(q**(state-outside)-1 if inside==0
                   else q**max(0,state-inside-outside))
            zeros=inside+outside
            total+=(comb(weight,inside)*comb(width-weight,outside)*count
                    *theta**(width-zeros)*(1-theta)**zeros)
    return min(theta**emission_minimum(weight,width,state),total/(q**state-1))


def local_transfer(weight,theta,width=64,state=16,q=65536,rounds=1,
                   syndrome_nonzero=None,uniform_lazy='density'):
    """Rows/columns: zero atom, arbitrary-live mass, uniform-live mass.

    Uniform mass U means density U/(q**s-1) on every nonzero state.
    Each entry is a positive domination, not a transition probability.
    None takes entrywise maxima of the feasible syndrome-zero/nonzero cases.
    """
    theta=F(theta)
    if (not 0<=weight<=width or not 0<theta<=1 or rounds<0
            or uniform_lazy not in ('density','mass')):
        raise ValueError('invalid local parameters')
    feasible=[False] if weight==0 else ([True] if weight<=state else [False,True])
    if syndrome_nonzero is None:
        choices=[local_transfer(weight,theta,width,state,q,rounds,c,uniform_lazy) for c in feasible]
        return tuple(tuple(max(m[i][j] for m in choices) for j in range(3)) for i in range(3))
    if syndrome_nonzero not in feasible:
        raise ValueError('syndrome flag impossible by MDS feedback distance')
    point,line,refresh=law_weights(q,state,rounds)
    lazy=point+line
    denominator=q**state-1
    low=theta**emission_minimum(weight,width,state,syndrome_nonzero)
    average=min(low,emission_upper(weight,theta,width,state,q))
    return_numerator=lazy*low+refresh*average
    uniform_live=(return_numerator if uniform_lazy=='density' else refresh*average)
    uniform_arbitrary=F(0) if uniform_lazy=='density' else lazy*average
    zero=F(0)
    if syndrome_nonzero:
        maximum_return=point+line/(q-1)+refresh/denominator
        return ((zero,theta**weight,zero),
                (maximum_return*low,lazy*low,refresh*low),
                (return_numerator/denominator,uniform_arbitrary,uniform_live))
    return ((theta**weight,zero,zero),
            (zero,lazy*low,refresh*low),
            (zero,uniform_arbitrary,uniform_live))


def apply_local(mass,matrix):
    return tuple(sum((mass[i]*matrix[i][j] for i in range(3)),F(0)) for j in range(3))


if __name__=='__main__':
    import json
    print(json.dumps(native_receipt(),indent=2))
    for weight in (0,1,16,17,32,48,49,64):
        matrix=local_transfer(weight,F(1,2))
        print('weight',weight,'emission',float(emission_upper(weight,F(1,2))),
              'uniform_return',float(matrix[2][0]),'uniform_live',float(matrix[2][2]))
