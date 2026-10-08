"""Blockwise q-ary IMT refresh laws; no distance or timing certificate.

The primary coefficient field is GF(2^16), not the 16-element field.
Exact small-field enumeration validates rank-one laws independently.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import product


def multiply(a,b,q):
    """Tiny characteristic-two fields for independent exhaustive tests."""
    modulus = {2:0b11,4:0b111,8:0b1011,16:0b10011}[q]
    answer = 0
    while b:
        if b&1:
            answer ^= a
        b >>= 1
        a <<= 1
        if a&q:
            a ^= modulus
    return answer


def dot(left,right,q):
    answer = 0
    for a,b in zip(left,right):
        answer ^= multiply(a,b,q)
    return answer


def add(left,right):
    return tuple(a^b for a,b in zip(left,right))


def scale(scalar,vector,q):
    return tuple(multiply(scalar,x,q) for x in vector)


def inverse(value,q):
    if not value:
        raise ZeroDivisionError
    return next(other for other in range(1,q) if multiply(value,other,q)==1)


def apply_rank_one(vector,u,v,q):
    return add(vector,scale(dot(v,vector,q),u,q))


def allowed_pairs(q,state,strict=False):
    """u uniform nonzero; v uniform vu!=1 (or vu=0 for strict transvection)."""
    vectors = tuple(product(range(q),repeat=state))
    for u in vectors:
        if any(u):
            for v in vectors:
                inner = dot(v,u,q)
                valid = inner==0 if strict else inner!=1
                if valid:
                    yield u,v


def law_weights(q,state,rounds=1,scalar_symmetrized=False):
    """Exact nonzero-state kernel: identity, same-line uniform, full uniform.

    Unscaled uses the invertible rank-one update vu!=1. Symmetrized can use
    either this update or a strict transvection vu=0, then a fresh scalar.
    """
    if q<2 or state<1 or rounds<0:
        raise ValueError("invalid field size, state size or round count")
    if rounds==0:
        return F(1),F(0),F(0)
    lazy = F(1,q)**rounds
    point = F(0) if scalar_symmetrized else F(q**(state-1)-1,q**state-1)**rounds
    return point,lazy-point,1-lazy


def theoretical_image(source,target,q,rounds=1,scalar_symmetrized=False):
    if len(source)!=len(target):
        raise ValueError("mismatched state dimension")
    if not any(source):
        return F(int(not any(target)))
    if not any(target):
        return F(0)
    point,line,refresh = law_weights(q,len(source),rounds,scalar_symmetrized)
    in_line = any(scale(c,source,q)==tuple(target) for c in range(1,q))
    return point*(tuple(source)==tuple(target))+line*in_line/(q-1)+refresh/(q**len(source)-1)


def enumerated_image(source,q,strict=False,symmetrize=False):
    counts,total = Counter(),0
    for u,v in allowed_pairs(q,len(source),strict):
        value = apply_rank_one(source,u,v,q)
        for scalar in (range(1,q) if symmetrize else (1,)):
            counts[scale(scalar,value,q)] += 1
            total += 1
    return {target:F(count,total) for target,count in counts.items()}


def step_law(state,syndrome,q,rounds=1,scalar_symmetrized=False):
    """Exact fixed-state successor law after ONE feedback addition."""
    size = len(state)
    if len(syndrome)!=size:
        raise ValueError("mismatched syndrome")
    result = {}
    for refreshed in product(range(q),repeat=size):
        probability = theoretical_image(state,refreshed,q,rounds,scalar_symmetrized)
        if probability:
            result[add(refreshed,syndrome)] = probability
    return result


def return_probability(state,syndrome,q,rounds=1,scalar_symmetrized=False):
    return theoretical_image(state,syndrome,q,rounds,scalar_symmetrized)


def matvec(matrix,vector,q):
    return tuple(dot(row,vector,q) for row in matrix)


def transpose(matrix):
    return tuple(zip(*matrix))


def encode(inputs,emission,feedback,schedule,q):
    """Literal SPIN recurrence: output before refresh, feedback uses raw X."""
    state = (0,)*len(feedback)
    outputs = []
    for x,pairs in zip(inputs,schedule):
        outputs.append(add(x,matvec(emission,state,q)))
        for u,v in pairs:
            state = apply_rank_one(state,u,v,q)
        state = add(state,matvec(feedback,x,q))
    return tuple(outputs)


def decode(outputs,emission,feedback,schedule,q):
    state = (0,)*len(feedback)
    inputs = []
    for y,pairs in zip(outputs,schedule):
        x = add(y,matvec(emission,state,q))
        inputs.append(x)
        for u,v in pairs:
            state = apply_rank_one(state,u,v,q)
        state = add(state,matvec(feedback,x,q))
    return tuple(inputs)


def adjoint(outputs,emission,feedback,schedule,q):
    """Reverse product order and transpose each rank-one map."""
    state = (0,)*len(feedback)
    inputs = []
    for y,pairs in reversed(tuple(zip(outputs,schedule))):
        inputs.append(add(y,matvec(transpose(feedback),state,q)))
        following = state
        for u,v in reversed(pairs):
            following = apply_rank_one(following,v,u,q)
        state = add(matvec(transpose(emission),y,q),following)
    return tuple(reversed(inputs))


def sample_fiber(u,q):
    """Enumerate the exact pivot sampler, used only as a small-field test.

    Sample all but one v coordinate uniformly and vu uniformly in F minus {1}.
    Every permitted v is generated exactly once for each fixed nonzero u.
    """
    if not any(u):
        raise ValueError("u must be nonzero")
    pivot = next(i for i,x in enumerate(u) if x)
    others = [i for i in range(len(u)) if i!=pivot]
    for values in product(range(q),repeat=len(others)):
        for pairing in range(q):
            if pairing==1:
                continue
            v = [0]*len(u)
            for i,value in zip(others,values):
                v[i] = value
            v[pivot] = multiply(inverse(u[pivot],q),pairing^dot(v,u,q),q)
            yield tuple(v)


def operation_receipts(q=65536):
    """Arithmetic counts only; no native instruction or elapsed-time model."""
    for width,state in ((16,16),(32,16),(64,16),(64,32),(64,64)):
        point,line,refresh = law_weights(q,state)
        yield dict(width=width,state=state,field_size=q,
                   rank_one_products_per_packet=2*state,
                   rank_one_products_per_raw_symbol=str(F(2*state,width)),
                   scaled_transvection_products_per_packet=3*state,
                   dense_GL_products_per_packet=state*state,
                   dense_fixed_map_products_per_packet=2*width*state,
                   MDS_emission_minimum=max(1,width-state+1),
                   MDS_feedback_kernel_minimum=state+1 if state<width else None,
                   lazy_mass=str(point+line),
                   extra_random_products_beyond_fixed_maps=True,
                   timing_claim=False,whole_certificate=False)


if __name__=='__main__':
    import json
    for receipt in operation_receipts():
        print(json.dumps(receipt),flush=True)
