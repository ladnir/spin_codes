"""Exact dyadic factorization screen for the sixteen-coordinate refresh.

This is a circuit-cost diagnostic, not a faster implementation or a new
distance certificate. It preserves H: interpolation from 0..15 to 16..31.
The operation count describes the explicit packed factor schedule only.
"""
from gf8_packet_maps import MUL, INV, product, xor_sum


def kernel():
    scale=MUL[product(16^j for j in range(16))][INV[product(range(1,16))]]
    return tuple(MUL[scale][INV[16^j]] for j in range(16))


def factors():
    # In characteristic two, epsilon_i=I+P_i satisfies epsilon_i^2=0,
    # where P_i shifts a coordinate index by XOR with 2^i. The upper
    # subset transform converts the dyadic kernel to this nilpotent basis.
    h=kernel()
    f=[xor_sum(h[t] for t in range(16) if t&s==s) for s in range(16)]
    assert f[0]==1
    result=[]
    for s in range(1,16):
        c=f[s]
        if not c:continue
        result.append((s,c))
        # Each factor is its own inverse. Earlier coefficients stay zero.
        old=f[:]
        for t in range(16):
            if not t&s:f[t|s]^=MUL[c][old[t]]
    assert f==[1]+[0]*15
    return tuple(result)


def evaluate(x):
    x=list(x)
    for s,c in factors():
        old=x[:]
        for i in range(16):
            x[i]^=MUL[c][xor_sum(old[i^t] for t in range(16) if t&s==t)]
    return tuple(x)


def fft_product_count(offset):
    # Match naturalFft16<Offset,*,true>: pair constants 0/1 need no
    # product; each nonzero packed small-stage constant vector needs one.
    def normalized(level,x):
        return MUL[product(x^j for j in range(1<<level))][
            INV[product((1<<level)^j for j in range(1<<level))]]
    count=2*int(normalized(3,offset) not in (0,1))
    count+=sum(normalized(2,offset^j) not in (0,1) for j in (0,8))
    for first in (0,4,8,12):
        count+=int(bool(normalized(1,offset^first)))
        count+=int(bool(normalized(0,offset^first) or normalized(0,offset^first^2)))
    return count


def check():
    h=kernel()
    for j in range(16):
        x=tuple(int(i==j) for i in range(16))
        value=evaluate(x)
        # Independent Lagrange formula, not the dyadic derivation.
        expected=tuple(MUL[product((16+i)^k for k in range(16) if k!=j)]
                       [INV[product(j^k for k in range(16) if k!=j)]]
                       for i in range(16))
        assert value==expected==tuple(h[i^j] for i in range(16))
        assert evaluate(value)==x
    # Four records occupy one ZMM. High index bits select ZMMs; low bits
    # select 128-bit lanes. A factor reduces high-bit parities before its
    # fixed product, but retains all low-bit lanes in each packed vector.
    count=sum(4//(1<<((s>>2).bit_count())) for s,c in factors() if c!=1)
    assert count==32
    fft=fft_product_count(0)+fft_product_count(16)
    assert fft==20
    return dict(basis_columns=16,exact=True,factors=factors(),
                factor_schedule_gfni=count,current_fft_schedule_gfni=fft,
                native_implementation=False,
                faster_circuit_claim=False)


if __name__=='__main__':
    import json
    print(json.dumps(check()))
