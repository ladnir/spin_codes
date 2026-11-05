"""Exact rational obstruction for sourcewise Cauchy quota comparisons.

This is a proof-comparison obstruction, not a codeword or physical-quota
failure witness. No numerical work occurs on import. Run --check through
the root's bounded file-free remote runner.
"""
from fractions import Fraction as F
from math import comb
from functools import lru_cache

RAW = (
    (255, 204, 0, 0, 0),
    (0, 816, 1530, 1020, 255),
    tuple(comb(4,j)*(256//(4**j)-1) for j in range(5)),
    tuple(comb(4,j)*(256-256//(4**j)) for j in range(5)),
    tuple(comb(4,j)*(4096//(4**j)-16//(2**j)) for j in range(5)),
    tuple(comb(4,j)*(4096//(2**j)-4096//(4**j)) for j in range(5)),
)
EDGES=((1,0,0),(0,0,1),(3,1,0),(2,1,2),(5,2,0),(4,2,2))

def mm(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

def mdadd(A,B):
    return [[A[i][j]+B[i][j] for j in range(3)] for i in range(3)]

def restricted_power(x):
    """Integer numerator of L_x^16 and its x*d/dx derivative."""
    a,b=x.numerator,x.denominator
    values=[sum(c*a**j*b**(4-j) for j,c in enumerate(row)) for row in RAW]
    derivs=[sum(j*c*a**j*b**(4-j) for j,c in enumerate(row)) for row in RAW]
    A=[[0]*3 for _ in range(3)]; DA=[[0]*3 for _ in range(3)]
    for k,i,j in EDGES:
        A[i][j]=values[k]**128
        DA[i][j]=128*values[k]**127*derivs[k]
    # The denominators are constant at this evaluation point. Derivative
    # numerators contain the derivative of the original polynomial.
    for _ in range(4):
        next_DA=mdadd(mm(DA,A),mm(A,DA))
        A=mm(A,A); DA=next_DA
    denominator=(255*b**4)**2048
    return A,DA,denominator

@lru_cache(maxsize=32)
def ln_unit_interval(x):
    """Rigorous rational bounds for ln(x),1<=x<=2,24 positive terms."""
    assert 1<=x<=2
    z=(x-1)/(x+1)
    total=F(0); power=z
    for k in range(24):
        total+=2*power/(2*k+1)
        power*=z*z
    tail=2*power/(49*(1-z*z))
    return total,total+tail

def log2_bounds(n,d=1):
    """Dyadic mantissa enclosure followed by rational log intervals."""
    assert isinstance(n,int) and isinstance(d,int) and n>0 and d>0
    e=n.bit_length()-d.bit_length()
    below=n<(d<<e) if e>=0 else (n<<(-e))<d
    if below: e-=1
    nn,dd=(n,d<<e) if e>=0 else (n<<(-e),d)
    assert dd<=nn<2*dd
    bits=64; mantissa=(nn<<bits)//dd
    ml=F(mantissa,1<<bits); mu=F(mantissa+1,1<<bits)
    lower,_=ln_unit_interval(ml); _,upper=ln_unit_interval(mu)
    log2lower,log2upper=ln_unit_interval(F(2))
    return e+lower/log2upper,e+upper/log2lower

def check():
    C=comb(8192,4096)
    assert (1<<8178)<C<(1<<8186)
    # Explicit dyadic vector is a subeigenvector only for the relaxed
    # comparison. It is zero on physical states outside O,P,H.
    v=(1,1<<314,1<<677)
    brackets=((F(1819,1000),F(1822,1000)),
              (F(2112,1000),F(2116,1000)),
              (F(2316,1000),F(2320,1000)))
    results=[]
    for row,(xl,xu) in enumerate(brackets):
        data=[]
        for x in (xl,xu):
            A,DA,den=restricted_power(x)
            P=sum(A[row][j]*v[j] for j in range(3))
            DP=sum(DA[row][j]*v[j] for j in range(3))
            slope=F(DP,P)-4096
            logratio=log2_bounds(P,den*C*v[row])
            logx=log2_bounds(x.numerator,x.denominator)
            lower_f=logratio[0]-4096*logx[1]
            data.append((slope,lower_f,logx))
        (sl,fl,tl),(su,fu,tu)=data
        assert sl<0<su, (row,float(sl),float(su))
        # Convexity locates every minimum in [tl,tu]. The tangent at the
        # left endpoint is decreasing, so its right-end value is a
        # uniform lower bound on that interval and on the global minimum.
        width_upper=tu[1]-tl[0]
        lower=fl+sl*width_upper
        assert lower>400, (row,float(lower))
        results.append(dict(row=row,bracket=(str(xl),str(xu)),
                            tangent_slope_left=float(sl),tangent_slope_right=float(su),
                            certified_minimum_log2_lower_display=float(lower),
                            exact_asserted_lower_bits=400))
    assert 8178-60518+255*400==49660
    print(dict(exact_rational_checks=results,all_positive_x_bulk_floor_bits=400,
               all_positive_potentials=True,complete_comparison_floor_bits=49660,
               physical_exact_quota_operator_refuted=False,encoder_refuted=False,
               geometry_tail_proved=False),flush=True)

def main():
    import argparse,os,resource,shutil,threading
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true",required=True)
    parser.parse_args()
    if shutil.disk_usage("/tmp").free<20*(1<<30):
        raise RuntimeError("20GiB free-space guard")
    resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    timer=threading.Timer(175,lambda:os._exit(124));timer.daemon=True;timer.start()
    try: check()
    finally: timer.cancel()

if __name__=="__main__":
    main()

