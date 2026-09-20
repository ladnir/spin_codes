"""Outward Q1-only bound for one transvection round; not a full certificate.

Recomputes retained map/outer spectra and single-input cancellations. No
floating-point diagnostic or neighbor table is trusted. Verify mode replays
at 512 bits and checks the saved upward dyadic bounds from 256 bits.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

from flint import arb, ctx
import transvection as tv

HERE=Path(__file__).resolve().parent


def identity(n):return tuple(arb(i==j) for i in range(n) for j in range(n))


def mul(a,b,n):
    return tuple(sum((a[i*n+k]*b[k*n+j] for k in range(n)),arb(0))
                 for i in range(n) for j in range(n))


def add(a,b):return tuple(x+y for x,y in zip(a,b))


def transfers(spectrum,columns,lam,rounds):
    weights=sorted(spectrum);n=len(weights)+2;t=len(columns);m=sum(spectrum.values())
    z=(-lam).exp();epsilon=arb(1)/(1<<rounds);eta=1-epsilon
    f0=[z**w for w in weights]
    f1=[(w*z**(w-1)+(t-w)*z**(w+1))/t for w in weights]
    d0=max(f.upper() for f in f0);d1=max(f.upper() for f in f1)
    cw=tv.cancellation_weights(columns)
    state_weights=[sum(tv.dot(col,q) for col in columns) for q in columns]
    zero=[arb(0) for _ in range(n*n)];one=zero.copy()
    zero[0]=arb(1);one[1]=z
    for matrix,f,d in ((zero,f0,d0),(one,f1,d1)):
        matrix[n+1]=epsilon*d
        for j,w in enumerate(weights):
            matrix[n+j+2]=eta*d*spectrum[w]/m
            for i,v in enumerate(weights):matrix[(i+2)*n+j+2]=eta*f[i]*spectrum[w]/m
    for i,v in enumerate(weights):
        zero[(i+2)*n+i+2]+=epsilon*f0[i]
        one[(i+2)*n+1]=epsilon*f1[i]
        lazy=sum((z**y for w,y in zip(state_weights,cw) if w==v),arb(0))/(t*spectrum[v])
        one[(i+2)*n]=epsilon*lazy+eta*f1[i]/m
    one[n]=epsilon*z**min(cw)/t+eta*d1/m
    return tuple(zero),tuple(one),n


def regions(zero,one,n,epochs):
    rz=identity(n);ra=(arb(0),)*(n*n);remaining=epochs
    while remaining:
        if remaining&1:ra=add(mul(ra,zero,n),mul(rz,one,n));rz=mul(rz,zero,n)
        remaining>>=1
        if remaining:one=add(mul(one,zero,n),mul(zero,one,n));zero=mul(zero,zero,n)
    return rz,tuple(x/epochs for x in ra)


def moments(zero,one,n,length):
    def vm(v,m):return tuple(sum((v[k]*m[k*n+j] for k in range(n)),arb(0)) for j in range(n))
    current=[(arb(1),)+(arb(0),)*(n-1)]
    for _ in range(length):
        updated=[(arb(0),)*n for _ in range(len(current)+1)]
        for w,v in enumerate(current):
            updated[w]=add(updated[w],vm(v,zero));updated[w+1]=add(updated[w+1],vm(v,one))
        current=updated
    return [sum(v,arb(0))/math.comb(length,w) for w,v in enumerate(current)]


def pack(value):
    assert value.is_finite()
    m,e=map(int,value.upper().man_exp());assert m>=0
    shift=max(0,m.bit_length()-160)
    if shift:m=(m+(1<<shift)-1)>>shift;e+=shift
    return dict(mantissa=str(m),exponent=e)


def unpack(row):return arb(int(row['mantissa']))*arb(2)**row['exponent']


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();ctx.prec=512 if args.verify else 256
    receipt=HERE/'TRANSVECTION_Q1_CERTIFICATE.json'
    saved=json.loads(receipt.read_text()) if args.verify else None
    if saved:
        assert saved['status']=='OUTWARD_Q1_ONLY_NOT_FULL_SPIN_CERTIFICATE'
        assert saved['outer']==[128,32,32] and saved['message_bits']==1<<20
        assert (saved['t'],saved['s'],saved['transvection_rounds'])==(128,19,1)
        assert saved['covered_occupations']==[1,1] and saved['full_occupation_coverage'] is False
        assert saved['precision_bits']==256 and len(saved['results'])==2
        for name,digest in saved['source_sha256'].items():
            assert hashlib.sha256((tv.ROOT/name).read_bytes()).hexdigest()==digest,name
    spectrum,_,sources=tv.fixed.load_inner()
    path=tv.ROOT/'workstreams/bare_bch_rm2sub/generated/manifest.json'
    columns=json.loads(path.read_text())['t128_s19']['columns']
    counts={w:n for w,n in enumerate(tv.smaller_outer.spectrum()) if w and n}
    construction=tv.smaller_outer.construction()
    assert (construction['length'],construction['dimension'],construction['distance_lower_bound'])==(128,32,32)
    results=[]
    for delta,log_lam,target in ((tv.Fraction(33,200),tv.Fraction(-939,100),40),
                                (tv.Fraction(19,100),tv.Fraction(-966,100),30)):
        lam=(arb(log_lam.numerator)/log_lam.denominator).exp()
        zero,one,n=transfers(spectrum,columns,lam,1)
        mm=moments(*regions(zero,one,n,256),n,128)
        cutoff=4194304*delta.numerator//delta.denominator
        union=32768*(lam*cutoff).exp()*sum((n*mm[w] for w,n in counts.items()),arb(0))
        bound=pack(union)
        assert unpack(bound)<arb(2)**(-target)
        row=dict(distance_target=str(delta),bad_weight=cutoff,log_surprisal=str(log_lam),
                 q1_failure_upper_dyadic=bound,q1_margin_bits=float(-unpack(bound).log()/arb(2).log()),
                 q1_target_bits=target)
        if saved:
            old=saved['results'][len(results)]
            for key in ('distance_target','bad_weight','log_surprisal','q1_target_bits'):assert old[key]==row[key]
            assert union.upper()<=unpack(old['q1_failure_upper_dyadic'])
        results.append(row);print(json.dumps(row),flush=True)
    sources += [Path(__file__),Path(tv.__file__),Path(tv.fixed.__file__),
                Path(tv.smaller_outer.__file__),Path(tv.fixed.maps.__file__),
                Path(tv.fixed.outer.__file__),Path(tv.fixed.outer.bch.__file__),path,
                tv.fixed.HERE/'sources/EBCH128_29.wd',tv.fixed.HERE/'sources/EBCH128_36.wd']
    payload=dict(status='OUTWARD_Q1_ONLY_NOT_FULL_SPIN_CERTIFICATE',precision_bits=ctx.prec,
                 outer=[128,32,32],message_bits=1<<20,t=128,s=19,transvection_rounds=1,
                 covered_occupations=[1,1],full_occupation_coverage=False,
                 source_sha256={p.relative_to(tv.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sources},results=results)
    if args.verify:print('512-bit replay is enclosed by the saved 256-bit Q1 bounds.',flush=True)
    else:receipt.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
