"""FLOAT joint matching/route moment for the structural rank deficit.

Keeps turnover and crossing endpoints in one local polynomial. Selected
row counts only; no whole-code or directed certificate. Serial <=30s.
"""
import math
import shutil
import tempfile
import time
from diagnose_random_pairing_capped import choose2,crossing_logs,minimize_convex,sum2


def profile(R,B,s,field=False):
    M=8192;w=256;L=16*R;T=L-1;W=460*(1<<21)//1000
    A=L*w;z=max(0,A+B-W);I=T*w-B;p=s/M
    assert 0<s<M and z<=I
    ks,law=crossing_logs(512,w)
    # A positive structural gap means C_target <= A-128M-z+128s-1.
    threshold=A-128*M-z+128*s-(0 if field else 1)
    if field:threshold+=(1+choose2(I,z))/32
    base=choose2(M,s)+choose2(512,w)+R*math.log2(M+1)

    def objective(values):
        lv,lt,lx,lu=values
        a0=sum2([math.log2(1-p),math.log2(p)+lt])
        a1=sum2([math.log2(1-p),math.log2(p)+lt+lv])
        lh=a1-a0
        # (x+u)^w (h+u)^k (1+u)^(256-k), coefficient u^w.
        first=w*sum2([lx,lu]);cross=sum2([lh,lu]);empty=sum2([0,lu])
        kernel=sum2([mass+first+k*cross+(256-k)*empty-w*lu
                     for k,mass in zip(ks,law)])
        return (base-threshold*lv-s*R*lt+R*M*a0+(A+B)*lh
                +T*kernel-B*lx)

    best=None
    # Positive-weight proposals only. Never claim a globally optimal bound.
    for initial in ((-5.,3.,-6.,-3.),(-16.,10.,-8.,-6.)):
        values=list(initial);previous=objective(values)
        for sweep in range(100):
            for j,bounds in enumerate((((-32,0) if field else (-64,0)),(-96,96),(-64,64),(-64,64))):
                def at(x):
                    trial=values.copy();trial[j]=x
                    return objective(trial)
                values[j],_=minimize_convex(at,*bounds,iterations=28)
            value=objective(values)
            if abs(value-previous)<.001:break
            previous=value
        record=dict(R=R,B=B,target_rows=s,field_count=field,upper_bits=round(value,3),
                    log_weights=[round(x,7) for x in values],sweeps=sweep+1,float_only=True)
        if best is None or value<best['upper_bits']:best=record
    return best


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic()
    for R,B,s in ((256,131072,1974),(256,32768,1200),(192,131072,4096)):
        for field in (False,True):
            print(profile(R,B,s,field),flush=True)
            if time.monotonic()-began>25:raise TimeoutError('joint geometry budget')
    print(dict(seconds=time.monotonic()-began,whole_code_certificate=False))


if __name__=='__main__':main()
