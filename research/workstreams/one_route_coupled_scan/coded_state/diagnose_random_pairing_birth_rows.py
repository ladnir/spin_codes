"""FLOAT screen: birth compatibility in the exact active-row moment.

Unchanged code distribution. Serial <=30s, no data files.
"""
import math
import shutil
import tempfile
import time
from diagnose_random_pairing_capped import (
    active_profile, choose2, crossing_logs, history_count, minimize_convex, sum2,
)


def birth_profile(R,B,optimistic=False):
    M=8192;T=16*R-1;W=460*(1<<21)//1000
    A=(T+1)*256+B+256;I=T*256-B;z=max(0,(T+1)*256+B-W)
    forced=(2 if optimistic else 1)*B
    allowance=0 if optimistic else .5
    p=A/(R*M);beta=forced/(R*M)
    assert z<=I and 0<=beta<=p<=1
    ks,law=crossing_logs(512,256)
    density=T*(256+max(a-choose2(256,k) for k,a in zip(ks,law)))
    overhead=history_count(T,256,B)+density+2*R*math.log2(M+1)
    coefficient=(1+choose2(I,z))/32-B-z
    baseline=active_profile(R,256,B)

    def objective(values):
        lr,lt,lv=values
        active=(p-beta)*2**(lr-lt)+beta*2**(lr-lt-lv)+1-p
        inactive=(p-beta)*2**(-lt)+allowance*beta*2**(-lt-lv)+(1-p)/2
        return (overhead+coefficient*lr+A*lt+forced*lv
                +M*sum2([R*math.log2(active)-128*lr,R*math.log2(inactive)]))

    values=[baseline['log2_r'],baseline['log2_t'],0.]
    unoptimized=objective(values)
    previous=unoptimized
    for sweep in range(120):
        for j,bounds in enumerate(((0,32),(-64,64),(-64,64))):
            def at(x):
                trial=values.copy();trial[j]=x
                return objective(trial)
            values[j],_=minimize_convex(at,*bounds,iterations=32)
        value=objective(values)
        if abs(value-previous)<0.001:break
        previous=value
    return dict(R=R,B=B,optimistic_not_code_bound=optimistic,baseline_bits=baseline['upper_bits'],
                birth_no_tuning_bits=round(unoptimized,3),upper_bits=round(value,3),
                log_weights=[round(x,7) for x in values],sweeps=sweep+1,
                float_only=True)


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--optimistic',action='store_true')
    args=parser.parse_args()
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic()
    for R,B in ((256,32768),(256,65536),(256,131072),(192,65536),(192,131072)):
        print(birth_profile(R,B,args.optimistic),flush=True)
        if time.monotonic()-began>25:raise TimeoutError('birth-row screen budget')
    print(dict(seconds=time.monotonic()-began,whole_code_certificate=False))


if __name__=='__main__':main()
