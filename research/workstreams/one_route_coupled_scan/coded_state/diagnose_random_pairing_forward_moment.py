"""FLOAT diagnostic: sum input supports and state paths before applying markers.

Constant half-full state after activation; one region-aligned suffix.
Unchanged construction. Not a full certificate. Serial <=30s, no files.
"""
import math
import shutil
import tempfile
import time
from diagnose_random_pairing_capped import choose2, crossing_logs, minimize_convex, sum2


def profile(R,s,B=None):
    M=8192;qlog=32;w=256;T=16*R-1;W=460*(1<<21)//1000;p=s/M
    ks,law=crossing_logs(512,w);lp=math.log2(p)
    base=1+choose2(M,s)-128*s*qlog+R*math.log2(M+1)+choose2(512,w)

    def evaluate(values):
        ly,lu,lt,lx=values;lz=lu-ly-qlog
        zero=sum2([math.log2(1-p),lp-lt]) if p<1 else -lt
        nonzero=lp+qlog-lt
        inputsum=sum2([zero,nonzero])
        # Relax q-1 and q-2 nonzero amplitude counts to q.
        active=sum2([ly+lx,lz,qlog+ly+lz])
        crossed=sum2([0,qlog+ly+lz])
        empty=sum2([zero,nonzero+ly+lz])
        terms=[mass+(w+k)*(inputsum-qlog)+w*active
               +k*crossed+(256-k)*empty-w*lz for k,mass in zip(ks,law)]
        kernel=sum2(terms)
        return (base+w*(nonzero+ly)+(512-w)*zero+T*kernel-W*ly+R*s*lt
                -(B*lx if B is not None else 0))

    values=[-16.,0.,32.,0.];previous=evaluate(values)
    bounds_list=((-64,0),(-96,96),(-64,96))+(((-96,96),) if B is not None else ())
    for sweep in range(120):
        for j,bounds in enumerate(bounds_list):
            def at(x):
                trial=values.copy();trial[j]=x
                return evaluate(trial)
            values[j],_=minimize_convex(at,*bounds,iterations=32)
        value=evaluate(values)
        if abs(value-previous)<0.001:break
        previous=value
    return dict(R=R,s=s,B=B,upper_bits=round(value,3),log_weights=[round(x,6) for x in values],
                sweeps=sweep+1,all_turnover=B is None,float_only=True)


def main():
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--turnover',action='store_true')
    args=parser.parse_args()
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic()
    if args.turnover:
        cases=((256,6144,32768),(256,6218,131072),(192,4096,131072))
    else:
        cases=((R,s,None) for R in (192,256) for s in (256,2048,4096,6144,8192))
    for R,s,B in cases:
        print(profile(R,s,B),flush=True)
        if time.monotonic()-began>25:raise TimeoutError('forward screen budget')
    print(dict(seconds=time.monotonic()-began,whole_code_certificate=False))


if __name__=='__main__':main()
