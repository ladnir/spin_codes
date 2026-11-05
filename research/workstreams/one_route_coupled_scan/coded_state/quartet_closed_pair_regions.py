"""Exact small regional composition of closed-set two-update kernels.

Carries seven pending types and extracts the original-inactive quota every
region; applies the terminal test only once. Explicit negative controls
detect replacing regional quotas by one total or inserting an extra terminal.
Serial, <=60s, small integer dictionaries, no files. Not a full-size cert.
"""
from collections import defaultdict
from functools import lru_cache
import math
import os
from pathlib import Path
import shutil
import threading
from quartet_closed_quotient import kernel,accepted
from quartet_closed_paid_sets import tests as closed_tests


@lru_cache(maxsize=210)
def pair(hist):
    current={(hist,0,0,0):1}
    for _ in range(2):
        after=defaultdict(int)
        for (src,inactive,paid,score),count in current.items():
            for (dest,u,v,a,d),ways in kernel(src):
                for c in range(v+1):
                    after[dest,inactive+u+c,paid+c,score+2*a-5*d]+=count*ways*math.comb(v,c)
        current=dict(after)
    return current


def compose(regional=True,extra_terminal=False):
    current={((0,0,0,0),0,0,0):1};peak=1
    for region in range(2):
        after=defaultdict(int)
        for (src,total,paid,score),count in current.items():
            for (dest,i,c,e),ways in pair(src).items():
                if regional and i!=2:continue
                if total+i>4:continue
                if extra_terminal and region==0 and not accepted(dest):continue
                after[dest,total+i,paid+c,score+e]+=count*ways
        current=dict(after);peak=max(peak,len(current))
    result=defaultdict(int)
    for (hist,total,paid,score),count in current.items():
        if total==4 and accepted(hist):result[paid,score]+=count
    return dict(result),peak


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    try:
        expected=closed_tests();actual,peak=compose();assert actual==expected
        global_only,_=compose(False);early,_=compose(True,True)
        assert global_only!=actual and early!=actual
        assert sum(v for (p,s),v in actual.items() if p==0)==math.comb(8,2)**2
        assert sum(v for (p,s),v in global_only.items() if p==0)==math.comb(16,4)
        print(dict(exact_pair_regional_coefficients=len(actual),peak_states=peak,
                   global_quota_negative_control=True,extra_terminal_negative_control=True,
                   exact_pending_boundary_types=True),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
