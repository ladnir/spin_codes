"""IID screening of an explicit equal-split flow, not a distance bound.

Each target gate assigns half its service to each incoming coordinate
timeline. Exact timeline min-cut recurrence, FLOAT IID pattern weights.
Actual regional conditioning and correlations are not discharged here.
Serial <=30s, standard library only, no files.
"""
import shutil
import tempfile
import time
from itertools import product


def advance(d,service,supply,cap):
    # d=F1-F0, while the accumulated value is min(F0,F1).
    normalization=max(0,-d)
    g0=normalization+min(0,d+cap)+supply
    g1=normalization+min(d,cap)+service
    return g1-g0,min(g0,g1)


def self_test():
    checks=0
    for supply,cap in ((1,2),(3,4)):
        for services in product((0,cap,2*cap),repeat=5):
            d=0;flow=0
            for service in services:
                d,increment=advance(d,service,supply,cap);flow+=increment
            exact=min(sum(service*((mask>>j)&1)+supply*(1-((mask>>j)&1))
                          for j,service in enumerate(services))
                      +cap*sum(((mask>>j)^(mask>>(j-1)))&1 for j in range(1,5))
                      for mask in range(32))
            assert flow==exact and flow>=0
            checks+=1
    return checks


def screen(supply,scale,p,length):
    assert scale%2==0
    cap=scale//2
    services=((0,(1-p)**2),(cap,2*p*(1-p)),(scale,p*p))
    law={0:1.};expected=0.
    for step in range(length):
        new={}
        # The last coordinate has no next gate inside this window.
        choices=services if step<length-1 else ((0,1.),)
        for d,mass in law.items():
            for service,probability in choices:
                dn,increment=advance(d,service,supply,cap)
                expected+=mass*probability*increment
                new[dn]=new.get(dn,0.)+mass*probability
        law=new
    assert abs(sum(law.values())-1)<1e-9
    lower_mean=(expected-supply*length/2)/scale
    return dict(supply=supply,scale=scale,length=length,mean_flow=expected/scale,
                mean_dual_lower_per_coordinate=lower_mean,
                extrapolated_full_cost=lower_mean/length*(1<<21),
                iid_screen_only=True)


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();p=1974/8192
    print(dict(exact_timeline_checks=self_test()),flush=True)
    for supply,scale in ((3,16),(4,16),(17,64)):
        for length in (32,256):print(screen(supply,scale,p,length),flush=True)
    assert time.monotonic()-began<25
    print(dict(seconds=time.monotonic()-began,no_setup_probability_or_distance_certificate=True))


if __name__=='__main__':main()
