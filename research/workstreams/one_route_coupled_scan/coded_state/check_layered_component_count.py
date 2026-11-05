"""Exact rational replay of ONE mixed-history class; no whole-code claim.

No floating point is used in the replay. Logarithms are bounded through
integer powers and bit lengths. Source retains only twelve rational markers.
"""
from fractions import Fraction as F
from itertools import product
import time
from layered_component_transfer import PAIR_TERMS,LAST_TERMS,STATES,coordinate_step
from layered_fixed_rank import path_nullity

MARKERS=[F(n,d) for n,d in (
    (1741164313309419,2147483648),
    (7118972087826913,70368744177664),
    (8117769220996979,2147483648),
    (6801569740989937,2417851639229258349412352),
    (3687855806834923,562949953421312),
    (5904255591893207,70368744177664),
    (8194147899841693,1125899906842624),
    (6332543051645819,70368744177664),
    (5354985777050427,562949953421312),
    (8573329971315407,70368744177664),
    (3928009700023673,2251799813685248),
    (5641738097505231,281474976710656))]


def rational_rows(markers,depth,growth=F(32),zero_factor=F(1,2**27)):
    size=1<<depth;contexts=size//2;mf,me,mb,mw=markers[:4]
    h=[[F(1),markers[4+2*z],markers[5+2*z]] for z in range(size)]
    rows=[[F(0)]*6 for _ in range(size)]
    last=[[F(0)]*6 for _ in range(contexts)]
    for src,dst,gain,full,exc,birth,weight in PAIR_TERMS:
        active=STATES[src//3][0]+STATES[src%3][0]
        zeros=active+birth-weight
        assert zeros>=0
        base=growth**gain*zero_factor**zeros*mf**full*me**exc*mb**birth*mw**weight
        for z in range(size):
            if not z&1 and (dst//3==0 or dst%3==0):continue
            context=z&(contexts-1)
            rows[z][src-3]+=base*h[z][dst//3]*h[z][dst%3]/h[2*context][src//3]/h[2*context+1][src%3]
    for src,gain,full,exc in LAST_TERMS:
        for context in range(contexts):
            last[context][src-3]+=growth**gain*mf**full*me**exc/h[2*context][src//3]/h[2*context+1][src%3]
    return h,rows,last


def log_upper_units(value,power):
    assert value>0
    # log2(n^p)<bit_length(n^p), log2(d^p)>=bit_length(d^p)-1.
    return (value.numerator**power).bit_length()-(value.denominator**power).bit_length()+1


def replay(power=512,profile=(589632,368520,49136),env=None):
    started=time.monotonic();R,M,b,N=192,8192,512,2**21
    f,e,B=profile;W=838860
    P=(R*M-b)//2;V=P+f;D=P+f-2*e;KR=R*M-N//2
    h,rows,last=rational_rows(MARKERS,2)
    # Recompute the prior simultaneous route envelope; do not assume dbar.
    from layered_mixed_gate import envelopes,mixed_dimension
    dbar=mixed_dimension(R,P,0,0,f,e,envelopes(R) if env is None else env)
    if profile==(589632,368520,49136):assert dbar==163923
    factors=[]
    for context in range(2):
        initial=h[2*context][1]*MARKERS[3]*(1+h[2*context+1][1]*MARKERS[3])
        factors.append((initial,b//4))
    factors.extend((max(row),(P-b//2)//4) for row in rows)
    factors.extend((max(row),b//4) for row in last)
    factors.extend((1/marker,count) for marker,count in zip(MARKERS[:4],(f,e,B,W)))
    exponent=1+27*(V+b+dbar)+5*KR-32*(D+e)
    upper=power*exponent+sum(count*log_upper_units(value,power) for value,count in factors)
    upper_ceiling=-((-upper)//power)
    if profile==(589632,368520,49136):assert upper_ceiling<-100
    assert time.monotonic()-started<20
    return dict(power=power,conditional_log2_upper=upper_ceiling,
                R=R,f=f,r=0,e=e,births=B,weight_threshold=W,
                class_closed=upper_ceiling<-100,whole_code_certificate=False)


def neighbors():
    from layered_mixed_gate import envelopes
    env=envelopes(192);worst=None;passed=0
    for df,de,db in product((-1024,0,1024),(-1024,0,1024),(-512,0,512)):
        value=replay(256,(589632+df,368520+de,49136+db),env)
        passed+=value['class_closed']
        if worst is None or value['conditional_log2_upper']>worst['conditional_log2_upper']:worst=value
    return dict(selected_neighbors=27,passed=passed,worst=worst,interval_coverage=False)


def supports(width,bit):
    return [mask for mask in range(1<<width)
            if all(mask>>j&1 for j in range(width) if not j>>bit&1)]


def small_checks():
    started=time.monotonic();width=4;T=2;cases=0;total=F(0)
    # Enumerate original masks directly, independently of the three-state
    # transfer. Two incoming states and a free terminal state.
    for s1,s2 in product(supports(width,1),supports(width,0)):
        supports_now=(s1,s2);full=[];fixed_base={j:0 for j in range(width)}
        for t,mask in enumerate(supports_now,1):
            bit=t%2
            for j in range(width):
                if j>>bit&1:continue
                k=j^(1<<bit)
                if mask>>k&1:full.append((t,j,k))
                else:fixed_base[t*width+j]=1
        inactive={t*width+j for t,mask in enumerate(supports_now) for j in range(width) if not mask>>j&1}
        common=[j for j in range(width) if (s1&s2)>>j&1]
        births=(s2&~s1).bit_count()
        for emask in range(1<<len(full)):
            fixed=fixed_base.copy()
            for k,(t,i,j) in enumerate(full):
                if emask>>k&1:fixed[t*width+i]=0;fixed[t*width+j]=1
            for zmask in range(1<<len(common)):
                ze={width+j for k,j in enumerate(common) if zmask>>k&1}
                kappa=path_nullity(width,T+1,inactive,ze,fixed)
                scanned=0
                for j in range(width):
                    src=STATES.index(((s1>>j)&1,0));sn=(s2>>j)&1
                    gamma=fixed.get(width+j,-1)
                    zero=(width+j in ze) or not ((s1|s2)>>j&1)
                    candidates=[item for item in coordinate_step(src,gamma)
                                if STATES[item[0]][0]==sn and item[3]==int(not zero)]
                    assert len(candidates)==1
                    dst,gain,_,_=candidates[0]
                    scanned+=gain+STATES[dst][1]+(2*width+j not in fixed)
                assert scanned==kappa
                zeros=len(ze);weight=s1.bit_count()+(s1|s2).bit_count()-zeros
                total+=F(2)**kappa*F(1,2)**zeros*F(2)**len(full)*F(1,3)**emask.bit_count()*F(2,3)**births*F(1,2)**weight
                cases+=1
    # Exact positive-potential bound for this complete toy partition sum.
    toy=[F(2),F(1,3),F(2,3),F(1,2),F(2),F(3),F(3),F(4)]
    h,rows,last=rational_rows(toy,1,F(2),F(1,2))
    initial=(h[0][1]*toy[3]*(1+h[1][1]*toy[3]))**2
    bound=initial*max(rows[0])*max(rows[1])*max(last[0])**2
    assert total<=bound
    assert time.monotonic()-started<15
    return dict(exhaustive_histories=cases,component_oracle_agrees=True,
                exact_partition_bound_passes=True,local_rows=12,terminal_rows=6)


if __name__=='__main__':
    print(small_checks())
    print(replay(256))
    print(replay(512))
    print(neighbors())
