"""Joint incoming/outgoing population diagnostic for closed paid sets.

Retains the two population vectors and original inactive count. Maximizes
run ages coefficientwise, so this remains an upper relaxation. FLOAT only.
Serial Peach180s, <256MiB, no files or certificate claims.
"""
import itertools
import math
import time
import numpy as np
from quartet_closed_transfer import TYPES,local
from quartet_closed_population import multinomial
from quartet_two_cohort import compositions


def caps(logb=2.5,logh=None):
    h=dict(zip(TYPES,np.zeros(len(TYPES)) if logh is None else logh))
    result={}
    for hist in itertools.combinations_with_replacement(TYPES,4):
        incoming=tuple(sum(t%3==i for t in hist) for i in range(3));row={}
        for (dest,unpaid,covered,activity,deaths),ways in local(hist):
            outgoing=(unpaid,4-unpaid-activity,activity)
            value=math.log2(ways)+logb*(707*activity-1905*deaths)/2560
            value+=sum(h[t] for t in dest)-sum(h[t] for t in hist)
            for paid in range(covered+1):
                key=(incoming,outgoing,unpaid+paid)
                term=value+math.log2(math.comb(covered,paid))-logb*6*paid/5
                row[key]=np.logaddexp2(row.get(key,-np.inf),term)
        for key,value in row.items():result[key]=max(result.get(key,-np.inf),value)
    return result


def blocked_caps(logb=2.5):
    """Keep pending coordinates that forbid a full next input quartet.

    Categories: U0; ordinary U1,S0; blocked U1,S0; S1.
    This is a projection of the same twelve types, not a new construction.
    """
    def category(t):
        if t%3==0:return 0
        if t%3==2:return 3
        return 2 if t in (7,16) else 1
    result={}
    for hist in itertools.combinations_with_replacement(TYPES,4):
        incoming=tuple(sum(category(t)==i for t in hist) for i in range(4));row={}
        for (dest,unpaid,covered,activity,deaths),ways in local(hist):
            outgoing=tuple(sum(category(t)==i for t in dest) for i in range(4))
            value=math.log2(ways)+logb*(707*activity-1905*deaths)/2560
            for paid in range(covered+1):
                key=(incoming,outgoing,unpaid+paid)
                term=value+math.log2(math.comb(covered,paid))-logb*6*paid/5
                row[key]=np.logaddexp2(row.get(key,-np.inf),term)
        for key,value in row.items():result[key]=max(result.get(key,-np.inf),value)
    return result


def bound(table,incoming,outgoing,inactive,iterations=120,replay=0.,return_fit=False):
    """Cauchy in both populations, positive extraction in original inactivity.

    Fixed zero population categories are eliminated, not approximated with
    large finite markers. Only the incoming multinomial is a denominator.
    """
    width=sum(incoming);assert sum(outgoing)==width and width%4==0
    groups=width//4
    oldactive=[i for i,n in enumerate(incoming) if n]
    newactive=[i for i,n in enumerate(outgoing) if n]
    oldkeep=[i for i in oldactive if i!=max(oldactive,key=lambda i:incoming[i])]
    newkeep=[i for i in newactive if i!=max(newactive,key=lambda i:outgoing[i])]
    interior=0<inactive<width
    base=[];powers=[];js=[];keys=[]
    for (old,new,j),value in table.items():
        if any(old[i] for i in range(len(incoming)) if i not in oldactive):continue
        if any(new[i] for i in range(len(outgoing)) if i not in newactive):continue
        if inactive==0 and j!=0 or inactive==width and j!=4:continue
        base.append(value+math.log2(multinomial(old)))
        powers.append([old[i] for i in oldkeep]+[new[i] for i in newkeep]+([j] if interior else []));js.append(j)
        keys.append((old,new,j))
    if not base:return -math.inf,-math.inf,0.
    base=np.array(base);powers=np.array(powers,dtype=float);js=np.array(js)
    target=np.array([incoming[i]/groups for i in oldkeep]+[outgoing[i]/groups for i in newkeep]+([inactive/groups] if interior else []))
    x=np.zeros(len(target))
    def evaluate(x):
        terms=base+powers@x;peak=max(terms);raw=np.exp2(terms-peak);p=raw/sum(raw)
        total=peak+math.log2(sum(raw));mean=p@powers;centered=powers-mean
        return total,p,mean,(centered.T*p)@centered
    for _ in range(iterations):
        total,p,mean,cov=evaluate(x)
        if not len(x) or max(abs(mean-target))<1e-10:break
        direction=np.linalg.solve(math.log(2)*cov+1e-10*np.eye(len(x)),target-mean)
        direction/=max(1.,max(abs(direction))/4);rate=1.;old=total-target@x
        for _ in range(25):
            proposal=np.clip(x+rate*direction,-256,256)
            if evaluate(proposal)[0]-target@proposal<=old+1e-11:break
            rate*=.5
        x=proposal
    if interior:x[-1]+=replay
    total,p,mean,cov=evaluate(x);q=np.bincount(js,weights=p,minlength=5)
    current=np.zeros(inactive+1);current[0]=1.
    for _ in range(groups):
        after=np.zeros_like(current)
        for j,weight in enumerate(q):
            if j<=inactive:after[j:]+=weight*current[:inactive+1-j]
        current=after
    cauchy=groups*(total-target@x)-math.log2(multinomial(incoming))
    if return_fit:
        return dict(keys=keys,probabilities=p,cauchy=cauchy,groups=groups,
                    powers=powers,target=target,tilts=x,base=base,
                    oldkeep=oldkeep,newkeep=newkeep,interior=interior,
                    maximum_residual=float(max(abs(mean-target),default=0.)))
    # Underflow never provides a purported improvement.
    mass=float(current[inactive]);logmass=math.log2(mass) if mass>0 else 0.
    return cauchy+logmass,cauchy,float(max(abs(mean-target),default=0.))


def tests():
    # Exact integer coefficients at width8, compared with FLOAT bounds.
    cats=list(compositions(4,3));table={}
    for old in cats:
        for new in cats:
            for j in range(5):table[old,new,j]=1+old[1]+2*new[2]+j
    logs={k:math.log2(v) for k,v in table.items()};checked=0
    pops=((8,0,0),(0,0,8),(4,0,4),(2,3,3),(0,4,4))
    for incoming,outgoing in itertools.product(pops,repeat=2):
        for inactive in (0,3,8):
            value=0
            for (old,new,j),a in table.items():
                otherold=tuple(n-m for n,m in zip(incoming,old))
                othernew=tuple(n-m for n,m in zip(outgoing,new))
                other=(otherold,othernew,inactive-j)
                if other not in table:continue
                value+=multinomial(old)*multinomial(otherold)*a*table[other]
            exact=value/multinomial(incoming)
            upper,cauchy,residual=bound(logs,incoming,outgoing,inactive)
            assert math.log2(exact)<=upper+1e-9
            if inactive==3:
                assert abs(upper-bound(logs,incoming,outgoing,inactive,replay=.2)[0])<1e-9
            checked+=1
    print(dict(joint_population_coefficients_checked=checked),flush=True)


def main():
    tests();start=time.monotonic();table=caps();print(dict(local_joint_coefficients=len(table)),flush=True)
    for a,s in ((381,381),(384,384),(408,408),(448,448),(480,480),(512,512),(384,128),(408,128)):
        incoming=(512-a,a-s,s);outgoing=incoming
        values=[bound(table,incoming,outgoing,i) for i in (131,132)]
        regional=12*values[0][0]+4*values[1][0]+math.log2(math.comb(16,4))-math.log2(math.comb(8192,2100))
        print(dict(allowed=a,selected=s,log2_conditional_coefficients=[v[0] for v in values],
                   conditional_region_proposal=regional,max_fit_residual=max(v[2] for v in values)),flush=True)
    print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)


if __name__=='__main__':main()
