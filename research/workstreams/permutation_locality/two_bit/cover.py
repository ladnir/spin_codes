"""Complete mixed-support covers for the 4096-pair ensemble.

Only numerical witness selection uses binary64. Every accepted cover is
replayed with exact counts, rational witnesses and outward arithmetic.
No four-bit driver, location constant, or module global is changed.
"""
from collections import Counter
from heapq import heappop,heappush
from math import comb,log

import numpy as np
from flint import arb,arb_mat,arb_poly
from scipy.optimize import minimize,minimize_scalar
from scipy.special import logsumexp

import model
from probe import aq
from shell_cover import IntervalFolds
from occupancy_adaptive import split,volume
from occupancy_cdf_cover import RetainedCover
from occupancy_screen import matrix_for_probabilities
from model import log_power_moment
from measure import DensityFolds
from mixture_cover import checked_threshold

DENOMINATOR=10**9


def outward_occupancy_masses(numerators,denominator=DENOMINATOR):
    """Positive Poisson-binomial coefficients, grouping equal probabilities.

    Arb polynomial multiplication replaces the Python quadratic recurrence.
    Coefficient upper endpoints retain entrywise domination even if the
    fast polynomial product produces unequal interval radii.
    """
    if denominator<1 or any(type(p) is not int or not 0<p<denominator for p in numerators):
        raise ValueError('interior rational probabilities required')
    polynomial=arb_poly([1])
    for numerator,copies in Counter(numerators).items():
        p=arb(numerator)/denominator
        polynomial*=arb_poly([1-p,p])**copies
    masses=[model.up(polynomial[i]) for i in range(len(numerators)+1)]
    if any(x<0 for x in masses):raise ArithmeticError('negative probability upper endpoint')
    return masses


def outward_box(regions,folds,tilt,box,multiplicity,numerators,*,groups=model.GROUPS,output_cutoff=model.THRESHOLD):
    output_cutoff=checked_threshold(output_cutoff)
    q=len(box);n=regions[0].nrows()
    if (not 1<=q<=groups or len(numerators)!=q or multiplicity<1
            or len(regions)<=q or n!=len(model.terminal())):
        raise ValueError('valid labeled support box and location count required')
    if any(type(p) is not int or not 0<p<DENOMINATOR for p in numerators):
        raise ValueError('strictly interior rational support witnesses required')
    masses=outward_occupancy_masses(numerators)
    matrix=sum((mass*region for mass,region in zip(masses,regions)),arb_mat(n,n))**256
    moment=sum((matrix[0,j] for j in range(n) if model.terminal()[j]),arb(0))
    value=model.up(moment*(aq(tilt)*output_cutoff).exp()*comb(groups,q)*multiplicity)
    # Repeated intervals have the same witness. Cache scalar folds locally.
    for (lo,hi,p),copies in Counter((lo,hi,p) for (lo,hi),p in zip(box,numerators)).items():
        value=model.up(value*folds.outward(lo,hi,arb(p)/DENOMINATOR)**copies)
    return value


def check_partition(leaves,q):
    """Replay the exact split tree; volume equality alone is insufficient."""
    if not 1<=q<=model.GROUPS or not leaves:raise ValueError('nonempty occupancy cover required')
    cells={}
    for leaf in leaves:
        box=tuple(tuple(x) for x in leaf['box']);mult=leaf['multiplicity']
        if (len(box)!=q or tuple(sorted(box))!=box or box in cells or type(mult) is not int or mult<1
                or any(len(x)!=2 or any(type(v) is not int for v in x) or not 38<=x[0]<=x[1]<=256 for x in box)):
            raise ValueError('invalid or duplicated support cell')
        cells[box]=mult
    pending=[(((38,256),)*q,1)];used=set();visited=0
    while pending:
        box,mult=pending.pop();visited+=1
        if visited>2*len(cells)-1:raise ValueError('incomplete support partition')
        if box in cells:
            if cells[box]!=mult:raise ValueError('incorrect label multiplicity')
            used.add(box)
        else:
            if all(lo==hi for lo,hi in box):raise ValueError('uncovered support point')
            pending.extend(split(box,mult))
    if used!=set(cells):raise ValueError('overlapping or unreachable support cells')


def replay_cover(operators,cdf,shells,q,leaves,*,density_fold=False,output_cutoff=model.THRESHOLD):
    output_cutoff=checked_threshold(output_cutoff)
    check_partition(leaves,q)
    folds=(DensityFolds if density_fold else IntervalFolds)(cdf,shells);total=arb(0)
    for leaf in leaves:
        if leaf['tilt'] not in operators:raise ValueError('missing local-family witness')
        total=model.up(total+outward_box(operators[leaf['tilt']][0],folds,leaf['tilt'],
                                        leaf['box'],leaf['multiplicity'],leaf['numerators'],output_cutoff=output_cutoff))
    return total


def complete_cover(operators,cdf,shells,q,*,target_bits=55,max_splits=100,joint=False,screen_only=False,certificate=None,density_fold=False,output_cutoff=model.THRESHOLD):
    """Return an outward bound only after covering {38,...,256}^q completely."""
    if not 1<=q<=model.GROUPS or target_bits<1 or max_splits<0:
        raise ValueError('bounded occupancy, positive budget and nonnegative split limit required')
    output_cutoff=checked_threshold(output_cutoff)
    assert all(len(exact)>q and len(floating)>q for exact,floating in operators.values())
    folds=(DensityFolds if density_fold else IntervalFolds)(cdf,shells)
    cache={};terminal=model.terminal()
    def proposal(tilt,interval):
        key=tilt,interval
        if key not in cache:
            region=operators[tilt][1];function=folds.function(*interval)
            def objective(z):
                p=1/(1+np.exp(-z))
                return log_power_moment(matrix_for_probabilities(region,[p]*q),256,terminal)+q*function(p)
            fit=minimize_scalar(objective,bounds=(-10.,18.),method='bounded')
            number=max(1,min(DENOMINATOR-1,round(DENOMINATOR/(1+np.exp(-fit.x)))))
            cache[key]=number,function(number/DENOMINATOR)
        return cache[key]

    def witness(box):
        best=None;candidates=[]
        for tilt,(_,region) in operators.items():
            items=[proposal(tilt,interval) for interval in box];nums=[p for p,_ in items]
            value=(log_power_moment(matrix_for_probabilities(region,[p/DENOMINATOR for p in nums]),256,terminal)
                   +float(tilt)*output_cutoff+sum(value for _,value in items))
            candidate=value,tilt,nums;candidates.append(candidate)
            if best is None or value<best[0]:best=candidate
        if joint:
            intervals=sorted(set(box));indices=[intervals.index(i) for i in box]
            functions=[folds.function(*interval) for interval in intervals]
            for _,tilt,nums in sorted(candidates)[:2]:
                region=operators[tilt][1]
                start=np.array([nums[box.index(interval)]/DENOMINATOR for interval in intervals])
                def objective(logits):
                    ps=1/(1+np.exp(-logits))
                    value=log_power_moment(matrix_for_probabilities(region,[ps[i] for i in indices]),256,terminal)
                    terms=[f(p) for f,p in zip(functions,ps)]
                    return value+float(tilt)*output_cutoff+sum(terms[i] for i in indices)
                fit=minimize(objective,np.log(start/(1-start)),method='L-BFGS-B',bounds=[(-10.,18.)]*len(start),
                             options={'maxiter':40,'ftol':1e-10})
                numbers=[max(1,min(DENOMINATOR-1,round(DENOMINATOR/(1+np.exp(-x))))) for x in fit.x]
                probabilities=np.array(numbers)/DENOMINATOR
                value=objective(np.log(probabilities/(1-probabilities)))
                if value<best[0]:best=value,tilt,[numbers[i] for i in indices]
        if not np.isfinite(best[0]):raise ArithmeticError('nonfinite floating proposal; no certificate')
        return best

    heap=[];tree=RetainedCover();serial=0;locations=comb(model.GROUPS,q)
    def push(box,mult,parent=None):
        nonlocal serial
        value,tilt,nums=witness(box);serial+=1
        item=(-value-log(locations*mult),serial,box,mult,tilt,nums)
        heappush(heap,item);tree.add(item,parent)
    push(((38,256),)*q,1)
    threshold=-(target_bits+2)*log(2);steps=0
    while heap and steps<max_splits and tree.nodes[1]['best']>threshold:
        _,identifier,box,mult,_,_=heappop(heap)
        if all(lo==hi for lo,hi in box):continue
        for child,child_mult in split(box,mult):push(child,child_mult,identifier)
        tree.update(identifier);steps+=1
        if steps%10==0:
            print('TWO-BIT COVER q',q,'splits',steps,'queue',len(heap),
                  'binary64 log2 upper',tree.nodes[1]['best']/log(2),flush=True)
    chosen=tree.selected()
    assert sum(volume(box,mult) for _,_,box,mult,_,_ in chosen)==219**q
    leaves=[dict(box=box,multiplicity=mult,tilt=tilt,numerators=nums) for _,_,box,mult,tilt,nums in chosen]
    check_partition(leaves,q)
    numerical=float(logsumexp([-item[0] for item in chosen]))/log(2)
    print('TWO-BIT COMPLETE DOMAIN q',q,'splits',steps,'selected leaves',len(chosen),
          'binary64 log2 upper',numerical,flush=True)
    if numerical>-(target_bits+2) or screen_only:
        print('No outward occupancy certificate from this cover.',flush=True)
        return None
    total=arb(0)
    for index,(_,_,box,mult,tilt,nums) in enumerate(chosen,1):
        term=outward_box(operators[tilt][0],folds,tilt,box,mult,nums,output_cutoff=output_cutoff)
        total=model.up(total+term)
        if index%100==0:print('TWO-BIT outward leaves',q,index,'/',len(chosen),flush=True)
    if not 0<total<arb(2)**-target_bits:
        raise ArithmeticError('outward occupancy replay failed its budget')
    print('VERIFIED TWO-BIT ALL SUPPORTS q',q,'locations',locations,'upper',total,
          'margin',-total.log()/arb(2).log(),flush=True)
    if certificate is not None:
        certificate.update(q=q,leaves=leaves,density_fold=density_fold,threshold=output_cutoff,upper_dyadic=[int(v) for v in total.upper().man_exp()])
    return total
