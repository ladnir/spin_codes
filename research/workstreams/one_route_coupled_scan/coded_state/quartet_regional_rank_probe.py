"""Regional global-rank diagnostics and an exact rational burst certificate.

The floating regional operators are NOT certificates. --burst-cert proves
only the full/empty-state class, not minimum distance for the whole code.

Keep the fixed quartet construction. The root runs this file-free on Peach,
single-threaded, with a180s watchdog and2GiB address-space bound. No dumps.
The kernel sums every destination; it carries only boundary state weight.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, log2, isqrt
import os
import resource
import shutil
import sys
import threading
import time
import numpy as np

try:
    reference_transfer
except NameError:
    from quartet_global_coverage_probe import transfer_row as reference_transfer

LD=np.longdouble
LN2=np.log(LD(2))
N=1<<21;K=1<<20;W=943718


def quartet_local(theta,logz,logx,logh,field_bits=32):
    z=np.exp2(LD(logz));a=np.exp2(LD(field_bits)*theta);m=a*z
    c=np.exp2(theta*np.log2(1+np.exp2(LD(field_bits)+LD(logz)/theta)))
    x=np.exp2(LD(logx));h=np.exp2(LD(logh))
    local=np.zeros((5,5),dtype=LD)
    for u in range(5):
        for death in range(u+1):
            survived=u-death
            for birth in range(5-u):
                v=survived+birth;covered=4 if u else birth
                local[u,v]+=LD(comb(u,death)*comb(4-u,birth))*z**death*m**birth*c**survived*x**covered*h**(v-u)
    return local


def update_kernel(width,theta,logz,logx,logh,field_bits=32):
    """T_uv*h^(v-u) using a hypergeometric block recurrence.

    Incoming-support normalization is applied during convolution, avoiding
    the large binomial coefficients of the raw bivariate polynomial.
    """
    assert width%4==0
    local=quartet_local(theta,logz,logx,logh,field_bits)
    old=np.ones((1,1),dtype=LD)
    for slots in range(4,width+1,4):
        size=len(old);new=np.zeros((slots+1,slots+1),dtype=LD)
        denominator=LD(comb(slots,4))
        for incoming in range(5):
            weights=np.array([LD(comb(u,incoming)*comb(slots-u,4-incoming))/denominator
                              for u in range(incoming,incoming+size)],dtype=LD)
            weighted=old*weights[:,None]
            for outgoing in range(5):
                new[incoming:incoming+size,outgoing:outgoing+size]+=weighted*local[incoming,outgoing]
        old=new
    assert abs(old[0,0]-1)<LD('1e-15')
    return old


def positive_power16(kernel):
    """Compute T^16-E00 without subtracting nearly equal numbers.

    Q_(2n)=Q_n^2+E00 Q_n+Q_n E00. Return Q scaled by a power of two.
    """
    positive=kernel.copy();positive[0,0]=0
    scale=LD(0)
    for _ in range(4):
        norm=positive.sum(axis=1).max()
        if norm:
            positive/=norm;scale+=np.log2(norm)
        # Choose a common scale for quadratic and linear terms.
        common=max(2*scale,scale)
        squared=(positive@positive)*np.exp2(2*scale-common)
        linear=np.exp2(scale-common)
        squared[0,:]+=linear*positive[0,:]
        squared[:,0]+=linear*positive[:,0]
        positive=squared;scale=common
    norm=positive.sum(axis=1).max()
    positive/=norm;scale+=np.log2(norm)
    return positive,scale


def occupied_tail(positive,scale,regions=256,minimum=129):
    """Exact scalar renewal identity evaluated in long double.

    Coefficient with r positive-footprint regions is
    [x^r] A(x)^(regions-r) B(x), A_k=(Q^k)00,B_k=e0 Q^k1.
    A separate positive scale contributes scale*r to its logarithm.
    """
    row=np.zeros(len(positive),dtype=LD);row[0]=1
    aa=np.empty(regions+1,dtype=LD);bb=aa.copy();row_scale=LD(0)
    for i in range(regions+1):
        aa[i]=(np.log2(row[0])+row_scale if row[0] else LD('-inf'))
        bb[i]=np.log2(row.sum())+row_scale
        if i<regions:
            row=row@positive;norm=row.sum()
            assert norm>0
            row/=norm;row_scale+=np.log2(norm)
    poly=bb.copy();logs=[]
    for zero_count in range(regions-minimum+1):
        occupied=regions-zero_count;logs.append(poly[occupied]+scale*occupied)
        # Log-domain convolution: ordinary coefficients can be <2^-100000.
        new=np.full(regions+1,LD('-inf'))
        for j in range(regions+1):
            new[j:]=np.logaddexp(new[j:]*LN2,(poly[:regions+1-j]+aa[j])*LN2)/LN2
        poly=new
    if not logs:return LD('-inf')
    biggest=max(logs)
    return biggest+np.log2(sum(np.exp2(v-biggest) for v in logs))


def exact_renewal_checks():
    # Compare the scalar identity with an explicit positive-region counter.
    checks=0
    for size in (2,3):
        q=[[F((i+1)*(j+2),11) for j in range(size)] for i in range(size)]
        a=[];b=[];row=[F(int(i==0)) for i in range(size)]
        regions=6
        for _ in range(regions+1):
            a.append(row[0]);b.append(sum(row))
            row=[sum(row[i]*q[i][j] for i in range(size)) for j in range(size)]
        paths={(0,0):F(1)}
        for _ in range(regions):
            new={}
            for (count,u),weight in paths.items():
                if u==0:new[count,0]=new.get((count,0),F(0))+weight
                for v in range(size):new[count+1,v]=new.get((count+1,v),F(0))+weight*q[u][v]
            paths=new
        poly=b.copy()
        for zero_count in range(regions+1):
            r=regions-zero_count
            assert poly[r]==sum(weight for (count,u),weight in paths.items() if count==r)
            checks+=1
            poly=[sum(poly[j]*a[k-j] for j in range(k+1)) for k in range(regions+1)]
        matrix=np.array([[LD(v.numerator)/LD(v.denominator) for v in r] for r in q])
        for minimum in (1,3,6):
            expected=sum(weight for (count,u),weight in paths.items() if count>=minimum)
            truth=np.log2(LD(expected.numerator)/LD(expected.denominator))
            assert abs(occupied_tail(matrix,0,regions,minimum)-truth)<LD('1e-14')
            checks+=1
    print(dict(exact_renewal_comparisons=checks),flush=True)


def small_kernel_checks():
    # Independent enumerator from the prior helper is injected by the loader.
    checks=0
    for width in (4,8):
        matrix=update_kernel(width,LD('0.5'),np.log2(LD(3)/16),np.log2(LD(5)/6),0,4)
        for incoming in range(width+1):
            expected=reference_transfer(width,incoming,F(4),F(3,16),F(5,6),F(5,4))
            for outgoing in range(width+1):
                truth=LD(expected[outgoing].numerator)/LD(expected[outgoing].denominator)
                assert abs(matrix[incoming,outgoing]-truth)<=LD('1e-15')*max(1,truth)
                checks+=1
        positive,scale=positive_power16(matrix)
        direct=np.linalg.matrix_power(matrix,16);direct[0,0]-=1
        assert np.max(abs(positive*np.exp2(scale)-direct))<LD('1e-12')*max(1,direct.max())
    print(dict(small_kernel_and_reference_comparisons=checks),flush=True)


def regional_diagnostic():
    theta=LD('0.25');logz=-8.;logh=-2.
    # Only a bounded first grid. Entrywise minima are valid upper envelopes.
    etas=[LD(0),LD(1)/128,LD(1)/64,LD(1)/32,LD(1)/16]
    matrices=[]
    for eta in etas:
        start=time.monotonic()
        update=update_kernel(512,theta,logz,-32*theta*eta,logh)
        positive,scale=positive_power16(update)
        zeros=int(np.count_nonzero(positive==0))
        # Upward floor makes underflow pessimistic; this remains undirected.
        positive=np.maximum(positive,np.exp2(LD(-15000)))
        matrices.append((eta,np.log2(positive)+scale))
        print(dict(kernel_eta=float(eta),matrix_log2_scale=float(scale),
                   floored_entries=zeros,seconds=round(time.monotonic()-start,3)),flush=True)
    for median in (512,2048,4096,6144):
        logs=np.minimum.reduce([logmatrix+32*theta*eta*median for eta,logmatrix in matrices])
        scale=logs.max();upper=np.exp2(logs-scale)
        norm=upper.sum(axis=1).max();upper/=norm;scale+=np.log2(norm)
        tail=occupied_tail(upper,scale)
        logbound=theta-32*theta-W*logz-32*theta*128*median+tail+max(0.,-512*logh)
        print(dict(median_coverage=median,all_history_log2_diagnostic=float(logbound),
                   minimum_occupied_regions=129,certificate=False),flush=True)


def adaptive_diagnostic():
    # Markers vary by the actual median-coverage class. No cert claim.
    theta=LD('0.25');logh=-2.
    for median in (512,2048,4096,6144):
        logz=-32*float(theta)*128*median/(N-W)-0.15
        logs=None;start=time.monotonic();floored=0
        for eta in (LD(0),LD(1)/4,LD(1)/2,LD(3)/4,LD(1)):
            update=update_kernel(512,theta,logz,-32*theta*eta,logh)
            positive,scale=positive_power16(update)
            floored+=int(np.count_nonzero(positive==0))
            entry=np.log2(np.maximum(positive,np.exp2(LD(-15000))))+scale+32*theta*eta*median
            logs=entry if logs is None else np.minimum(logs,entry)
        scale=logs.max();upper=np.exp2(logs-scale)
        norm=upper.sum(axis=1).max();upper/=norm;scale+=np.log2(norm)
        tail=occupied_tail(upper,scale)
        logbound=theta-32*theta-W*logz-32*theta*128*median+tail+max(0.,-512*logh)
        print(dict(adaptive_median=median,logz=logz,log2_diagnostic=float(logbound),
                   floored_entries=floored,seconds=round(time.monotonic()-start,3),
                   certificate=False),flush=True)


def burst_rational_certificate():
    """All empty/full-state histories, all row counts. Exact rational checks.

    Proof uses t=1/10, grouped zero constraints, row negative association,
    and the number2r+1 of possible partial regions. No empirical spectrum.
    """
    regions=256;rows=8192;steps=4096;t=F(1,10);mu0=F(576,5)
    c=F(13,2500);exp_c=F(7,5)
    ln2_upper=F(7,10)
    assert 1+ln2_upper+ln2_upper**2/2+ln2_upper**3/6+ln2_upper**4/24>2
    assert K-W-N//32==39322
    assert t*t/(2*(1-t/3))<c
    v=64*c
    assert 1+v+v*v/(2*(1-v/3))<exp_c
    assert 1+t*F(64,5)+(t*F(64,5))**2/2>3
    low_c=t/(4*F(64,5))+c*exp_c/3
    end_c=t/16+c*exp_c
    worst_low=(F(-10**9),None);worst_high=(F(-10**9),None)
    checked=0
    for deaths in range(2049):
        count=comb(steps+1,2*deaths+1)
        count_bits=(count-1).bit_length()
        entropy=F(7,10)*count_bits  # ln2<7/10.
        variance=min(F(64),F(deaths+1,2))
        if deaths>=3:
            endpoint=entropy-t*(512*deaths+F(11,32)*rows*(124-mu0))+rows*variance*end_c
            if deaths>=7:
                start=entropy-t*512*deaths+rows*variance*low_c
            else:
                mu=128-2*deaths;gap=128-mu
                start=(entropy-t*(512*deaths+F(11,32)*rows*(mu-mu0))
                       +rows*variance*(t/(4*gap)+c*exp_c))
            value=max(start,endpoint)
            assert value<F(-112)
            if value>worst_low[0]:worst_low=(value,deaths)
            checked+=1
        # For coverage >=124M, remove the deterministic excess when mean>128.
        # The selected-zero entropy costs at most I/32 field dimensions.
        # The remaining field deficit is at least512r+39322.
        scaled=variance*256
        root=isqrt(scaled.numerator//scaled.denominator)
        if F(root*root)<scaled:root+=1
        sqrt_upper=F(root,16)
        assert sqrt_upper*sqrt_upper>=variance
        value=entropy-t*(512*deaths+39322)+rows*(variance*c+t*sqrt_upper/2)
        assert value<F(-112)
        if value>worst_high[0]:worst_high=(value,deaths)
        checked+=1
    assert checked==4095 and checked<1<<12
    # Each branch is <exp(-112)<2^-160; their union is <2^-148.
    print(dict(exact_all_burst_branches=checked,all_burst_failure_below_2_to=-148,
               worst_low_r=worst_low[1],worst_low_natural_log_upper=str(worst_low[0]),
               worst_high_r=worst_high[1],worst_high_natural_log_upper=str(worst_high[0]),
               full_code_certificate=False),flush=True)


def footprint_burst_certificate():
    """Extend the same scalar certificate to binary covered footprints.

    Intermediate states can be arbitrary. Run counts are grouped before
    the fractional moment; see QUARTET_GLOBAL_COVERAGE_RUNS.md.
    """
    checks=0
    for q in (2,16):
        for x in (F(1,32),F(1,16),F(1,2),F(1)):
            lam=1+(q+1)*x
            # First covered update is a forced nonzero state.
            row=[F(0),q*x]
            for length in range(1,17):
                assert row[0]<=x*lam**(length-1)
                assert sum(row)<=(q+1)*x*lam**(length-1)
                checks+=2
                row=[row[0]+x*row[1],q*x*row[0]+(1+q*x)*row[1]]
    q=1<<32
    # e<3 by the positive geometric Taylor tail, hence gamma*e^2<18<q.
    assert 1+1+F(1,2)/(1-F(1,3))<3
    assert 18<q
    assert F(2*N,q)<2 and F(N,N)+F(N,q)<2
    assert F(N,q)<2
    # Preserve 2^F in the below-K high-coverage branch.
    assert F(31,32)*(124*8192)-W==40346
    assert K-W-N//32==39322
    burst_rational_certificate()
    print(dict(exact_run_generating_checks=checks,
               covered_footprint_class='each update empty or full; arbitrary intermediate states',
               all_footprint_burst_failure_below_2_to=-148,
               full_code_certificate=False),flush=True)


def main():
    if shutil.disk_usage('/tmp').free<20*(1<<30):raise RuntimeError('20GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(2*(1<<30),2*(1<<30)))
    timer=threading.Timer(175,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    if '--burst-cert' in sys.argv:
        burst_rational_certificate();timer.cancel();return
    if '--footprint-burst-cert' in sys.argv:
        footprint_burst_certificate();timer.cancel();return
    exact_renewal_checks();small_kernel_checks()
    if '--adaptive' in sys.argv:adaptive_diagnostic()
    else:regional_diagnostic()
    print(dict(seconds=round(time.monotonic()-start,3),full_certificate=False),flush=True)
    timer.cancel()


if __name__=='__main__':main()
