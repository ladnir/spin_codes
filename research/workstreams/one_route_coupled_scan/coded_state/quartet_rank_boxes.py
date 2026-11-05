"""Joint canonical-profile boxes; injected lowrank helper supplies transfers.

All large results are FLOAT diagnostics, never a full certificate.
The geometry quantity counts footprints, not the probability of a geometry event.
"""


def conditional_field_log(rank, available, forced, weight=W, terminal=512, q=1<<32):
    """Optimize the forced-start run bound exactly in x, then cap at one.

    gamma<2 is enlarged to 2. Profile endpoints must be canonical integers.
    The endpoint formula assumes weight<q (true for the production target).
    """
    assert weight<q and 0<=forced<=min(available,weight)
    if 2*forced-terminal>weight:return LD('-inf')
    d=LD(q)+1; ld=np.log2(d)
    candidate=forced*ld+(available-forced)*np.log2(LD(q)+2)
    if forced<=terminal:
        n=available-forced;k=weight-forced;constant=weight
    else:
        n=available-terminal;k=weight-forced;constant=k+terminal
    if n>k>=0:
        def term(v):return LD(v)*np.log2(LD(v)) if v else LD(0)
        candidate=min(candidate,constant*ld+term(n)-term(k)-term(n-k))
    return min(LD(0),1-np.log2(LD(q))-rank*np.log2(LD(q))+candidate)


def geometry_logmoment(h, loga, logb, logu, logcoverage=0., tilts=9):
    """E sum_V a^-R b^L u^f chi^F, with one choice of eligible H.

    The existing positive relaxation also adds formal annotations. Canonical
    coefficients have the stated meaning; added ones are not code statistics.
    """
    rows=8192;actual=LD(h)/rows;a=LD(2)**loga;b=LD(2)**logb;m=b*LD(2)**logu
    choose=(lgamma(rows+1)-lgamma(h+1)-lgamma(rows-h+1))/float(np.log(LD(2)))
    low=actual/(a-(a-1)*actual);logs=None;endlogs=None
    for fraction in np.linspace(0.,1.,tilts):
        odds=(low/(1-low))**(1-fraction)*(actual/(1-actual))**fraction;p=odds/(1+odds)
        logpi=choose+h*np.log2(p)+(rows-h)*np.log2(1-p)
        C,terminal=complete_force_matrix(512,m,b,p,p+(1-p)/a,LD(2)**logcoverage)
        Q,scale,end,end_scale=complete_force_region(C,terminal)
        if not np.all(Q>0) or not np.all(end>0):raise ArithmeticError('Positive entry underflow')
        candidate=np.log2(Q)+scale-logpi;candidate_end=np.log2(end)+end_scale-logpi
        logs=candidate if logs is None else np.minimum(logs,candidate)
        endlogs=candidate_end if endlogs is None else np.minimum(endlogs,candidate_end)
    scale=logs.max();Q=np.exp2(logs-scale);norm=Q.sum(axis=1).max();Q/=norm;scale+=np.log2(norm)
    end_scale=endlogs.max();end=np.exp2(endlogs-end_scale)
    return choose-loga*128*h+complete_force_tail(Q,scale,end,end_scale)


def box_diagnostic(moment, markers, box):
    """Box order R,L,f,F. Geometry must not be capped before field multiplication."""
    loga,logb,logu,logchi=markers
    minimum=sum(coef*(bounds[0] if coef>=0 else bounds[1])
                for coef,bounds in zip((-loga,logb,logu,logchi),box))
    count=moment-minimum
    field=conditional_field_log(box[0][0],box[1][1],box[2][0])
    return dict(log2_geometric_count=float(count),log2_field=float(field),
                log2_product=float(count+field),certificate=False)


def box_checks():
    # Compare the algebraic optimizer to its minimizing x and a dense x grid.
    checks=0
    for available in range(1,13):
        for forced in range(min(available,5)+1):
            if 2*forced-2>5:continue
            n=available-min(forced,2);k=5-forced
            candidates=[LD(1)]
            if n>k>0:candidates.append(LD(k)/(17*(n-k)))
            if k==0 and n>0:candidates.append(LD('1e-60'))
            def direct(x):
                lam=1+17*x;e=min(LD(1),lam/17)
                return min(LD(0),1-4-4*9-5*np.log2(x)
                           +forced*np.log2(17*x)+(available-forced)*np.log2(lam)
                           +max(0,forced-2)*np.log2(e))
            exact=conditional_field_log(9,available,forced,5,2,16)
            assert abs(exact-min(map(direct,candidates)))<LD('1e-16')
            assert exact<=min(direct(LD(2)**(-j/4)) for j in range(161))+LD('1e-16')
            checks+=1
    # Recover the same unboxed bound from independent moment and field markers.
    h=2048;theta=LD(1)/16;logz=LD(-7)/16;x=LD(2)**(logz/theta)
    logb=theta*np.log2(1+((1<<32)+1)*x);logu=logz
    moment=geometry_logmoment(h,32*theta,logb,logu,tilts=3)
    run_log=logb-theta*np.log2(LD((1<<32)+1))
    prior=complete_force_value(h,float(theta),float(logz),tilts=3,birth_run_bonus=True)
    reproduced=theta-32*theta-W*logz-512*run_log+moment
    assert abs(prior-reproduced)<LD('1e-7')
    print(dict(optimized_field_checks=checks,moment_recomposition=True,certificate=False),flush=True)


def closure_class_certificate():
    """Exact rational constants for the all-footprint 2^-76 class theorem."""
    q=1<<32;root=1<<16;x=F(1,root);lam=1+(q+1)*x
    e=min(F(1),lam/(q+1));d=max(e,(1+q*x*x)/lam)
    assert d<F(1,1<<15)
    # euler number < sum_{j=0}^5 1/j! + (1/6!)/(1-1/7).
    e_upper=F(163,60)+F(7,4320)
    assert e_upper<F(87,32)
    assert F(87,32)*F(65536,65535)<F(11,4)
    assert 11**32<2**111   # (11/4)^32 < 2^47.
    assert N*(F(1,root)+F(1,q))==32+F(1,2048)
    threshold=16*W+N+92
    assert N+16+16*W-threshold==-76
    print(dict(closure_class_threshold=threshold,remainder_cap=threshold-1,
               strict_failure_upper='2^-76',exact_rational_constant_checks=True,
               full_code_certificate=False),flush=True)


def box_probe():
    start=time.monotonic();h=2048;theta=LD(1)/16;logz=LD(-7)/16
    x=LD(2)**(logz/theta);logb=theta*np.log2(1+((1<<32)+1)*x)
    markers=(32*theta,logb,logz,LD(0))
    moment=geometry_logmoment(h,*markers,tilts=9)
    profiles=((903775,891035,71426,1132763), (900000,1000000,40000,1200000),
              (900000,1100000,20000,1200000), (900000,900000,20000,1200000))
    for center in profiles:
        box=tuple((v-1024,v+1024) for v in center)
        print(dict(h=h,center_R_L_f_F=center,halfwidth=1024,**box_diagnostic(moment,markers,box)),flush=True)
    print(dict(seconds=round(time.monotonic()-start,3),certificate=False),flush=True)


if __name__=='__main__':
    if shutil.disk_usage('/tmp').free<20*(1<<30):raise RuntimeError('20GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(2*(1<<30),2*(1<<30)))
    timer=threading.Timer(175,lambda:os._exit(124));timer.daemon=True;timer.start()
    if '--checks' in sys.argv:box_checks();closure_class_certificate()
    else:box_probe()
    timer.cancel()
