"""Positive exponential-potential action of the exact global two-count kernel.

Requires injected lowrank/direct helpers. FLOAT proposals are not certificates.
"""


def joint_potential_action(width,m,b,p,r,alpha,lt,ls,coverage=1.):
    groups=width//4;t=LD(2)**lt;s=LD(2)**ls
    a=LD(m)*p*coverage;r0=LD(r)*coverage
    hit=np.array([LD(0)]+[LD(comb(4,i))*r0**4*(s+b*t)**(4-i)*(alpha*s+b*t)**i for i in range(1,5)])
    partial=sum(LD(comb(4,i))*(a*t)**i for i in range(4))
    empty=np.array([(1+a*t)**4]+[LD(comb(4,i))*partial for i in range(1,5)])
    apowers=np.zeros((width+1,groups+1),dtype=LD);bpowers=apowers.copy()
    left=np.ones(1,dtype=LD);right=left.copy()
    for k in range(groups+1):
        apowers[:len(left),k]=left*comb(groups,k)
        bpowers[:len(right),groups-k]=right
        if k<groups:left=np.convolve(left,hit);right=np.convolve(right,empty)
    coefficients=apowers@bpowers.T
    result=np.zeros((width+1,width+1),dtype=LD)
    for u in range(width+1):
        denominator=np.array([LD(comb(width,u))*LD(comb(width-u,e)) for e in range(width-u+1)])
        result[u,:width-u+1]=coefficients[u,:width-u+1]/denominator
    if not np.all(np.isfinite(result)):raise ArithmeticError('nonfinite potential action')
    return result


def joint_potential_checks():
    # Enumerate masks and every width4/8 quartet partition, independently.
    entries=0
    for width in (4,8):
        groups=width//4;full=(1<<width)-1
        partitions=[(full,)] if width==4 else []
        if width==8:
            for rest in combinations(range(1,8),3):
                first=1+sum(1<<i for i in rest);partitions.append((first,full^first))
        m,b,p,r,alpha=map(F,(1,2,F(1,3),F(2,5),F(3,4)))
        t,s=F(1,2),F(2)
        actual=joint_potential_action(width,LD(1),LD(2),LD(1)/3,LD(2)/5,LD(3)/4,-1,1)
        for u in range(width+1):
            oldA=(1<<u)-1
            for e in range(width-u+1):
                oldE=((1<<e)-1)<<u;truth=F(0)
                for partition in partitions:
                    hit=sum(quad for quad in partition if quad&oldA)
                    if hit&oldE:continue
                    for newA in range(full+1):
                        newV=hit|newA
                        if any((quad&oldE) and ((quad&newV)==quad) for quad in partition):continue
                        newE=hit&~newA;birth=newA&~hit
                        truth+=(r**hit.bit_count()*(m*p)**birth.bit_count()
                                *b**(newA&hit).bit_count()*alpha**(oldA&~newA).bit_count()
                                *t**newA.bit_count()*s**newE.bit_count())/len(partitions)
                target=LD(truth.numerator)/LD(truth.denominator)
                assert abs(actual[u,e]-target)<LD('1e-15')*max(LD(1),target)
                entries+=1
    print(dict(exact_joint_potential_entries=entries,full_code_certificate=False),flush=True)


def joint_potential_probe():
    h=2048;theta=LD(1)/16;logz=-LD(7)/16
    q=LD(2)**32;x=LD(2)**(logz/theta);c=(q+1)*x;lam=1+c
    d=max(min(LD(1),lam/(q+1)),(1+q*x*x)/lam)
    m=c**theta;b=lam**theta;alpha=d**theta;field=q**theta
    actual=LD(h)/8192;low=actual/(field-(field-1)*actual)
    choose=(lgamma(8193)-lgamma(h+1)-lgamma(8193-h))/float(np.log(LD(2)))
    started=time.monotonic()
    for lt in (-4,-3,-2,-1):
        for ls in (0,1):
            best=None
            for frac in (0.,.5,1.):
                odds=(low/(1-low))**(1-frac)*(actual/(1-actual))**frac;p=odds/(1+odds)
                logpi=choose+h*np.log2(p)+(8192-h)*np.log2(1-p)
                action=joint_potential_action(512,m,b,p,p+(1-p)/field,alpha,lt,ls)
                worst=LD('-inf');at=None
                for u in range(513):
                    row=action[u,:513-u]
                    feasible=(u+3)//4+(np.arange(513-u)+3)//4<=128
                    if not np.all(row[feasible]>0):raise ArithmeticError('positive action underflow')
                    with np.errstate(divide='ignore'):
                        ratios=np.log2(row)-u*lt-np.arange(513-u)*ls
                    e=int(np.argmax(ratios))
                    if ratios[e]>worst:worst=ratios[e];at=(u,e)
                logregion=16*worst-logpi
                if best is None or logregion<best[0]:best=(logregion,float(p),at,float(worst))
            terminal=max(0,-512*lt)
            bound=theta-32*theta-W*logz+choose-32*theta*128*h+256*best[0]+terminal
            print(dict(logt=lt,logs=ls,log2_bound=float(bound),p=best[1],worst_state=best[2],
                       log2_lambda=best[3],certificate=False),flush=True)
    print(dict(seconds=round(time.monotonic()-started,3)),flush=True)


if __name__=='__main__':
    if shutil.disk_usage('/tmp').free<20*(1<<30):raise RuntimeError('20GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(2*(1<<30),2*(1<<30)))
    timer=threading.Timer(175,lambda:os._exit(124));timer.daemon=True;timer.start()
    if '--checks' in sys.argv:joint_potential_checks()
    else:joint_potential_probe()
    timer.cancel()
