"""Exact finite census checks for genuine quartet availability histories.

No encoder change. Load with the lowrank helper namespace for guarded Peach
execution. New class exclusions do not establish the full certificate.
"""
from itertools import product


def census_choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def availability_history_counts(width, steps):
    """Independent enumeration of one-coordinate words, then exact product."""
    coordinate = {}
    for word in range(1 << steps):
        activity = word.bit_count()
        closed = sum(bool(word & (1 << t)) and not bool(word & (1 << (t+1)))
                     for t in range(steps-1))
        endpoint = (word >> (steps-1)) & 1
        key = activity, closed, endpoint
        coordinate[key] = coordinate.get(key, 0) + 1
    for (activity, closed, endpoint), count in coordinate.items():
        formula = (census_choose(activity-1, closed+endpoint-1)
                   * census_choose(steps-activity, closed)) if activity else 1
        assert count == formula
    result = {(0,0,0):1}
    for _ in range(width):
        after = {}
        for (l,d,e), count in result.items():
            for (ll,dd,ee), number in coordinate.items():
                key = l+ll,d+dd,e+ee
                after[key] = after.get(key,0) + count*number
        result = after
    return result


def census_checks():
    counted = weighted = recovered = 0
    for width,steps in ((1,4),(2,4),(4,4),(3,6)):
        size=width*steps;counts=availability_history_counts(width,steps)
        by_pair={};by_death={}
        for (l,d,e),count in counts.items():
            if l:
                padded=census_choose(l-1,d+e-1)*census_choose(size-l+width,d+e)
                assert count<=padded
                assert d<=min(l,size-l)
            by_pair[l,d]=by_pair.get((l,d),0)+count
            by_death[d]=by_death.get(d,0)+count
        for (l,d),count in by_pair.items():
            if not l:continue
            padded=sum(census_choose(l-1,d+e-1)*census_choose(size-l+width,d+e)
                       for e in range(width+1))
            assert count<=min(comb(size,l),padded)
            counted+=1
        for d,count in by_death.items():
            assert count<=min(1<<size,census_choose(size+width,2*d+width))
        for j in range(5):
            exact=sum(F(count,1<<(2*j*d)) for d,count in by_death.items())
            pf=(1<<(j*width))*(1+F(1,1<<j))**size
            assert exact<=pf
            weighted+=1
    # Exhaust actual state masks and verify V is recovered from availability.
    for width,steps,partitions in ((4,3,((15,),)),(8,2,((15,240),(51,204),(85,170)))):
        full=(1<<width)-1
        for partition in partitions:
            gamma=[sum(quad for quad in partition if quad&mask) for mask in range(full+1)]
            seen={}
            for states in product(range(full+1),repeat=steps):
                previous=0;footprint=[]
                for state in states:
                    footprint.append(gamma[previous]|state);previous=state
                available=[footprint[t]&sum(quad for quad in partition if footprint[t+1]&quad==quad)
                           for t in range(steps-1)]+[footprint[-1]]
                previous=0;rebuilt=[]
                births=0
                for allowed in available:
                    births+=(allowed&~gamma[previous]).bit_count()
                    rebuilt.append(gamma[previous]|allowed);previous=allowed
                assert rebuilt==footprint
                covered=sum(mask.bit_count() for mask in footprint)
                activity=sum(mask.bit_count() for mask in available)
                assert 5*covered-activity<=4*width*steps
                assert 2*covered-activity-births<=width*(steps-1)
                key=tuple(available);value=tuple(footprint)
                assert seen.setdefault(key,value)==value
                recovered+=1
    # Small rational comparisons certify the weighted full-size constants.
    assert 9**64<2**203
    assert F(1024,989)*F(1024,1023)<2
    # ln2>2/3 follows by integrating 1/t>4/(t+1)^2 for 1<t<2.
    # log2(1+u)<3u/2 then gives the general integer census bound below.
    constants=[]
    for j in range(1,17):
        numerator=3*N;denominator=1<<(j+1)
        exponent=512*j+(numerator+denominator-1)//denominator
        if j==3:exponent=361984
        constants.append((j,exponent,exponent+45,16*W+exponent+92))
    print(dict(exact_profile_checks=counted,weighted_census_checks=weighted,
               realizable_footprint_recovery_checks=recovered,
               weighted_census_constants_j_c_rank_threshold_closure_threshold=constants,
               each_class_strict_failure_upper='2^-76',full_code_certificate=False),flush=True)


if __name__=='__main__':
    if shutil.disk_usage('/tmp').free<20*(1<<30):raise RuntimeError('20GiB guard')
    resource.setrlimit(resource.RLIMIT_AS,(2*(1<<30),2*(1<<30)))
    timer=threading.Timer(175,lambda:os._exit(124));timer.daemon=True;timer.start()
    census_checks()
    timer.cancel()
