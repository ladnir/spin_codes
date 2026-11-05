"""Exact scalar replay for the all-row-count low-state-activity bound.

Uses the established conditional pin-and-shortening inequality. This verifies
its final combinatorial union, not the lemma's algebra or whole-code distance.
No floating-point arithmetic, generated files, or external dependencies.
"""
from fractions import Fraction as F
from pathlib import Path
import shutil

CERTIFIED_BANDS=((1,1,84),(2,3,86),(4,7,88),(8,15,90),(16,31,92),
                 (32,63,95),(64,127,97),(128,255,100),(256,511,102),
                 (512,1023,105),(1024,2047,108),(2048,4095,111),
                 (4096,8191,114),(8192,8192,118))


def remaining_activity_floor(active_rows):
    """Smallest activity not covered by the certified low-activity event."""
    assert 1<=active_rows<=8192
    return next(b*active_rows+1 for lo,hi,b in CERTIFIED_BANDS if lo<=active_rows<=hi)


def certificate():
    q=1<<32;rows=8192;row_dimension=128;N=1<<21;activity_per_row=84
    assert rows*256==N
    density=1/(1-(N//4)*F(280,q))
    assert density==F(1024,989)<2
    assert activity_per_row*rows<=N//2
    # e < 11/4 follows, for example, by bounding the tail starting at 1/4!
    # with a geometric series whose subsequent term ratio is at most 1/5.
    e_upper=F(1)+1+F(1,2)+F(1,6)+F(1,24)/(1-F(1,5))
    assert e_upper<F(11,4)
    # choose(rows,n)<=rows^n; choose(N,84n)<=(11N/(4*84n))^(84n).
    # Replace n in the denominator by one to obtain a geometric row union.
    ratio=rows*F(11*N,4*activity_per_row)**activity_per_row/q**(row_dimension-activity_per_row)
    assert ratio==F(1<<33)*F(11,21)**84
    assert 5*(1<<78)*11**84 < 4*21**84
    assert ratio<F(4,5*(1<<45))<1
    # The factor q/(q-1) sums all activities <=84n. The factor two is the
    # already-established global density inflation for accepted field maps.
    upper=2*F(q,q-1)*ratio/(1-ratio)
    assert upper<F(1,1<<44)
    # Sharper thresholds for larger row sets. Bound each dyadic row-count
    # range by a geometric series; its ratio is largest at the left end.
    bands=[(1,1,84)];tail=F(0)
    for j in range(1,14):
        lower=1<<j;upper_rows=min(rows,2*lower-1)
        price=(60+lower-1)//lower
        accepted=[]
        for b in range(84,128):
            candidate=rows*F(11*N,4*b*lower)**b/q**(row_dimension-b)
            if candidate<=F(1,1<<price):accepted.append(b)
        assert accepted
        b=max(accepted);bands.append((lower,upper_rows,b))
        assert b*upper_rows<=N//2
        assert price*lower>=60
        # ratio<=1/2 and ratio^lower<=2^-60.
        tail+=F(1,1<<59)
    adaptive_upper=2*F(q,q-1)*(ratio+tail)
    assert adaptive_upper<F(1,1<<44)
    assert bands[0][0]==1 and bands[-1][1]==rows
    assert tuple(bands)==CERTIFIED_BANDS
    assert all(left[1]+1==right[0] and left[2]<=right[2] for left,right in zip(bands,bands[1:]))
    print(dict(certified_failure_less_than='2^-44',all_active_row_counts=True,
               arbitrary_support_histories=True,total_state_activity_per_active_row=84,
               free_terminal=True,base_field_bits=32,whole_code_distance=False,
               arithmetic='exact rational and integer',files_written=0))
    print(dict(adaptive_activity_thresholds=bands,adaptive_failure_less_than='2^-44',
               inactive_rows_2100_threshold=next(b for lo,hi,b in bands if lo<=6092<=hi)*6092))
    return adaptive_upper


if __name__=='__main__':
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    certificate()
