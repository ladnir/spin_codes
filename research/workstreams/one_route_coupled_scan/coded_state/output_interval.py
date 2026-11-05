"""Integer certificate for a dense mixed class, not full minimum distance."""
import time

N, K, n, M = 1 << 21, 1 << 20, 16384, 8192
W = 491*N//1000


def touched_regions(length, per_region=64):
    return (length+2*per_region-2)//per_region if length else 0


def dense_mixed_certificate():
    started = time.monotonic()
    worst = (-10**20,None)
    cells = 0
    for v in range(1,n//2+1):
        live = v-1
        shortening = M*max(0,touched_regions(live)-128)
        for empty in range(33):
            j = n-v-empty
            if j<0 or 65*j>W:
                continue
            allowed = 128*(v+j)
            log2_bound = allowed+128*(W-64*j+shortening-65*live)
            if log2_bound>worst[0]:
                worst = log2_bound,(v,j,empty,shortening)
            cells +=1
        if time.monotonic()-started>10:
            raise RuntimeError('Integer cover cap')
    # At most n^2 interval placements,33 empty counts,n^32 choices of empties.
    union_bits = 28+6+448
    margin = -worst[0]-union_bits
    assert cells==252978 and margin>=40000,(cells,worst,margin)
    print(f'Dense one-episode certificate: {cells} cells, worst={worst}, '
          f'union cost <{union_bits} bits, margin >={margin} bits at >49.1%.',flush=True)
    # The former forward obstruction: a terminal episode in the final half.
    v,j = 8192,8192
    assert 128*(v-1)<=128*M
    exact_aligned = N+128*(W-64*j-65*(v-1))
    assert exact_aligned==-1359104
    print(f'Aligned former obstruction: log2 failure <= {exact_aligned}. '
          'Other mixed classes remain open.',flush=True)
    return cells,worst,margin


if __name__=='__main__':
    dense_mixed_certificate()
