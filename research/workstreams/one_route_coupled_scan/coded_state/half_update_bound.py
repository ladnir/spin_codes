"""Integer consequence of the one-half sparse-update root bound."""


def second_holes_bound(s,h0,r,i,w):
    assert 1<=w<=s and 0<=h0<s
    assert 0<=r<=s-w and 0<=i<=w and 0<=h0-r-i<=s
    return min(s-1-r,2*(s-1-h0+r+i+w))


def pair_weight_bound(s,w):
    worst = (-1,None)
    for h0 in range(s):
        for r in range(min(h0,s-w)+1):
            for i in range(min(w,h0-r)+1):
                h1 = second_holes_bound(s,h0,r,i,w)
                if h0+h1>worst[0]:
                    worst = h0+h1,(h0,r,i,h1)
    return 4*s-worst[0],worst[1]


def live_nonzero_update_count(q,dimension,old_zero,new_zero,update_zero):
    """Exact kernel inclusion-exclusion when B restricted to W is injective."""
    return q**dimension-q**old_zero-q**new_zero-q**update_zero+2


def dimension_cap(s,w,h0,h1,r):
    return min(max(0,s-h1-r)+max(0,w-h0+r),
               s-h0+max(0,w-h1))


if __name__=='__main__':
    for w in (1,2,4,8,12,16,32,64):
        lower,profile = pair_weight_bound(64,w)
        print(f'Half-supported update weight {w}: feedback pair weight >= {lower}/256; '
              f'worst relaxed profile (h0,r,i,h1)={profile}',flush=True)
    print('Local all-field bounds, not a full distance certificate or measured performance.')
