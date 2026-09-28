"""Outward rank-one distance bound for a single active 16-row group.

Includes all 512 possible groups and all nonempty row subsets. This is not a
full distance certificate: higher-rank tuples and multiple groups are omitted.
No files are written. Inputs are authenticated/reconstructed on every run.
"""

import argparse
from math import comb
from flint import arb, arb_mat, ctx

import group_moment as diagnostic
import bch_joint_support as joint

TILTS = ('0.00016','0.00025','0.0004','0.00064','0.001','0.0016','0.0025')


def up(value):
    assert value.is_finite()
    return arb(value.upper())


def transfers(data,tilt,a):
    spectrum,windows,zero,maximum,cancellation = data
    levels = sorted(spectrum)
    n,m = len(levels)+2,(1 << 19)-1
    lam = arb(tilt)
    powers = [up((-lam*w).exp()) for w in range(145)]
    empty = [[arb(0) for _ in range(n)] for _ in range(n)]
    active = [[arb(0) for _ in range(n)] for _ in range(n)]
    empty[0][0] = arb(1)
    empty[1][1] = powers[min(levels)]/2
    for k,w in enumerate(levels):
        empty[1][k+2] = powers[min(levels)]*spectrum[w]/(2*m)
    for i,v in enumerate(levels):
        empty[i+2][i+2] += powers[v]/2
        for k,w in enumerate(levels):
            empty[i+2][k+2] += powers[v]*spectrum[w]/(2*m)
    choices = 8*comb(16,a)
    active[0][0] = powers[a]*zero[a]/choices
    active[0][1] = powers[a]*(choices-zero[a])/choices
    d = powers[max(0,min(levels)-a)]
    min_cancel = min(w for row in cancellation[a].values() for w in row)
    cd = min(up(d),up(powers[min_cancel]*maximum[a]/choices))
    active[1][0] = cd/2+d/(2*m)
    active[1][1] = d/2
    for k,w in enumerate(levels):
        active[1][k+2] = d*spectrum[w]/(2*m)
    for i,v in enumerate(levels):
        moment = arb(0)
        for local,count in windows[v].items():
            for overlap in range(max(0,a+local-16),min(local,a)+1):
                moment += count*comb(local,overlap)*comb(16-local,a-overlap)*powers[v+a-2*overlap]
        moment = up(moment/(choices*spectrum[v]))
        cancel = up(sum((count*powers[w] for w,count in cancellation[a][v].items()),arb(0))/(choices*spectrum[v]))
        active[i+2][0] = cancel/2+moment/(2*m)
        active[i+2][1] = moment/2
        for k,w in enumerate(levels):
            active[i+2][k+2] = moment*spectrum[w]/(2*m)
    return arb_mat([[up(v) for v in row] for row in empty]),arb_mat([[up(v) for v in row] for row in active])


def regions(empty,active):
    n = empty.nrows()
    rz = arb_mat([[int(i==j) for j in range(n)] for i in range(n)])
    ra = arb_mat(n,n)
    for _ in range(64):
        ra = ra*empty+rz*active
        rz = rz*empty
    return (arb_mat([[up(rz[i,j]) for j in range(n)] for i in range(n)]),
            arb_mat([[up(ra[i,j]/64) for j in range(n)] for i in range(n)]))


def update_best(best,weights,empty,active,tilt):
    n = empty.nrows()
    z,a = empty.transpose(),active.transpose()
    first = arb_mat([[active[0,j]] for j in range(n)])
    current = arb_mat([[1]*n])
    factor = (arb(tilt)*209715).exp()
    for remaining in range(256):
        conditional = current*first
        for w in weights:
            if w <= remaining+1:
                value = up(conditional[w-1,0]*factor/comb(remaining,w-1))
                best[remaining][w] = min(best[remaining][w],value)
        if remaining == 255:
            break
        ez,ea = current*z,current*a
        count = current.nrows()
        # Nonnegative coefficient recurrence; upper endpoints prevent interval
        # radii from accumulating across polynomial dynamic-programming steps.
        current = arb_mat([[up((ez[w,j] if w<count else arb(0))
                              +(ea[w-1,j] if w else arb(0)))
                            for j in range(n)] for w in range(count+1)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--rows',type=int,nargs='+',default=list(range(1,17)))
    args = parser.parse_args()
    assert args.precision >= 128 and all(1 <= a <= 16 for a in args.rows)
    assert len(set(args.rows)) == len(args.rows)
    spectrum = joint.authenticated_caps()
    weights = [w for w in range(1,257) if spectrum[w]]
    data = diagnostic.census()
    # Historical verifier imports set their own precision at module scope.
    # Apply this run's setting after those imports/authentication have ended.
    ctx.prec = args.precision
    assert ctx.prec == args.precision
    total = [[arb(0)]*257 for _ in range(256)]
    for a in args.rows:
        best = [[arb(1)]*257 for _ in range(256)]
        for tilt in TILTS:
            update_best(best,weights,*regions(*transfers(data,tilt,a)),tilt)
        for l in range(256):
            for w in weights:
                if w <= l+1:
                    total[l][w] = up(total[l][w]+comb(16,a)*best[l][w])
        print(f'Outward replay {ctx.prec} bits: completed active-row count {a}',flush=True)
    upper = arb(0)
    for w in weights:
        probability = sum((comb(l,w-1)*min(arb(1),total[l][w]) for l in range(w-1,256)),arb(0))/comb(256,w)
        upper += 512*spectrum[w]*probability
    upper = up(upper)
    print('Outward union upper:',upper,flush=True)
    print('Margin (display only):',-upper.log()/arb(2).log(),flush=True)
    if set(args.rows) == set(range(1,17)):
        assert 0 < upper < arb(2)**-41
        print('VERIFIED: rank-one, one-active-group distance failure union < 2^-41',flush=True)
    print('Full SPIN distance certificate: NO (higher ranks and multiple active groups missing).')


if __name__ == '__main__':
    main()
