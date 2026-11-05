"""FLOAT-only whole-class diagnostics for small independent reset groups.

Uniform full feedback; batch/state b split into b/d globally-MDS groups.
Every nonzero closed segment has weight >=d+1 and reset rank >=d, not b.
Forward reset probability is still q^-b. No files or certificate claims.
"""
import argparse
import math
from pathlib import Path
import sys
import time
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'theory'))
import persistent_mds_certificate as old


def box(N, B, b, d, e, vl, vh, cl, ch, density_bits=0, field_bits=128, sum_rank=False):
    if b*vh*B+2*ch*N <= (B//2)*N:
        return -math.inf, 'outer_excludes', (5, 4)
    n = N//b
    n1, n2 = vh-cl, n-vl+ch
    count = (old.entropy_logchoose(n1, min(ch-1, n1//2))
             + old.entropy_logchoose(n2, min(ch, n2//2)))
    free = b*vh
    volume = min(float(free), e+free/field_bits)
    if e < free:
        volume = min(volume, e+old.entropy_logchoose(free, e)/(field_bits*old.LN2))
    volume -= d*(cl-1)
    if sum_rank:
        assert d == 1
        r = cl-1
        # Ambient reset rank b*r. On a support of weight w, active-coordinate
        # rank R satisfies w <= 2*R+b*(v-2*r), with a free terminal episode.
        volume = min(volume, free-b*r, (e+b*(vh-2*r))/2+free/field_bits)
    al, ah = max(0, b*(vl-ch)), b*(vh-cl)
    ul, uh = b*vl, b*vh
    closest_a = min(ah, N//2) if al <= N//2 else N-al
    closest_u = min(uh, N//2) if ul <= N//2 else N-ul
    candidates = []
    for marker in old.MARKERS:
        theta = math.log(marker[0]/marker[1])
        for kind, selected, coefficient in (
            ('inverse', closest_a, volume-min(al, N//2)+density_bits*(vh-cl)/field_bits),
            ('forward', closest_u, max(0, uh-N//2)-b*(cl-1)+density_bits*(ch-1)/field_bits)):
            excess = old.excess_log(B, B//2, ch, selected*B/N, theta)
            softplus = max(0, excess)+math.log1p(math.exp(-abs(excess)))
            routed = (N//B)*softplus
            if kind == 'forward' and uh <= N//2:
                corrected = -math.inf if excess == -math.inf else (
                    math.log(N//B)+excess+(N//B)*math.exp(excess)-math.log(math.expm1(theta)))
                routed = min(routed, corrected)
            candidates.append(((count+theta*coefficient+routed)/old.LN2, kind, marker))
    return min(candidates)


def cover(B, b, d, seconds, distance_permille=491, density_bits=0, field_bits=128, sum_rank=False):
    N = 2**21
    assert 0 < distance_permille < 500
    assert field_bits >= 16 and all(a/b <= 2**field_bits for a,b in old.MARKERS)
    e = distance_permille*N//1000
    n = N//b
    cmax = min((n+1)//2, (e-1)//(d+1)+1)
    term_bits = 40+(n*cmax-1).bit_length()
    todo = [(1, n, 1, cmax)]
    accepted = []
    start = time.monotonic()
    checks = 0
    while todo:
        if time.monotonic()-start > seconds or checks >= 50000:
            return dict(status='TIMEOUT_NOT_CERTIFICATE', checks=checks)
        vl, vh, cl, ch = todo.pop()
        vl = max(vl, 2*cl-1)
        ch = min(ch, (vh+1)//2, cmax)
        if vl > vh or cl > ch:
            continue
        checks += 1
        value, kind, marker = box(N, B, b, d, e, vl, vh, cl, ch, density_bits, field_bits, sum_rank)
        if value <= -term_bits:
            accepted.append((vl, vh, cl, ch, kind, marker, value))
        elif vl == vh and cl == ch:
            return dict(status='FAILED_BOUND_NOT_COUNTEREXAMPLE', v=vl, c=cl,
                        log2_upper=value, method=kind, marker=marker, checks=checks)
        elif vh-vl >= 2*(ch-cl) and vl < vh:
            mid = (vl+vh)//2
            todo.extend(((vl, mid, cl, ch), (mid+1, vh, cl, ch)))
        else:
            mid = (cl+ch)//2
            todo.extend(((vl, vh, cl, mid), (vl, vh, mid+1, ch)))
    return dict(status='FLOAT_COMPLETE_COVER_NOT_CERTIFICATE', boxes=accepted,
                checks=checks, term_bits=term_bits, cmax=cmax)


def activity_point(B, b, d, v, c, active):
    """Condition on active reset groups in closed episodes; charge their choices.

    Each closed active group removes d dimensions AND restricts the first
    input batch to those groups. These first-batch holes stay in endpoint
    regions, so the same <=2c partial-region moment applies.
    """
    N, W = 2**21, 491*2**21//1000
    n, r, groups = N//b, c-1, b//d
    assert r <= active <= groups*r and (d+1)*active+1 <= W
    count = (old.entropy_logchoose(v-c, c-1)+old.entropy_logchoose(n-v+c, c)
             + old.entropy_logchoose(groups*r, active))
    rank = d*active
    # Every closed episode spans >=2 batches. Inactive groups are zero in
    # both. Additional length can only increase this removal; allow all of
    # it to the unconstrained final episode for a safe capacity upper bound.
    free = b*(v-2*r)+2*rank
    volume = min(float(free), W+free/128)
    if W < free:
        volume = min(volume, W+old.entropy_logchoose(free, W)/(128*old.LN2))
    a = b*(v-c)
    # Noninitial batches, chosen closed-episode activation groups, free last start.
    selected_forward = a+rank+b
    answers = []
    for marker in old.MARKERS:
        theta = math.log(marker[0]/marker[1])
        for kind, selected, coefficient in (
                ('inverse', a, volume-rank-a),
                ('activation_forward', selected_forward, -b*r)):
            small = min(selected, N-selected)
            excess = old.excess_log(B, B//2, c, small*B/N, theta)
            route = (N//B)*(max(0, excess)+math.log1p(math.exp(-abs(excess))))
            if kind == 'activation_forward' and selected <= N//2:
                corrected = -math.inf if excess == -math.inf else (
                    math.log(N//B)+excess+(N//B)*math.exp(excess)-math.log(math.expm1(theta)))
                route = min(route, corrected)
            coefficient += max(0, selected-N//2)
            answers.append(((count+theta*coefficient+route)/old.LN2, kind, marker))
    return min(answers)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outer', type=int, choices=(256, 512), required=True)
    ap.add_argument('--group', type=int, choices=(1, 4, 8, 16, 32, 64, 128, 256, 512), required=True)
    ap.add_argument('--cover', action='store_true')
    ap.add_argument('--activity', action='store_true')
    ap.add_argument('--screen', action='store_true')
    ap.add_argument('--distance-permille', type=int, default=491,
                    help='Cover-only distance threshold; 450 means 45 percent.')
    args = ap.parse_args()
    if not 0 < args.distance_permille < 500:
        ap.error('distance must be strictly between 0 and 500 permille')
    if args.distance_permille != 491 and (not args.cover or args.screen or args.activity):
        ap.error('a changed distance is supported only with --cover')
    b = 512 if args.outer == 256 else 256
    assert args.group <= b
    start = time.monotonic()
    if args.screen:
        worst = (-math.inf, None)
        checked = 0
        n = 2**21//b
        for v in (n//2-128, n//2, n//2+128):
            for c in (16, 32, 64, 128, 256, 512, 1024):
                if 2*c-1 > v:
                    continue
                slots = (b//args.group)*(c-1)
                for active in sorted({c-1, max(c-1, slots//4), max(c-1, slots//2), max(c-1, 3*slots//4), slots}):
                    if (args.group+1)*active+1 > 491*2**21//1000:
                        continue
                    result = activity_point(args.outer, b, args.group, v, c, active)
                    checked += 1
                    if result[0] > worst[0]:
                        worst = result[0], (v, c, active, result[1:])
                    assert time.monotonic()-start < 30, 'Bounded screen expired; no certificate'
        print(dict(kind='FLOAT_SELECTED_SCREEN_NOT_CERTIFICATE', outer=args.outer,
                   group=args.group, checked=checked, worst=worst), flush=True)
    elif args.activity:
        n = 2**21//b
        c = 32 if b == 512 else 64
        for v in (n//2, n//2+128):
            slots = (b//args.group)*(c-1)
            for active in sorted({c-1, max(c-1, slots//4), max(c-1, slots//2), max(c-1, 3*slots//4), slots}):
                if (args.group+1)*active+1 > 491*2**21//1000:
                    continue
                print(dict(kind='FLOAT_ACTIVITY_CLASS_NOT_CERTIFICATE', outer=args.outer,
                           group=args.group, v=v, c=c, active=active,
                           bound=activity_point(args.outer, b, args.group, v, c, active)), flush=True)
                assert time.monotonic()-start < 30
    elif args.cover:
        result = cover(args.outer, b, args.group, 45, args.distance_permille)
        result['distance_permille'] = args.distance_permille
        if 'boxes' in result:
            result['boxes'] = len(result['boxes'])
        print(result, flush=True)
    else:
        n = 2**21//b
        for v, c in ((n//2, 16), (n//2, 32), (n//2, 64), (n//2, 128),
                     (n//2, 256), (n//2+128, 64), (n//2-128, 64)):
            print(dict(outer=args.outer, batch=b, group=args.group, v=v, c=c,
                       bound=box(2**21, args.outer, b, args.group, 491*2**21//1000, v, v, c, c)), flush=True)
            assert time.monotonic()-start < 30


if __name__ == '__main__':
    main()
