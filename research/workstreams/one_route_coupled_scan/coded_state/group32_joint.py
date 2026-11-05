"""Bounded FLOAT diagnostics for grouped resets; never a certificate.

Charge inactive interior output coordinates under uniform feedback, and
sum the associated excursion lengths and active-group masks jointly.
No files are generated. Frozen certificate helpers are read-only.
"""
import math
import time
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'theory'))
import persistent_mds_certificate as old
from striped_reset_probe import activity_point

FINE_MARKERS = tuple(sorted(set(old.MARKERS) | {(128+j, 128) for j in range(1, 33)},
                            key=lambda pair: pair[0]/pair[1]))
MARKERS = old.MARKERS


def logsum(values):
    top = max(values)
    if top == -math.inf:
        return top
    return top + math.log(sum(math.exp(x-top) for x in values))


def weighted_lengths(groups, d, r, active, extra, theta, parameters=False):
    """Upper log coefficient of x^active z^extra in

    (1-z)^-1 [sum_{a=1}^G C(G,a)x^a/(1-z*tau^(-d(G-a)))]^r.
    Extra = v-2r-1. All coefficients are nonnegative.
    Any x>0, 0<z<1 is valid; optimization affects tightness only.
    """
    if not r:
        assert active == 0
        return (0.0, 0.0, -math.log(extra+2)) if parameters else 0.0
    logs = [math.log(math.comb(groups, a)) for a in range(1, groups+1)]

    def objective(lx, lminuslz):
        # z=exp(-exp(lminuslz)), avoiding z rounding to one.
        lz = -math.exp(lminuslz)
        terms = [logs[a-1]+a*lx-math.log(-math.expm1(lz-theta*d*(groups-a)))
                 for a in range(1, groups+1)]
        return -active*lx-extra*lz-math.log(-math.expm1(lz))+r*logsum(terms)

    fraction = min(.999, max(.001, active/(groups*r)))
    lx = math.log(fraction/(1-fraction))
    lminuslz = -math.log(extra+2)
    best = objective(lx, lminuslz)
    step = 2.0
    for _ in range(80):
        choices = [(objective(lx+dx, lminuslz+dz), lx+dx, lminuslz+dz)
                   for dx, dz in ((step, 0), (-step, 0), (0, step), (0, -step))
                   if -50 < lminuslz+dz < 10 and -60 < lx+dx < 60]
        candidate = min(choices)
        if candidate[0] < best:
            best, lx, lminuslz = candidate
        else:
            step *= .5
            if step < 1e-7:
                break
    # The unweighted combinatorial bound is another valid choice.
    simple = (old.entropy_logchoose(extra+r, r)
              + old.entropy_logchoose(groups*r, active))
    return (best, lx, lminuslz) if parameters else min(best, simple)


def point(B, b, d, v, c, active):
    N, n, r = 2**21, 2**21//b, c-1
    selected = b*(v-c+1)+d*active
    gaps = old.entropy_logchoose(n-v+c, c)
    answers = [(*activity_point(B, b, d, v, c, active), 'old')]
    rank, a = d*active, b*(v-c)
    E = b*(v-2*r)+2*rank
    W = 491*N//1000
    # Concave output-volume tangents retain the SAME inactive-coordinate
    # cost H in both counting directions. Each (intercept,slope) bounds
    # log_q(volume(E-H,W)) by intercept-slope*H for all feasible H.
    tangents = [(float(E), 1.0), (W+E/128, 1/128)]
    if E > W:
        for fraction in (1, .5, .25, .125):
            anchor = W+(E-W)*fraction
            slope = math.log(anchor/(anchor-W))/(128*old.LN2)
            intercept = W+old.entropy_logchoose(anchor, W)/(128*old.LN2)+slope*(E-anchor)
            tangents.append((intercept, slope))
    for marker in MARKERS:
        theta = math.log(marker[0]/marker[1])
        count = gaps+weighted_lengths(b//d, d, r, active, v-2*r-1, theta)
        excess = old.excess_log(B, B//2, c, min(selected, N-selected)*B/N, theta)
        route = (N//B)*(max(0, excess)+math.log1p(math.exp(-abs(excess))))
        if selected <= N//2:
            corrected = -math.inf if excess == -math.inf else (
                math.log(N//B)+excess+(N//B)*math.exp(excess)-math.log(math.expm1(theta)))
            route = min(route, corrected)
        result = (count+theta*(max(0, selected-N//2)-b*r)+route)/old.LN2
        answers.append((result, 'joint_forward', marker, 'new'))
        inverse_excess = old.excess_log(B, B//2, c, min(a, N-a)*B/N, theta)
        inverse_route = (N//B)*(max(0, inverse_excess)+math.log1p(math.exp(-abs(inverse_excess))))
        for intercept, slope in tangents:
            count = gaps+weighted_lengths(b//d, d, r, active, v-2*r-1, theta*slope)
            result = (count+theta*(intercept-rank-min(a, N//2))+inverse_route)/old.LN2
            answers.append((result, 'joint_inverse', marker, 'new'))
    return min(answers)


def bennett_excess(mu, variance, k, theta):
    """Upper log E[(exp(theta*(T-k)_+)-1)] from a variance cap.

    T is a sum of independent Bernoullis for ONE outer block. The actual
    blocks remain negatively associated, as in the retained route proof.
    Bennett's centered mgf plus an integer tail sum gives this envelope.
    """
    if variance == 0:
        return math.log(math.expm1(theta*(mu-k))) if mu > k else -math.inf

    def value(lam):
        return (math.log(math.expm1(theta))-theta-math.log(math.expm1(lam-theta))
                -lam*(k-mu)+variance*(math.expm1(lam)-lam))

    lo, hi = theta+1e-9, max(2.0, theta+1)
    for _ in range(55):
        left, right = (2*lo+hi)/3, (lo+2*hi)/3
        if value(left) < value(right):
            hi = right
        else:
            lo = left
    return value((lo+hi)/2)


def variance_cap(B, b, d, c, active, holes, forward):
    """Joint activation geometry; holes=sum (b-da_i)(length_i-2)."""
    r, G, m = c-1, b//d, (2**21//B)//b
    full = min(r, (active-r)//(G-1))
    base = (r+active/G)/m if forward else r/m
    extra_full = (.5-2/m) if forward else (.5-1/m)
    return min(B/4, c/2, .5+base+full*extra_full+holes/(m*d))


def joint_hole_box(B, b, d, v, c, active, hlo, hhi):
    N, W, n, r = 2**21, 491*2**21//1000, 2**21//b, c-1
    gap = old.entropy_logchoose(n-v+c, c)
    simple = (old.entropy_logchoose(v-c, c-1)
              + old.entropy_logchoose((b//d)*r, active))
    rank, a = d*active, b*(v-c)
    free = b*(v-2*r)+2*rank-hlo
    volume = min(float(free), W+free/128)
    if W < free:
        volume = min(volume, W+old.entropy_logchoose(free, W)/(128*old.LN2))
    answers = []
    for marker in old.MARKERS:
        theta = math.log(marker[0]/marker[1])
        length = min(simple-theta*hlo,
                     weighted_lengths(b//d, d, r, active, v-2*r-1, theta))
        for forward, selected, coeff, count in (
                (True, a+rank+b, -b*r, gap+length),
                (False, a, volume-rank-a, gap+simple)):
            mu = min(selected, N-selected)*B/N
            variance = variance_cap(B, b, d, c, active, hhi, forward)
            excess = min(old.excess_log(B, B//2, c, mu, theta),
                         bennett_excess(mu, variance, B//2, theta))
            routed = (N//B)*(max(0, excess)+math.log1p(math.exp(-abs(excess))))
            if forward and selected <= N//2:
                corrected = -math.inf if excess == -math.inf else (
                    math.log(N//B)+excess+(N//B)*math.exp(excess)-math.log(math.expm1(theta)))
                routed = min(routed, corrected)
            value = (count+theta*(coeff+max(0, selected-N//2))+routed)/old.LN2
            answers.append((value, 'forward' if forward else 'inverse', marker))
    return min(answers)


def joint_point(B, b, d, v, c, active):
    """Partition H into bounded bins; a union of all bins for this class."""
    maximum = (b-d)*(v-2*(c-1)-1)
    bounds = sorted({0, maximum+1} | {min(maximum+1, x) for x in
                     (256, 1024, 4096, 16384, 65536, 262144, 1048576)})
    parts = [joint_hole_box(B, b, d, v, c, active, lo, hi-1)
             for lo, hi in zip(bounds, bounds[1:])]
    value = logsum([x[0]*old.LN2 for x in parts])/old.LN2
    return min((value, 'joint_holes', tuple(x[0] for x in parts)),
               point(B, b, d, v, c, active))


def main():
    global MARKERS
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--holes', action='store_true')
    ap.add_argument('--selected', action='store_true')
    ap.add_argument('--fine', action='store_true')
    ap.add_argument('--group', type=int, choices=(32, 64, 128), default=32)
    args = ap.parse_args()
    if args.fine:
        MARKERS = FINE_MARKERS
    start = time.monotonic()
    for B, b in ((512, 256), (256, 512)):
        worst = (-math.inf, None)
        checked = 0
        n = 2**21//b
        for v in (n//2-128, n//2, n//2+128):
            for c in (16, 32, 64, 128, 256, 512, 1024):
                if 2*c-1 > v:
                    continue
                slots = (b//args.group)*(c-1)
                for active in sorted({c-1, max(c-1, slots//4), max(c-1, slots//2),
                                      max(c-1, 3*slots//4), slots}):
                    if (args.group+1)*active+1 > 491*2**21//1000:
                        continue
                    if args.selected and (v, c) != (n//2+128, 128):
                        continue
                    result = (joint_point if args.holes else point)(B, b, args.group, v, c, active)
                    if args.selected:
                        print(dict(B=B, v=v, c=c, active=active, result=result), flush=True)
                    if result[0] > worst[0]:
                        worst = result[0], (v, c, active, result[1:])
                    checked += 1
                    assert time.monotonic()-start < 50, 'Bounded diagnostic expired'
        print(dict(kind='FLOAT_SELECTED_NOT_CERTIFICATE', B=B, checked=checked,
                   worst=worst, seconds=round(time.monotonic()-start, 3)), flush=True)


if __name__ == '__main__':
    main()
