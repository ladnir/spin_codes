# Mixed state histories: transitions, coverage, and the remaining gap

The construction remains the layered pair-shear inner with RS[256,128],
one regional route, and GF(2^32) keys. The GF32 sandwich is unchanged.
This note extends the [one-sided result](LAYERED_ONE_SIDED.md), then explains
why a broader mixed-state count still fails. No whole-code claim or timing
is attached to the new bounds.

## A completed extension: occasional empty and full pairs

Keep the prior notation: M=8192, b=512, N=2^21, q=2^32. For the suffix
of R regions put u=RM and P=(u-b)/2. The first incoming state is zero.
After that batch, call a pair defective if its type is 00 or 11 rather
than 10 or 01. Locations may be arbitrary and may recur in every batch.

The one-sided restricted theorem extends to histories with at most H=2048
defective pair locations: except with setup probability below 2^-52,
no such nonzero output has weight at most W=838860=floor(.40N).
The statement unions all messages, defect locations, amplitudes, orientations,
and regional suffixes. The terminal state is free. It is over J=GF(2^32),
not a full scalar-extension theorem for arbitrary GF(2^128) histories.

Here is the proof modification. Let z,l,r,f count the four pair types
00,10,01,11, respectively, and let e count exceptional full pairs. Then

    V=l+r+2f,       D=l+2r+2f-2e,       V-(D+e)=e-r.

V is the number of nonzero incoming-state variables, and D counts free
inverse-input coordinates. Enlarge the free set to all live heads and
r+f tails. The previous bound d_R therefore gives d_pi(F)<=d_R(min(P,r+H)).

Replace every defective pair by type10 to obtain a one-sided reference
history. This preserves r. Adding a coordinate to repair a 00 pair changes
the potential output support by at most two positions; removing one from
a 11 pair cannot increase that reference support. Thus the actual potential
output support is at least 3u/4-b-min(r,P-r)-2H. Also it is at least V>=P-H.
The number of independent temporal zero equations is consequently at least

    C_H(r)=max(0,P-H-W,3u/4-b-min(r,P-r)-2H-W).

Encode a set of at most H defect locations by an ordered padded list with
H entries in {0,...,P}. There are at most (P+1)^H lists. Given locations,
each pair has two possible support types. There are at most 2^H exceptional
masks and at most 2^(P+H) temporal zero-edge sets. Put
A=H*ceil(log2(P+1)); the verifier uses P.bit_length() for the ceiling.
The inverse log2 bound is

    2P+A+2H+32*(b+H-r-C_H(r)+d_R(min(P,r+H)))+1.       (I_H)

For the forward argument, next-state defects remove at most one of the
b/4 required second-coordinate zeros per defect. Current 00 pairs remove
at most one random constraint each. An exceptional full pair can remove
that constraint, but pays its own fresh root instead. Summing over internal
transitions loses at most 2H roots. The forward bound becomes

    P+A+H+32*(max(0,u-K)-u/4+b/2+2H)+1.               (F_H)

These are complete bounds for each r-class. The location count already
includes every defect count up to H; there is no missing union over H.
The same route event as before suffices.

At H=2048, exact integer interval coverage checks all 100892800 (R,r)
classes. The forward bound handles 38 suffix lengths; 1019 inverse intervals
cover the rest. The largest accepted interval upper bound is -1663.
All classes meet -80, so the prior union gives failure below 2^-52.
H=4096 does not pass: at R=206,r=416587 the per-class bound is -63,
missing the -80 target. This is not a counterexample.

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/layered_one_sided_gate.py --defects 2048
```

## Track new activations instead of treating transitions as free

For the broader diagnostic, allow arbitrary mixtures of 10,01,11 but no00
after the first batch. This is still a restricted family. Its state supports
S_1,...,S_(L-1), where L=u/b, each have size at least b/2. Define

    B=sum_(t=1)^(L-2) |S_(t+1) minus S_t|.

B counts newly active coordinates at internal transitions. The number of
deactivations is B+|S_1|-|S_(L-1)|. Initial activation plus both kinds of
transition force output weight at least 2B+2|S_1|-|S_(L-1)|>=2B.
Thus B<=W/2. The temporal zero rank is at least

    C=max(0,V-W,V+B-b/2-W,3u/4-b-min(r,l)-W).

The last term follows by comparing the two half-support references for
consecutive pairings, now applied to the zero-coordinate supports.

For the inverse free set, every full exceptional pair removes two positions.
The simultaneous head and internal-edge envelopes give

    d <= min(floor(D/2), max(0,u-K),
             max(0,max_s[min(D,P-e+E_R(s),H_R(s)+r+f-e,R*s)-128s])).

Let A_pat bound the log2 number of state-support and exceptional masks with
these counts. We use the smaller of two bounds:

    log2 binom(P,f)+log2 binom(P-f,r)+log2 binom(f,e),

and a temporal-change count. The latter chooses an initial support in at
most 2^b ways, at most 2B+b/2 support changes among u-2b positions, and
an exceptional mask in binom(f,e) ways. Binomial logarithms are bounded
upward using exact integer dyadic Chernoff bounds, not floating-point logs.

There are at most V-B-b/2 internal temporal edges. The state-amplitude
bound, before union over count classes, is therefore

    A_pat+V-B-b/2+32*(b+e-r-C+d)+1.                  (I_state)

The ordinary output enumeration also remains available:

    u+P+32*(W-D-e+d)+1.                             (I_output)

For forward bounds, each one-sided pair has one free coordinate outside
its current state support. Unless a new activation occupies that coordinate,
the next-state zero forces a fresh root. Including exceptional roots and
discarding at most b/2 last-batch one-sided pairs gives

    roots >= max(e,l+r+e-B-b/2).

Use A_pat+32*(max(0,u-K)-roots)+1 for its log2 bound.
We can also charge positions whose next state OR output is zero. For a
fixed support/zero-edge mask, these impose at most one constraint on each
free coordinate. They give the additional lower bound

    roots >= D+e+B-W-b.

This joint bound must include the extra V-B-b/2 bits for choosing temporal
zero edges. It must not reuse the cheaper support-only pattern count.
Local tests explicitly check this distinction and the free terminal boundary.

## The nine-step wiring cycle rules out some count profiles

For no00 histories, each zero coordinate contributes to r exactly when
its current pairing bit is zero. Let T=L-1 and H0=P-f be the total number
of zero coordinates across these states. A fixed coordinate label v has
d_v zero bits among its nine bits. For any contiguous temporal interval,

    |9*(number of zero pairing bits)-d_v*(interval length)| <=20.

At most B+b/2 zero-coordinate episodes occur: each ends with a new activation
or survives into the last state. Minimize sum_v d_v*time_v subject to
0<=time_v<=T and sum time_v=H0. Fill capacities T*binom(9,j) in increasing
j=0,...,9; call the resulting integer value G_T(H0). Necessarily

    9*min(l,r) >= G_T(H0)-20*(B+b/2).

For example, a static support with r=0 can omit only the all-ones label
if it persists through a complete wiring cycle. The old aggregate count
incorrectly allowed many omitted labels while retaining r=0.

A stronger necessary test uses exact two-state temporal paths for each
coordinate. For integer weights a,h,m, minimize

    a*r-h*H0+m*B

over independent zero/nonzero indicator paths. Each period has nine steps.
Min-plus powering computes the exact relaxed minimum L_T(a,h,m). This
gives a*r+m*B>=h*H0+L_T. Replacing every label by its bitwise complement
gives the same inequality for l. Ignoring spatial pairing restrictions
makes this a necessary test, not a claim that passing profiles are feasible.
The source uses twenty fixed integer choices of (a,h,m), without a search.

## A sharper uniform routing bound from ordered head degrees

The earlier head envelope unions over all vertex subsets. Instead, let D_i
be block i's number of head positions in R regions. For every threshold k,

    max_(|S|=s) sum_(i in S) D_i
       <= k*s + Z_k,      Z_k=sum_i(D_i-k)_+.

Each D_i individually is Bin(R,1/2), although the block degrees are dependent.
Thus E Z_k=M*E(Bin(R,1/2)-k)_+, computed exactly from binomial integers.

Expose the R independent uniform M/2-head subsets without replacement,
using M/2 draws each. Changing one revealed choice admits a coupling that
swaps at most one head and one tail in the final subset. Such a swap changes
Z_k by at most one. The Doob martingale therefore has conditional range at
most one per reveal, and

    Pr[Z_k >= E Z_k+t] <= exp(-2*t^2/(R*M/2)).

The integer slack ceil(sqrt(7*(R*M/2)*80/20)) makes each tail <2^-80.
Union over at most 2^17 threshold/suffix pairs gives failure <2^-63.
Minimize k*s+ceil(E Z_k)+slack over k, and retain deterministic caps.
Intersecting this event with the old head/edge events is valid; independence
between these events is not required. This improvement is used only in
the mixed diagnostic, not silently substituted into the earlier certificates.

## Results, checks, and next discriminating step

The broad screen considers 55161 selected count profiles. The wiring tests
exclude 9194 as impossible. Of the remaining 45967, 40 still fail the
per-class -100 target. This is not complete coverage, regardless of how
many selected profiles pass. The weakest retained point is

    R=192, f=589632, r=0, e=368520, B=49136,
    d<=163923, C>=585828,
    best log2 upper bound = +734997.

The positive bound is uninformative. It is neither a witness nor a failure
probability lower bound. The state has mostly full pairs, many exceptional
updates, and a slowly moving set of zero coordinates. The new ordered-head
bound lowers its shortened-dimension bound from172026 to163923, but does
not resolve it. Merely spending more time on integer coverage will not
repair this selected failure.

Exact checks pass4096 defect-support cases,729 small-cycle orientation cases,
400 local defect/root cases,87552 nine-bit discrepancy checks,4000 joint
birth/rank/root fixtures,5049 binomial upper bounds,18 exhaustive temporal
dual comparisons,1080 ordered-degree identities, and576 swap-range checks.
These test the stated inequalities; they are not random-instance evidence
for minimum distance. Runs are serial, bounded, and generate no data files.

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/test_layered_mixed.py
python -B research/workstreams/one_route_coupled_scan/coded_state/layered_mixed_gate.py
```

Next: couple the head bias and internal-edge count for the SAME outer-block
subset. Current envelopes bound them separately even though the route samples
a head subset and then a matching between heads and tails. Also retain the
fixed-coordinate consistency equations at exceptional full pairs if that
joint routing improvement is insufficient. Preserve the completed restricted
results and GF32 sandwich; no construction change is currently justified.
