# Count activation support, routing, and resets jointly

The fixed three-shear model now excludes the five-clean-crossing-episode
family with 131 touched regions except with probability at most 2^-89.
This includes all episode lengths, active-row counts, and interval placements.
It remains a restricted history result, not whole-code distance closure.

The family has five clean proper closed episodes of lengths one through
eight. Each input interval crosses one regional boundary. All touched
regions of different intervals are disjoint. A final interval remains live
through the last internal state; its internal activity is unrestricted.
There are exactly R=131 touched regions, so k=R-128=3.

## The constraint discarded by the earlier count

In a clean episode, every input lies in the coordinate support U of its
activation state. That state is the activation input, not an independently
chosen vector. Its support is a union of pairs in the first matching, and
|U|>=2^l for an episode of length l. Thus routing must simultaneously:

- place the activation input on complete matched pairs;
- place every later input on coordinates from that same U;
- provide a state which the subsequent maps return to zero.

Counting regional eligibility alone forgets the first two requirements.
This loss explains the previous positive bound for five length-three
episodes. The new count does not add a randomization step to the encoder.

## One episode with fixed nonzero input counts

Fix a nonzero message before exposing the route. The outer rows use the
same fixed RS[256,128] encoder over J=GF(2^32). Fix the interval endpoints.
An episode's l+1 input batches split into u and v batches across two
regions, with u+v=l+1. The activation batch is the first of the u batches.
Let n_u and n_v be the message's nonzero coordinate counts in these regions.

Write b=512 and M=8192. A fixed set of n nonzero rows has a uniformly
distributed n-subset of positions in its region. For an even activation
weight w, there are binom(b/2,w/2) possible pair-closed supports U.
Once U is fixed, the remaining left-region inputs have at most (u-1)w
allowed positions. Right-region inputs have at most vw allowed positions.

The probability of this containment event, weighted by the reset bound,
is at most

    H_l(u,v,n_u,n_v)
      = sum_(even w >= 2^l)
          binom(b/2,w/2) binom((u-1)w,n_u-w) binom(vw,n_v)
          / [binom(M,n_u) binom(M,n_v)] * (l/q)^w,

where q=2^32 and invalid binomial coefficients are zero. The sum also has
w<=min(n_u,b). It counts the route and support choices together. In
particular, no independent pair-closure factor is multiplied into an
already conditioned route distribution.

For each fixed route, the activation support is known before the future
maps. The reset argument bounds the map contribution by (l/q)^w in the
independent uniform pair-image comparison. The accepted pair maps incur
one total density factor two over the entire encoding. We do not pay
that factor separately for each episode.

Different episodes use disjoint regions and disjoint sets of future maps.
For the fixed message and endpoints, multiply their H factors. All input
coordinates in the tail are given probability one, an upper bound which
does not restrict its internal state activity.

## RS zeros couple the episode counts

Outside the R touched regions every input must be zero. Each active RS
row therefore lies in a shortened space of dimension k, giving at most
q^(ks) messages on s specified active rows before route exposure. Every
nonzero RS word has weight at least 129. Among its R positions it has at
most k-1 zeros. Consequently the ten early-region counts satisfy

    sum_i (2s-n_(u_i)-n_(v_i)) <= (k-1)s.

Zeros used in the tail only reduce the available early-region deficit.
A small dynamic program maximizes the product of H factors over all early
counts satisfying this budget. These counts and splits are determined by
a fixed message and interval geometry. Taking a maximum therefore needs
no union over them. Let J_s denote the resulting upper bound.

For fixed endpoints the contribution of s active rows is at most

    2 * binom(M,s) * q^(ks) * J_s.

The eligible-row route bound from the preceding note limits s to a finite
cap C, except with probability at most 2^(-156-100) per geometry. Also
s>=2^(max_i l_i). The expression above still chooses rows from M, not C:
route rarity is not counted twice. Summing over s and at most 2^156
endpoint choices gives the good-route contribution. The union of route
exceptions contributes at most 2^-100.

## Checks and scope

The verifier `diagnose_pair_activation_route.py` computes each H with
integer numerators and a common q^s denominator. Its upper logarithms
are directed integer bounds. It runs all 792 multisets of five lengths;
the endpoint union already specifies their ordering and values.

For the previously worst profile (3,3,3,3,3), the route cap is 28 and all
21 source counts from 8 through 28 are included. The good-route upper
bound is 2^-429, after the endpoint union. With the route exception it
is at most 2^-99.

Across all 792 multisets, 707 are excluded by the eligible-row cap and
minimum activation size. The worst remaining profile is (1,1,1,1,1),
with good-route upper 2^-90. Adding the route exception yields the full
restricted-family upper **2^-89**. No output-weight restriction was used.

The toy test exhaustively enumerates two regional images for width four.
It checks 80 count/split cases, including the reset weighting, against
the closed-form H expression. The selected profile took under a tenth
of a second on Peach; the full family took about 2.3 seconds. Jobs used
256 MiB address-space limits, 10/45-second process limits, and Python -B.
Only source and these compact results are retained.

## Audit repair and next decision

During this pass, an audit found that the previous eligibility note had
identified available shortened dimension with actual nonzero word weight.
An actual word can have additional zeros inside its available positions.
The revised note sums those possibilities. Each additional zero loses
32 message bits, has at most 256 location choices, and saves at most
16 reset bits. Its correction is at most 256/255 per active row. The
repaired verifier was replayed for all 129/130-region length profiles;
their reported family bounds are unchanged. The H argument above avoids
that dependence entirely by using q^(ks) before routing and maximizing
over the allowed actual zero deficits.

No structural need for redesign has emerged from these clean episodes.
Continue with shared regions and internal births/one-sided pairs, which
are not covered here. The next discriminating extension is to retain the
same activation-support count when two episodes share a region. Their
routes are then coupled by sampling without replacement; multiplying
the present independent-region factors would be invalid. Preserve this
conditional dependence rather than silently relaxing it away.

The certified control, old helpers, construction, library, and paper remain
unchanged. These restricted calculations do not establish whole-code
distance over either J or the application's larger alphabet.
