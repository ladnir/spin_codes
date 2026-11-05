# One defective batch in a closed episode

For the fixed rejected-uniform SL2 construction, the complete h<=6,a=1,c=6
no-full-reset branch now has failure probability at most 2^-71 at 46%
distance. Together with the previous branches, this closes all histories
with at most six defective incoming batches and no full-state reset,
with failure at most 2^-69. These probabilities are over setup and quantify
over all nonzero J-valued messages, for J=GF(2^32). They are not whole-code
or extension-field certificates.

The full-reset partition concerns incoming states Q_1,...,Q_(L-1).
The final state Q_L remains unconstrained, as in the retained verifier.

The construction, RS outer, regional routing, and state size are unchanged.
The proof uses the [episode envelope](UNIFORM_PAIR_DEFECT_ENVELOPE.md)
and the [unaligned frontier](UNIFORM_PAIR_UNALIGNED_FRONTIER.md).

## What the remaining defect budget implies

The branch has six live intervals and an active first state. Five subsequent
activations consume five defective batches. There is at most one internal
defective batch. The previous result covers the case where every closed
early episode is clean; it allows the internal defect in the terminal
live interval. It remains to put that defect in one closed early episode.

Fix this exceptional episode, of length l, with incoming supports
S_0,...,S_(l-1) and S_l empty. Every incoming support is nonempty and proper:
a full state in this closed episode would be a forbidden full-state reset.
Fix the exceptional update j, indexed from zero, and put r=l-j-1.
Let U=S_0, w=|U|, and e=H+P, where H counts empty-pair births and P counts
one-sided incoming pairs. All these defects occur at update j.

Before j, supports decrease and each incoming support is closed under its
current matching. Tracing one surviving coordinate backward through j
distinct matching bits forces at least 2^j initial coordinates. After j,
the clean closed suffix of length r forces at least 2^r coordinates in
S_(j+1), when r>0. Properness therefore gives

    0<=j<=8,  0<=r<=8,  l=j+r+1<=17,
    w>=2^j,  w+e>=2^r.

The last inequality also holds for r=0 because w>=1. The earlier envelope
identity supplies at most w+e possible input coordinates and exactly w+e
prescribed zero next-state coordinates in live pairs. Its probability
cost is q^(-w-e) in the independent-image comparison, q=2^32.

## Count births at one update, not at arbitrary times

For fixed U and j, every birth is at j. A birth has a coordinate label
and at most r possible death times. Every initial component has at most
l possible death times. Since the number of births is at most e, histories
with charge at most h are bounded by

    D(l,w,j,h) = min(2^(b(l-1)),
                     l^w * sum_(v=0)^min(h,b) binom(b,v) r^v),
    b=512.

The first bound counts all intermediate support masks. The second counts
birth labels and death times. Both may include invalid histories; neither
counts field-valued amplitudes. For r=0 the sum has only its v=0 term.

If j>0, the initial support is closed under the first matching, giving
binom(b/2,w/2) choices for even w and zero for odd w. If j=0, use the
previous A_b(w,h), which counts initial supports with at most h one-sided
pairs. Restrict j to w>=2^j and w+h>=2^r, then sum these initial-support
counts times D over j. Denote the resulting coefficient by C(l,w,h).

For a charge band [lo,hi], this coefficient bounds all histories in the
band together. Use probability q^(-w-lo) and envelope size
v_w=min(b,w+hi). No additional factor for the number of charges is needed.
The verifier covers e=1,...,512 with contiguous bands

    [1,1], [2,2], [3,4], [5,8], [9,16], [17,32],
    [33,64], [65,128], [129,256], [257,512].

One batch has e<=512, so this covers every possible positive charge.
The accepted-map density correction is paid once for the whole count,
not once per position, band, or episode.

## A simple joint route factor

For one regional fragment with n nonzero inputs, the contribution of
weight w in the band is at most

    C(l,w,hi) * binom(l*v_w,n-w)
      / [binom(M,n) q^(w+lo)],     M=8192.

For two fragments with n total nonzero inputs, let n_1 be the count in
the activation fragment. With s active RS rows,

    max(w,n-s)<=n_1<=min(s,n).

For a split u+v=l+1, Vandermonde's identity gives

    binom((u-1)*v_w,n_1-w) binom(v*v_w,n-n_1)
      <= binom(l*v_w,n-w).

The denominator binom(M,n_1)binom(M,n-n_1) is log-concave in n_1, so its
minimum on the allowed interval is attained at an endpoint. Divide by
that minimum. This bound applies to every regional split; the endpoint
tuple already fixes the actual split. It deliberately overcounts possible
input allocations, avoiding a separate expensive optimization over them.

Sum over w and the charge bands with directed integer logarithmic bounds.
Combine this exceptional factor with four clean-episode factors. Retain
actual RS zeros Z and virtual deficit Z+gs, where g is the number of
duplicate regional copies. The existing multifragment correction pays
for shared-region dependence and allocations.

## Geometry and full coverage of the residual branch

Let d be the maximum length among the four clean episodes, E_c their total,
and l the exceptional length. Then d+3<=E_c<=4d. The terminal interval
uses at most ceil((T-E_c-l+1)/16) regions. There are two cases:

- If the exceptional episode is first, it is aligned and uses
  ceil((l+1)/16) fragments. All four clean episodes may cross a boundary.
- Otherwise, the first clean episode is aligned, so at most three clean
  episodes cross. The exceptional episode can use one or two fragments;
  three fragments are possible only for l=17.

For m total early fragments, the mean fragment fraction is
(E_c+l+5)/(16m). Use the retained factorial-moment route bound with
threshold 129 minus the terminal-region count. Request 259 bits in its
eligible-row cap. Six endpoint pairs and five possible exceptional-episode
designations cost at most 156+3 bits, leaving a route exception at most
2^-100. The number of active rows is at least

    max(2^d, 2^max(0,l-9)).

Every three-fragment cell is excluded by this eligible-row cap throughout
the actual residual range. Thus the verifier needs no unproved
three-fragment exceptional factor. It asserts exclusion before using
the one- or two-fragment factor.

At T=1941, the screen checks 6162 count cells and excludes 4256 by the
route cap. All remaining lengths, charges, positions, sharing counts,
and active-row counts give good-route exponent at most -103. The worst
is an aligned exceptional episode of length one, four clean length-one
episodes, no sharing, and eligible cap 21.

The complete integration recomputes all 34 residual lengths T=1908,...,1941.
It checks 64128 exceptional count cells. The worst exceptional exponent
is -103 at T=1941. The clean-early subclass is replayed and has worst
exponent -74. The retained verifier covers other lengths with exponent
at most -6244. Including the old route exception and both new exceptions,

    2^-6244 + 2^-80 + 2^-103 + 2^-100 + 2^-74 + 2^-100
      <= 2^-71.

The old verifier uses W=floor(.46N); the new episode exclusions do not
need an output-weight restriction.

## All six-defect histories without a full reset

A fresh audit of the ten other branches (a=0 or1, c=1,...,5) finds no
residual. Their worst old good-route exponent is -1267, at a=1,c=4.
Each has the retained route exception at most 2^-80. Even paying that
exception ten times, their union is below 2^-75.

The previously proved a=0,c=6 and a=1,c=7 branches each have failure
at most 2^-78. Adding the new a=1,c=6 branch gives

    2^-75 + 2^-78 + 2^-78 + 2^-71 <= 2^-69.

These branches cover c>=1. A nonzero message with no nonzero incoming
state would have input supported only in the final batch, contradicting
the RS outer's minimum support of 129 regions. Thus the no-full-reset
family with at most six defective batches is fully covered.

## Tests and remaining work

The retained source is `diagnose_pair_single_defect.py`. Its `--test`
passes 24232 support-history checks, 6588 cumulative single-batch history-count checks,
1006 denominator-endpoint checks, and 30180 Vandermonde allocation checks.
These validate the implementation; the arguments above supply the general
claims. The support tests enumerate proper histories at widths2,4,8,
including every length up to five at width8. Cumulative checks validate
the charge-band use of the history cap, not just separate exact charges.
`--integrate` checks every residual length and both subclasses;
`--audit-six` rechecks all ten other branches from the old verifier.

The integration took about24.7s locally; the other-branch audit took7.1s.
Jobs were serial, used bounded caches, Python-B, and explicit internal
time limits. They produced no numerical files. These are verifier times,
not encoder timings. Numerical work began with more than22GiB free.

The frozen full-reset verifier was also replayed at h=6. It still leaves
25 residuals, with worst exponent +9147 at a=1,c=7,T=1953. Its first
residual is a=1,c=5,T=1969. This is the next target, not an established
code failure. In particular, examine the duration of the final decay
from a full state and the route cost of reaching full support. The old
full-reset bound permits each coordinate's death time anywhere in the
whole stream; that relaxation may discard the small defect budget.

Higher defect counts and extension-field closure also remain open.
No structural obstruction requiring a redesign has been established.
