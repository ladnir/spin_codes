# Independent audit of the regional subfield proof

## Main conclusion

A projective-source-word first moment cannot certify the whole current
regional construction at 10% distance. A rare routing event creates a
large subspace of bad words. Its contribution alone makes that expected
count exceed 2^80. This is an obstruction to this proof method, not a
comparable lower bound on the probability that a sampled code is bad.
The local two-state transfer remains exact. Increasing the state size
does not remove this first-moment obstruction in the checked examples.

The following arguments concern the GF256 ensemble with a repeated
GRS[16,8] outer and independently uniform regional permutations of
individual symbols. They do not transfer a uniform global permutation
bound onto this route.

## A window event with many bad source lines

Write q=256, b=16, k=8, and L=8192 for the field size, outer length,
outer dimension, and number of repeated outer rows. There are K=kL
message symbols and N=bL output symbols. Each region is divided into
t=16-symbol steps. The causal linear inner has an s-symbol state,
starts at zero, and discards its final state.

Fix a set of a outer rows. In every region choose the suffix containing
W=816 symbols, or 51 whole steps. Consider the routing event that all
a coordinates from these rows land in that suffix. Its probability is

    p_event = (binom(W,a) / binom(L,a))^b.

The regional permutations are independent, so the product is exact.
The event is specified before choosing a message. Independently sampled
nonzero diagonal coefficients do not change it.

Restrict the message to the selected rows, leaving ka field variables.
For a realized setup satisfying the event, impose zero state at the
end of each of the first b-1 regions. These impose at most (b-1)s
homogeneous linear equations. Their solution space has dimension at
least

    D = ka - (b-1)s.

Every such message produces zero output outside the suffix windows.
Indeed, the first region starts in zero state, and its prefix input is
zero. Each subsequent region starts in zero state by the imposed
equation. All possibly nonzero output lies inside its suffix. The last
region needs no end-state equation because the stream ends at that
suffix and discards its state. Thus every output has weight at most

    bW = 13056 < 13108 = ceil(N/10).

There are at least (q^D-1)/(q-1) nonzero projective source lines in this
subspace. Consequently the expected number of bad source lines is at
least

    (binom(W,a) / binom(L,a))^b * (q^D-1)/(q-1).

Only one fixed row set is used; no union over row sets or independence
between source-line events is needed. Exact integer comparisons give:

| State symbols s | Selected rows a | Guaranteed dimension D | Expected bad-line lower bound |
|---|---|---|---|
| 2 | 32 | 226 | greater than 2^80 |
| 4 | 64 | 452 | greater than 2^80 |
| 8 | 128 | 904 | greater than 2^80 |

The dimension argument does not assert that every resulting message
uses all a rows. It counts distinct source lines supported inside that
row set. This is enough to obstruct the complete source-word first
moment, including one that sums separately over exact active-row counts.

The routing event itself is very rare. Large expected bad-word counts
arise because many bad messages occur together. Multiplying its
probability by the dimension factor is valid for the expectation;
it must not be presented as a code-failure probability lower bound.

## Sound partial bounds that remain useful

Let R_h(z) be the exact two-state transfer for a region containing h
nonzero input symbols, normalized by binom(L,h), with 0<=z<=1. State
categories are zero and nonzero. For a specified nonzero row support S
of weight w, the MDS shortening dimension gives at most q^(w-8) words
when w>=9. Assigning the same positive envelope to impossible supports
w<9 only enlarges the count.

For exactly a selected active rows, sum the relaxed binary row-support
patterns columnwise and define

    B_a(z) = sum_(h=0)^a binom(a,h) q^h R_h(z).

The expected bad-source-line count in this class is at most

    binom(L,a)/(q-1) * q^(-8a) * z^(-e)
        * (1,0) B_a(z)^16 (1,1)^T,

where e=13107. This preserves state across regions. R_h already
averages the diagonal-randomized values, so no extra (q-1)^h belongs
in B_a. Only one global projective division by q-1 is valid.
This bound is sound but, by the window obstruction, cannot close all
active-row classes at the desired threshold.

A simple dense-tail bound is independent of the routing and state
dimension. Put f=z+(1-z)/(q-1). Condition on the state after its
refresh, before exposing this step's independent nonzero input
scalars. Each active lane has output moment z if its state contribution
is zero, and f otherwise. Each inactive lane has moment at most one.
Consequently a step with w active inputs has conditional moment at
most f^w. Repeated conditioning bounds a fixed outer word's output
moment by f raised to its outer weight. This argument requires fresh
independent input scalars, but not a particular state-refresh distribution.

Let P be the exact nonzero GRS[16,8] word enumerator. The a-row class
is then bounded by binom(L,a)*P(f)^a*z^(-e)/(q-1). At fixed z,
successive binomial terms have decreasing ratio. Starting at a0, the
entire tail is bounded by its first term divided by
1-(L-a0)*P(f)/(a0+1), when this denominator is positive.

At the rational tilt z=7/2500, an exact integer comparison proves that
the combined contribution of all a=2879,...,8192 is below 2^-64.
This partial coverage is consistent with the obstruction, which occurs
among much smaller row supports. It is not a whole-code certificate.

## A possible proof pivot: minimal block-support witnesses

Fix an output support T. Its bad-message space is the kernel of the
generator restricted to output coordinates outside T. Suppose A is
an inclusion-minimal nonempty set of outer message blocks on which
this kernel is nonzero. Every nonzero vector in the restricted kernel
then uses every block of A. Projection onto any one k-symbol block is
injective; otherwise a nonzero vector would avoid that block. Hence
the restricted kernel has dimension at most k, regardless of |A|.

This lemma explains the overcount above. The large window subspace
contains smaller-support witnesses: the same end-state equations
already leave a nonzero solution on floor((b-1)s/k)+1 rows. Counting
all combinations of those witnesses multiplies the first moment
without multiplying the underlying bad-setup event.

A useful next proof should track independent constraints or minimal
support witnesses, not all field-valued messages. The lemma alone is
not a probability bound: one must still control routing, output
supports, and the rank of the resulting equations. The two-state
fixed-message transfer does not retain enough information for that
task by itself.

## Audit of the local rank bounds

The local lemmas in [SUBFIELD_PROOF.md](SUBFIELD_PROOF.md) pass a separate
review. Four independent GF4 tests enumerate dependent input forms, every
fixed two-by-two return target, and the empty-input rank operator. A
counterexample checks that dropping required-nonzero conditions can invalidate
the diagonal lemma.

The active-input bounds also preserve the exposure order. Diagonal-only
return constraints are imposed before the fresh SL2 matrix is exposed.
The state rank is recomputed after those constraints. In particular, an
old rank-two state cannot automatically receive both the diagonal penalty
and the rank-two return penalty. The minimum of the two separately proved
upper bounds avoids that mistake.

For nonreturns, canceled inputs lie in the refreshed-state span. Only
uncanceled input forms can restore rank. Maximizing the remaining-computation
bound over feasible next ranks is therefore sound. This validates local
inequalities, not their yet-unevaluated composition across regional routing.

## Reproduction

The independent audit module is `theory/subfield_bound_audit.py`.
Its exact checks contain no floating-point logarithms. A separately
labeled floating-point bulk screen is only a search aid.

```sh
python -B research/workstreams/char2_field_spin/theory/subfield_bound_audit.py --verify-bulk
python -B research/workstreams/char2_field_spin/theory/subfield_bound_audit.py --verify-window-obstruction
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_subfield_bound_audit.py -v
python -B -m unittest discover -s research/workstreams/char2_field_spin/lane_mixing -p test_rank_lemmas_audit.py -v
```

Six tests check the MDS envelope, compare the coupled bound against
literal GF4 support enumeration, check the scalar-potential alternative,
verify the dense row-sum inequality, and verify the three state-size
window witnesses. No production kernel or benchmark was changed.

## GF16, four-state audit

The GF16 candidate uses the same regional route with GRS[16,8] outer
blocks. Each inner step has sixteen lanes and four state symbols.
Write A_j=(1,j,j^2,j^3), for all j in GF16, and C=A^T. A fresh uniform
GL4 matrix refreshes the state before emission. Each input coordinate
has its own independent uniform nonzero GF16 scalar.

The two-state fixed-word moment calculation records whether the state
is zero. Uniform GL4 refresh makes every nonzero refreshed state
uniform on GF16^4 minus zero, even after conditioning on prior outputs.
Fresh input scalars preserve the required independence. The local
envelope in `state_screen.py` follows by counting zero positions in
the [16,4] MDS code generated by A and counting the nonzero solutions
of C_S X=0. Uniform permutations within each region justify division
by binom(L,h), separately for that region. The linear matrix products
retain the state between regions. They do not use the unresolved
nonlinear maximum from the separate quotient-rank investigation.

For one active lane j, return requires the refreshed state to equal
C_j x. At an inactive lane l, the output divided by x is
(1+lj)^3. The active output divided by x is j^2(1+j+j^2)^2.
Among the sixteen possible j, thirteen give output weight fifteen,
one gives weight sixteen, and two give weight fourteen. Thus the
averaged marked return probability is exactly

    J_1(z) = (13 z^15 + z^16 + 2 z^14)/(16(16^4-1)).

The monomial maps also satisfy CA=0, so B=I+AC is invertible with
B^2=I. On a return, Y=BX. For a fixed active input set S of size w
and output-zero set Z, rank(B[Z,S]) is at least |Z|-(16-w).
Exposing pivot coordinates bounds its zero probability by
15^(-max(0,|Z|-16+w)). Expanding the marked weight as a sum over Z
gives the entrywise return upper bound

    J_w(z) <= 1/(16^4-1) * sum_{h=0}^{16} binom(16,h)
              * z^(16-h) (1-z)^h * 15^(-max(0,h-16+w)).

This bound may be combined with the earlier return upper bound by
taking their minimum. It does not justify replacing the lower return
bound used in the nonreturn entry.

The independent `test_state_screen.py` suite enumerates all 65,535
nonzero GF16 states, all singleton lanes and nonzero values, and every
local GF4 support for state dimensions one through four. It checks
the local envelopes, singleton polynomial, self-orthogonality, and
regional normalization. These exact local checks do not certify the
floating-point whole-code screens.

For the GF16 candidate at K=2^16 and target distance 10%, an exact
integer comparison certifies the combined classes a=5998,...,8192
below 2^-64 at z=1/52. Reproduce this partial dense-tail certificate
with `theory/subfield_bound_audit.py --verify-bulk-gf16`.

### Complete coverage checkpoint

The local audit and the three coverage mechanisms now give a whole-code
bound for this ideal ensemble. The code has K=65536 and N=131072
field symbols; the bad event is a nonzero codeword of weight at most
13107. These ranges partition every possible nonzero outer-row count:

| Active rows | Checked bound |
| --- | --- |
| 1 through 127 | Directed local-coefficient calculation; exact outer supports for rows 1 through 4 |
| 128 through 5997 | Eight directed coefficient-Chernoff intervals |
| 5998 through 8192 | Exact-integer geometric tail, below 2^-64 |

The sparse receipt sums to at most 1.7634592448654532e-13. Its recorded
sum was checked against an exact rational sum of all 127 entries.
Fresh directed replays of classes 1, 4, 5, 23, and 127 match the
receipt. The middle calculation independently replays to at most
3.740908268729343e-22. Exact rational addition of the unrounded stored
bounds and the conservative dense-tail bound gives a total below
1.763459790707448e-13, which is strictly below 2^-42.

For completeness, the middle coefficient identity retains the regional
conditioning. If T_w(z) is the local two-state envelope, put

    H(x,z) = sum_w binom(16,w) (16x)^w (1+x)^(16-w) T_w(z).

For a active rows, the coupled regional mixture equals
[x^a]H(x,z)^512/binom(8192,a). Its nonnegative coefficients permit
evaluation at a positive x as an upper bound. Composing sixteen
regions yields the factor binom(8192,a)^(-15), one global division
by 15, and the source factor 16^(-8a). For fixed x,z, the resulting
class bound is log-convex in a. The two interval endpoints therefore
bound every intermediate class. The verifier checks all interval
endpoints and rejects gaps or overlaps.

This certificate concerns independent uniform regional permutations,
independent uniform nonzero input scalars, and independent uniform GL4
refreshes. A deterministic seeded implementation needs a separate
statement about its sampling distribution. No such sampling claim is
part of this audit.

For each realized generator over GF16, extension of scalars preserves
its minimum symbol distance. A nonzero extension-field codeword has a
nonzero GF16 coordinate component, whose support is contained in the
extension-field support. Conversely, embedded GF16 codewords attain
the base-field minimum. Thus the same setup event controls all the
requested fields GF(2^m), for m=8,32,64,128, when their coefficients
use a genuine GF16 embedding. No additional union bound over fields
or coordinate components is needed.

### K18/six-state extension

An independent audit and complete directed replay also passed for the
[K18/six-state profile](GF16_SIZES.md). The audit checked CA=0, the actual
six-state singleton return weights, the rational source-spectrum envelope,
and continuous coverage of all 32768 possible active-row counts.

The eight middle intervals cover 224--23991; the exact 80-bit dense tail
covers 23992--32768. Combined with the 223 sparse classes, the failure bound
is at most 8.272540769181098e-14 < 2^-43. The largest sparse contribution
is near a=57, so this result is not inferred from singleton behavior.

The audit found and corrected a receipt-integrity issue: an earlier verifier
checked interval contributions but ignored the displayed endpoint bounds.
It now replays both endpoints, rejects nonfinite or understated values, and
checks that the displayed contribution covers the entire interval. Regression
tests include zero, missing, nonfinite, and inconsistently inflated endpoints.
The final full replay passed after this correction. No K16 numeric witness
or certificate engine changed.
