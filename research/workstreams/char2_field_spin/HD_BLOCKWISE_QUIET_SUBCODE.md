# A deterministic quiet-subcode distance cap

The width-64 identity-feedthrough construction is rejected at state dimension 4
or 8, independently of the outer length, route, and state mixers. At K=2^20 and
rate 1/2, every setup has a nonzero codeword of relative weight at most 12.5% or 25%,
respectively. Thus neither parameter choice can meet 40% distance with any setup
failure margin. State dimension 16 has a 50% upper cap; this obstruction does not
reject it at 40%, but does not certify it.

This result supersedes the weaker state4/state8 conclusions in
[HD_BLOCKWISE_PARTIAL.md](HD_BLOCKWISE_PARTIAL.md). Unlike that note's window
events, the present obstruction holds for every routing permutation.

## Model and cap

Work over a field F. Encode R independent message rows by an injective linear
[b,k] outer code, then apply an arbitrary permutation of its N=Rb coordinates.
The message dimension is K=Rk. Partition the routed input into P=N/t complete
packets X_i in F^t. The inner map is

    Y_i = X_i + A_i Q_i,
    Q_(i+1) = M_i Q_i + C_i X_i,     Q_0 = 0,

where Q_i lies in F^s and the matrices have the indicated compatible dimensions.
The matrices may vary with i. The encoder emits Y_0 through Y_(P-1), and does
not emit or constrain the final state Q_P. There is no additional source
transform between the independent outer rows and the support-preserving route.

For every integer a with 1<=a<=R and

    ak > s(P-1),

the code has a nonzero word of weight at most ab. In particular, if

    a_* = floor(s(P-1)/k)+1 <= R,

then its relative minimum distance satisfies

    d_min/N <= a_*/R.

If a_*>R, this particular argument makes no nonzero-word claim. Equality of the
upper cap with a requested distance does not reject an "at least" distance
requirement; strict inequality does.

### Proof

Choose any a source rows, and set all other message rows to zero. Their message
space V has dimension ak. Impose the homogeneous linear equations

    C_i X_i = 0,     i=0,...,P-2.

There are at most s(P-1) independent equations on V. Their common kernel
therefore has dimension at least ak-s(P-1)>0. Choose a nonzero message from
this kernel. Starting from Q_0=0, induction gives Q_i=0 for every emitted packet
i=0,...,P-1. Consequently Y_i=X_i throughout the emitted word. The last feedback
syndrome may make Q_P nonzero; it does not affect any emitted packet.

Only the selected a outer rows can contain nonzero coordinates, so the routed
input has weight at most ab. It is nonzero because the outer map is injective
and the route is a permutation. Since Y=X, the emitted word is also nonzero.
This proves the claim without assuming that the full inner map is injective.

The proof uses neither the MDS property nor CA=0. It does not require M_i to be
invertible or random. Arbitrary choices of A_i, C_i, and M_i cannot eliminate
this subcode. If each record comprises several field coordinates with shared
coefficients, using the witness in one coordinate gives the same record-weight
cap.

## Exact K20 parameters

Here K=1048576, N=2097152, t=64, and P=32768. For the two outer codes under
consideration, the exact bounds are:

| Outer [b,k] | R | s | Selected rows a_* | Guaranteed dimension | Weight cap | Relative cap |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| [128,64] | 16384 | 4 | 2048 | 4 | 262144 | 12.5% |
| [128,64] | 16384 | 8 | 4096 | 8 | 524288 | 25% |
| [128,64] | 16384 | 16 | 8192 | 16 | 1048576 | 50% |
| [256,128] | 8192 | 4 | 1024 | 4 | 262144 | 12.5% |
| [256,128] | 8192 | 8 | 2048 | 8 | 524288 | 25% |
| [256,128] | 8192 | 16 | 4096 | 16 | 1048576 | 50% |

For either outer code at t=64, the smallest integer state dimension not rejected
by this cap is 13 for 40% distance, 15 for 45%, and 16 for 50%. These are necessary
thresholds, not sufficient state dimensions. Increasing b at the same outer
rate does not repair the obstruction. Up to the exact rounding above, an outer
rate r=k/b gives the cap s/(rt), which is 2s/t at rate 1/2.

No setup conditioning, refreshed state mixer, or sharper probabilistic proof
can establish 40% distance for the rejected width-64/state-4 or state-8 geometry:
its distance failure probability is 1. A global quiet-subcode dimension count
alone did not previously give this conclusion. Restricting the message to a
small set of outer rows is the step that also controls the witness's support.

## Two smaller-packet alternatives

Reducing packet width together with state dimension changes the necessary
state-to-packet ratio. The following two alternatives both have an exact 50%
quiet-subcode cap:

| Outer | t | s | P | a_* | Quiet dimension lower bound | Minimum s not rejected at 40% |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| [128,64] | 32 | 8 | 65536 | 8192 | 8 | 7 |
| [256,128] | 16 | 4 | 131072 | 4096 | 4 | 4 |

As a second necessary screen, apply the exact balanced-window law from
[HD_BLOCKWISE_GEOMETRY.md](HD_BLOCKWISE_GEOMETRY.md). Choose a suffix in each
even-indexed region and a prefix in the following region. Adjacent windows
merge, leaving b/2 connected components, each with one state-exit constraint.
The windows use a total of 838848 symbols, strictly below 40% of N, and all
endpoints are packet-aligned. Independent uniform regional permutations give
the following exact-rational calculations, with displayed means rounded:

| Outer, t, s | Expected shortened dimension | Exit constraints | Failure lower bound | Minimum s not rejected by this window gate at 40 bits |
| --- | ---: | ---: | ---: | ---: |
| [128,64],32,8 | 344.6784748603134 | 512 | 0 | 6 |
| [256,128],16,4 | 10.49452346472164 | 512 | 0 | 1 |

The window probability bound is the earlier lower bound
max(0,(E[S]-c)/(T/2-c)), where c is the exit-constraint count and T is the
window weight. Its zero value is not a distance certificate. It only says
that these two necessary screens do not reject the proposed alternatives.

The next useful step is to examine these changed packet geometries with the
local emission/rank bounds and then the full source-aware composition. No
native performance prediction follows from either screen. In particular,
outer length 256 alone was not the repair: increasing s/t from 1/8 or 1/16 to 1/4 is what
removes the current deterministic rejection.

## Boundaries and possible changes

Imposing Q_P=0 as well adds at most s equations and replaces s(P-1) by sP.
The helper reports this optional stronger restriction separately. An encoder
that appends flush outputs requires a separate count of its output length and
weight; the free-terminal formula must not be applied to those outputs.

Replacing X_i in the output by an invertible diagonal or monomial image does
not repair the support bound, since such a map preserves packet weight when
Q_i=0. A dense input-to-output map can change that step, but is not thereby
proved to meet the target distance. Other possible changes must likewise
address the localized quiet subcode rather than only improve state mixing.

## Reproducible checks

The helper is [theory/hd_blockwise_quiet_subcode.py](theory/hd_blockwise_quiet_subcode.py),
with tests in [theory/hd_blockwise_quiet_subcode_tests.py](theory/hd_blockwise_quiet_subcode_tests.py).
The tests cover exact K20 arithmetic, target inequalities, rounding, free versus
constrained terminal state, the one-packet case, zero state dimension, exhausted
source dimension, both alternative window gates, and a finite-field witness.

The GF4 witness uses eight outer [4,2] rows, selects four of them, routes the
coordinates by fixed regional permutations, and uses t=4 and s=1. Seven
feedback equations on eight message variables have nullity 1. The computed
nonzero message emits a word of weight 12<=16 with all eight emitted states
zero and terminal state 2. The time-varying maps have CA nonzero in every packet,
explicitly checking that CA=0 and a zero terminal state are unnecessary.

From the theory directory, run:

    python -B -m unittest hd_blockwise_quiet_subcode_tests
    python -B hd_blockwise_quiet_subcode.py

The 11 tests pass. This work changes no native encoder or benchmark files.
