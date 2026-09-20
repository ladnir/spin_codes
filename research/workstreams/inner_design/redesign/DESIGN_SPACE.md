# Larger redesigns: interfaces, costs, and rejection tests

The target is the quarter-rate [128,32,32] outer at K=2^20 and N=2^22.
The existing asymmetric t=128,s=19 implementation is the timing control.
Its distance margins are numerical screens; the balanced implementation is
the separate certified control. All new files here are isolated experiments.

## 1. Factored macroblocks

Fix the certified 128-by-19 expansion A0. Initially choose h=2 or h=4 and
set A to a vertical stack of h identical copies of A0. Let t=128h. Choose
a full-rank 19-by-t feedback B with distinct, nonzero weight-three columns.
For each epoch, use r independently sampled transvections from the balanced
ensemble. Write M=M_(r-1)...M_0 for their product.

The forward map starts at state zero and computes Y=X+Aq, q'=Mq+BX.
The transpose starts with zero terminal state and computes V=U+B^T rho,
rho'=A^T U+M^T rho. Apply the transposed factors in reverse order.
The final forward state is not emitted.

If U is partitioned into h chunks, A^T U=A0^T(U_0+...+U_(h-1)). This removes
h-1 applications of the dense feedback circuit per macroblock, but adds
128(h-1) block XORs. The state mixer runs only once per t positions, with r
rounds. The sparse transpose emission still touches every input position.
The implementation keeps the 128-position feedback transform and fully
unrolls the macroblock's fixed sparse emissions.

The A image spectrum is exact: multiply every A0 image weight by h, retaining
its multiplicity. B's image spectrum is separately enumerated over 2^19
states. Repeating the old B columns is deliberately excluded: it would
create weight-two kernel inputs across chunks. Distinct triple columns avoid
weight-one/two/three kernel words but do not prevent many weight-four words.

For r mixer rounds, a fixed nonzero state has marginal law
2^(-r) delta_q + (1-2^(-r)) Uniform(nonzero states). The existing weight-class
transfer therefore supports a first Q=1 screen with the enlarged t and fewer
epochs per region. `MACRO_Q1.json` also tests exact refresh as a diagnostic
ceiling within this transfer family, not as an implemented cheap mixer.

Failure checks include large B-syndrome fibers, cancellation within a long
epoch, persistence of a poorly emitting state, and the unflushed last epoch.
Q=1 passes are not all-occupancy certificates. The retained cover witnesses
cannot be reused as certified numbers after changing t.

## 2. Accumulator cascades

Let C be the length-N binary prefix accumulator: (Cx)_i=sum_(j<=i) x_j.
Its transpose is the suffix accumulator. Let P_j be coordinate permutations,
sampled independently and uniformly for the mathematical ensemble. For
r=2 or 3, the alternative encoder is

    C P_(r-1) ... C P_0 O,

where O is the unchanged block-diagonal outer encoder. Its transpose applies
suffix scans and inverse permutations in reverse order, followed by O^T.
This replaces the SPIN inner and its particular permutation ensemble; it is
not merely a faster implementation of the same code.

Each scan needs one XOR per position, but intervening permutations may require
additional full-buffer traffic. Precomputed schedules remove online permutation
generation at the cost of setup storage. The existing binary RAA comparison
uses a different outer and on-demand permutations, so its timings do not
settle the performance of this proposed BCH-accumulator cascade.

For a fixed input weight a>0, the exact number of accumulator inputs that
produce output weight b is

    T_N(a,b) = binom(N-b,floor(a/2)) binom(b-1,ceil(a/2)-1).

Use zero for out-of-range binomial arguments and T_N(0,b)=1_{b=0}. To derive
the expression, decompose the output into runs of ones. If a is even, the
output ends in zero; if a is odd, it ends in one. The two binomial factors
count positive run lengths and permissible intervening zero lengths.

After an independent uniform permutation, the weight transition probability
is T_N(a,b)/binom(N,a). Thus a rigorous first-moment route exists: form the
outer weight enumerator, compose these transitions, then sum expected output
counts through the bad-weight cutoff. This does not require codewords to
have independent outputs under one shared setup.

The unresolved issue is efficient, outward evaluation at N=2^22, together
with a performant implementation of the required permutations. Dense N-by-N
tables are not a practical certificate strategy. Near-canceling pairs and
terminal runs must be included, not dismissed through typical-weight arguments.
`test_designs.py` checks the transition counts exhaustively through N=10 and
tests the two/three-pass adjoints on all bilinear basis pairs at N=12.

## 3. Routing and locality

Start by distinguishing an exact schedule change from a new code. Changing
the order of reads, writes, and outer transforms can preserve the existing
linear map and permutation distribution. Restricting which coordinates can
cross cache chunks changes the ensemble and requires a new distance proof.

The first necessary rejection test is a dimension argument. Suppose a chunk
contains k_local independent outer-message coordinates and n_local output
coordinates. Suppose every influence of those messages outside the chunk
factors through a linear interface of dimension at most r. If k_local>r,
some nonzero local message annihilates that interface. Its codeword is then
confined to the chunk, so d_min<=n_local (or the encoder is not injective).

This also rules out a tempting design with a single global s-bit state,
local outer rows, and coordinate permutations confined to each chunk.
With zero incoming state, the outgoing state has rank at most s. If the
chunk's message dimension exceeds s, a nonzero local message returns the
state to zero and cannot affect later chunks. Carrying the state rather
than resetting it does not remove this obstruction.

For example, eight equal chunks at K=2^20 have n_local/N=1/8. This design
cannot attain 16.5% or 19% relative distance. Multiple local passes linked
only by a few short state vectors must pass the same total-interface-rank
test. Actual coordinate exchange between chunks can avoid this obstruction,
but its communication and memory cost must be measured.

Candidate directions are therefore exact scheduling changes first, or a
multistage block permutation that moves enough independent coordinates between
chunks. For the latter, specify the exact distribution before borrowing any
SPIN occupancy argument. Graph connectivity alone is not a distance proof.

## Performance gates

`generate_headroom.py` produces explicitly invalid-distance controls: remove
the inner while retaining scalar or unrolled routing, or retain only the
contiguous outer transform. They preserve the measured N-to-K buffer geometry
but do not preserve the code. They estimate headroom and expose scheduling
effects; they are not rigorous runtime lower bounds or additive phase timings.

All timed processes use the shared benchmark lock and run serially. A new
construction advances only after correctness tests and structural screening.
The aspiration of a 20–30% improvement is a selection criterion, not a forecast.
