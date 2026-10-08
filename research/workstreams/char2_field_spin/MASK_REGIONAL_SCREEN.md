# Actual-mask regional screen with a pooled fallback

This screen does not close the sparse gap for the higher-rate GRS[32,22]
construction. At K20, a bounded optimization with exact zero-, one-, and
two-input mask kernels gives margins **70.870, -49.025, and -249.904 bits**
for 4, 8, and 16 active source rows. The latter two bounds are vacuous.
All reported margins use ordinary floating arithmetic, not directed bounds.

The construction is unchanged: GF256, width 16, four residue-state coordinates,
the four fixed Cauchy phases, eta=2, fresh nonzero input and state diagonals,
and one regional routing stage. This note concerns a proof calculation,
not a native implementation or a performance measurement. It neither
re-proves the existing half-rate result nor claims a whole-code certificate.

## State interface and conditioning

The first sixteen states are the actual coordinate-support masks
H subset {0,1,2,3}. Mask 0 means actual zero state. The seventeenth state,
P, records only that the actual state is nonzero.

Fix the source support and the common equation history before a packet.
The current packet diagonals remain fresh. For input weight w at most the
chosen cutoff, known masks use `mask_state_rank.local_mask_transfer`.
That kernel includes the equations for the next mask's zero coordinates
and requires its other coordinates to remain nonzero. Its entries are
rank-potential bounds, not fixed-word transition probabilities.

For larger w, a known mask uses the uniform zero/nonzero kernel and forgets
its nonzero mask. The zero destination remains known. Every other destination
is P. A packet starting in P uses this uniform bound at every weight;
P can leave for mask 0 but never for a known nonzero mask. Empty packets
preserve their incoming mask, and an empty packet starting in P stays in P.

The uniform return entry is the minimum of two independently valid
whole-packet bounds: the old return entry and, for w=1,2,3, the audited
`rank_return_moment` with input placements averaged. The code uses the
selected-zero **rank** API from `cauchy_joint_return.py`. It does not reuse
current diagonals in a fixed-word return probability and a rank bound.

Each packet matrix is averaged over its C(16,w) input placements. For each
w, the screen then takes an entrywise maximum over the four phase matrices.
This replaces the actual phase schedule conservatively; it does not average
phases. Denote the resulting 17 by 17 matrix by T_w.

The exact-mask bound can split one conditional history among several
possible successor masks. Its sum can exceed the old uniform bound.
Consequently, retaining more mask information need not improve this linear
upper bound. An entrywise minimum with an old nonzero marginal would not
justify distributing that marginal among masks.

## Exact regional support averaging

For a block of n slots, define the matrix F_h by summing its packet products
over all h-element slot subsets and dividing by C(n,h). A single packet has
F_w=T_w; its coefficient must not receive another factor C(16,w).

For successive blocks of lengths n1 and n2, the normalized product is

\[
 F^{(12)}_h=
 \sum_i
 \frac{\binom{n_1}{i}\binom{n_2}{h-i}}
      {\binom{n_1+n_2}{h}}
 F^{(1)}_i F^{(2)}_{h-i}.
\]

The matrix order follows the packet order. Binary powering computes the
regional coefficients, truncated at h=A. Such truncation preserves every
coefficient subsequently used. Python integer binomials and exact rational
ratios avoid overflowing raw support counts; the screen converts each
ratio to float only for matrix arithmetic.

Fix A selected source rows out of L. Independent regional permutations send
these rows to a uniform A-subset in each region. Allowing an arbitrary subset
of those selected rows to be active gives the regional matrix

\[
 R=\sum_{h=0}^{A}\binom Ah\tau^h F_h.
\]

This is sampling without replacement, not a Bernoulli approximation.
The relaxation includes infeasible outer supports; adding them is
conservative. Actual supports using exactly A source rows have shortened
outer dimension at most W-10A, where W is their total outer weight.
The common-rank source potential is therefore bounded by

\[
 \binom LA\tau^{-10A-1}y^{-e}
 e_0^T R^{32}{\bf1},
 \qquad 0<y\le1,\quad 1\le\tau\le255.
\]

Here e is the largest forbidden output weight, e0 selects the initial zero
state, and the terminal sum permits every final state. The implementation
rescales R and each terminal-vector multiplication. It rejects nonfinite
packet values and terminal underflow; its screen guard also limits regional
attenuation. These guards do not turn floating evaluation into a certificate.

## Bounded numerical outcome

K20 uses L=47,680 rows, N=1,525,760 outer symbols, and e=152,575. The frozen
witnesses are parameterized by
`y=exp(-exp(logtilt)*A/N)` and `tau=exp(logtau)`.

| A | Frozen logtilt | Frozen logtau | Cutoff 1 margin | Cutoff 2 margin |
|---:|---:|---:|---:|---:|
| 4 | 5.32479700865 | 4.96215263668 | 69.957897 | 70.727203 |
| 8 | 5.13883585678 | 2.88863774319 | -53.709356 | -53.601286 |
| 16 | 5.07867455120 | 2.33096175675 | -265.752691 | -267.242819 |

One L-BFGS-B start per A, capped at 30 iterations, gave:

| A | Optimized logtilt | Optimized logtau | Cutoff 2 margin | Iterations |
|---:|---:|---:|---:|---:|
| 4 | 5.365146462364348 | 5.073732175001389 | 70.869910 | 5 |
| 8 | 5.306013970826702 | 3.065941017050765 | -49.025464 | 8 |
| 16 | 5.315550552465184 | 2.561780602656156 | -249.903705 | 5 |

All three local searches reported convergence. They do not establish global
optimality. The preceding two-state screen gave approximately 71.98, -19.294,
and -143.748 bits, respectively. Thus this mask composition is not a useful
replacement for that screen, even though its local kernels retain more
information. No additional parameter search was performed.

## Replay and audit

From `research/workstreams/char2_field_spin/theory`:

```text
python -B -m unittest test_mask_regional_screen.py
python -B mask_regional_screen.py --active 4 8 16 --cutoff 1
python -B mask_regional_screen.py --active 4 8 16 --cutoff 2 --optimize --maxiter 30
```

Seven tests pass. They compare exact Fraction convolutions with independent
enumeration of tiny slot subsets, test noncommuting matrix order, verify the
selected-row regional factor, and check pooled-state routing, phase maxima,
and invalid inputs. The existing local mask tests remain separate.

An independent agent audited the composition and source factor. It also
checked 144 small-field successor-mask cases with dependent state/input
forms, enumerating all current diagonals and source vectors. Those checks
support the conditional mask interface; they do not certify the floating
whole-stream values.

## Limited next step

The alternatives below were proposed before the [realized-route
obstruction](WINDOW_DIMENSION_SEARCH.md) was found. They are preserved
proof ideas, not the active next step for the unchanged four-state
GRS[32,22] construction. Its geometry must be reconsidered first.

A fixed policy can avoid the mask-splitting loss without redistributing P.
Retain the known full mask after a singleton birth and through subsequent
empty packets. Pool it at the next active packet, using the old uniform
row. All other births enter P. Only states 0, full, and P are reachable.
After collapsing full and P to nonzero, active transitions match the old
two-state bound; the exact full-mask empty contraction is no larger than
the old uniform empty contraction. This policy therefore cannot worsen
that bound. Its expected gain is limited to empty stretches after births;
it has not been screened here.

A stronger continuation-dependent minimum would compare two entire covers,
not their individual destination entries. That is a separate nonlinear
proof calculation. The present negative result alone does not justify
expanding the search or implementing the construction.
