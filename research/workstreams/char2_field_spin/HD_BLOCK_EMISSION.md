# Cheaper emission maps: proof and cost screen

The alternating-feedback checkpoint still spends about 3.17 ms in the
routed inner, versus 2.92 ms after routing is removed. This experiment
reduces fixed-map arithmetic on both sides of the recurrence. It does not
yet provide a whole-code certificate.

## Construction and conditional bounds

Keep the outer, regional permutations, and shared-diagonal packet refresh
from [the certified reference](HD_SHARED_DIAGONAL.md). Partition each
64-symbol packet into b=16/h blocks of length L=4h. Set A to the direct
sum of b RS[L,h] generators and set C=A^T. Each local generator evaluates
the normalized degree-below-h polynomial basis at nodes 0,...,L-1.
The state still contains sixteen GF65536 coordinates. The update remains

    Y = X + A Q,
    Q_next = Q + D2 (H D1 C X + 1 (v^T Q)).

The native candidate uses h=8: two RS[32,8] maps. Unlike global RS[64,16]
emission, its minimum nonzero emitted weight is 25, not 49. This weaker
floor must be included in the proof; the old certificate does not transfer.

Fix a potential input set U with block occupancies l_1,...,l_b and let
j=sum l_i. In the outer comparison experiment X is uniform on F^U,
including zero coordinate values. For a mark 0<z<1, define

    m_h(n,z) = sum_(c=0)^n binom(n,c) z^(n-c)(1-z)^c q^-min(h,c),
    k_h(l,z) = sum_(c=0)^l binom(l,c) z^(l-c)(1-z)^c q^(-c-min(h,l-c)).

Independent state blocks and their MDS ranks give

    K_l = product_i k_h(l_i,z),
    T_j = (z+(1-z)/q)^j,
    V_l = T_j [q^16 product_i m_h(L-l_i,z) - 1] / (q^16-1).

Here V_l averages over a uniform nonzero incoming state. The subtracted
one removes the all-zero state; its inactive output coordinates all have
zero weight. These formulas follow by expanding z^wt as a positive sum
over coordinate-zero constraints.

For a fixed nonzero incoming state, at least one state block is nonzero.
Its contribution to the output has weight at least
max(0,L-h+1-max_i l_i). With rho_l=sum_i min(h,l_i), the quiet cap is

    Gamma_l = min(K_l, q^-rho_l z^max(0,L-h+1-max_i l_i), (q^16-1)V_l).

The homogeneous/affine zero-constraint argument proves the first bound.
The second combines the syndrome rank with the output floor. The third
keeps one term from the nonnegative sum defining V_l. Substitution into
the [block-feedback operator](HD_BLOCK_FEEDBACK.md) gives a valid local
comparison, with the same shared-refresh density loss. Composition and
complete outer-row coverage remain separate obligations.

The screen retains every occupancy profile. Feedback-only grouping of
small occupancies is invalid for emission: equal quiet moments need not
give equal V_l. The implementation also expands the subtraction in V_l
as positive terms to avoid numerical cancellation. Exhaustive GF4 tests
check every input support and nonzero state. A separate exact-rational
comparison checks the FLOAT paired coefficient builder.

## Initial screens and next checks

Paired-regional FLOAT screens at active rows 128/512/2048 give about
-220/+648/+1850 for h=4, and -254/-210/-4748 for h=8. The positive
values are failed upper bounds, not counterexamples. These runs use
the actual weaker emission map and local feedback. They are not a
certificate or a reason to change the selected profile.

The h=8 native mode is `bucket8-both`. It removes the highest emission
FFT level and keeps the explicitly unrolled striped kernel. The outer,
route, and random refresh are unchanged. Its dense-map digest is
9075152494480514616; correctness checks and timing must precede selection.

Release and ASAN/UBSAN checks now pass for all 33 native modes. The h=8
mode takes 5.717/5.786 ms around 5.717/5.714-ms frozen controls, so it is
not selected. Further h=8 screens give -2370/+80/-1695/-8096 at active
rows16/256/1024/3072. In particular, it is not a complete certificate.

The h=4 mode, `bucket4-both`, takes 5.484/5.506 ms around a 5.719-ms
control. Its checksum is 8134999575494699104 and its dense-map digest is
5069833586872545664. This is a validated cost candidate, not a certified
replacement. It removes two emission FFT levels and retains the four-check
feedback circuit.

## Tracking sparse states

The next screen retains the number of nonzero state blocks. Component u
has density at most B_u/d_u on each state with exactly u nonzero blocks,
where d_u=binom(b,u)(q^h-1)^u. Component zero is a point mass. For fixed
input occupancies, V_u averages the marked output over that state class.
It is the degree-u elementary symmetric polynomial in the local nonzero-
state moments, divided by binom(b,u), and multiplied by T_j.

Let f_u be the sum of the u smallest values max(0,L-h+1-l_i). Then
Gamma_u=min(K_l,q^-rho_l z^f_u,d_u V_u) is a quiet cap for every state
in class u. Before the common density loss, the transition bounds are

    P[0,0] = K_l,
    P[0,v] = (d_v/d)(T_j-K_l),
    P[u,0] = 1[j>0] V_u/d,
    P[u,v] = (d_v/d)(1[j>0]+1[j=0](1-1/q)) V_u
             + 1[u=v] Gamma_u/q,              u,v > 0.

Here d=q^16-1. The dense terms use the existing pointwise refresh bound;
the last term retains the incoming state and hence its class. Exact GF4
enumeration checks the class moments and quiet caps. The FLOAT module
offers packetwise fractional powers and a Jensen bound after multiplying
two profile-averaged packet matrices. Tests check the latter bound and
agreement at exponent one.

Neither screen currently closes the sparse middle. At rows256/512, the
packetwise bounds are about +1989/+2973 and the averaged-pair bounds are
+1113/+4481. These are failed bounds. Retaining classes introduces more
state paths before the fractional-power inequality; it is not a free
tightening. Longer or profile-preserving composition remains unexplored.

## Rank-two refresh

A different tradeoff is to spend one additional random dot product while
keeping the cheaper emission. Let E have columns 1 and x, where x lists
the sixteen distinct output nodes of H. Set U=D2 E. Sample a two-row
matrix V uniformly, rejecting when det(I_2+V U)=0. Use

    Q_next = Q + D2 (H D1 C X + E V Q).

The complete certificate and native follow-up are recorded below. The
determinant condition makes I_s+U V invertible. Since
U has rank two, V U is uniform over two-by-two matrices before rejection.
The acceptance probability is (1-1/q)(1-1/q^2).

For a fixed nonzero Q before rejection, VQ is uniform on F^2. If CX=0,
the identity branch has probability q^-2 instead of q^-1. A nonzero
linear polynomial on the output nodes has at most one zero. Combining
that fact with D2's independent nonzero scalars bounds every nonzero
target by (1-q^-2)(q-1)^-s. For CX!=0, the Cauchy columns together with
1,x have every required row rank: multiplication by their pole product
reduces a dependence to a polynomial of degree at most w+1. The same
pivot argument therefore bounds every target by (q-1)^-s. Conditioning
costs the reciprocal acceptance probability.

These local arguments suggest replacing q^-1 by q^-2 in the ideal quiet
term, with only a tiny additional density loss. `--lazy-rank 2` explores
that operator. Exact GF8 enumeration checks every feedback support and
the complete unconditioned pointwise law. Exact GF4 enumeration checks
the acceptance probability for every nonzero diagonal. Selected concrete
augmented minors also pass. These tests supplement the local argument;
they do not replace a complete directed certificate or native timing.

The initial rank-two paired-regional screens give -133/-261/+515 at
rows256/512/2048. The extra mixing repairs the sampled sparse-middle
obstruction but does not close the denser four-check feedback class.
Do not report this variant as certified or assume its extra arithmetic
preserves the 5.5-ms rank-one timing.

The [larger-feedback follow-up](HD_RETURN_DESIGN_SCREEN.md) now implements
and tests the rank-two construction. Native mode38 retains four RS[16,4]
emission blocks but uses two RS[32,8] feedback blocks. It takes about5.6 ms.
The follow-up records a correction to its grouped transition bound; earlier
partial receipts are withdrawn. Corrected complete 256-bit coverage now
gives more than 50.9966 bits at 40% distance. The 133 rational witnesses
are pinned, and a separate 384-bit replay confirms all 8192 classes and
the same margin. Streaming output
reduces the measured cost to about 5.4--5.5 ms.
The earlier pure four-check screen and its failure are preserved above.

## Further structural direction: transform the carried state

Another unproved design could avoid repeatedly carrying an exceptional
state unchanged. The expression

    Y = X + A Q,
    Q_next = D2 H (Q + D1 C X)

would reuse the existing H circuit for both carried state and feedback,
while removing the separate random dot product. The incoming-state map
D2 H is invertible, but this is not the rank-one construction and none
of its complete certificates transfers automatically.

There are two useful structural facts to investigate. First, [I;H] is MDS,
so wt(Q)+wt(HQ)>=17 for nonzero Q. For four-coordinate state blocks,
at least five blocks are active across those two states. Without input
interference, two block-emission packets therefore have combined weight
at least 5*13=65. Second, C A=0 for these additive-subspace RS blocks,
so I+A D1 C is invertible for every D1. This may help analyze a return
to the zero state. Neither fact controls an entire noisy state history;
the weighted transition with correlated state and feedback is the missing
proof obligation. The native mode `bucket4-carry` now passes Release and
ASAN/UBSAN checks and measures 5.367/5.373 ms. It remains uncertified.

A [common-refresh variant](HD_COMMON_REFRESH.md) instead adds Q before
D1. It has a simpler local proof and measures 5.318/5.328 ms, but its
current whole-code screens still fail. Neither replaces the certified
5.7--5.8-ms alternating-feedback checkpoint.
