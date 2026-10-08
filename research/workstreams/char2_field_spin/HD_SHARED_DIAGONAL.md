# Sharing a diagonal between state mixing and feedback

The systematic RS inner spends two field multiplications on terms that
can share one final diagonal. Correlating their setup changes the
construction, so the independent-map proof does not apply directly.
The following joint bound recovers its operator with a small extra loss.

## Construction

Keep the outer, route, field F=GF65536, and fixed maps A,C from
[the systematic RS construction](HD_SYSTEMATIC_MDS.md). Write q=65536,
s=16, d=q^s-1, and let e be the all-one vector in F^s. The fixed H still
interpolates from points 0,...,15 to points 16,...,31 in the AES subfield.

For each packet, sample D1 and D2 with independent uniform nonzero diagonal
entries. Put u=D2 e. Sample v uniformly from {v in F^s: v^T u != 1}.
The joint tuple (D1,D2,v) is independent across packets and independent
of the outer and route. Fix all setup choices before encoding any message.
Starting with Q=0, compute

    Y = X + A Q,
    Q_next = Q + D2 (H D1 C X + e (v^T Q)).

Equivalently M=I+u v^T and T=D2 H D1, but M and T share D2.
M is invertible because v^T u != 1. The encoder is F-linear and its
identity feedthrough makes it invertible as a map on the inner input.
There is no restriction on the terminal state.

The fused update saves s GF65536 multiplications per packet compared with
separately computing u(v^T Q) and the final feedback scaling. It leaves
the expensive fixed-map circuits unchanged. Native setup must explicitly
tie u to D2; applying the fused kernel to independent maps would be wrong.

## A bound for the joint setup

We first analyze unconditioned v sampled from all F^s. For each fixed D2,
the acceptance event v^T D2 e != 1 has probability 1-1/q. Consequently
conditioning increases any nonnegative expectation by at most
a=q/(q-1). It leaves the marginal distribution of D1,D2 unchanged.

Fix a nonzero incoming state Q and any feedback vector R=CX. In the
unconditioned experiment, tau=v^T Q is uniform in F and independent of
D1,D2. We bound the distribution of

    W = D2 (H D1 R + e tau),       Q_next = Q + W.

Every square minor of H is nonzero. Also every restriction to w+1 rows
of w columns of H together with e has full rank, for w<s. To see the
second assertion, H[i,j]=c/(16+i+j) for one nonzero constant c. A
linear combination of those w columns and e is a rational function
with numerator degree at most w. If it vanishes at w+1 distinct output
points, its numerator is zero. Its residues and constant coefficient
then force every coefficient of the combination to vanish. The constant
c is common to all entries because the input and output sets are additive
cosets: the subspace polynomial and its derivative give constant factors.

Suppose R != 0 and its support has size w. Then D1 R is uniform over
nonzero entries on that support. Let L=H D1 R+e tau. For a specified
set of z<=w zero coordinates of L, first fix tau and all but z pivot
coordinates of D1 R. The corresponding H minor fixes the pivots, so
the probability is at most (q-1)^(-z). If z>w, the augmented rank
argument makes that event impossible. In particular L=0 is impossible
when w<s; when w=s, its probability is at most (q-1)^(-s).

Fix a nonzero target h with z zero coordinates. Conditional on L having
exactly its support, D2 L=h with probability (q-1)^(-(s-z)). Other
supports cannot reach h. Combining the two bounds, including the
separate target h=0, gives the pointwise bound

    Pr[W=h] <= (q-1)^(-s) = b/d,    b=d/(q-1)^s.

This holds for every h, including zero. It is a density upper bound,
not a statement that W is uniform or never zero.

When R=0, tau=0 gives W=0 with probability 1/q. Otherwise W is uniform
over full-support vectors. Thus for every nonzero h,

    Pr[W=h] = (1-1/q) (q-1)^(-s) = (1-1/q)b/d.

After the acceptance conditioning, R=0 cannot send a nonzero Q to zero,
because the resulting M is invertible.

## Transfer to the existing complete certificate

Use K_j,T_j,V_j and the fixed-state quiet cap Gamma_j from the systematic
RS proof. Let the incoming weighted measure be bounded by mass Z at zero
and density B/d on each nonzero state. For the zero-state row, the joint
setup gives bounds (K_j,b(T_j-K_j)): tau is zero and the old diagonal-
Cauchy lemma applies. Acceptance does not change these diagonal marginals.

For nonzero incoming states and j>0, the preceding pointwise bounds give
at most a b V_j/d for the outgoing zero mass. The outgoing nonzero
density coefficient is at most a(b V_j+Gamma_j/q). The first term sums
the dense part against the incoming flat envelope; the second term is
the quiet event R=0 and tau=0, which leaves that particular state fixed.
For j=0, replace the two bounds by zero and
a(b(1-1/q)V_0+Gamma_0/q), respectively.

Let P_j be the original ideal two-component operator, before its
diagonal-Cauchy density transfer. Each entry of the new operator is at
most a b times the corresponding entry of P_j, including the exact
zero-to-nonzero and nonzero-to-zero entries when j=0. Positive operator
products, regional coefficient extraction, and the outer union bound
preserve this comparison. For a witness exponent alpha<=1, the loss is
at most n_p log2(a b), with n_p=32768 packets. The dense inverse-output
bound is unchanged.

The existing complete certificate already pays n_p log2(b). Sharing the
diagonal therefore costs at most an additional n_p log2(q/(q-1)), less
than 0.722 bits. The conservatively rounded source bound -50.06 implies
more than 49.338 bits of whole-code margin at the same 40% distance.
[The directed helper](theory/hd_shared_diagonal.py) can either use that
rounded source bound or replay every underlying rational witness first.

The distance statement is over the entire sampled setup and concerns the
existence of any nonzero message after that setup is known. Eight parallel
field lanes in a 128-bit record reuse the same scalar-code guarantee.
The deterministic native benchmark seed is not itself a sampled-setup
certificate. Performance and native validation are separate obligations.

## Checked checkpoint

The independent 384-bit replay covers all 8192 active-row counts and gives
log2(failure) < -49.33892764811194. Exact enumeration over GF4 tests every
incoming state and feedback vector, including zero, against the pointwise
bounds before and after conditioning. Separate exact-field checks confirm
the concrete H formula and selected augmented minors. The universal rank
claim follows from the rational-function argument, not these sampled tests.

The isolated `spin_char2_hd_sumrank_forward` target provides modes `shared`
and `shared-striped`. Both agree with the independent dense scalar oracle
on all 80 packet basis directions, both field components, random nonzero
states, routed offsets, full encoders, and exact aliasing. The striped
variant takes approximately 6.1 ms at K20 on the established Peach CPU2
4.5-GHz harness, with precomputed setup and 128-bit records. A 501-call
comparison gives 6.094/6.168/6.042 ms for shared-striped/independent/
shared-striped; the independent-map control is `mds-affine`.
The shared variants' checksum is 4678095192827670896.
No production-library profile or default is changed.

All 25 native research variants pass Release and ASAN/UBSAN checks. The
high-distance regression suite passes 397 tests; all four focused joint-
refresh tests pass after adding a complete weighted-operator check on the
small GF4 packet. The validated executable SHA256 is
`46cf36539f86e8d004fe88df78a998088d6429f4599c707b341e9f0bf067a288`.
Reproduce timing with `--bench 501 full 20 shared-striped`.
