# Why independent expansion and feedback remain tractable

This note derives the changed local bound. The accompanying outward numerical
certificate is described in `CERTIFICATE_RESULT.md`. The instance retains the quarter-rate [128,32,32] outer,
K=2^20, N=2^22, and the baseline SPIN permutation distribution.

## Construction and remaining randomness

Fix full-rank binary maps A:F_2^s -> F_2^t and B:F_2^t -> F_2^s,
with t=128 and s=19. The forward inner has state q_0=0 and recurrence

    Y_i = X_i + A q_i,
    q_(i+1) = M_i q_i + B X_i.

For each epoch, sample u from the nonzero s-bit vectors. Conditional on u,
sample v from its orthogonal hyperplane, including zero. Set M=I+uv^T.
The pairs are independent across epochs and independent of the permutations.
These are distributions of a single shared linear setup, not per-message
fresh coins. The distance argument unions over nonzero outer messages.

For each fixed nonzero q, this mixer has the exact marginal

    Mq ~ (1/2) delta_q + (1/2) Uniform(F_2^s \ {0}).

The baseline transvection argument proves this identity. It concerns M alone
and imposes no relation between A and B. When q=0, Mq=0 deterministically.
Conditioning on earlier epochs leaves the current mixer fresh.

The adjoint, initialized with zero reverse terminal state, is

    V_i = U_i + B^T r_(i+1),
    r_i = A^T U_i + M_i^T r_(i+1).

In particular A^T acts on U_i, not V_i. Sparse B therefore accelerates
transpose emission without replacing the dense, balanced state expansion A.
The forward inner is invertible because it is block triangular with identity
diagonal blocks. Neither this fact nor the adjoint identity requires BA=0.

## Separate the two spectra

Write a_v for the number of nonzero states q with wt(Aq)=v. Write b_w for
the number of nonzero dual states a with wt(B^T a)=w. Both sum to
m=2^s-1. The retained A has image weights 48,56,64,72,80.

The transfer representation upper-bounds a weighted state measure using:

- Z: mass at zero;
- D: arbitrary mass supported on nonzero states;
- C_v: mass dominated pointwise by a uniform measure on the a_v states of
  image weight v.

More precisely, a nonnegative coordinate vector (Z,D,(C_v)_v) represents an
upper bound of the form Z delta_0 + nu + sum_v C_v U_v, where nu has mass
at most D and U_v is the uniform probability measure on that shell.
Terminal total mass is at most the sum of the coordinates.

For uniform X of fixed input weight j, the emitted moment at state q depends
only on wt(Aq). In contrast, the counts of X with BX=0, and the sizes of
the other syndrome fibers, depend on B. `general_occupancies.epochs` remains
applicable when given the A spectrum, B kernel spectrum, and B fiber caps
separately. MacWilliams inversion of the B^T spectrum gives the kernel counts.

Lazy cancellation occurs when BX=q. For one or two input positions, the
new screen explicitly enumerates wt(X+A BX), grouped by BX and wt(A BX).
This replaces the old symmetric-map cancellation data. No injectivity of
the pair-syndrome map is assumed.

For larger j, the generic transfer uses integer fiber caps and the inequality
wt(X+Aq) >= |j-wt(Aq)|. Packing, Fourier triangle, and Parseval bounds use
the B kernel and B^T spectrum. If every weight-j input has zero syndrome,
the lazy branch preserves its A shell. Otherwise its remaining mass may be
assigned to D. This includes j=t; no special symmetry shortcut is needed.

## A dense bound without the old intersection identity

For the dense calculation, temporarily take independent Bernoulli(theta)
input coordinates, with 0<theta<1. Let z=exp(-lambda), lambda>0. For a
fixed nonzero q of image weight v, define

    g0 = 1-theta+theta*z,       g1 = theta+(1-theta)*z,
    F_v = g0^(t-v) g1^v,
    rho0 = |1-2*theta*z/g0|,
    rho1 = |1-2*(1-theta)*z/g1|.

F_v is exactly E[z^wt(X+Aq)]. Normalize the weighted input law by F_v.
Its coordinates remain independent. Their absolute Fourier biases are rho0
outside the support of Aq and rho1 inside it.

For a nonzero dual state a, put w=wt(B^T a) and let h be the intersection
size of the supports of B^T a and Aq. The Fourier coefficient has magnitude
rho0^(w-h) rho1^h. Without relating A to B, the available restriction is

    max(0,v+w-t) <= h <= min(v,w).

Fourier inversion and the triangle inequality therefore bound every syndrome
point probability of the tilted input law by

    H_v = 2^(-s) [1 + sum_w b_w
            max_h rho0^(w-h) rho1^h].

The maximum is at an endpoint because its logarithm is affine in h.
`screen_dense.bernoulli` uses the upper endpoint when rho1>=rho0 and the
lower endpoint otherwise. For an outward implementation, evaluating both
endpoints and taking an outward maximum also avoids an uncertain comparison.

The old symmetric calculation restricted h through wt(A(q+a)). That
restriction is unavailable for independent maps. The new bound does not use
it and does not single out a=q.

Fix an arbitrary target state y. On the lazy branch, the weighted probability
of q+BX=y is at most F_v H_v. On the fresh branch, conditional on X,
the probability of Mq+BX=y is at most 1/m. Thus

    E[z^wt(X+Aq) 1_{Mq+BX=y}] <= F_v (H_v/2 + 1/(2m)).

This holds for every target, including zero. Assign the right-hand side to
Z and a_w times that bound to each output C_w. For an incoming D component,
maximize over v. For incoming C_v, the same bound applies to every source
in the shell and hence to its uniform average.

At q=0, use the exact weighted B-kernel enumerator for the Z-to-Z entry;
assign the remaining emitted mass to D. These choices give a nonnegative
transfer that bounds the weighted measure after every epoch.

A second valid dense transfer is the binomial mixture of the fixed-j
transfers. The screen evaluates both complete moment bounds and takes their
minimum. It does not take an unjustified entrywise minimum of representations.
The existing outer coefficient bounds and occupancy-cover geometry can then
be reused, with witnesses retuned for the new transfer.

## Evidence and certificate scope

Exact integer audits cover all seven B maps. Rational toy-state tests compare
both transfer families to every fixed-weight input law and to selected
Bernoulli laws, including zero, arbitrary nonzero atoms, and uniform shells.
The tests also establish the adjoint identity on all bilinear basis pairs
of a three-epoch example with BA nonzero. These tests supplement the argument;
they do not prove the actual instance's numerical margin.

Both shortlisted maps have numerically passing covers of every occupancy
Q=1..32768. No unexamined occupancy interval remains in that screen. The
selected `greedy3_2` map has subsequently passed the following checks:

1. The outward producer reconstructs B's exact spectra, fiber caps, and A/B
   cancellation histograms, keeping A and B separate throughout.
2. It uses the independent-map Q=1 transfer and the dense H_v bound above,
   without the symmetric Q=1 or intersection shortcuts.
3. It evaluates every retained witness with 256-bit outward arithmetic and
   packs upward dyadic bounds. A 512-bit replay checks every bound; a separate
   verifier checks the exact union and complete occupancy partition.
4. The verifier binds the certificate maps to the tested generated encoder
   and checks that its outer generator spans the certified outer code.

The old balanced certificate does not certify a changed B. The new certificate
is specific to `greedy3_2`, this outer, and K=2^20. Higher-precision replay uses
the same producer, not an independently implemented proof. Rational toy tests
exercise the transfer formulas but do not replace review of their derivation.
Other finite parameters and the asymptotic theorem remain separate migration
gates; see `MIGRATION_PLAN.md`. The supported default is unchanged.
