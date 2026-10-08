# Randomizing the feedback direction

The fixed feedback map can create a concentrated state distribution when
the recursion first leaves zero. The old proof keeps this distribution as
arbitrary residual mass. This candidate instead randomizes the feedback
direction. It changes the construction, not just its analysis.

Use the outer, route, and fixed maps from [the sum-rank candidate](HD_SUMRANK_INNER.md).
For state dimension s, identify F^s with L=GF(q^s), where F=GF65536.
Fix this F-linear identification before setup. For each packet t, sample
rho_t independently from L\{0}. Let R_t be multiplication by rho_t.
Independently sample the old rank-one M_t. The recursion is

    Y_t=F X_t+A Q_t,
    Q_(t+1)=M_t Q_t+R_t C X_t,       Q_0=0.

All sampled maps are fixed in setup and reused for every message.
For each fixed nonzero h in F^s, R_t h is uniform nonzero in F^s.
For h=0, R_t h=0. This statement is pointwise in h; it does not assert
independence between different messages under the same R_t.

The s=4 version requires a GF(2^64) multiplication per scalar lane and
packet. This is not a native timing claim. The s=8 and s=12 screens are
separate, larger-state candidates, not changes to that cost.

## A two-component density bound

Condition on a packet's potential positions U. Under the outer comparison,
X_U is uniform and independent of the preceding weighted state measure.
Put D=q^s-1, epsilon=1/q, and c=1-epsilon. Let U_* be uniform nonzero state.
Represent the weighted state measure by mu<=Z delta_0+B U_*.
For output mark 0<z<1, define

    K=E_X[1[CX=0] z^wt(FX)],
    T=E_X[z^wt(FX)],
    V=E_(h uniform nonzero,X)[z^wt(FX+Ah)],
    V0=E_(h uniform nonzero,X)[1[CX=0] z^wt(FX+Ah)].

The zero-state row is exactly (K,T-K). If feedback is nonzero, R_t makes
the successor uniform nonzero independently of the current output weight.
There is no arbitrary residual component at birth.

For a flat incoming state, the feedback-active branch has weighted mass
V-V0. Adding the independent uniform nonzero feedback to M_t h puts at
most (V-V0)/D mass on any one target. Its zero-target mass is exactly
(V-V0)/D, because M_t h is nonzero.

On the feedback-zero branch the rank-one kernel is a mixture of uniform
refresh, of mass c, and maps preserving scalar lines, of mass epsilon.
The uniform part gives flat coefficient c V0. For each nonzero target
on a line-preserving branch, the previous state is fixed up to a scalar;
X remains constrained by CX=0. Let Gamma bound this fixed-state quiet
emission moment. The nonnegative Fourier coefficients of the Hamming tilt
imply Gamma<=K: an affine coset cannot exceed its zero coset.

Thus a positive upper operator is

    [ K       T-K             ]
    [ V/D     V+epsilon Gamma ].

At empty input its lower-left entry is zero and its lower-right entry is
c V+epsilon Gamma. Nonempty input uses the displayed conservative operator.
It does not subtract an upper bound for V0.

The sum-rank formulas imply a pointwise live-output floor d(U).
If r is the rank of C_U, then Gamma<=q^-r z^d(U). The implementation takes
the minimum of this bound, K, and D V. With an input mark v>=1 it uses
v^|U| times the unmarked quiet bounds. It does not apply the coset argument
to an input tilt whose Fourier coefficients have the wrong sign.

The [small-field regression](theory/hd_feedback_refresh_tests.py) checks
every input subset of a four-coordinate GF4 example. It checks the
zero-state row exactly and the flat-state bound pointwise at every target.

## Current result

The [s=16 refinement](HD_SUMRANK_FEEDBACK_CERTIFICATE.md) is fully certified:
K=2^20, rate 1/2, 40% symbol distance, and more than 59.87 bits of margin.
It retains the profile-conditioned intermediate-state sum across pairs
and extracts the actual regional occupancy coefficients. Its native cost
has not been measured. The smaller-state diagnostics below remain open.

The feedback refresh improves several middle and dense numerical bounds.
Packetwise fractional bounds still leave intermediate active-row counts
open. The initial s=8 screen at r=1024 was positive, about +6336 bits.
A floating projective-state calculation reduced a related witness to about
+3755 bits. Neither is a certificate or a counterexample.

Exact regional coefficients and short-block grouping supplied that closure.
An initial FLOAT regional screen underflowed and printed negative infinity;
those outputs are invalid and discarded. The helper now forms short-block
products in log space and rejects zero regional coefficients. Directed
arithmetic is required before accepting any negative result.

No native implementation, timing, or promotion has been made. All current
encoders and earlier proofs remain available.
