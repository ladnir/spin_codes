# Removing the local Cauchy feedthrough

The current encoder spends work transforming each input tile before adding
state-dependent output. This candidate instead uses a systematic packet:

    Y = X + A Q,
    Q_next = M Q + T A^T X.

Here X,Y belong to F^64, Q belongs to F^16, and F=GF65536. Initialize Q=0
and leave the final state unrestricted. Retain the randomized RS[256,128]
outer, its independent output scalars, and the single regional permutation.
Each packet independently samples the same rank-one mixer M and feedback
randomizer T as the preceding candidate. These are setup coins, fixed for
all messages. Setup remains independent of the outer and route.

This is a different construction. The old Cauchy-feedthrough certificate
does not transfer by a basis change: its local moments must be replaced.
The original certified encoder remains available as the control.

## Fixed maps

Use B=GF16, its ordered eight-point additive subspace H, and the first four
normalized additive-basis polynomials n_0,...,n_3 from the
[tower-basis construction](HD_TOWER_BASIS_INNER.md). Let A0 evaluate these
polynomials on H. Its image is an MDS [8,4,5] code over B. Moreover,

    A0^T A0 = 0.

To see this, the derivative of the subspace polynomial of H is a nonzero
constant. Lagrange interpolation therefore gives sum_(h in H) p(h)=0 for
every polynomial p of degree below seven. Apply this to n_i n_j, whose
degree is at most six.

Let R_i be the sixteen-by-four LRS blocks from the
[FFT state-map refinement](HD_FFT_STATE_MAPS.md). Define

    A_i = A0 R_i^T,
    A = vertical_stack_i A_i,     C = A^T.

Then C A=0. The matrices remain F-linear. The encoder is invertible after
the outer: given the previous state and Y, recover X=Y+A Q, then compute
the next state. Thus the dense-tail bound based only on injectivity still
applies.

## Exact local moments

Fix a potential input support U and an output-zero mask Z. In tile i write
a_i=|U_i|, c_i=|U_i intersect Z_i|, and d_i=|Z_i minus U_i|. The identity
feedthrough fixes the c_i intersecting input variables directly. The
remaining constraints have ranks

    r0 = sum_i c_i,
    rho = min(16, sum_i min(4,a_i-c_i)),
    eta = min(16, sum_i min(4,d_i)).

Indeed, after eliminating the intersecting variables, the quiet constraint
uses C0 columns on U_i minus Z_i. Their B-rank is min(4,a_i-c_i). The live
output constraint uses A0 rows on Z_i minus U_i, of B-rank min(4,d_i).
The LRS subspace-rank theorem adds these residual dimensions, up to sixteen.

The probabilities for the quiet, zero-state, and uniform-nonzero-state
moments K,T,V are consequently the same functions of r0,rho,eta as before:

    q^(-r0-rho),  q^(-r0),
    q^(-r0) (q^(16-eta)-1)/(q^16-1),       q=65536.

For each tile sum over 0<=c_i<=a_i and 0<=d_i<=8-a_i with positive weight

    binom(a_i,c_i) binom(8-a_i,d_i)
      z^(8-c_i-d_i) (1-z)^(c_i+d_i).

The resulting moments depend only on the potential occupancies a_i, not
their positions. In particular T=(z+(1-z)/q)^|U|. Maximizing the number of
zeros subject to eta<16 gives the same live-output floor as previously:
for a specified zero count b_i, its least constraint cost is
min(4,max(0,b_i-a_i)).

The positive zero/live-state operator therefore has the same form, now
using these new moments. Its quiet live-state cap is still the minimum of
K, q^(-rank C_U) z^d(U), and (q^16-1)V. The unchanged route permits the
same paired-region coefficient extraction and projective message count.
The diagonal-Cauchy feedback density loss remains 11.541649 bits at K20.

## Current evidence

The exact rank and local-moment tests pass. Native packet basis tests and
full encodings agree with an independently constructed dense-map oracle.
The Python/native F,A,C digest is 1274355414671835662. The first full
native comparison is 7.184 ms versus 7.405/7.371-ms certified controls.
This timing alone does not establish the new distance guarantee.

The adaptive run covered 7519 active-row counts but stopped at the first
uncovered count 2290 when its paired-screen proposal lacked margin. This
is an incomplete certificate, not a counterexample. The fixed-coefficient
GFNI-affine refinement measures about 7.1 ms with the same output. It is
superseded by the simpler [global RS inner](HD_SYSTEMATIC_MDS.md), which
has a complete certificate and measures about 6.2 ms. Preserve this
intermediate construction as a separate proof/implementation experiment.

Sources: [exact local tests](theory/hd_identity_feedthrough_tests.py),
[FLOAT screening](theory/hd_identity_feedthrough_screen.py), and
[directed coverage](theory/hd_identity_feedthrough_certificate.py).
