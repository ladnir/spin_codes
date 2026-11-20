# The BCH × RM2Sub t64_s20 distance target is closed

The fixed BCH [256,128] outer with the selected RM2Sub t64_s20 inner meets
the finite distance target at K=2^20 message bits and rate1/2. The
computer-assisted failure bound covers every nonzero message and satisfies

\[
U<2^{-50}<2^{-40}.
\]

This is a first-moment proof under the fresh-independent-multiplier setup
defined below. It does not establish seeded-randomness equivalence,
implementation performance, or the full SPIN protocol's security.

## Fixed encoder and setup

Use the fixed code C defined in the **Fixed outer and random setup**
section of `../../bch_spectrum_work/bch_spectrum_codex_bundle/RANDOM_INNER_M22_CLOSURE.md`.
Specifically, C is the five-dimensional p37-syndrome restriction of the
extended designed-distance-37 BCH code over GF(256), with modulus0x14D
and syndrome subspace{0,...,31}. Its dimension is128 and minimum distance
is at least38. The argument uses certified constraints on this code's
spectrum, not a random-code spectrum model.

Fix a linear isomorphism G from F2^128 to C. Split each K-bit message into
8192 rows, and encode each row with G. Independently permute the256
coordinates of every encoded row. Transpose into256 regions of8192 bits,
then independently permute every region. Serialize the regions in order.
The resulting length is N=2^21 bits.

Fix the map A:F2^20->F2^64 whose20 generators appear in
`generated/larger_state_inputs_v1/t64_s20_selection.json`, selected seed3390111757.
The file SHA256 is
`23219aa244e3fd6b426f5e5474e3824dff2570459853ed8d3a257a53f408d0f5`.
Set B=A^T. Exact audits give rank(A)=20, BA=0, minimum nonzero A weight16,
and minimum nonzero B-kernel weight8. All generators are evaluations of
binary polynomials of degree at most2.

Fix an identification of F2^20 with GF(2^20). Divide the serialized input
into32768 epochs of64 bits. For epoch input X_i, use the recurrence

\[
q_0:=0,\qquad Y_i:=X_i+Aq_i,\qquad
q_{i+1}:=\alpha_iq_i+BX_i.
\]

Each alpha_i is uniform on the nonzero field elements. All multipliers
and permutations are mutually independent. The state carries across region
boundaries, output precedes update, and there is no final flush.

Let theta denote the complete setup, sampled once and shared by every
message. Write E_theta:F2^K->F2^N for the resulting encoder. The inner
transform is invertible: recover X_i=Y_i+Aq_i before updating q_i. Thus
E_theta is an injective linear map for every setup and has rate1/2.

## Finite distance theorem

**Theorem (computer-assisted, fresh-multiplier RM2Sub t64_s20).**
For the fixed encoder and setup distribution above,

\[
\Pr_\theta\!\left[
  \exists x\in\mathbb F_2^{2^{20}}\setminus\{0\}:
  \operatorname{wt}(E_\theta(x))\le209716
\right]\le U<2^{-50}.
\]

Consequently, outside a setup event of probability less than2^-50, the
image code has minimum distance at least209717 and relative distance
greater than0.1. The exact rational U is the `covered_union_upper` field
of `generated/larger_coverage_full.json`.

The final audit receipt is `generated/larger_t64_s20_full_closure_audit.json`.
It records diagnostic margin50.43906854330953bits, reconstructs the1163-row
outer joint LP, replays46 exact shell caps and the OA29 caps, and passes
nine tests in six modules. The exact rational inequality, not this decimal
margin, is the claimed bound.

## Proof and coverage

For Q=1,...,8192, let Z_Q count bad messages with exactly Q nonzero BCH
rows. Every nonzero message belongs to exactly one such class. The
certificates supply deterministic bounds E_theta[Z_Q]<=U_Q. Hence

\[
\Pr_\theta\!\left[\sum_Q Z_Q>0\right]
\le\mathbb E_\theta\!\left[\sum_QZ_Q\right]
\le\sum_{Q=1}^{8192}U_Q=U<2^{-50}.
\]

No independence between different messages is assumed. The one-row bound
uses the activation-aware transfer and the existing exact BCH weighted
shell certificate. Higher occupancies use the exact kernel spectrum,
fixed-weight region coefficients, and certified BCH shell caps.
`LARGER_STATE_GAP_CLOSURE.md` records the adaptive counting inequality.

| Occupancies | Evidence selected by the final ledger |
|---|---|
| 1 | Larger-state Q1 certificate and independent coefficient replay |
| 2..511 |40 short gap witnesses and512-bit replay |
| 512..1024 | Two overlapping range certificates and512-bit replays |
| 1025..2047 | Seven gap witnesses and512-bit replay |
| 2048..4096 | Three overlapping range certificates and512-bit replays |
| 4097..8025 | Checked extraction from upper-half shards and512-bit replay |
| 8026..8192 | Tilt6 endpoint-range certificate and512-bit replay |

Every integer occupancy is evaluated. The ledger takes the smallest
eligible certified bound where ranges overlap and sums with rational
arithmetic. It checks that the covered set equals{1,...,8192}.

The producer uses256-bit Arb for region coefficients and final powers,
with directed binary64 arithmetic in the positive scaled recurrence.
Range replay uses512-bit Arb and repeats the directed recurrence without
running an optimizer. Exact toy tests check the recurrence separately.
The Q1 replay additionally uses a different coefficient implementation.

The initial upper-half search failed at Q8026 because its tilt grid was
too coarse. The failed bound is retained. A new witness at tilt6 covers
every Q8026..8192. No failed row or floating screen contributes to U.
The prepared constant-row box fallback was not needed for closure.

## Handoff

The theorem concerns this selected map and the ideal setup distribution.
It does not certify t128_s15, every RM2Sub map, a PRG-generated multiplier
stream, or another message length. The complete BCH spectrum remains
unnecessary for this distance result.

Next: independently review and package the proof, then benchmark this
fixed baseline. Only after retaining this proved baseline should state
size, epoch size, or implementation cost be optimized.
