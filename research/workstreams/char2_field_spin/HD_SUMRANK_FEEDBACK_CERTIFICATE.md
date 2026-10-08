# Larger-state feedback candidate: complete distance certificate

This is a proof-side refinement, not the timed four-coordinate encoder.
It keeps the randomized RS outer, one regional route, and 64-symbol
packets. It uses sixteen GF65536 state coordinates and an independent
random multiplier on each feedback vector. Its whole-code certificate is
closed at K=2^20, rate 1/2, and 40% GF65536-symbol distance, with more than
59.87 bits of setup-failure margin. The directed replay covers every
positive active-row count, not a sample of counts or messages.

The larger state and feedback multipliers have no native timing yet.
In particular, this note makes no sub-5-ms claim. The earlier four-state
proofs and implementations remain separate and unchanged.

The subsequent [encoding work](HD_ENCODING_OPTIMIZATION.md) preserves this
reference and implements a different feedback distribution. Its faster
diagonal-Cauchy/tower-basis refinement transfers this certificate with an
explicit loss, leaving more than 48.32 bits of margin. The further
[FFT state-map refinement](HD_FFT_STATE_MAPS.md) preserves its local
moments; the current implementation takes about 7.4 ms.

## Fixed maps and sampled setup

Put F=GF65536, q=65536, R=8192, K=128R, N=256R, and s=16. Hamming weight
counts nonzero F symbols. The outer applies RS[256,128,129] independently
to R message rows. Its coefficients lie in GF256, embedded in F. Setup
multiplies each outer output coordinate by an independent uniform element
of F excluding zero. Each of the 256 coordinate regions independently
permutes its R symbols uniformly. The route is sampled once, not resampled
for each message; a region contains 128 successive 64-symbol packets.

Use the local F0, A0, C0 maps and the F/GF16 basis from the
[sum-rank construction](HD_SUMRANK_INNER.md). Extend its generator blocks
to sixteen rows:

    R_i[j,l] = a_i^j beta_l^(16^j),   0<=j<16, 0<=l<4, 0<=i<8.

Their concatenation is an LRS[32,16] generator. For arbitrary GF16-defined
local subspaces W_i contained in GF16^4, its restriction has F-rank
min(16, sum_i dim(W_i)). The distinct norm classes and the sum-rank theorem
are the same as in the linked construction. Define

    F_packet = diag(F0,...,F0),
    C = [R_0 C0 | ... | R_7 C0],
    A = vertical_stack(A0 R_0^T,...,A0 R_7^T).

Fix an F-linear identification of F^16 with L=GF(q^16). For every packet t,
independently sample rho_t uniformly from L excluding zero, and let T_t
be multiplication by rho_t. Also sample u_t uniformly nonzero in F^16,
then v_t uniformly subject to v_t^T u_t != 1, and put M_t=I+u_t v_t^T.
All these coins are independent of the outer scalars and route. Encode by

    Y_t = F_packet X_t + A Q_t,
    Q_(t+1) = M_t Q_t + T_t C X_t,       Q_0=0.

There is no constraint on the final successor state. Every sampled map is
fixed and reused for all messages. The setup-failure event is the existence
of any nonzero message whose output has weight at most D=838860.
The multiplier T_t is a GF(2^256) operation on the internal state; it does
not change the symbol alphabet. Its implementation cost remains open.

The encoder is injective. Each output packet determines its input packet
once the preceding state is known, because F_packet is invertible; the
outer is also injective. Thus its rate is exactly K/N=1/2. Applying the
same sampled scalar code to eight F lanes in each 128-bit record preserves
the scalar symbol-distance guarantee. No union over eight independent
setups is needed in that shared-setup representation.

**Theorem.** For the setup distribution above,

    Pr_setup[there exists m in F^K excluding zero:
                 wt(Encode_setup(m)) <= 838860] < 2^(-59.87).

Consequently, except with that setup-failure probability, every nonzero
message has at least 838861 nonzero output symbols. Its relative distance
is at least 40%. The message may be chosen after seeing the entire setup.

## A positive two-component state bound

Fix the route and a set of r potential outer rows. The
[outer counting-measure comparison](HD_RANDOMIZED_OUTER_PROOF.md) gives
domination by beta^r times independent uniform F values in the 256r
potential positions, where beta=q^256/(q-1)^128. It includes artificial
zero values; a potential position is not asserted to be nonzero.

Condition on one packet's potential set U. For an output mark 0<z<1,
define K=E_X[1[CX=0] z^wt(F_packet X)], T=E_X[z^wt(F_packet X)], and
V=E_(h uniform nonzero,X)[z^wt(F_packet X+Ah)]. These expectations use
uniform X on F^U. Write D_s=q^s-1, c=1-1/q, and epsilon=1/q.

Bound the output-weighted state measure pointwise by Z delta_0+B U_*,
where U_* is uniform on the D_s nonzero states. The density coefficient B
does not assert that the actual weighted state distribution is uniform.
The [feedback-refresh argument](HD_FEEDBACK_REFRESH.md) gives

    (Z,B)_next <= (Z,B) P_U,
    P_U = [ K       T-K                ]
          [ V/D_s   V+epsilon Gamma_U  ].

At U empty the lower-left entry is zero and the lower-right entry is
c V+epsilon Gamma_U. Here Gamma_U bounds
E_X[1[CX=0] z^wt(F_packet X+Ah)] for every fixed nonzero h. This uses
the actual rank-one law M=a I+b Line+c Uniform on nonzero states,
where a=(q^(s-1)-1)/D_s and b=1/q-a. Both residual terms preserve a
pointwise density cap on each scalar line. Independently, T_t sends each
fixed nonzero feedback vector uniformly to the nonzero state space.

The verifier takes the minimum of three valid upper bounds for Gamma_U:
K; q^(-rank C_U) z^d(U); and D_s V. The first follows by expanding the
Hamming mark into positive output-zero constraints: an affine coset either
has no solutions to a selected constraint or has the same number as the
zero coset. The second also uses the minimum live-output weight d(U).
The third keeps the chosen nonzero state as one term of the sum defining V.

## Exact local moments

Let a_i count the potential inputs in tile i. For a selected output-zero
mask let b_i count its positions in that tile. Set

    r0 = sum_i min(a_i,b_i),
    rho = min(s, sum_i min(4,max(a_i-b_i,0))),
    eta = min(s, sum_i min(4,max(b_i-a_i,0))).

Sum-rank separation gives ranks r0+rho for the quiet constraints
[F_(zeros,U); C_U], and r0+eta for [F_(zeros,U), A_zeros]. Therefore the
mask's probabilities in K, T, and V are respectively

    q^(-r0-rho),   q^(-r0),
    q^(-r0) (q^(s-eta)-1)/(q^s-1).

Multiply by product_i binom(8,b_i) z^(8-b_i)(1-z)^b_i and sum over all
b_i. These positive sums are exact rational numbers at rational z. They
depend only on the eight occupancies, not their positions within a tile.
For d(U), maximize sum_i b_i subject to eta<s, and subtract this maximum
from 64. If eta=s, the selected zeros would force the state to be zero.

## Retaining state sums and the actual route

Pair consecutive packets within each region. For profiles u and v, first
form P_u P_v; only then raise its entries to a power 0<alpha<=1. This
retains the sum over the intermediate state. Subadditivity of x^alpha
justifies splitting histories at pair boundaries.

An occupancy profile a has multiplicity

    m(a) = 8! / product_j multiplicity_a(j)! * product_i binom(8,a_i).

Define a positive matrix polynomial H whose degree-j coefficient is
the sum of m(u)m(v) (P_u P_v)^[alpha] over pairs with |u|+|v|=j. The
power on a matrix in brackets is entrywise, not a matrix power. Then

    B_r = [x^r] H(x)^64 / binom(R,r)

bounds one entire region. The coefficient enforces exactly r potential
positions in that region. The 256 regions are independent, so a bound for
all low-weight messages with r active outer rows is

    binom(R,r) * (beta^r / ((q-1) z^D))^alpha
        * (1,0) B_r^256 (1,1)^T.

The factor binom(R,r) is outside alpha. The divisor q-1 counts projective
message classes: every realized encoder is F-linear. Conditional Markov
and min(1,x)<=x^alpha are applied before averaging the route. No uniform
global permutation or independent sampling with replacement is substituted.

For large r, a cheaper packetwise coefficient-marker bound suffices. Its
logarithm at a fixed witness is affine in r minus 255 log2 binom(R,r).
This is convex on integer row counts, so directed endpoint comparisons
certify the entire returned interval.

## Arithmetic and reproducibility

The [paired verifier](theory/hd_sumrank_paired_certificate.py) computes
local moments as exact rationals. It rounds operator entries upward to
dyadic rationals before grouping identical profiles. Above a declared
occupancy cutoff it uses the entrywise maximum over every profile at that
count, with their full binomial multiplicity. No tail is omitted.

The pair products and fractional powers use Arb intervals. Matrix
polynomial powers retain every coefficient needed for the requested rows;
truncating higher degrees cannot affect them. Positive upper endpoints
are used when a coefficient ball crosses zero. A positive marker rescales
the polynomial without changing the coefficient ratio. Acceptance and
the final union use directed comparisons, never floating signs.

[Pinned rational witnesses](theory/hd_sumrank_certificate_witnesses.py)
are the replay inputs. Run `python hd_sumrank_certificate_witnesses.py`
from the theory directory. The replay reports uncovered row counts and
sets `whole_code_certificate` only when that list is empty and the
combined logarithmic failure bound is below -40. A failed numerical
witness is not a low-weight-codeword counterexample.

The clean 256-bit replay returned 53 receipts, no uncovered row counts,
and the interval

    log2(failure) <=
    [-59.87186558535774360100320290125852863375538893712309295746130786021714
       +/- 7.90e-69].

Forty paired-profile witnesses cover rows 1 through 962. Twelve
packetwise witnesses cover rows 962 through 8192; an independently
proved inverse-output tail is included conservatively in the union.
Every integer in a paired interval is evaluated. Packetwise intervals
use the convexity argument above. Higher-precision 384-bit checks of
selected weak paired and packetwise witnesses agree. All 375 high-distance
regression tests pass, including the actual sixteen-state map ranks.

The sparse/quiet certificates for the earlier designs remain preserved.
No library implementation or benchmark has been changed by this proof.
Next, reduce the state or feedback-refresh cost against this certified
reference before optimizing a native implementation. The original
four-coordinate encoder's nonquiet classes are still open.
