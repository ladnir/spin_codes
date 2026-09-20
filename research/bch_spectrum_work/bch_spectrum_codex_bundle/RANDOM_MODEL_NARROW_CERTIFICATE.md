# Low-occupation bounds and BCH shell tests for the random inner

Updated: 2026-09-04. Both full higher-shell experiments and the final audit passed.

Subsequent completion: RANDOM_MODEL_FULL_DISTANCE_BOUND.md closes the entire
q>=4 tail below 2^-83. The all-occupation random-inner bound is now below 2^-40,
conditional on the same three tested shell caps. This note retains the detailed
low-occupation calculations and historical receipts.

The target is the fixed BCH-derived outer C=[256,128,38], repeated in 8192
rows, followed by the RandomStepConv-M22 inner. There are 2^20 message bits
and 2^21 output bits. The distance screen uses cutoff 209716 and target 2^-40.
This note does not transfer the calculation to RM2Sub.

For each q, let F_q be the expected number of nonzero messages with exactly q
active outer rows whose output weight is at most the cutoff. The expectation
is over independent row-coordinate permutations, region permutations, and the
RandomStepConv setup. Linearity of expectation does not require independence
between the outputs of different messages. To bound the probability of a bad
setup by 2^-40, the first-moment argument makes it sufficient to prove

\[
\sum_{q=1}^{8192}F_q\le2^{-40}.
\]

The setup is sampled once and shared by all messages. The calculations below
bound each message's probability before summing. They do not estimate a variance.

Assuming the three shell caps below, outward arithmetic gives the following
bounds. The displayed bit values approximate minus the base-two logarithm of
each certified upper bound; the threshold comparisons use exact fractions.

| Contribution | Certified upper bound, in bits (approximate) | Shell evidence used |
|---|---:|---|
| F_1 | 40.0476888687 | Three statistical caps plus deterministic caps |
| F_2 | 73.4876764115 | Deterministic caps only |
| F_3 | 114.5938571911 | Deterministic caps only |
| F_1+F_2+F_3 | 40.0476888685 | Same assumptions as F_1 |

The combined upper bound is strictly below 2^-40. Its exact residual exceeds
2^-45 and is approximately 3.2515045185% of the target budget. Consequently,
proving sum_(q=4)^8192 F_q<=2^-45 would suffice to complete the first-moment
bound under these shell assumptions. The subsequent full-distance note proves
that tail bound, while retaining the statistical shell hypotheses.

## Shell evidence and its probability space

The predeclared shell targets are

\[
A_{38}\le3{,}827{,}351{,}840{,}403,\quad
A_{40}\le50{,}000{,}000{,}000{,}000,\quad
A_{42}\le5{,}000{,}000{,}000{,}000{,}000.
\]

All three tests completed with zero hits and matching full independent replay.
Every retained audit record passed independent Python reconstruction.

| Weight | Fixed trials | Endpoint hits | Singular queries | Python audit records |
|---:|---:|---:|---:|---:|
| 38 | 150,361,035 | 0 | Not applicable to this sampler | 407 |
| 40 | 436,628,033 | 0 | 1,707,138 | 723 |
| 42 | 139,431,904 | 0 | 544,641 | 427 |

Singular queries at weights 40 and 42 counted as non-hits and were never
resampled. The higher-shell tests ran sequentially with fresh retained tapes,
not preflight prefixes. Neither test stopped early or retried its outcome.
Both received conditional statistical acceptance. Their receipt is
`generated/higher_endpoint_certificate_20260904/certificate.json`.

The endpoint-incidence identities concern a fixed code. Their statistical
soundness concerns erroneous acceptance of a false shell-count claim under
independent sampling. It is not a posterior probability that a fixed claim is
true. The implemented byte source is Windows CNG; the information-theoretic
false-accept bounds assume ideal IID bytes. Computational replacement requires
the corresponding distinguishing-advantage qualification.

The three error allocations are 2^-40, 2^-41, and 2^-41. Their familywise union
bound is 2^-39, separate from the desired 2^-40 bound on bad SPIN setups.
Neither this statistical evidence nor the arithmetic below is a deterministic
BCH spectrum theorem.

For the **joint cap claim**, acceptance requires all three tests to pass, as
they did in this execution. Fix a code violating at least one of the three caps,
and fix an index i of a violated cap. Joint acceptance is a subset of acceptance
by test i, so its probability is at most alpha_i<=2^-40. This argument does not
require independence between the tests. It concerns the experiment before
conditioning on the observed outcomes, not a posterior probability.
The 2^-39 bound remains valid for familywise error when reporting separate
claims; it is not the tight bound for this all-tests-must-pass decision.

## Independent occupation-one arithmetic

Let s>0 be a fixed Chernoff witness and z=exp(-s). Write b=(1+z)/2 and
epsilon=2^-22. The two-state moment matrices are

\[
Z=\begin{pmatrix}1&0\\\epsilon b&(1-\epsilon)b\end{pmatrix},\qquad
T=\begin{pmatrix}\epsilon b&(1-\epsilon)b\\\epsilon b&(1-\epsilon)b\end{pmatrix}.
\]

The states are inactive and active; the process starts inactive. Z handles an
input zero and T handles an input one. Matrix entries include the output-bit
generating-function factor.

For a length-L region, define U_j(L) as the sum of ordered matrix products
over all supports of size j. Then

\[
U_j(n+1)=U_j(n)Z+U_{j-1}(n)T,\qquad U_0(0)=I,
\]

with all other initial entries zero. The uniform-support region matrix is
R_j=U_j(L)/binom(L,j), where L=8192. The independent Arb computation uses
these unnormalized sums and exact binomial denominators.

For occupation one, an outer support of weight w selects w of the 256 regions
to use R_1; the other regions use R_0. A second positive recurrence sums these
products, starting with the row vector (1,0), before dividing by binom(256,w).
The terminal row sum is a moment bound M_w(s). Define

\[
c_w(s)=8192\min\{1,\exp(209716s)M_w(s)\}.
\]

Thus F_1<=sum_w A_w c_w(s_w), with a separate fixed witness s_w allowed for
each shell. The full Arb computation encloses every possible nonzero shell:
even weights 38 through 218, and the all-one shell 256.

All 91 interior coefficients passed comparison against the earlier directed
coefficient bounds. The separate rational witness s=1/128 bounds the all-one
contribution below 2^-642. No output region is deleted in this computation.
The receipt is `generated/bch256_q1_full_arb_transfer.json`.

The retained LP witnesses were rechecked without overwriting their certificates.
Exact rational dual inequalities reproduce all seven low-shell caps; the
Johnson certificates for weights 52 and 54 also pass fresh primal and dual
checks. The receipt is `generated/bch256_deterministic_caps_reaudit.json`.
This checks the arithmetic and cap mapping against the retained BCH-sandwich
model, rather than rederiving every algebraic constraint in that model.

The final evidence-aware sum substituted the three accepted caps and retained
the existing deterministic LP, Johnson, and packing caps for other shells.
Its receipt is `generated/bch256_q1_full_arb_evidence_audit.json`.
That frozen Q1-only receipt retains its earlier Q2/Q3 status text; the joint
receipt and the completed certificates below supersede those status fields.

## Occupation two is outward-certified

For two distinct active row positions, fixed codewords of weights a and b give
independent uniform supports after their independent row-coordinate permutations.
This remains true when their local messages coincide. Their counting factor is
A_a A_b; no pair enumerator of C is needed.

At each outer coordinate, the two support indicators have four possible values.
Appending that coordinate multiplies the current row by R_0, R_1, R_1, or R_2.
The positive two-index recurrence sums all ordered support pairs. At the end,
divide coefficient (a,b) by binom(256,a)binom(256,b), apply the Chernoff bound,
and sum with factor binom(8192,2)A_a A_b.

Arb encloses the region matrices. The support recurrence rounds every positive
binary64 operation upward and repairs positive underflow. The final sum is exact
rational arithmetic. Exact small examples check all 25 weight pairs and the
underflow branch.

Using only existing deterministic shell caps, the certified Q2 upper bound has
73.4876764 bits, as a decimal display approximation. Its exact rational value is
below 2^-60. The receipt is `generated/bch256_q2_positive_outward.json`.

## Higher occupations

The first cumulative-reference diagnostic does not establish an all-occupation
bound. Its failure is a failure of that upper bound, not evidence of a bad code.
The receipt is `generated/cumulative_reference_q3_diagnostic.json`.

A direct Q3 calculation is now outward-certified. It retains exact low weight classes
38,40,...,68 and replaces every weight at least 70 by weight 70, with total
class count at most 2^128-1. Deleting input ones increases the RandomStepConv
moment in expectation through an order-preserving coupling of its two-state
chain. This statement is not pointwise monotonicity for one fixed inner code.
Uniform row supports can be coupled by a shared random ordering; subsequent
region permutations preserve input-set inclusion.

To verify the monotonicity claim, draw independent coins B_i~Bernoulli(1/2)
and K_i~Bernoulli(2^-22) at each inner position. For input x_i and current
activity a_i, put u_i=a_i OR x_i, output y_i=u_i B_i, and next activity
a_(i+1)=u_i(1-K_i). This chain generates Z and T above. Coupling two chains
with the same coins preserves activity and output order when input ones are
added. Thus deleting ones increases the moment for 0<z<=1. The coupling is
applied separately to each message's marginal law; it does not assert
pointwise monotonicity of one fixed sampled linear encoder on all messages.

The calculation uses a three-index positive support-count recurrence and divides
by the product of the three exact binomial counts. Truncating indices at 70
does not affect any retained coefficient. Exact rational small examples check
all 64 weight triples and verify that truncation preserves the retained values.

The resulting Q3 upper bound has 114.5938572 bits, as a decimal display
approximation. The exact rational value is below 2^-60. It uses no statistical
shell caps. The receipt is `generated/bch256_q3_capped_outward.json`.

These low-occupation calculations alone do not cover all nonzero messages.
The subsequent full-distance note supplies the required q>=4 bound, giving
an all-message theorem conditional on the tested shell inequalities.

## Low-occupation verification and subsequent tail completion

`generated/bch256_q123_joint_evidence.json` records the exact aggregate,
residual, accepted caps, joint-test soundness, and source hashes. A separate
read-only checker reaggregates all three contributions and checks 23 JSON
dependencies, including retained numerical arrays and LP witnesses:

    python -B code/verify_bch_q123_evidence.py

This checker passed. The final Q1 audit also reran both endpoint-certificate
verifiers, including tape hashes and Python audit records. The long replays
had already completed; they were not rerun a second time.

The remaining proof target at this milestone was

\[
\sum_{q=4}^{8192}F_q\le2^{-45}.
\]

The full-distance note now establishes this inequality using cumulative moment
envelopes and certified transfer bounds, not extrapolation. A deterministic BCH
spectrum theorem and transfer to RM2Sub remain separate tasks. The present
calculation uses RandomStepConv-M22.
