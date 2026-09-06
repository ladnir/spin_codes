# A rigorous obstruction to the unconditional first-moment route

For the fixed BCH [256,128] and RM2Sub t128_s15 construction, the expected
number of bad messages at occupancy Q=2620 exceeds **2^19600**. This is a
lower bound for the actual construction, not a random-code surrogate.

Consequently, no valid unconditional first-moment upper bound can close
the full range at these parameters. The result does **not** establish a
large setup-failure probability. A large mean can arise from rare setups
with many bad messages. A probability conclusion needs further evidence,
such as a second-moment bound.

All previously verified partial ranges remain valid. The original goal
requires setup failure below 2^-40; this note neither proves nor refutes
that probability statement. It identifies why refining the current
unconditional first-moment calculation cannot finish the goal.

## The branch that leaves the input unchanged

Use the setup and state update from `README.md`, including q_0=0:

\[
Y_i=X_i+Aq_i,\qquad q_{i+1}=\alpha_iq_i+BX_i.
\]

If BX_i=0 in every epoch, induction gives q_i=0 and Y_i=X_i throughout.
This event is independent of the multiplier values. Every input of weight
at most H=209716 on this branch is a bad output.

Let T be the set of BCH codewords with weights from 38 through 80, and
write a=|T|. The retained exact BCH constraint model gives

\[
a\ge a_0:=443405397513809550622719915749.
\tag{1}
\]

All words in T have even weight. Fix Q=2620 row positions and consider
all messages whose selected rows lie in T, with zero rows elsewhere.
Their input weights are at most 80Q=209600<H.

For counting only, sample the Q row codewords independently from T, then
apply the construction's independent row permutations. This auxiliary
sampling averages over messages; it is not added setup randomness.
Let J_r be the number of ones in region r, for r=1,...,256.

## Even column counts occur often enough

Let V be one sampled and permuted row in F_2^256, and let mu be its law.
For u in F_2^256, define its real Fourier coefficient

\[
\widehat\mu(u):=\mathbb E[(-1)^{u\cdot V}].
\]

The vector of region parities is the sum of Q independent copies of V.
Fourier inversion gives

\[
\Pr[J_r\text{ even for every }r]
 =2^{-256}\sum_u\widehat\mu(u)^Q\ge2^{-255}. \tag{2}
\]

Indeed, Q is even, so every summand is nonnegative. The characters u=0
and u=1 both contribute one because every V has even weight.

For each fixed region, J_r has distribution Bin(Q,p), where
38/256<=p<=80/256. Define

\[
D(x\Vert p):=x\log(x/p)+(1-x)\log((1-x)/(1-p)),
\]

with natural logarithms. Chernoff bounds and a union bound give

\[
\Pr[\exists r:J_r\notin[64,1536]]\le\varepsilon,
\]
\[
\varepsilon:=256\left(
 e^{-Q D(63/Q\Vert38/256)}+
 e^{-Q D(1537/Q\Vert80/256)}\right).
\]

The probability that all region counts are even and in this interval is
therefore at least gamma:=2^-255-epsilon. The outward calculation gives
gamma>2^-256. This subtraction requires no independence between region
counts and their parities.

## Exact probability of remaining in the kernel within a region

Let k_j count weight-j words in ker(B), and define
K(u):=sum_(j=0)^128 k_j u^j. The exact audited kernel spectrum determines K.
A region contains 64 epochs. Given its weight j, its independent uniform
permutation lies in (ker B)^64 with probability

\[
\beta_j:=\frac{[u^j]K(u)^{64}}{\binom{8192}{j}}.
\]

The outward calculation verifies

\[
\min_{j\in\{64,66,\ldots,1536\}}\beta_j>2^{-959}. \tag{3}
\]

Conditional on the row outcomes, region permutations are independent.
Combining (2), the interval restriction, and (3), the average probability
that all epochs have zero B-syndrome is at least

\[
\gamma\left(\min_j\beta_j\right)^{256}
 >2^{-256}\,2^{-959\cdot256}=2^{-245760}. \tag{4}
\]

The exponent 245760 equals 15 times the 16384 epochs. Equation (4) is
proved from the fixed kernel and row permutations; it does not assume
independent random parity checks.

## Counting the bad messages

Let Z_0 count messages in this row family that remain on the zero-state
branch. Distinct choices of occupied positions and nonzero BCH row words
give distinct messages. The auxiliary averaging identity and (4) imply

\[
\mathbb E_{\rm setup}[Z_0]
 \ge \binom{8192}{2620}\,a^{2620}\,2^{-245760}
 \ge \binom{8192}{2620}\,a_0^{2620}\,2^{-245760}
 >2^{19600}. \tag{5}
\]

Every message counted by Z_0 has output weight at most 209600. Hence the
same lower bound applies to the bad-message first moment at Q2620.
The simplified bound (5) has diagnostic logarithm 19669.8982893 bits.
The less-rounded stored lower bound has logarithm 19670.8994179 bits.

## Verification and next obligation

`exact_tail_lower_bound.py` maximizes the negative of the tail count in
the retained BCH model. It checks every primal constraint, dual sign,
dual domination inequality, and primal-dual equality with rational arithmetic.
Negating the resulting upper bound and rounding upward proves (1).
The lower bound is stored separately from all shell upper caps.

`zero_state_lower_bound.py` computes the kernel polynomial, the two tails,
and the resulting lower factor using Arb. A separate 512-bit replay gives
a lower factor no smaller than the stored one. Exact toy tests check the
even-convolution Fourier identity and the zero-state kernel event.

All replays below passed:

```powershell
python -B workstreams/bch_rm2sub_bridge/exact_tail_lower_bound.py --verify
python -B workstreams/bch_rm2sub_bridge/verify_zero_state_obstruction.py
python -B workstreams/bch_rm2sub_bridge/test_zero_state_parity.py
```

The next substantive question is how concentrated Z_0 is across setups.
A second-moment upper bound would turn (5) into a failure-probability lower
bound through E[Z_0]^2/E[Z_0^2]. That second moment is not yet available.
Do not silently replace s=15 with a larger state: any parameter change
requires a separate user decision and a new full-range objective.
