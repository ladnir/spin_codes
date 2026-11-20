# A region-type table and a central-count bad-message family

The second-moment calculation needs actual shared-permutation probabilities
for many region types. The earlier certificate covered 125 types. We now
certify 73,600 types with one outward Fourier transform. Separately, we can
restrict the counted messages to a wider central weight window while retaining
a large first moment. Both results concern the unchanged fixed BCH/RM2Sub
construction. They do not establish its setup-failure probability.

## Extracting all Fourier aliases at once

Use the pair kernel polynomial K_2 from `ZERO_STATE_SECOND_MOMENT.md`.
Let m=(m_00,m_01,m_10,m_11) be the counts of the four coordinate pairs in
one region. The exact probability of two zero-syndrome region vectors under
their shared uniform permutation is

\[
\beta_2(m)=\frac{[x^m]K_2(x)^{64}}{\binom{8192}{m}}.
\]

Fix a positive tilt p whose entries sum to one. A discrete Fourier transform
on grid sizes M=(M_1,M_2,M_3), with x_00 fixed at p_00, returns

\[
A_r(p)=\sum_{m':\ (m'_{01},m'_{10},m'_{11})\equiv r\pmod M}
 [x^{m'}]K_2(x)^{64}p^{m'}.
\]

All coefficients are nonnegative. Therefore A_{m mod M}(p) bounds the
target coefficient times p^m from above. Unlike the earlier L1 bound,
the transform retains the complex phases separately for each residue class.
One table consequently bounds many types. Extra aliased coefficients remain
an upper-bound cost; they are not discarded.

The 128x128x64 diagnostic grid was insufficient for the new box: its largest
ratio bound was about 15.94. A 256x256x128 grid reduced the screened maximum
to about 1.02446. Both screens remain retained. Screen results are not used
as outward certificates.

## Explicit outward transform

`certify_pair_alias_table.py` encloses every stage of the larger calculation.
The positive tilt is the exact dyadic tilt from the retained central receipt.
`pair_fourier_error.py` supplies the previously audited error for evaluating
K_2 at each grid point. The code multiplies by 2^30 before six complex
squarings; the resulting polynomial is scaled by 2^1920.

Let z be a computed complex value with absolute error e, and let u:=2^-53.
An upward bound for the error of its computed square is

\[
2|z|e+e^2+8u|z|^2+2^{-1000}. \tag{1}
\]

The first two terms propagate the input error. The third covers normal
binary64 complex multiplication; the last covers subnormal errors.
Magnitude calculations and the nonnegative operations in (1) round upward.

`outward_radix2_fft.py` then performs an explicit radix-two transform.
Arb encloses every twiddle's real and imaginary components within 2^-52 of
their binary64 values. Thus each complex twiddle error is at most 4u.

For a butterfly with computed inputs a,b and errors e_a,e_b, let t be its
computed twiddle product with b. Both outputs a+t and a-t have error at most

\[
e_a+e_b+16u|b|+2u(|a|+|t|)+2^{-1000}. \tag{2}
\]

The exact twiddle has unit magnitude, so propagating b's input error costs
e_b. Its approximation error and the complex product's normal roundoff fit
inside 16u|b|. The following addition or subtraction fits inside the next
term. The implementation rounds the bound upward at every operation.
No statistical independence of arithmetic errors is assumed.

The final normalization divides by a power of two. Its radius rounds upward,
including subnormal cases. The certificate assumes ordinary IEEE binary64
round-to-nearest arithmetic and square root, as did the earlier evaluator.
It checks finite outputs and verifies that every imaginary residual is
inside its radius; exact coefficient aliases are real.

Arb evaluates the final multinomial, tilt, and single-region normalization
in the log domain. Per-overlap maxima are rounded up to dyadics with
denominator 2^30. The 256-bit producer and 512-bit root/final replay passed.
Direct 512-bit DFT tests cover one-, two-, and three-dimensional transforms,
both with zero input error and with perturbed inputs. Separate tests cover
the six squarings. These tests supplement the error bounds (1)--(2).

## Certified box

For every x,y in {780,782,...,858} and every k in {60,61,...,105}, define

\[
m=(8192-x-y+k,y-k,x-k,k).
\]

All 40*40*46=73,600 types satisfy

\[
\beta_2(m)<\frac{41}{40}\,\beta_x\beta_y.
\]

The largest stored upper ratio is approximately 1.0244557494297624.
The largest diagnostic relative transform radius among these coefficients
is about 3.29e-7. The compact receipt retains 46 per-overlap maxima and
their ranges; every integer type is checked by the producer and replay.
Large Fourier arrays are temporary and are not written to disk.

The box is not a full second-moment bound. In particular, multiplying
(41/40)^256 does not cover messages whose types leave the box, nor does
it evaluate the dependence between their single-region factors.

## Restricting the counted messages without changing the construction

For a setup theta, let Z_G(theta) count messages whose Q=2620 occupied rows
belong to the actual weight-80 shell T_80, whose every region count lies in
I={740,742,...,900}, and whose every epoch has zero B-syndrome.
The region-count restriction depends on theta through its row permutations.
It restricts only this auxiliary count; no encoder or setup randomness changes.

Every message counted by Z_G has output weight 209600 and is bad for cutoff
209716. Thus proving a sufficiently large probability that Z_G>0 would
refute the desired probability bound for the same construction. That final
probability implication has not yet been obtained.

We now lower-bound E Z_G without paying a large normalization loss from the
window restriction. After row permutation, each occupied row is an independent
uniform weight-80 subset. Set n:=256, p:=5/16, mu:=Qp, and v:=p(1-p).
For the vector J of region counts, define its centered squared norm

\[
S:=\sum_{r=1}^n(J_r-\mu)^2.
\]

Its exact mean and variance are

\[
M:=\mathbb E S=Qnv=144100,
\qquad
V:=\operatorname{Var}S=2Q(Q-1)\frac{n^2}{n-1}v^2. \tag{3}
\]

For completeness, write each centered row as W_i. Its squared norm is nv
deterministically. Then S=Qnv+2 sum_{i<j} W_i dot W_j. Distinct summands
are uncorrelated, including pairs sharing one index, because each row has
mean zero. The row covariance has diagonal v and off-diagonal -v/(n-1),
so E(W_i dot W_j)^2=n^2v^2/(n-1). This proves (3).
`test_count_energy.py` verifies both formulas by exhaustive toy slices.

Let E_1 be the all-even column parity event. The single-message version
of the marked-row lemma in `WEIGHTED_PARITY_AND_JENSEN.md` implies

\[
|\mathbb E[f1_{E_1}]-2^{-255}\mathbb E f|
 \le 2^{-255}C_f\delta_d,
\quad \delta_d:=(2^{256}-2)(3/8)^{Q-d}, \tag{4}
\]

for degree-d polynomials in normalized row averages. This bound is
conservative by a factor of two in its nontrivial-character count.
For S, the coefficient norm is at most C:=n(Q+mu)^2. For S^2 it is at
most C^2. Since Q is even, Pr[E_1]>=2^-255. Equations (3)--(4) give

\[
m_E:=\frac{M-C\delta_2}{1+\delta_0}
 \le\mathbb E[S\mid E_1],
\]
\[
\operatorname{Var}(S\mid E_1)
 \le M^2+V+C^2\delta_4-m_E^2=:v_E.
\]

The conditional one-column parity bound from the previous note also gives
Pr[G]>=2^-255 g, where G means all region counts lie in I and

\[
g:=1-512(1+\eta)\Pr[\operatorname{Bin}(Q,p)\notin I
                     \text{ and even}],\qquad \eta<2^{-3300}.
\]

In particular, Pr[G|E_1]>=gamma:=g/(1+delta_0)>0.8365269.
Here G is a subset of E_1 because I contains only even integers.
Conditional Cauchy--Schwarz, applied to S and the indicator of G, yields

\[
\mathbb E[S\mid G]\ge m_E-
 \sqrt{v_E\frac{1-\gamma}{\gamma}}=:s_G. \tag{5}
\]

Finally, verify a quadratic lower envelope at every j in I:

\[
\log\beta_j\ge a+b(j-\mu)+c(j-\mu)^2,\qquad c>0.
\]

The producer uses an LP only to propose b,c. It reconstructs a by taking
the minimum of outward lower residuals over every eligible integer, so
optimizer accuracy is not a proof premise. The total count is exactly
80Q=n mu, which cancels the linear term. Jensen and (5) therefore give

\[
\Pr[G\text{ and zero-state everywhere}]
 \ge2^{-255}g\exp(na+cs_G). \tag{6}
\]

`zero_state_central_witness.py` certifies (3)--(6), with c approximately
0.00015656925675252602. The 512-bit replay passed without rerunning the LP.
Multiplying (6) by binom(8192,Q)a_80^Q gives E Z_G>2^19624, with diagnostic
log_2 lower 19624.688243732788. This loses only 10.405423268675321 bits
against the unrestricted Jensen lower. All region weights in this witness's
second moment are now inside 740..900 by definition.

## Wider types and remaining work

We can also retain parity when examining an arbitrary function of one
region's pair type. Fix the two occupied supports. Condition on that region's
complete input bit patterns, whose two marginal counts x,y are even.
At a common occupied row, the conditioned bits e,f each have probability
p_e,p_f, where p_1=5/16 and p_0=11/16. Coordinate transitivity makes these
probabilities independent of the shared permutation, so its conditional
law remains uniform. Given that permutation, the two row words remain
independent under their respective one-bit conditionings.

For either bit e, the average squared nontrivial character of its conditioned
row law is at most

\[
b_e:=\frac{91}{255}+
 \left(1-\frac{91}{255}\right)\frac{p_{\rm close}}{p_e^2}.
\]

Indeed, the two words in that square share the same conditioned bit, so
their difference has zero in that coordinate. Its close-pair probability
is at most p_close/p_e^2 by Bayes' rule. For every far even distance82..160
on the remaining255 positions, all nontrivial normalized Krawtchouk values
have absolute value at most91/255. Cauchy--Schwarz therefore bounds the
mixed e,f character by b_*:=max(b_0,b_1)<0.373748<3/8. Single conditioned
characters are bounded by a_*:=97/255, and a_*^2<=b_*.

The remaining255-column joint Fourier inversion has four trivial terms,
all positive because x,y are even. Repeating the support-overlap argument
gives, for every compatible conditioned type m,

\[
\left|\frac{\Pr[E_1\text{ for both messages}\mid m]}{2^{-508}}-1\right|
 \le (2^{255}-2)a_*^Q+
 \frac{(2^{255}-2)^2}{4}b_*^Q<2^{-3200}.
\]

The estimate is uniform over bit patterns and thus survives conditioning
only on their type. `certify_pair_type_conditioned_parity.py` checks its
constants exactly. `test_conditioned_pair_characters.py` verifies the
conditional character-square identity and mixed Cauchy--Schwarz calculation
on a transitive, nonuniform slice subset. Averaging the displayed estimate
against a nonnegative function of m is valid. Multiplying such functions
over different regions still requires a further argument.

A retained three-tilt screen covers weights 740..900 and overlaps 40..126,
using pair-swap symmetry. It is not an outward certificate. Its extreme
entries have noticeably larger imaginary residuals, and some suffer from
coefficient aliases. Targeted retilting confirms this matters: at type
(6514,860,778,40), the old atlas ratio about 1.317 drops to a screened
Fourier-L1 ratio about 1.0115. At (6838,614,614,126), the targeted ratio
remains around 1.302, so large-overlap growth is not solely numerical slack.

Next, extend the checked table to the central witness's complete weight
window and bound the overlap dependence, with explicit treatment of its
remaining tails. The shared-message quadratic dependence and weighted parity
lemmas remain available. The central witness removes outside-window weights
from its own second moment, but does not remove overlap tails or cross-region
dependence. The full positive-distance goal remains unchanged and open.

```powershell
python -B workstreams/bch_rm2sub_bridge/certify_pair_alias_table.py --verify
python -B workstreams/bch_rm2sub_bridge/test_outward_radix2_fft.py
python -B workstreams/bch_rm2sub_bridge/zero_state_central_witness.py --verify
python -B workstreams/bch_rm2sub_bridge/test_count_energy.py
python -B workstreams/bch_rm2sub_bridge/certify_pair_type_conditioned_parity.py
python -B workstreams/bch_rm2sub_bridge/test_conditioned_pair_characters.py
```
