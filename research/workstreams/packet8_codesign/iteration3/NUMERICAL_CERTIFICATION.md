# Outward evaluation of the fine-grouped bound

The expensive step is the G8 macro sum: 9^8=43,046,721 ordered occupancy
tuples, each with a 10-by-10 matrix product and an entrywise fractional power.
The current NumPy evaluation is a proposal. Its agreement checks do not make
`numpy.power`, BLAS, or the final sums outward-rounded.

A practical certificate can separate exact local arithmetic, a tightly
specified native macro evaluator, and a small outward global verifier.
This note proposes that design; it does not certify an existing receipt.

## Exact local inputs

Choose the weight variable z=u/v in (0,1) directly, preferably with v a power
of two. Choose alpha=a/2^b in (0,1]. These are proof parameters; replacing
nearby floating choices requires a fresh evaluation. An optional route marker
can use a rational r>=1 in place of exp(nu), avoiding that transcendental too.

For a character t of the 16-bit state, let r_h(t) be the Hamming weight of
the binary character restricted to packet h. Define

```
D = 255 v^8,
Q_r = (v+u)^(8-r) (v-u)^r - v^8,
N_j(t) = [x^j] product_(h=0..7) (1+x Q_(r_h(t))).
```

An integer Walsh transform gives the exact weighted syndrome law:

```
W_j(s) = sum_t (-1)^(t·s) N_j(t) / (2^16 C(8,j) D^j).
```

All cancellations occur in integer arithmetic. Verify nonnegativity and the
total-mass identity exactly; never clip negative floating coefficients.

For a state a, write w_h=wt((Aa)_h), I_h=u^w_h v^(8-w_h), and
B_h=(v+u)^8-I_h. Its emission moment is

```
M_j(a) = [x^j] product_h (I_h+x B_h) / (v^64 255^j C(8,j)).
```

The normalized birth means are rational sums of W_i(a)M_j(a), divided by
the exact positive nonzero birth mass. Aggregate by the existing expansion
profiles to avoid constructing a rational object for every state pair.
The uniform means and the potential-slot thinning with probability 255/256
are rational as well. This produces the nine exact potential matrices U_j.

Convert each U_j entry to a dyadic upper endpoint. Verify the endpoint by
integer cross-multiplication against the exact rational entry. Preserve exact
zeros. This makes conversion correctness independent of the conversion library.

## Streaming the macro sum

For tuple J=(j_1,...,j_8), define

```
P_J = U_(j_1) ... U_(j_8),
n_J = product_i C(8,j_i),
ell = sum_i j_i.
```

The macro matrix for ell is the sum of n_J times the entrywise alpha power
of P_J, divided by C(64,ell). An optional route marker r^(sum_i min(j_i,c))
multiplies each summand **after** the fractional power. Neither n_J nor the
route marker is raised to alpha.

Keep n_J as an exact integer. Its maximum is 70^8<2^53, so every such
multiplicity is exactly representable in binary64. Postpone division by
C(64,ell) until after summation. This avoids billions of inexact weight
divisions and leaves only 65-by-100 final endpoint divisions.

### Fractional powers without a general power function

For alpha=11/32, compute five successive square roots and multiply the
1/4,1/16,1/32 roots. This uses five square roots and two multiplications.
For alpha=3/8, three square roots and one multiplication suffice. A rounding
analysis requires correctly rounded square roots, not merely a function
whose name is `sqrt`.

A faster portable alternative is a verified tangent table. For a positive
binary64 input x, write x=2^e m with 1<=m<2. For rational bin center t,
concavity gives

```
m^alpha <= alpha t^(alpha-1) m + (1-alpha) t^alpha.
```

Store dyadic upper endpoints for both positive coefficients and for
2^(e alpha). For dyadic alpha, their inequalities can be checked once by
integer powers and cross-multiplication. Evaluate the selected tangent using
positive multiplication and addition. A bin need only choose a good tangent:
every verified tangent bounds the function globally. Finer bins reduce
slack, but no unproved interpolation error enters correctness. The tangent
table needs a fresh proposal check before choosing its resolution.

## Native rounding assumptions and error budget

The fast evaluator must specify its arithmetic circuit. One sufficient
contract is IEEE binary64, round-to-nearest, correctly rounded basic
operations, no reassociation by fast-math, and no overflow or underflow on
positive intermediates. Fused multiply-add may be allowed under an explicitly
counted circuit. A square-root variant also requires correctly rounded roots.
The tangent-table variant does not need a square-root assumption.

Underflow before the fractional power is especially dangerous: x can
underflow even when x^alpha is representable. Structural zero must not be
confused with underflow. Before the run, use the checked dyadic input matrices
to bound the positive range. For example, if m is their smallest positive
entry, every positive eight-step product contains a positive term at least
m^8. A normal-range proof must also cover intermediate products, weighted
summands, and accumulators. If that proof fails, use an exponent-tracking or
arbitrary-precision path; silently flushing values is not permitted.

Let eps=2^-53. Without exceptional values, every positive basic operation
has a relative error bounded by eps. A conservative 10-dimensional dot-product
count is 20 factors per matrix multiplication, hence 140 for an eight-step
product. Fractional powers attenuate that input error by alpha. The root
chain adds a small constant number of factors. Positive summation over at
most 9^8 terms adds at most 9^8-1 factors, irrespective of cancellation,
because there is none.

For an audited circuit with total error count K and K eps<1, inflate each
computed sum by the exact rational factor

```
1/(1-K eps) = 2^53/(2^53-K).
```

Bernoulli's inequality justifies this replacement for (1-eps)^(-K).
For the specified simple circuit, K=9^8+256 leaves ample allowance for the
small per-term overhead; the final implementation must check that count
against its actual operation graph. The resulting relative inflation is
about 4.8e-9 per macro, not an empirical error estimate. Even repeated across
256 macros, this correction costs only about 1.8e-6 bits. Tangent slack is a
separate, rigorously upper-bounding approximation and must be measured separately.

After inflation, divide each dyadic sum by the exact integer C(64,ell) using
integer arithmetic or an outward interval library. Save the resulting upper
macro matrices with all map, parameter, table, circuit, and source identities.

## Verifier and remaining work

Use outward interval arithmetic for the smaller regional recurrence, the
continuous-state product across 32 regions, and the final occupancy sum.
The prefactor C(512,q) stays outside alpha. The factor
(beta^q z^(-cutoff))^alpha is evaluated outward. Any route-failure terms are
added using their exact rational witnesses. Select proof parameters per q
only through an explicit receipt; cover every required q before claiming a
whole-code margin.

The portable reference verifier can evaluate the macro circuit with directed
MPFR arithmetic or a specified software binary64 implementation. It is slower
but does not trust a vendor BLAS or power function. The fast native verifier
is valid only on platforms satisfying the stated arithmetic contract. Source
hashes and regression tests authenticate a run; they do not replace that
contract or prove arithmetic correctness.

Recommended order: first certify one selected G4 point with exact local
inputs. Then implement the checked tangent or square-root G8 evaluator and
replay the promising occupancy range. No new code construction or encoder
benchmark is needed for this numerical work.
