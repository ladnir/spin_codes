# Clipping the message count before averaging routes

The route-conditioned bound is valid for the ideal setup distribution.
Its good event concerns regional permutations only. It does not concern
outer symbol maps or inner update maps. The proof removes bad routes first;
it does not assert independence after conditioning on the good event.

For a fixed subset S of q outer groups, let J record every physical step's
number of potential slots. These slots exist even when their byte labels are zero.
Conditional on J, positions within each step are uniform. They are independent
between steps. Outer setup remains independent of these regional permutations.
The fixed-route outer measure bound therefore still applies.

The existing H2 argument discards a rare set of routes and bounds the expected
number of low-weight words on the remainder. A stronger alternative clips that
conditional expectation at one before averaging routes. One unfavorable route
cannot contribute more than one to the conditional failure probability.

Fix an output-weight variable z in (0,1]. Write B=beta^q z^(-D), where D is the
low-weight cutoff. Let Ttilde_j(z) be the existing positive local operator,
after averaging the zero/nonzero labels at j potential slots. Define

```
F(J) = e0 product_t Ttilde_{J_t}(z) 1.
```

This is an upper bound after averaging positions conditional on J. The
probability that some message supported on S has low-weight output is at most
min(1, B F(J)). The group-subset union consequently satisfies

```
Pr[failure at occupancy q] <= C(512,q) E_J[min(1, B F(J))].
```

For 0<alpha<=1, min(1,x)<=x^alpha. Expand the matrix product defining F(J)
as a sum of nonnegative path weights. Subadditivity of x^alpha gives

```
F(J)^alpha <= e0 product_t (Ttilde_{J_t}(z))^(entrywise alpha) 1.
```

The resulting bound uses the existing ordered placement recurrence, with
entrywise-powered local operators. Its prefactor is C(512,q) B^alpha.
The factor C(512,q) is outside the power: geometry is averaged separately for
each group subset. Potential-label thinning must also precede the entrywise power.
Powering the nonzero-packet operators first would concern a different random
variable and does not establish this bound.

This fractional argument can retain a routing condition H_c>=h. Add its
bad-route probability once, then multiply each powered local operator by
exp(nu min(j,c)), for nu>=0. Multiply the final expression by exp(-nu h).
This uses a deterministic indicator inequality. No conditional independence
claim about the universal good-route event is needed.

The main possible loss is the pathwise subadditivity bound. It separates
hidden inner-state paths before taking the fractional power. A bounded
numerical gate should test that loss before introducing a larger state space
or a more expensive encoder. Alpha=1 must reproduce the previous exact-q
potential-slot comparison. Every floating result remains a proposal.

## First bounded numerical gate

The baseline 16-bit state was tested at q=119, weight tilts 0.3, 0.4, and
0.5, and fractional powers 0.25, 0.4, 0.55, 0.7, 0.85, and 1. No routing
threshold was imposed in this gate. The strongest sampled margin was
-297.463566 bits at tilt 0.5 and power 0.4. At that same tilt, power one
reproduced the ordinary comparison margin of -4015.972720 bits. All three
power-one comparisons passed the regression check. The source-pinned
receipt is `fractional_v1.json`; it records floating proposals only.

This improves the previous H2-conditioned baseline proposal of roughly
-617.49 bits, but it does not close the q=119 bound. A bounded next test is
to average two or four consecutive steps conditional on their total
potential occupancy before taking the entrywise fractional power. This
removes some hidden path splitting at the cost of conditioning on less
routing information. The following test evaluates that tradeoff.

## Conditioning on short runs of steps

Fix a group size g dividing 64. Partition each region's 64 physical steps
into consecutive macros of g steps. For a fixed outer-group subset S, let
Jbar record the number of potential slots in each macro. Conditional on
Jbar, supports in distinct macros are independent. A macro with count j
has a uniform j-subset of its 8g positions. Define its comparison operator

```
Q_j = sum_{j1+...+jg=j} [product_i C(8,ji) / C(8g,j)]
                             Ttilde_j1 ... Ttilde_jg.
```

The product retains chronological order. Conditional averaging of the
positive physical comparison gives

```
E[F(J) | Jbar] = e0 product_macros Q_Jbar 1.
```

The function min(1,Bx) is concave for nonnegative x. Consequently, clipping
after this additional averaging remains an upper bound. Applying the same
positive-path argument then gives entrywise powers Q_j^alpha. A region's
macro counts have weights product_macros C(8g,j)/C(512,q), with total q.
The outer-group union remains C(512,q), outside the fractional power.
At alpha=1, the two conditional averages recover the physical-step formula.

The bounded gate tested g=2 and g=4, weight tilts 0.4, 0.5, and 0.6, and
powers 0.3, 0.4, 0.5, and 0.6. The best g=2 margin was -983.467496 bits;
the best g=4 margin was -1498.634375 bits. Both occurred at tilt 0.5.
They are weaker than the physical-step result. The alpha=1 regressions
passed, as did five tiny exhaustive tests. `grouped_v1.json` records the
source-pinned floating calculations. These results do not change the code.

A possible refinement retains each physical occupancy while grouping
hidden state paths. For two steps, its macro family would be

```
R_l = sum_{j+k=l} [C(8,j) C(8,k) / C(16,l)]
                         (Ttilde_j Ttilde_k)^(entrywise alpha).
```

This is not the coarse family Q_l^alpha. Subadditivity bounds R_l by the
corresponding product of powered physical operators. Concavity also bounds
R_l by Q_l^alpha. Thus this refinement avoids the coarser-conditioning loss.

The same formula extends to g physical steps. Fix the complete fine-count
sequence J before applying the fractional bound. For each consecutive
g-tuple, multiply its physical operators and take that product's entrywise
power. This sums hidden state paths inside the tuple before subadditivity.
Only then average the tuple over its fine counts, conditional on their sum.
The conditional weights are product_i C(8,ji)/C(8g,l); they are not powered.
The resulting regional recurrence sums all original fine-count sequences
with their original probabilities. It does not discard their occupancy
information before the nonlinear step.

An H2 condition can be combined with this refinement. Multiply each tuple
by exp(nu sum_i min(ji,2)) after its matrix power. The final prefactor contains
exp(-nu h), and the universal bad-route probability is added once. The
marker is not raised to alpha. Neither this argument nor its implementation
adds unconditional zero-state terms to the existing physical operators.

The fine-count gate used the same three tilts and powers 0.3, 0.4, and 0.5.
Its best g=2 result was -180.318024 bits; g=4 improved to -59.554191 bits.
Both used tilt 0.5 and power 0.4. At the best g=4 parameters, the tested
H2 markers nu=1,2,3 gave approximately -192.34, -891.73, and -2044.00 bits.
Thus none improved the unconditioned fractional bound. The exact route
threshold was h=2559 with witness x=9/35; its union bound exceeds 60 bits.

Five additional tests verify the explicit fine-count average, both ordering
inequalities, marked alpha=1 regression, the pointwise clipping inequality,
and preservation of zero entries. The source-pinned receipts are
`fine_g2_v1.json`, `fine_g4_v1.json`, and `fine_g4_h2_v1.json`. They concern
only q=119 and remain floating proposals, not a whole-code certificate.
Refining the g=4 grid to powers 0.325 through 0.425 and tilts 0.45 through
0.55 improved its selected-occupancy margin to -39.486504 bits, at tilt
0.475 and power 0.35. The receipt is `fine_g4_refined_v1.json`.

At those same parameters, grouping eight physical steps improved q=119 to
44.998077 bits. The encoder and setup distribution remained unchanged.
The calculation merged two four-step tuple lists in bounded-memory chunks;
it did not enumerate the complete eight-step list in memory. A simultaneous
alpha=1 sum checked the conditional weights and chronological products.
The configured 64 MiB budget estimates live working arrays, not peak memory.
The first implementation retains the previous product chunk while allocating
the next, so its transient usage can exceed that estimate. This allocation
issue does not affect the numerical formulas or results.
Four tests cover agreement with the previous g=4 implementation, an explicit
tiny g=8 census, occupancy normalization, and parameter validation.

The resulting powered operator does not depend on q. Evaluating every
q from 2 through 512 with that single operator gives positive margins for
q=113,...,159 and at least 50 bits for q=120,...,151. The receipt
`long_g8_allq_v1.json` stores the operator, source hashes, regression error,
and complete occupancy screen. The selected g=4 operator is also retained
in `long_g4_allq_v1.json`. These remain floating proposals; other occupancies
and outward rounding must be handled before asserting a whole-code bound.
The next bounded step is to tune the g=8 weight tilt and fractional power,
then select a small collection of operators that covers the needed range.

## A bounded occupancy curve

A subsequent g=4 sweep varied the weight tilt with q and tested powers
0.2, 0.4, 0.6, and 1. Each selected winner was reevaluated using a fresh
all-logarithmic placement recurrence. The largest difference from its
initial evaluation was less than 1.4e-12 bits. No tilt failed its numerical
checks. These checks compare floating evaluations; they do not supply
outward rounding.

| q | Best sampled margin (bits) | Weight tilt | Fractional power |
| --- | ---: | ---: | ---: |
| 2 | 49.01 | 0.005 | 1 |
| 4 | -12.42 | 0.01 | 0.6 |
| 8 | -72.21 | 0.02 | 0.4 |
| 16 | -112.71 | 0.06 | 0.4 |
| 32 | -178.95 | 0.12 | 0.4 |
| 64 | -215.04 | 0.25 | 0.4 |
| 96 | -163.44 | 0.4 | 0.4 |
| 160 | 245.72 | 0.7 | 0.4 |
| 256 | 1607.15 | 1.2 | 0.6 |
| 384 | 4117.28 | 1.6 | 0.6 |
| 512 | 4064.24 | 2.2 | 1 |

The dense tail is not the sampled obstruction. The lower-middle range,
especially q=32,...,96, still needs a stronger comparison. The next targeted
g=8 computation should therefore use parameters near the q=64 winner,
not repeat the successful q=119 point. The receipts `curve_16_64_v1.json`
and `curve_remaining_v1.json` retain the selected powered operators and
their source hashes. No encoder or setup distribution changed in this sweep.

The targeted g=8 computation at tilt 0.25 and power 0.4 improved q=64 from
-215.04 to approximately -173.45 bits. It produced no positive occupancy
interval. Its receipt, `long_g8_q64_allq_v1.json`, retains the macro family.
The alpha=1 relative regression error was below 7e-14. This result does not
justify extrapolating the q=119 success to smaller supports.

## Birth families after potential-label thinning

The original comparison has one birth family for each nonzero-byte count.
A fixed potential count j therefore enters a mixture of these families.
Fractional path bounds split that mixture across several coordinates.
Changing the comparison basis can remove this artificial split without
changing the encoder or the underlying positive comparison moment.

Let Told_j be the already-thinned comparison operator for potential count j.
Its coordinates are delta0, the uniform nonzero distribution U, and eight
normalized birth families P_k. Define c_jk=(Told_j)_(0,k+1) and
b_j=sum_k c_jk. For j=1,...,8, use the normalized mixture

```
Ptilde_j = sum_k (c_jk/b_j) P_k.
```

Let Q fix delta0 and U. Row j+1 of Q contains the coefficients c_jk/b_j
on the old birth coordinates. Thus Q1=1 and e0Q=e0. Define the new zero
row to be (Told_j)_(0,0)e0+b_j e_(j+1). Every nonzero-source row is the
corresponding row of Q Told_j, restricted to columns delta0 and U.
At j=0 there is no birth. The old comparison sends every nonzero source
only into delta0 and U, so these definitions give the exact identity

```
Tnew_j Q = Q Told_j.
```

Telescoping this identity proves equality of the terminal comparison moment
for every chronological potential-count sequence. This is an identity
between positive comparison closures. It does not assert that either
closure is the exact physical state transition.

The first three g=4 tests retained their previous weight tilts and powers:

| q | Before rebasing (bits) | After rebasing (bits) | Gain (bits) |
| --- | ---: | ---: | ---: |
| 16 | -112.707212 | -111.731661 | 0.975550 |
| 64 | -215.043032 | -203.932667 | 11.110365 |
| 119 | -39.486504 | 1.262165 | 40.748669 |

All-j intertwining, row stochasticity, global alpha=1 identity, and fresh
all-log evaluations passed. Four tiny tests also check every length-five
sequence over three potential counts. `potential_birth_g4_v1.json` retains
Q and the new macro operators. This removes a genuine representation loss,
but the tested small-occupancy gaps remain open.

The larger gain at q=119 has a direct explanation. A potential count j=1
can produce only the old family P_1, so its birth needs no rebasing.
Counts j>=2 can mix several old families and occur more often at larger q.
For normalized mixture weights p_k=c_jk/b_j, the old entrywise power assigns
total birth mass b_j^alpha sum_k p_k^alpha. The new coordinate assigns only
b_j^alpha. Stronger weight tilts also increase the relative weight of zero
labels and can broaden this mixture. This explains a source of the observed
gain; it is not a claim that every occupancy benefits by a fixed amount.

## The outer density is already tight on allowed symbol supports

The four parallel GF16 RS[16,8] rows form an MDS code over Q=65,536 symbols.
For a fixed support of w symbols, with 9<=w<=16, its number of words is

```
a_w = sum_{i=0}^{w-9} (-1)^i C(w,i) (Q^(w-8-i)-1).
```

Independent transitive symbol maps give density a_w/(Q-1)^w at each word
with that support. The current uniform majorant gives density (Q-1)^(-8).
Their ratio is a_w/(Q-1)^(w-8). It equals one at w=9 and approximately
0.99987793 for w=10,...,16. Across 119 groups, replacing these coefficients
would save at most approximately 0.021 bits.

The substantial distinction is the exact absence of supports w<=8.
An outer refinement could help if the tilted comparison places appreciable
mass on those impossible supports. Otherwise, refining the allowed-support
coefficients cannot explain a gap of hundreds of bits.
