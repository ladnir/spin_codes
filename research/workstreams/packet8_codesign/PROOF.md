# Byte-native proof gate

The baseline passes the sparse occupancy-one screen but fails the tested
middle-occupancy bound by thousands of bits. Retaining full birth-state
shapes does not close this gap. Scaling the expansion to minimum weight
17 changes the middle bound by only a few bits.

These are **floating proposals**, not probability endpoints or distance
certificates. A negative margin means that this first-moment upper bound
exceeds one. It does not exhibit a low-distance codeword.

## Literal construction and geometry

Use the AES polynomial basis of GF256, with modulus `0x11b`. State integer
`a + 256*b` represents the ordered field pair `(a,b)`. For `h=0,...,7`,
the baseline forward maps are

\[
 (A(a,b))_h=a+h b,\qquad
 Cx=\left(\sum_h x_h,\sum_h h x_h\right).
\]

Each physical step emits `y=x+A(a,b)` and updates the state to `M(a,b)+Cx`.
Every `M` is a fresh independent uniform element of `GL(2,GF256)`.
This group acts transitively on the 65,535 nonzero binary states. Thus,
for every fixed nonzero entering state, its refreshed image is uniform
nonzero. The proof does not require uniform sampling from `GL(16,2)`.

The expansion and feedback use literal field multiplication. In the AES
polynomial basis, `C` is not generally the binary transpose of `A`.
`maps.py` builds their binary columns separately. It verifies `CA=0`,
single-byte restriction rank eight, and two-byte restriction rank sixteen.
The complete expansion census has minimum weight eight, attained eight times.

The distinct candidate `byte_native_A_scaled` multiplies output byte `h`
of `A` by `3^h`. The eight scales are
`[1,3,5,15,17,51,85,255]`. It leaves `C` unchanged. Its minimum expansion
weight is 17, attained three times, but `CA` is nonzero. The original-order
encoder remains invertible: its output at each step has identity diagonal
in the current input. Neither that fact nor the local analysis below
requires `CA=0`. This is not the inverse-scaled-feedback candidate.

Both candidates use the small outer: four GF16 RS[16,8] rows encode each
128-bit group into 256 bits. Independent nonzero-transitive 16-bit symbol
maps randomize the sixteen aligned symbols. Each symbol supplies two
eight-bit packets. Independent group shuffles assign the 32 packets to
32 regions; independent regional shuffles order the 512 groups.

Consequently `K=65536`, `N=131072`, and each region contains 64 physical
steps of eight byte slots. State starts at zero and persists through all
2,048 physical steps, including region boundaries. There is no final flush.
The low-weight cutoff is `floor(N/10)=13107`.

The separate `byte_native_t32` candidate restricts the literal evaluation
points to `0,1,2,3`. Its physical maps have 32 output bits and four byte
slots, not 64-bit maps with changed metadata. It keeps the 16-bit state
and uses 128 physical steps per region. Every byte pair still has rank
sixteen. There are 4,096 physical steps in the complete inner, and six
abstract coordinates in its finite-birth operator.

## Exact finite local formulas

Let `b` be the packet width, `W` the number of packet slots, and
`S=2^s` the state-space size. The candidate has `b=W=8` and `s=16`.
For occupancy `j`, let `X_j` be uniform on inputs with exactly `j`
nonzero packets. Their positions are uniform, and their nonzero labels
are independent. Fix `0<z<=1`.

Define the weighted feedback law and entering-state emission moment by

\[
 W_j(v)=\mathbb E[z^{\operatorname{wt}X_j}\mathbf1\{CX_j=v\}],
 \qquad
 M_j(a)=\mathbb E[z^{\operatorname{wt}(X_j+Aa)}].
\]

For binary character `t`, write `r_h(t)=wt(C_h^T t)`. Its weighted
character value is

\[
 \widehat W_j(t)=\binom Wj^{-1}[u^j]
 \prod_{h=0}^{W-1}
 \left[1+u\frac{(1+z)^{b-r_h(t)}(1-z)^{r_h(t)}-1}{2^b-1}\right].
\]

Binary Walsh inversion gives
`W_j(v)=S^-1 sum_t (-1)^(t dot v) widehat W_j(t)`.
The `-1` removes the zero input label. The character profiles come from
`C`, not from the expansion images.

For `r_h=wt((Aa)_h)`, the emission polynomial is

\[
 M_j(a)=\binom Wj^{-1}[u^j]
 \prod_{h=0}^{W-1}
 \left[z^{r_h}+u\frac{(1+z)^b-z^{r_h}}{2^b-1}\right].
\]

An active packet emits every byte except its unshifted expansion byte.
The implementation evaluates both formulas for all 65,536 states. It
also recomputes the single-packet weighted law by a positive census of
all `8*255` labels, and checks the Walsh result against that census.

## Ten-coordinate positive operator

Let `L=S-1`, and let `U` assign mass `1/L` to each nonzero state.
For `i=1,...,8`, restrict `W_i` to nonzero states, call its mass `m_i`,
and normalize that restriction to a probability measure `P_i`.
The coordinate measures are `delta_0, U, P_1,...,P_8`.
All have terminal mass one.

From zero, occupancy `j` produces exactly

\[
 W_j(0)\delta_0+m_jP_j
\]

when `j>0`; occupancy zero retains `delta_0`.
For a nonzero-source coordinate measure `P`, define
`T=sum_a P(a) M_j(a)`. Conditional on `a!=0` and `X_j=x`, the next
state is uniform over all states except `Cx`. The emitted weight does
not depend on the fresh `M`. The weighted output measure is therefore
dominated pointwise by

\[
 T U+\frac{T}{L}\delta_0.
\]

For `j=0`, return to zero is impossible, so omit the second term.
For `j>0`, the return term drops the necessary condition `Cx!=0`;
it is an upper bound, not an equality in general. The nonzero part
also fills excluded destinations. The resulting operator is a positive
measure bound, not an exact Markov transition.

Every entering measure dominated by a nonnegative combination of these
coordinates remains dominated after a step. Birth shapes need persist
for only one step: every nonzero source then enters the common uniform
density bound. No maximization over expansion-weight classes is used.

## Ordered placement and the outer comparison

Write `T_j` for the local operator. After `e` physical steps, condition
on `r` active packet slots. The regional recurrence is

\[
 R_{e,r}=\sum_k
 \frac{\binom Wk\binom{(e-1)W}{r-k}}{\binom{eW}r}
 R_{e-1,r-k}T_k,\qquad R_{0,0}=I.
\]

Given the occupancies of disjoint step intervals, their position subsets
and labels are independent. The factors multiply in chronological order.
The recurrence thus averages without-replacement placements and does not
reset state. The implementation uses 64 physical steps directly; no old
four-bit macro wrapper is used.

The expected nonzero-message measure of one outer group satisfies

\[
 \mu\le\beta\,\mathrm{Uniform}(\{0,1\}^{256}),\qquad
 \beta=\frac{2^{256}}{(2^{16}-1)^8}.
\]

The comparison includes artificial zero input mass. Under that comparison,
byte activity is `p=255/256`. For `q` active groups, put

\[
 Q_q=\sum_{r=0}^q\binom qr p^r(1-p)^{q-r}R_{64,r}.
\]

The nonnegative first-moment expression evaluated by the tail screen is

\[
 \binom{512}q\beta^q z^{-13107} e_0 Q_q^{32}\mathbf1.
\]

The matrix power retains state between regions. Replacing it by a scalar
power of a zero-start regional moment would not be the same bound.

For `q=1`, the screen instead uses exact expected outer byte-shell counts.
It computes the coefficient of `u^v` in
`e_0(R_0+u R_1)^32 1`, divides by `binom(32,v)`, and applies the output
tilt. A separate minimum over tilts is allowed for each shell before
summing its exact expected count. No old four-bit endpoint is reused.

## Numerical scope and results

Local polynomial sums and Walsh inversion use NumPy extended precision
where available, then convert to binary64. Small negative inversion
residuals are clipped and recorded. This is not outward rounding.
Regional placement uses separate scaling at every occupancy. Detected
underflow triggers a complete log-domain recomputation. The log path
retains positive contributions without an absolute cutoff, but it is
still floating arithmetic. Neither backend produces a certified endpoint.

The baseline and A-only scaling were screened at every `q=10,...,256`
over fourteen tilts from `.00256` through `1.6`. The old rank-corrected
map was also evaluated with the same birth-family operator. These coarse
development results precede the final source additions described below.

| Candidate | Fresh q1 margin, bits | Weakest tested middle margin, bits |
|---|---:|---:|
| Byte-native baseline | 55.443909 | -4150.934639 at q119 |
| A-only powers-of-three scaling | 55.444557 | -4145.154906 at q120 |
| Historical rank-corrected map | Not checked at sparse tilts in this run | -4149.300948 at q119 |
| Actual t32 byte-native maps | 55.465386 | -4098.960381 at q119 |

For the baseline, selected middle values are approximately -2979.54 bits
at q64, -4039.55 at q113, -4052.88 at q128, and -807.38 at q256.
At common q113, A-only scaling recovers about 6.28 bits. Its stronger
minimum expansion weight does not address the dominant gap in this bound.
Likewise, full two-byte feedback rank and exact birth shapes do not make
the byte-native bound substantially stronger than the historical control.

The t32 result uses the same fourteen tilts. At its selected q119 tilt
`.6`, the all-zero-state path contributes an expression with margin
`+3041.09` bits, while the full operator gives `-4098.96` bits. The zero
path is a lower contribution to the same comparison expression, not a
stronger failure upper bound or a lower bound on actual code failure.

A finer baseline grid covered every q80 through q160, with tilts `.35`
through `.70` in increments of `.025`. The weakest value improved to
`-4023.604887` bits at q126. At q119, the best value was `-4015.96980`
bits at tilt `.5`. The coarse grid hid about 135 bits, not the
four-thousand-bit gap. This completed fine-grid receipt authenticates
against the current sources.

### The omitted return condition is too small at the tested tilts

The transition to zero requires `CX!=0`; the local operator drops that
condition. This is genuine upper-bound slack. It cannot be called an
actual return event without checking the joint emitted-weight condition.

There is a quantitative check that avoids a large joint census. Fix an
entering state and an occupied packet support. After weighting by
`z^wt(X+Aa)`, the active packet labels remain independent. For a fixed
expansion byte `v`, the tilted label normalizer is

\[
 Z_v=(1+z)^8-z^{\operatorname{wt}v}\ \ge\ Z_{\min}:=(1+z)^8-1.
\]

Every label has numerator at most one, so its tilted probability is at
most `1/Z_min`. Fixing all but two active labels determines at most one
remaining pair with `CX=0`, because every two-byte restriction is invertible.
Thus, for occupancy at least two,

\[
 \mathbb E[z^{\operatorname{wt}(X+Aa)}\mathbf1\{CX=0\}]
 \le Z_{\min}^{-2} M_j(a).
\]

The zero-syndrome term is exactly zero at occupancies one and two.
The inequality holds for each entering state and support, and therefore
for every normalized source family. Correcting only the zero-return entry
can improve a `T`-step expression by at most
`-T log2(1-Z_min^-2)`, when `Z_min>1`.
For t64, this limit is about `0.8322` bits at tilt `.4` and `2.8667`
bits at `.6`. The missing condition does not explain the middle gap there.

A related bound controls the whole nonzero-source closure. Fixing all
but one injective active byte bounds every tilted syndrome atom by
`eta=1/Z_min`. For each coordinate probability measure, its exact weighted
output dominates `(1-eta)` times the represented upper row, pointwise
in the destination state. Zero-source birth rows and occupancy-zero refresh
rows are exact. Positivity therefore gives the same lower comparison for
each chronological product, with factor `(1-eta)^T`; placement averages
preserve it. The entire local closure can overestimate the exact inner
moment under the **same outer comparison measure** by at most
`-T log2(1-eta)`: about 50 bits at `.4` and 93.47 bits at `.6` for t64.

These limits are tilt-specific. They neither lower-bound actual code
failure nor remove the outer pointwise majorant. Their purpose is to
distinguish local-envelope slack from a weakness of the complete outer
comparison. Small-field exhaustive tests check the tilted-atom inequalities.

### Counterfactual return-denominator sensitivity

`return_sensitivity.py` holds the literal t64/s16 emission moments,
birth measures, uniform-density transitions, outer comparison, and routing
geometry fixed. It changes **only** the nonzero-to-zero entries, replacing
their denominator 65,535 by `2^s-1`. The altered operator is not asserted
to arise from any code. In particular, `s=24` does not describe a
24-bit-state construction.

Nine tilts from `.10` through `1.6` gave these best floating margins:

| Counterfactual denominator bits | q64 | q119 | q128 | q256 |
|---|---:|---:|---:|---:|
| 20 | -2889.93 | -3893.17 | -3924.46 | -668.76 |
| 24 | -2821.57 | -3817.93 | -3848.01 | -591.37 |
| 28 | -2757.15 | -3752.46 | -3782.32 | -525.60 |
| 32 | -2693.12 | -3688.35 | -3718.19 | -461.46 |

Extra cancellation bits alone do not close this expression on the tested
grid. This does not rule out actual larger-state codes: changing the
state space also changes expansion and birth laws, which this diagnostic
deliberately freezes. A real 24-bit proposal needs its actual local maps
and moments before drawing a conclusion.

These observations do not identify an actual low-distance event, nor do
they prove that every stronger analysis must fail. They do rule out
promoting these maps from q1 alone. Exact q2 and a complete outward replay
were not attempted. The next proof task should compare actual constructions,
not extrapolate from the counterfactual. Two distinct directions are a
real larger-state byte-native inner and a wider outer group with 64
byte-routing regions. The latter changes packet spreading and outer
randomization cost. Neither direction has a timing or proof prediction
from this work. A full outward replay of the current negative expression
is not useful.

## Reproduction and checks

From the repository root, use fresh output filenames:

```text
python -B -m unittest discover -s research/workstreams/packet8_codesign -p "test_*.py" -v
python -B research/workstreams/packet8_codesign/screen.py --map byte_native --tilts .00256 .00384 .00512 .00768 .01 .0256 .0512 .10 .20 .40 .60 .80 1.2 1.6 --output <fresh-baseline.json>
python -B research/workstreams/packet8_codesign/screen.py --map byte_native_A_scaled --tilts .00256 .00384 .00512 .00768 .01 .0256 .0512 .10 .20 .40 .60 .80 1.2 1.6 --output <fresh-scaled.json>
python -B research/workstreams/packet8_codesign/screen.py --map old_rankfix --tilts .00256 .00384 .00512 .10 .20 .40 .60 .80 --output <fresh-control.json>
python -B research/workstreams/packet8_codesign/screen.py --map byte_native_t32 --tilts .00256 .00384 .00512 .00768 .01 .0256 .0512 .10 .20 .40 .60 .80 1.2 1.6 --output <fresh-t32.json>
python -B research/workstreams/packet8_codesign/screen.py --map byte_native --tilts .35 .375 .40 .425 .45 .475 .50 .525 .55 .575 .60 .625 .65 .675 .70 --min-q 80 --max-q 160 --output <fresh-fine-grid.json>
python -B research/workstreams/packet8_codesign/return_sensitivity.py --output <fresh-counterfactual.json>
```

`test_screen.py` independently enumerates small input spaces to check
weighted syndromes, per-state emissions, and full pointwise transition
domination. Noncommuting fixtures verify ordered placement and the
regional support recurrence. Further tests check the literal candidate
maps, complete spectrum totals, outer counts, and the distinct A-only
scaling. The separate implementation audit is in `test_independent.py`.

Each receipt pins the numerical driver, map generator, imported placement
backend, and outer-count/envelope sources. A source change requires a
fresh receipt. The current complete, source-authenticated receipts are
`t32_middle_v1.json`, `byte_native_fine_q80_160.json`,
`byte_native_q1_current_v1.json`, and `return_sensitivity_v1.json`.

The earlier `byte_native_middle_v1/v2/v3.json`, `A_scaled_middle_v1.json`,
and `old_rankfix_middle_v1.json` are development receipts. Their numerical
results are recorded above, but their saved source hashes do not match
the current files after the added t32 geometry, exact structural zeros,
and zero-sector diagnostic. They are not current authentication receipts.
The reproduction commands create fresh outputs from the final sources.
All these ignored JSON files are proposal outputs, not certificates.

## Separating routing failures from message counting

The larger-state diagnostic identifies a concentration effect that the first
moment counts once for every message. At q=119 and tilt .5, its normalized
comparison expression has approximately 549 occupied steps out of 2,048.
Those steps carry about 3,568 active bytes. Approximately 1,499 steps have
both zero input and zero entering state. These are statistics of a tilted
upper-bound expression, not observed encoder failures.

A routing event is shared by every assignment of nonzero messages to the
same outer groups. Its probability can therefore be bounded before counting
those assignments. The following split changes the analysis, not the code.

### A property of the sampled route

Let L=512 be the number of outer groups, R=32 the number of regions, and
W=8 the number of byte slots per physical step. Fix a subset S of q groups.
In every region, mark the q slots assigned to S, whether or not their eventual
byte values are zero. Let H(S) count physical steps containing at least one
marked slot, summed over all regions.

For a fixed S, each region contains a uniform q-subset of its 512 slots.
The regional permutations are independent. If B is the number of occupied
steps in one region, then

\[
 \Pr[B=b]=\frac{\binom{64}{b}}{\binom{512}{q}}
 [x^q]\big((1+x)^8-1\big)^b.
\]

The coefficient counts choices of slots that hit each of the b selected
steps. Thus H(S) is a sum of 32 independent copies of B. Write
g_q(v)=E[v^B]. For an integer threshold h_q and 0<v<1,

\[
 \Pr_{\rm route}[\exists S\subseteq[L],\ |S|=q:\ H(S)<h_q]
 \le \binom Lq v^{-(h_q-1)}g_q(v)^R.
\]

This union counts group subsets, not messages. In particular, it has no
factor beta^q. The strict threshold H<h_q explains the exponent h_q-1.

### Counting messages only on well-spread routes

Choose a set of occupancies to which the routing condition will apply.
Let G be the event that H(S)>=h_q for every selected q and every q-subset S.
For any nonnegative nu and any such S,

\[
 \mathbf1_G\le\mathbf1_{\{H(S)\ge h_q\}}
 \le \exp\big(\nu(H(S)-h_q)\big).
\]

After removing the bad-route event once, use the outer measure bound on
the remaining expected number of low-weight words. For each fixed S,
the comparison measure assigns independent uniform bits to all 256q
potential input coordinates. This comparison remains valid for a fixed regional route;
the outer setup is independent of those permutations.

The local operator must now distinguish potential from nonzero slots.
Let T_k be the original local upper operator for k nonzero byte packets.
For j potential slots in a physical step, define

\[
 \widetilde T_j(z)=\sum_{k=0}^j\binom jk
  (255/256)^k(1/256)^{j-k}T_k(z),\qquad
 T^{(\nu)}_j(z)=e^{\nu\mathbf1_{j>0}}\widetilde T_j(z).
\]

Uniform placement makes the j potential positions uniform within the step.
Conditional on k nonzero bytes, their positions are uniform among its eight
slots. This justifies using T_k in the mixture.

Apply the ordered placement recurrence to T^(nu), with exactly q potential
slots among 512 positions. Denote its regional operator by R_(64,q)(z,nu).
There is no additional binomial mixture over q at the regional boundary.
The contribution for q groups on G is at most

\[
 \binom Lq\beta^q z^{-13107}e^{-\nu h_q}
 e_0 R_{64,q}(z,\nu)^R\mathbf1.
\]

State remains continuous across all regions. The full failure bound adds
the bad-route union and these restricted message contributions. Other
occupancies may retain their original bounds. No division by Pr[G] is
needed: this is a bound on an expectation with an indicator, not a claim
that conditioning preserves independent regional permutations.

The finite formulas above are suitable for a numerical proof gate.
Any floating evaluation remains a proposal until a complete outward replay
covers the routing union and every message-occupancy contribution.

### Counting the first two potential packets per step

Counting occupied steps alone does not distinguish a singly occupied step
from a step containing several potential packets. A second statistic retains
some of that information. Let j_i(S) be the number of potential packets from
S in physical step i. For an integer c in 1,...,8, define

\[
 H_c(S):=\sum_i\min(j_i(S),c).
\]

The previous statistic is H_1. For one region, the probability generating
function for its contribution to H_c is

\[
 g_{q,c}(v)=\frac{[x^q]
    \left(\sum_{j=0}^8\binom8j x^j v^{\min(j,c)}\right)^{64}}
    {\binom{512}{q}}.
\]

Use this function in the same bad-route Chernoff bound. In the message bound,
replace the local marker by exp(nu min(j,c)). All group-subset unions,
potential-slot thinning, and chronological state products remain unchanged.
The case c=2 constrains excess packet multiplicity that H_1 does not record.
This is a stronger description of routing geometry, not additional inner mixing.

Both routing conditions can be required simultaneously. Their failure
probabilities add by a union bound; their statistics need not be independent.
On H_1>=h_1 and H_2>=h_2, use the single indicator majorant

\[
 \exp\big(\nu_1(H_1-h_1)+\nu_2(H_2-h_2)\big),
 \qquad \nu_1,\nu_2\ge0.
\]

The joint local marker is exp(nu_1 1_{j>0}+nu_2 min(j,2)). Evaluate its
ordered matrix product on the same route. Multiplying two separately
averaged regional moments would not preserve that dependence.
