# How BCH, RM, and random constituents change the engineering surface

Fix relative distance at 10%. The inputs are message dimension K, outer
block size B, epoch size t, and state dimension s. The response is the
failure-bound margin in bits. This study compares continuous Q1 margins,
where Q1 counts messages with one nonzero outer row. Its purpose is to
explain parameter sensitivity before selecting a threshold.

The primary exact spectra are BCH lengths 8, 32, 64, 128 and RM lengths
8, 32, 128, 512. BCH and RM share the complete spectra at lengths 8 and 32.
Their divergence at length 128 is therefore a useful controlled comparison.
BCH-256 belongs to the other workstream and is excluded here.

Random constituents require a different evidence label. The main random
curve averages over one uniformly sampled rate-half binary subspace, reused
across all rows. It is not the spectrum or certificate of a particular
sampled code. A second calculation uses the existing simultaneous spectrum
caps with a 60-bit setup-event budget. Both meanings remain explicit in
the outputs.

All finite calculations use the four-state uniform-refresh transfer from
`BCH_GROWTH_ANALYSIS.md`, fixed nested inner maps, complete epochs, output
before update, no flush, and cutoff floor(BL/10), where L=K/(B/2).
The existing full-occupation ledger is unchanged. These binary64 Q1
diagnostics do not establish a complete distance certificate.

## Separate the row count from the constituent contribution

Let M1 be the Q1 margin and define C=M1+log2(L). The explicit row-count
factor accounts for the L possible nonzero rows. If the remaining transfer
changes little with K, then

    M1(K,t,s,B) approximately C(B,t,s) - log2(K/(B/2)).

This suggests three distinct questions. Does C stabilize as K grows?
How rapidly does C improve with s? Which spectral feature limits that
improvement for each family? Plotting only a passing state at one target
margin would obscure all three questions.

The study measures t64/s20 at every message exponent 12 through 26,
all available t/s pairs at K=2^20, and additional K/s interaction points.
The epoch sizes are 64, 128, and 256. Each uses states log2(t)+1 through 20.
The interaction points use t=64, s in {10,12,16,20}, and message exponents
16,18,20,22,24. Non-native configurations with fewer than one epoch are
excluded. Random block sizes are 8,32,64,128,256,512.
The completed study has 387 geometries and 1,292 family evaluations.

For identical geometry and inner maps, every family uses the same evaluated
coefficient vectors. Each spectrum then weights those coefficients. The
tilt search includes a grid in lambda L and a local refinement around
the shells contributing within 24 bits of each family's largest term.
This controls numerical grid noise without fitting away real changes.

The generated figures are:

- `rm_growth_curves.png`: RM margin and row-adjusted margin against K.
- `rm_state_tradeoff.png`: RM contribution against s, comparing all three t.
- `random_growth_curves.png` and `random_state_tradeoff.png`: the same
  engineering slices for the random ensemble through length 512.
- `constituent_engineering_comparison.png`: contribution against B and
  state-size tradeoffs across the three families.
- `engineering_surface_slices.png`: measured K/s interactions, with numeric
  margins and a separate color scale for each panel.

At K=2^20, t=64, s=20, the finite reference points are:

| Family | B | Q1 margin | Row-adjusted contribution C | Dominant weight |
| --- | ---: | ---: | ---: | ---: |
| BCH | 64 | 10.387 | 25.387 | 12 |
| BCH | 128 | 32.770 | 46.770 | 22 |
| RM | 128 | 10.619 | 24.619 | 16 |
| RM | 512 | 39.201 | 51.201 | 32 |
| Random ensemble | 128 | 21.196 | 35.196 | 13 |
| Random ensemble | 256 | 60.087 | 73.087 | 26 |
| Random ensemble | 512 | 137.351 | 149.351 | 52 |

Thus the ordering changes with the available block size. BCH-128 exceeds
the random ensemble at the same B, while RM-128 is below both. Increasing
RM from 128 to 512 buys about 28.6 margin bits; the same increase buys
about 116.2 bits for the random ensemble. These comparisons concern Q1
margin at matched t/s/K, not encoding speed or a fixed random realization.

RM's weight-32 shell supplies more than 99.9999999% of its length-512 Q1
bound at the reference point. By contrast, the random ensemble's largest
single shell contributes only about 7.8%. This supports a minimum-shell
explanation for the RM anchor and a spectrum-wide explanation for random.

For RM-512 at t64 and K=2^20, increasing s through 12,14,16,20 gives
37.695,38.746,39.085,39.201 bits. The last four state bits buy only 0.116
bits. The surface makes this diminishing return visible without selecting
a target margin first.

A practical knee can be defined relative to the largest tested state.
At K=2^20 and t=64, the smallest tested s within the indicated loss of
the s=20 margin is:

| Family | B | Within 1 bit | Within 0.25 bits |
| --- | ---: | ---: | ---: |
| BCH | 64 | 13 | 16 |
| BCH | 128 | 14 | 17 |
| RM | 128 | 12 | 14 |
| RM | 512 | 13 | 15 |
| Random ensemble | 128 | 11 | 13 |
| Random ensemble | 256 | 13 | 15 |
| Random ensemble | 512 | 15 | 17 |

This measures diminishing returns within the tested maps; it is not a
minimum state for a full-distance certificate. The random quarter-bit
knee increases by two state bits per doubling of B over these three sizes.
That is a local engineering observation, not an asserted asymptotic law.

## Why the families need different size models

The persistent-state model in `BCH_GROWTH_ANALYSIS.md` isolates the geometry
of first activation. It assumes that the state never cancels after becoming
nonzero and that subsequent output has rate one half. An active outer word
of weight w chooses w of the B regions uniformly. Its first active position
within the earliest such region is continuous and uniform.

Let p_w be the probability that this modeled output has relative weight
at most delta=0.1. Put m=floor(2 delta B) and f=2 delta B-m. Then

    p_w = [binomial(m,w) + f binomial(m,w-1)] / binomial(B,w).

The constituent contribution depends on the sum of N_w p_w, where N_w is
the number of outer words of weight w. The finite transfer uses a Chernoff
bound rather than this exact model probability, so the report also computes
the corresponding Chernoff model. The model is explanatory; it omits
cancellations and finite-epoch fluctuations.

For BCH, the complete low-weight tail controls this contribution. At length
128, minimum weight is 22. The relevant onset shells are 22,24,26. At length
64, weights 12 and 13 both matter. The measured high-state contribution
closely follows the Chernoff model at these exact-spectrum anchors.
Predicting a larger BCH therefore requires information about its low-weight
counts, not just its length or a line through four observed margins.

For rate-half RM(r,2r+1), the minimum weight and its count are known:

    B = 2^(2r+1),
    d = 2^(r+1) = sqrt(2B),
    N_d = 2^r [2r+1 choose r]_2.

The Gaussian binomial coefficient counts r-dimensional linear subspaces.
The factor 2^r counts the translates of a codimension-r subspace, whose
indicator is a minimum-weight RM word. Each step in the rate-half family
quadruples B and doubles d. The exact spectrum at B=128 starts at weight 16,
whereas the exact BCH spectrum at that length starts at weight 22.

The model therefore attributes a larger late-activation contribution to
RM's lower-weight words. Increasing s can suppress state-related losses;
it cannot remove those words or change their placement probabilities.
This is the reason to expect a different plateau at the same B.

The minimum-shell model also suggests a scale beyond the observed points.
When d/B is small, the leading placement cost is approximately
d log2(1/(2 delta)). Subtracting log2(N_d) gives the approximate intercept

    C_RM,min approximately sqrt(2B) log2(5) - log2(N_d)

at delta=0.1, with a finite-population correction. Here log2(N_d) grows
quadratically in log2(B). This is a minimum-shell model, not a theorem
about the complete RM spectrum or encoder. Other shells can only increase
the modeled union and reduce its margin. The exact anchors let us measure
how much those other shells contribute before using this approximation.

For the random ensemble, the entire mean spectrum is available. A fixed
nonzero binary vector belongs to a uniform D-dimensional subspace with
probability q=(2^D-1)/(2^B-1), so

    E[N_w] = binomial(B,w) q.

For Q1 the outer counts enter linearly. Consequently this mean spectrum
can be used directly while averaging over one reused outer constituent and
the encoder setup. No independence between the different words is required.
For higher occupations, powers and products of spectrum counts arise;
replacing those quantities by powers of their means would be invalid.

The binomial factors cancel in the onset sum, giving an exact identity
within this model:

    sum_w E[N_w] p_w = q [(1+f) 2^m - 1].

At rate one half and delta=0.1, its negative base-two logarithm is
approximately 0.3B, up to a bounded rounding correction. This gives a
linear-in-B contribution for the random ensemble's onset model, compared
with the square-root leading term in RM's minimum-shell model.
The finite transfer, state-size effects, and model Chernoff loss must still
be checked before turning either relationship into an engineering forecast.

These explanations concern competing effects: low-weight placement cost,
the number of low-weight words, and state dynamics. Relative constituent
distance alone does not determine the surface. Neither does one passing
failure-margin threshold.

## A common model for the state-size knee

The first-activation model explains the plateau, but deliberately removes
the cancellations that s controls. A two-state long-region model restores
those cancellations and gives a candidate surface C(B,s) directly.

Use states Z (zero) and V (live). Every active input starts a zero state.
For a live state, a fresh nonzero multiplier cancels the nonzero input
syndrome with probability p_s=1/(2^s-1). Between active inputs, the finite
encoder refreshes its live state. Since every inner coordinate is a nonzero
linear functional, its mean output rate under a uniform nonzero state is
mu_s=2^(s-1)/(2^s-1).

The model replaces the many epoch outputs inside each live interval by
their deterministic mean rate mu_s. Let a=lambda L and v=a mu_s. Averaging
the active input location uniformly within a region gives
h(v)=(1-exp(-v))/v. With rows and columns ordered Z,V, the region transfers
are

    R0 = [[1, 0],
          [0, exp(-v)]],

    R1 = [[0,       h(v)],
          [p_s h(v), (1-p_s) exp(-v)]].

R0 represents a region without an active input; R1 represents a region
with one active input. A zero-to-live route emits only after activation,
which gives h(v). A cancelling live route emits only before activation,
which gives p_s h(v). A surviving live route emits throughout the region.
Positive polynomial composition of these matrices averages the w active
regions over their binomial(B,w) placements, just as in the finite transfer.

The model contains B, s, and the outer spectrum, but neither K nor t.
K appears afterward through the explicit row count. The ratio t/L governs
how many epochs the finite encoder has per region, so finite-epoch effects
can prevent the approximation from working at short lengths. This gives
a testable explanation for near-parallel K curves and nearly coincident
matched-state t curves. Setting p_s=0 and mu_s=1/2 recovers the earlier
persistent-state moment exactly.

`check_continuous_state_model.py` compares this approximation at s=10,12,16,20
with the finite surfaces at K=2^16,2^20,2^24. It uses complete exact spectra
for BCH and RM, and mean spectra for the random ensemble. The comparison
measures model error; it does not assert an outward finite-size error bound.

Across 44 family/state cases at each message length, the maximum absolute
error in the adjusted contribution is:

| K | Maximum model error |
| --- | ---: |
| 2^16 | 0.652 bits |
| 2^20 | 0.105 bits |
| 2^24 | 0.020 bits |

These comparisons use t=64 and s=10,12,16,20, with BCH and RM at lengths
8,32,128, RM at 512, and random means at 8,32,128,512. The model contains
no fitted coefficients. Its agreement supports the decomposition into a
constituent/state contribution and the row-count cost in this regime.
The displayed errors are measured discrepancies, not confidence intervals
or guaranteed errors for untested parameters.

The finite interaction slices provide a separate check. Across all three
primary families, C changes by at most 0.090 bits over K=2^20 through 2^24
at the sampled states 10,12,16,20. At K=2^20 and common states 9 through 20,
the largest spread among t=64,128,256 is 0.147 bits. For exact RM alone,
the spread is at most 0.048 bits. These measurements delimit the range
where treating t as a small Q1 correction is useful.

## How to use the surface

In a region where C is stable in K, doubling K consumes approximately one
bit. Below a state-size knee, increasing s can compensate. Near a plateau,
the useful knob becomes the constituent spectrum or block size instead.
The K/s interaction plot checks the range where this separation is useful;
short regions can introduce a visible finite-epoch correction.

Comparing t at matched s isolates its effect on this Q1 bound. It does not
measure implementation cost. Larger t reduces the frequency of state
updates, while changing inner-map arithmetic and possibly higher-occupation
bounds. No performance benchmarks are part of this study.

The random cap calculation should not be confused with the ensemble curve.
Its conditional margin uses deterministic simultaneous bounds on the
realized spectrum. Adding the event-failure term gives the separate field
`setup_charged_margin_bits`. That field cannot exceed 60 bits, because the
chosen setup budget is 2^-60. This ceiling is a consequence of the proof
budget, not an inherent property of random codes. Loose spectrum caps can
also cause a large gap below the ensemble curve.

At B=512, K=2^20, t64/s20, these three quantities are 137.351 bits for
the ensemble mean, 77.576 bits under the conditional spectrum caps, and
59.999993 bits after charging the setup event. They answer different
questions. Mean shell counts can be below one because the average includes
rare outer-code draws; they are not literal counts of a typical realization.

The next refinement should target a few points near the state-size knees
and check their higher occupations. That would test whether the Q1 surface
also predicts the complete union. It is more informative than immediately
launching another full grid or fitting a formula across distinct mechanisms.

## Reproduction

Run these commands sequentially in `landscape_db`, after generating or
restoring the local pilot inputs described in the README:

```text
powershell -ExecutionPolicy Bypass -File build_activation_refresh.ps1
python -m unittest -v test_activation_refresh_native.py test_engineering_surfaces.py
python study_engineering_surfaces.py
python report_engineering_surfaces.py
python verify_engineering_surfaces.py
python check_continuous_state_model.py
```

The native implementation preserves the four-state log recurrence and
explicit fixed-width contractions. It is checked against rational
arithmetic and the Python implementation through length 512. It does not
modify the frozen three-state producers. The study checkpoints each
geometry and rejects cached results if its source fingerprint changes.
Generated checkpoints, tables, summaries, and figures stay local and ignored.

The tests also enumerate every two-dimensional subspace of F_2^4 to check
the ensemble mean spectrum, verify total mean spectrum mass, and check the
closed-form random onset identity against the literal shell sum.
An independent 90-digit replay uses positive arithmetic and linear epoch
iteration at selected RM and random witnesses. Its scope is numerical
validation of those coefficient vectors, not outward certification.

All 85 local tests pass. The 130 shared BCH/RM anchor pairs agree exactly.
All 212 earlier BCH points reproduce within 0.000138 bits under the revised
witness grid, and no dominant witness reaches a grid edge. The independent
90-digit replay checks 1,155 coefficients at three RM/random witnesses,
with maximum absolute log error 8.868e-12. It also authenticates all 77
recorded study dependencies. These checks support the numerical comparisons;
they do not promote Q1 into a full-distance claim.
