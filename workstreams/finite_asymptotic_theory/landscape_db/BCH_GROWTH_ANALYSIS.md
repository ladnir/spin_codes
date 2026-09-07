# BCH growth: separate message length, state size, and outer spectrum

The follow-up `BCH_DOMINANCE_ANALYSIS.md` tests whether these Q1 curves
represent the complete occupation bound. Until that audit controls the
remaining sum, the growth slopes and state-size knees below describe
Q1 only and do not establish full-margin parameter recommendations.

The question is how large the BCH constituent must become as the message
length grows, and which state sizes give useful tradeoffs. A 40-bit threshold
does not explain that relationship. We first study the continuous bound for
one nonzero outer row, called Q1, using the four complete exact spectra
already authenticated in the local landscape: lengths 8, 32, 64, and 128.

The main finding is a separation of effects. Doubling the message length
costs nearly one bit of Q1 margin in the tested long-message range.
Increasing the state dimension initially improves the margin. For BCH-64
and BCH-128, that improvement approaches a plateau explained by the position
of the first active outer region. The smaller BCH codes remain sensitive to
state size. A single size fit through all four points mixes these mechanisms.

These are Q1 diagnostics and an explanatory model. They do not establish a
full-distance certificate or an asymptotic law for unknown BCH spectra.
The existing higher-occupation bounds have not been changed by this study.

## The measured quantity

Let the outer constituent have length B, dimension D=B/2, and nonzero weight
counts N_w. Let k be the message dimension and L=k/D the number of outer
rows. The output has BL bits, divided into B regions of L bits each.
An epoch contains t bits; its state has s bits. All runs use complete epochs,
the recorded nested RM2Sub maps, and cutoff floor(BL/10).

For outer weight w, let G_w(lambda) be the transfer upper bound on
E[exp(-lambda Y)], where Y is output weight. The expectation is over the
encoder's permutations and fresh nonzero field multipliers, for fixed maps.
The code computes

    U1 = L sum_w N_w min(1, min_lambda exp(lambda floor(BL/10)) G_w(lambda)),
    M1 = -log2(U1),
    C_B = M1 + log2(L).

Each shell may use its own positive tilt. The numerical minimum is over
evaluated witnesses; it need not be the global minimum to give an upper
bound in exact arithmetic. Binary64 evaluation is not outward certification.
A negative M1 means this bound is vacuous, not that the code fails.

The row factor L is explicit. Plotting C_B removes that counting cost and
exposes the contribution of the constituent and state dynamics.

The refined study has 212 tuples: t=64, s=20 at every message exponent
12 through 26, plus all available t/s choices at k=2^20. It uses t in
{64,128,256} and s from log2(t)+1 through 20. Overlapping tuples run once.
The witness search includes a grid in lambda L and refines relevant shell
optima at spacing 0.002 in log(lambda). No dominant witness reaches the
outer edge of the evaluated grid. This removes visible grid-phase noise;
it is not a proof of global optimization.

At k=2^20, t=64, s=20:

| Exact constituent | Minimum weight | Q1 margin M1 | Adjusted contribution C_B |
| --- | ---: | ---: | ---: |
| BCH [8,4] | 4 | 2.545 | 20.545 |
| BCH [32,16] | 8 | 5.023 | 21.023 |
| BCH [64,32] | 12 | 10.387 | 25.387 |
| BCH [128,64] | 22 | 32.770 | 46.770 |

The length-64 input is the recorded shortened extended BCH spectrum. It has
odd weights, including weight 13; parity filtering would change the answer.

Fitting M1 against log2(k) over exponents 16 through 24 gives slopes
-0.993, -0.980, -0.988, and -0.978, respectively. Across that range, C_B
varies by at most 0.219 bits for any one constituent. Shorter messages show
a finite-epoch transition, so the local near-unit slope should not be
extended to arbitrary lengths without checking that regime.

The generated `bch_growth_curves.png` shows both M1 and C_B against log2(k).
The generated `bch_state_tradeoff.png` shows C_B against s, with separate
curves for each t. These continuous curves replace threshold staircases.

## A tighter transfer removes avoidable state-density loss

The old Q1 envelope carries a near-uniform density penalty through every
zero-input epoch. Fresh multiplication actually restores exact uniformity
after such an epoch. The new evaluator retains that information in a fourth
state class. It also separates the zero and surviving masses after an
active epoch. Existing producer sources and their receipts remain untouched.

Here is the encoder interface needed for the argument. Fix an injective
linear inner map A from F_(2^s) to binary t-vectors. Fix the linear syndrome
map J from binary t-vectors to F_(2^s), with J(e_j) nonzero for every unit
vector e_j. The selected RM2Sub maps have this property. At each epoch with
state x and input u, output u+A(x), then set

    x' = alpha x + J(u),

where alpha is freshly sampled from the nonzero field elements. The state
starts at zero and no final flush is appended. For Q1, each epoch input is
zero or one unit vector. Within an active region, the active epoch and its
coordinate are uniform; different regions use independent placements.

Put M=2^s-1 and kappa=M/(M-1). For a fixed positive tilt, classify the
normalized distribution of the current state after weighting past outputs:

- Z: the state is zero.
- D: an arbitrary distribution on nonzero states.
- U: exactly uniform on nonzero states.
- L: each nonzero state has probability at most 1/(M-1).

Each class also carries its unnormalized weighted mass. The normalization
is only for expressing the distribution invariant.

Fix z=exp(-lambda). Let a_v count nonzero inner words of weight v and let
d_A be their minimum weight. Define

    m0 = sum_v a_v z^v / M,
    m1 = sum_v a_v [(v/t) z^(v-1) + (1-v/t) z^(v+1)] / M,
    r0 = z^d_A,
    r1 = z^(d_A-1).

For zero input, the nonzero weighted transitions are

    Z -> Z : 1
    D -> U : r0
    U -> U : m0
    L -> U : min(r0, kappa m0).

For one unit input, Z goes to D with weight z. From D, U, or L, take r equal
to r1, m1, or min(r1,kappa m1), respectively. The transitions are

    live -> Z : r/M
    live -> L : r(M-1)/M.

To justify the refresh, fix any nonzero prestate x and current input u.
The current output weight is determined before alpha is sampled.
Multiplication by alpha sends x uniformly over the nonzero field elements,
independently of that output weight. If u=0, the destination is therefore
exactly uniform, even after weighting the past and current outputs.
If u=e_j, exactly one multiplier cancels J(e_j). Conditional on survival,
the destination is uniform on the nonzero states except J(e_j).
Its density is at most 1/(M-1). Averaging over prestates and coordinates
preserves this density bound, including after output weighting.

The factors m0 and m1 are the exact output moments for class U. A density
bound gives kappa times those moments for class L; minimum distance gives
the alternative pointwise bounds. These facts establish the displayed
transitions as weighted-mass upper bounds. Positive composition preserves
the invariant across epochs and regions.

Region matrices average the one active epoch over its L/t positions.
A polynomial recurrence then averages the w active regions over their
binomial(B,w) placements. The implementation batches tilts, uses binary
matrix powering for epochs, and keeps the four-state contractions explicit.

The correction is substantial at small s. At BCH-128, k=2^20, t=64, s=7,
the new margin is 18.171 bits, an improvement of 64.623 bits over the old
envelope evaluated on the same witnesses. At s=20 the improvement is only
0.009 bits. Thus the earlier small-state curves contained significant
avoidable proof loss; the high-state BCH-128 level was already stable.

## Why the larger exact BCH curves plateau

Consider a simplified persistent-state model. After the first activation,
the state never cancels, active output has rate exactly one half, and the
activation position inside its region is continuous and uniform.
This idealization omits finite-epoch output fluctuations and cancellations.

Fix an outer weight w and sample its support uniformly among the w-subsets
of the B regions. Let H count the regions from its first active region
through the end, including that first region. Then

    P(H=h) = binomial(h-1,w-1) / binomial(B,w),  w <= h <= B.

Let V be uniform on [0,1], independently of H. The modeled output satisfies
Y/L=(H-V)/2. With delta=0.1, put m=floor(2 delta B) and f=2 delta B-m.
The exact failure probability within this model is

    p_w = [binomial(m,w) + f binomial(m,w-1)] / binomial(B,w),

where invalid binomial coefficients are zero. This identity follows by
conditioning on H: values at most m contribute one, value m+1 contributes
f, and larger values contribute zero.

For scaled tilt a=lambda L, its log moment is

    log g_w(a) = log((1-exp(-a/2))/(a/2))
                 + log sum_h P(H=h) exp(-a(h-1)/2).

Optimize exp(a delta B) g_w(a) separately for each shell, then sum with N_w.
The following intercepts omit the outer row factor, just as C_B does:

| Constituent | Exact probability in onset model | Chernoff bound in onset model | Actual finite transfer C_B at k=2^20, t64/s20 |
| --- | ---: | ---: | ---: |
| BCH-64 | 28.443 bits | 25.408 bits | 25.387 bits |
| BCH-128 | 50.099 bits | 46.804 bits | 46.770 bits |

The finite transfer lies within 0.035 bits of this model's Chernoff level.
The contributing onset shells are weights 12 and 13 for BCH-64, and 22,
24, and 26 for BCH-128. Their counts and placement probabilities explain
the plateau more directly than a fit in constituent length alone.

For BCH-8 and BCH-32, every nonzero shell satisfies
(w-1)/2 >= delta B. Their modeled onset failure probability is zero.
The observed finite margin therefore depends on effects this idealization
omits. The roughly one-bit gain per added state bit is consistent with a
cancellation contribution proportional to 2^-s. Establishing its coefficient
and uniform error term remains open.

The roughly three-bit gap between the two model columns is Chernoff loss
inside the idealization. It suggests investigating a distributional bound;
it does not prove that the actual encoder's margin can improve by that amount.

## Consequences for parameter selection and extrapolation

At k=2^20 and common state sizes s=9 through 20, the three t curves differ
by less than 0.030 bits for every constituent. BCH-128 at t64 has Q1 margins
30.197, 32.495, 32.710, and 32.770 bits for s=12,16,18,20. The diminishing
returns are visible without choosing a security threshold.

This Q1 evidence supports considering t=64,128,256 on implementation cost
and higher-occupation behavior. It does not establish equal full-distance
bounds or equal runtime. No performance benchmarks were run here.

In the observed plateau regime, a useful planning relation is

    M1 approximately C_B - log2(k/D).

For a desired Q1 margin beta, it suggests looking for a constituent with
C_B at least beta+log2(k/D). This relation describes the measured local
message-length trend; beta can be any benchmark. Moving from BCH-64 to
BCH-128 raises C_B by about 21.4 bits at the reference point. Extrapolating
that gain to arbitrarily large k or to a full-distance claim is unsupported.

The next size model should predict low-weight shell counts and their onset
cost, while separately accounting for cancellation and finite-epoch terms.
Two exact constituents in the nonzero-onset regime do not yet determine a
robust law for that spectrum tail. The earlier fit through all four BCH
sizes should not be used as a reliable forecast.

BCH-256 remains a secondary calibration point with bounded spectral input.
The other workstream reports a completed outward closure. This study neither
replays that closure nor treats its spectrum bounds as exact counts.
A useful comparison must match the map, cutoff, t/s, and transfer, then
identify how much of the discrepancy comes from the outer caps.

The next smallest useful task is to transfer the uniform-refresh invariant
to selected higher occupations for BCH-64 and BCH-128. Inputs with multiple
bits may have zero syndrome, so the one-bit transition must not be copied
unchanged. Compare continuous contributions at a few matched t/s settings
before launching another whole grid. This will determine whether the Q1
size trend also governs the complete union or hides a different bottleneck.

## Reproduction and checks

Run these commands sequentially in this directory, after restoring or
regenerating the local exact-spectrum pilot inputs described in the README:

```text
python study_bch_growth.py
python verify_bch_refresh.py
python plot_bch_growth.py
python -m unittest -v test_bch_growth.py
```

The study authenticates its source receipts and writes source hashes and
the CSV hash to `bch_growth.json`. Generated CSV, JSON, logs, and figures
remain local and ignored. The new study is separate from the old complete
grid ledger; no existing certificate or producer receipt is rewritten.

Five tests check exact rational domination of an actual GF(4) encoder,
agreement of the log recurrence with rational arithmetic, removal of repeated
zero-epoch density loss, exhaustive onset support counts, numerical
quadrature of the model moment, and Chernoff domination within that model.
An independent 90-digit positive-arithmetic replay checks all 236 weight
coefficients at the four t64/s20, k=2^16 reference points. Its maximum
absolute log error is 1.081e-12. That replay uses linear epoch iteration,
separate matrices, and positive polynomial arithmetic. It validates numerical
evaluation at those points; it is not an outward proof of every plotted row.
