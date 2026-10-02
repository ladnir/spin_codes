# Retaining Variance in the Shuffle Comparison

The current dense proof pays a worst-case count-conditioning factor in
every region. That bound permits deterministic packet counts, even when
the comparison inputs still contain substantial Bernoulli randomness.
This note gives a bound that retains that randomness. The optional
`--variance-shuffle` path applies it within the full composition cover;
the original comparison remains available as a fallback.

## Local Comparison

Fix n independent Bernoulli variables with parameters p_1,...,p_n.
Let Y be their sum, mu=sum_i p_i, and V=sum_i p_i(1-p_i). Suppose
0<mu<n. Write b(k)=Pr(Y=k), and let h(k) be the mass function of
Bin(n,mu/n). After uniformly shuffling the Bernoulli variables, each
binary word of weight k has probability b(k)/binom(n,k). Its ratio to
the iid Bernoulli(mu/n) probability is therefore b(k)/h(k).

We first locate a maximum of that ratio: it suffices to check
k=floor(mu) and k=ceil(mu). This is a finite-binomial comparison,
not a Poisson approximation. Related density-ratio methods for a
Poisson reference appear in [Dümbgen and Wellner (2020)](https://arxiv.org/abs/1910.03444).
That paper's Poisson-reference bound is not imported into a SPIN
certificate.

## Where the Ratio Is Largest

First assume 0<p_i<1 and put q_i=p_i/(1-p_i). For an m-subset L,
assign probability proportional to product_{i in L} q_i. Denote this
subset distribution by W_m. Its inclusion probabilities are ordered
in the same way as the p_i: swapping i and j gives

    Pr_Wm(i in L)-Pr_Wm(j in L)
      = (q_i-q_j) e_{m-1}(q excluding i,j) / e_m(q),

where e_j denotes the elementary symmetric polynomial. Consequently,
the rearrangement inequality gives

    E_Wm[(1/m) sum_{i in L} p_i] >= mu/n.

Fix an integer k with mu<=k<n. Expanding the elementary symmetric sums
and then grouping terms by their (k+1)-subset gives the identity

    b(k)/((k+1)b(k+1))
      = E_W{k+1}[(1/(k+1)) sum_{i in L}
                    (1-p_i)/(p_i+sum_{j notin L} p_j)].

For fixed L, convexity in p_i bounds the inner average below by
g(a_L), where a_L is the average of p_i in L and

    g(a)=(1-a)/(mu-k a).

All denominators here are positive. On the interval containing these
averages, g is increasing and convex because k>=mu. A second application
of Jensen's inequality and the previous inclusion-probability bound give

    b(k)/((k+1)b(k+1)) >= g(mu/n)
                           = (n-mu)/(mu(n-k)).

The expression on the right is h(k)/((k+1)h(k+1)). Thus b(k)/h(k)
is nonincreasing for integer k>=mu. Apply the same argument to the
failure indicators 1-X_i to obtain the reverse monotonicity below mu.
The maximum is attained at floor(mu) or ceil(mu).

Parameters equal to zero or one follow by continuity. Replace every
p_i by (1-epsilon)p_i+epsilon mu/n; this preserves the mean and puts
all parameters strictly inside (0,1). The comparison binomial remains
fixed as epsilon tends to zero. If mu is zero or n, both distributions
are the same point mass and their density ratio is one.

## A Variance-Sensitive Upper Bound

Fourier inversion bounds every atom of Y by the integral of the modulus
of its characteristic function. For real t,

    |1-p_i+p_i exp(it)|^2 = 1-4 p_i(1-p_i) sin^2(t/2).

The inequality sqrt(1-x)<=exp(-x/2) therefore implies

    b(k) <= (1/(2 pi)) integral_{-pi}^{pi}
                exp(-2 V sin^2(t/2)) dt
          = exp(-V) I_0(V),

where I_0 is the modified Bessel function of order zero. The integral
decreases with V, so any proved lower bound V_0<=V can replace V.
Combining this atom bound with the location of the maximum gives

    max_k b(k)/h(k)
      <= max_{k in {floor(mu),ceil(mu)}} exp(-V_0) I_0(V_0) / h(k).

At V_0=0 the numerator is one. No central-limit approximation enters
this bound. `shuffle_variance.py` evaluates it with upward rounding.
Tests enumerate small rational Bernoulli families, check the exact
location of the maximum, and compare every ratio with the outward bound.

For a numerical illustration, take 2048 slots: 60 probabilities equal
to 1/2, 70 equal to 4/5, and the remaining 1918 equal to zero. Then
mu=86 and V=26.2. At 256-bit precision, the new density factor is below
1.783685; the existing fixed-count comparison with active-count cap
130 gives a factor above 27.67582. Repeated across 256 independent
comparison regions, this saves more than 1012 bits in the logarithm
of that comparison factor. This example is not a SPIN composition
cover or a distance claim.

## Uniform Bounds Within a Composition Cell

Fix the scalar cover's base activity tilt. Component type i has tilted
activity probability f_i and active-group indicator a_i in {0,1}.
Write v_i=f_i(1-f_i). A composition assigns proportions r_i to these
types, with r_i>=0 and sum_i r_i=1. A cell [l,h] and its occupancy
threshold q_min impose

    l <= sum_i r_i f_i <= h,
    sum_i r_i a_i >= q_min/G.

Each comparison region contains G independent Bernoulli inputs before
shuffling. Its count has mean mu=G sum_i r_i f_i and variance
V=G sum_i r_i v_i. These are properties of the positive comparison
measure, not independence assertions about actual BCH rows.

Choose rational slopes b and c with c>=0, and set

    a = min_i (v_i-b f_i-c a_i).

Then v_i>=a+b f_i+c a_i for every component, including the deterministic
zero and one types. Every composition in the cell therefore satisfies

    V >= V_0 := G max(0, a+min(b l,b h)+c q_min/G).

A floating linear program proposes b and c. The verifier reconstructs
a and V_0 with rational arithmetic; it does not trust the numerical
objective or assume the proposed slopes are optimal. The witness
stores the two rational slopes. A negative c is rejected.

The mean lies in [G l,G h]. For every integer k from floor(G l) through
ceil(G h), intersect that interval with [k-1,k+1]. Only means in this
intersection can make k one of the two central indices. At fixed k,
the binomial mass h(k) is unimodal as a function of its mean, attaining
its maximum at mean k. Its minimum on each intersection is therefore
at an endpoint. The verifier checks both endpoints with interval
arithmetic and takes the largest resulting density bound.

The cell uses the smaller of this factor and the original universal
factor. Cells touching mean zero or G retain the original factor.
Alternate activity tilts also retain the original factor: the cell
constraints control the base-tilted mean only. In particular, a variance
bound from one tilt is never silently reused for another tilt.

The search, saved metadata, dense replay, and full assembler support
this option. Older covers default to the original comparison. Tests
cover rational mean intervals across integer boundaries, exact small
compositions including deterministic types, invalid dual slopes, the
alternate-tilt fallback, and agreement between proposals and outward
evaluation.

## Current Evidence

At distance 8.5% and q>=129, coordinate 0.04 improves from a proposed
log2 bound of -119.0705 to -1041.2467. At coordinate 0.0405, a fresh
256-bit evaluation gives -1020.3369136198469, matching the proposal.
Selected points also remain negative with a lower dense handoff q=97.
The full 8.5% q>=97 cover has closed with 155 accepted cells and no
unresolved cells. No whole-code certificate yet uses this refinement;
independent replay verified the dense partition and all sparse occupancies
except q=4, where the original tilt grid missed its stopping budget.
A fresh 9% assembly uses a finer grid that closes this case.
The strongest complete result is 8.2% with
52.4502533222 bits of setup-failure margin, using an earlier dense bound.

The independent 256-bit dense replay gives aggregate log2 upper
-72.9420481758. This number excludes q=1--96 and is not the margin
of a whole-code certificate.

## Diagnostic: Charge the Count of Low-Variance Compositions

The cell minorant deliberately permits every composition in the cell.
Its worst variance need not occur in the compositions with the largest
outer count. `variance_split_probe.py` tests whether separating these
two cases helps. That diagnostic fixes the mean x to a single point.
The separate `variance_partition.py` module now evaluates the same
argument outward on whole mean intervals, as described below.

Partition the average variance v=sum_i r_i f_i(1-f_i) into intervals
[v_l,v_h]. Write c'_i for the positive, base-tilted component coefficients.
For any real eta and gamma, and any mu>=0, the total outer comparison
mass of compositions in one interval is at most

    (sum_i c'_i exp(eta f_i + mu a_i + gamma f_i(1-f_i)))^G
      exp(-G eta x - mu q_min - G min(gamma v_l,gamma v_h)).

This is the multinomial theorem with nonnegative terms, followed by the
mean, occupancy, and variance constraints. No optimizer optimality is
needed for the inequality. Each interval can then use variance lower
bound G v_l in the shuffle comparison. Adding its contribution to those
of the other intervals bounds their union; shared endpoints only cause
conservative double counting. The current diagnostic chooses rational
dual parameters numerically and evaluates scores in floating arithmetic.
Small exact composition counts test the witness formula, but those tests
do not substitute for an outward evaluation of the real instance.

At 9%, q>=97, and x=0.04, the existing birth-density and cell-variance
bounds give score +1806.99. A 16-bin variance partition reduces this to
+1200.38. Combining the birth-weight-class diagnostic with a 32-bin
partition gives +758.49. At x=0.06 that combined score is -239.32, but
the positive scores at lower means remain unresolved. All are selected
points with the same encoder and fixed inner/output witnesses; none
establishes a full 9% certificate or an actual distance obstruction.

## Outward Variance Partitions on Whole Mean Cells

The `--variance-shuffle --variance-bins 8` path splits the variance range
only when the base-tilt mean cell still needs sharpening. For a mean
cell [l,h], every composition has average variance in

    [0, min(h,1-l,1/4)].

This follows from f(1-f)<=f, f(1-f)<=1-f, and f(1-f)<=1/4 for each
component. The saved witness lists a contiguous partition of this
entire interval, together with rational eta, mu, and gamma for every
part. Replay rejects a gap, overlap, missing endpoint, or negative mu.
The number or placement of parts need not be numerically optimal.

For each part, replace eta x in the preceding counting inequality by
min(eta l,eta h). The verifier recomputes the component generating sum,
the three dual penalties, and the regional density bound using Arb.
Its variance lower bound is the part's lower endpoint times G. It may
also retain the previously checked whole-cell density bound when that
bound is smaller. The same inner comparison applies to all variance
parts; the verifier adds their outer-count-times-density contributions
before multiplying by that inner bound. Alternate activity tilts do
not use this partition.

The option and witnesses pass through dense search, saved metadata,
replay, and full assembly. Tests exhaust small compositions, including
deterministic component types, and compare against their exact shuffled
Bernoulli density ratios. Other tests reject damaged partitions and
check proposal/outward agreement. Existing covers omit the option and
retain their original interpretation.

At 8.6%, q>=97, the mean interval [0.029999,0.030001] improves from
proposed log2 bound +64.9955 to an outward bound of -256.7044267586677.
The interval [0.049999,0.050001] improves from -44.4356 to
-1285.4995392012914. Independent evaluations at 256- and 384-bit precision
agree to the displayed digits. The full dense cover has closed with
154 accepted cells and no unresolved cells; its independent full replay
was stopped after finding the same q=4 grid gap. Neither
these interval checks nor the selected sparse trials certify the whole
8.6% construction yet. The birth-weight-class refinement now has a
separate outward evaluator, but this 8.6% path does not use it.

## Retaining the Regional Packet Count

A single density factor charges every packet count for the worst ratio
b(j)/h(j). The inner can put most of its low-output moment on quite
different counts. We therefore retain a separate upper bound R_j for
each ratio and combine it with the inner before taking a regional
matrix power. This changes the proof, not the encoder or setup
distribution.

Fix one comparison composition. Region size is n=2048; there are 256
regions, each containing 64 consecutive inner steps of 32 packets.
Type i has tilted activity f_i and active-group indicator a_i. Its
frequency is r_i. The mean cell and variance part constrain

    x=sum_i r_i f_i in [l,h],
    v=sum_i r_i f_i(1-f_i) in [v_l,v_h],
    n sum_i r_i a_i >= q_min.

These independent Bernoulli activities belong to the positive outer
comparison measure. They are not an assertion that the actual BCH
coordinates are independent. After the independent regional shuffle,
conditional on packet count J=j, support is a uniform j-subset. The
nonzero GF(16) labels are independent and uniform.

### A Checked Upper Bound for Each Count

For a real MGF tilt t define

    g_i(t)=log(1-f_i+f_i exp(t)).

Choose rational b,c,d with d<=0, and recompute

    a=max_i (g_i(t)-b f_i-c f_i(1-f_i)-d a_i).

Summing this affine majorant over all n slots proves

    log E[exp(t J)] <= K(t)
      := n(a+max(b l,b h)+max(c v_l,c v_h))+d q_min.

Since the expectation includes the nonnegative contribution of J=j,
Pr(J=j)<=exp(K(t)-t j), for either sign of t. Linear programs propose
the rational slopes; the verifier checks the majorant on every finite
component type using outward arithmetic. Their numerical objectives
are not proof inputs.

Let h_x(j)=binom(n,j)x^j(1-x)^(n-j). For 0<l<=h<1, its logarithm is
concave in x, so its minimum on the mean cell is attained at an
endpoint. If C is the already verified uniform density cap for this
variance part, then

    R_j = min(C, min_t exp(K(t)-t j)/min(h_l(j),h_h(j)))

is a uniform upper bound on Pr(J=j)/h_x(j). Any finite list of checked
tilts is valid, including an empty list that retains C. The current
proposal uses t in {0,+/-1/4,+/-1/2,+/-1,+/-2,+/-4,+/-8,+/-16}.
Cells touching mean zero or one retain the existing proof path.

### Retaining Spread After Exponential Tilting

The bound exp(K(t)-t j) can be improved without changing its dual.
For the fixed composition, exponential tilting preserves independence
and replaces activity f_i by

    f_i(t)=f_i exp(t)/(1-f_i+f_i exp(t)).

If M(t)=E[exp(t J)] under the original comparison measure, the tilted
count J_t satisfies the exact identity

    Pr(J=j)=M(t) exp(-t j) Pr(J_t=j).

The preceding Chernoff bound used only Pr(J_t=j)<=1. Instead, retain
its variance. For every type with 0<f_i<1,

    f_i(t)(1-f_i(t))
      = f_i(1-f_i) exp(t)/(1-f_i+f_i exp(t))^2.

Let kappa(t) be a lower bound on the minimum of the ratio on the right
over the finite interior types. Deterministic types contribute zero
variance before and after tilting. Hence Var(J_t)>=n v_l kappa(t).
The universal kappa(t)=exp(-|t|) is also valid. Applying the Fourier
atom bound from above with V_t=n v_l kappa(t) gives

    Pr(J=j) <= exp(K(t)-t j) exp(-V_t) I_0(V_t).

For v_l=0 the extra factor is one. The implementation computes a lower
endpoint for V_t and an upper endpoint for the decreasing atom-bound
function. A saved `tilted_atom` flag selects this refinement for each
MGF witness; omission retains the original interpretation. The
selected-interval tool exposes it through `--tilted-atom`. This does
not change packet randomness or add inner updates.
The optional `--fine-tilts` adds a spacing-1/8 grid between -2 and 2;
these extra witnesses are checked by the same inequality.

The optional `--tilted-variance` also uses the finite component types
when bounding V_t. Put y_i(t)=f_i(t)(1-f_i(t)). Choose rational slopes
b,c,d with d>=0, and reconstruct the minorant

    a=min_i (y_i(t)-b f_i-c f_i(1-f_i)-d a_i).

Every composition in the cell and variance part then satisfies

    Var(J_t) >= n(a+min(b l,b h)+min(c v_l,c v_h))+d q_min.

The evaluator takes the larger of zero, this lower bound, and the
previous ratio-based lower bound. Linear programming proposes the
slopes, but outward evaluation checks the minorant on every type.
The saved MGF witness contains its three rational slopes. A negative
d is rejected; this sign is opposite to the upper-MGF occupancy slope.

### Ordered Inner Operators Conditional on Count

Let K_k be the fixed-occupancy birth-class comparison matrix for one
32-packet step with exactly k active packets. Its coordinates and
nonzero-label averaging are defined in the fixed-occupancy note. The
whole-region conditional matrix is bounded by

    P_j = [u^j](sum_{k=0}^{32} binom(32,k) u^k K_k)^64 / binom(n,j).

Here [u^j] selects a polynomial coefficient. The binomial factors count
the supports in each consecutive step; the denominator counts all
j-subsets of the region. Matrix products retain chronological order.
No commutation of distinct K_k is assumed. Outward polynomial matrix
squaring computes all coefficients together.

The base activity tilt tau is removed by multiplying each active input
by tau^-1. Uniformly over the mean cell, the resulting reference input
weights are bounded by w_0=1-l and w_1=h/tau. Set S=w_0+w_1 and p=w_1/S.
For this variance part define

    R = sum_{j=0}^n Pr(Bin(n,p)=j) R_j P_j.

The weighted regional comparison is bounded entrywise by S^n R.
Regions are independent conditional on the fixed comparison composition;
all entries and terminal coordinates are nonnegative. Thus their total
weighted inner moment is bounded by S^(256n) e_0 R^256 1. The vector e_0
starts the recursion in its zero state, and 1 sums all terminal comparison
coordinates. No state flush is required.

Multiply this quantity by the independently checked outer-count bound
for the variance part and by exp(lambda D), where D is the output-weight
cutoff and lambda>0 the output Chernoff tilt. Sum over the entire
variance partition, then over the complete mean-cell partition. Finally
add the separately regenerated sparse prefix. Only this final aggregate
can establish a whole-code certificate.

### Bound Count Masses Directly

The interval argument above separately enlarges inactive and active
weights to 1-l and h/tau. This can inflate a wide cell even when its
pointwise bounds are good. A direct count bound avoids that step.

For each j, let M_j upper-bound Pr(J=j) throughout the mean and variance
cell. The checked ratio cap C supplies C max_{x in [l,h]} h_x(j).
The maximum occurs at x equal to j/n clamped into [l,h]. The checked
MGF bounds supply additional caps exp(K(t)-tj), optionally multiplied
by the tilted atom bound. Take M_j to be the minimum of these caps
and one.

Undoing the base activity tilt contributes exactly tau^-j at count j.
Consequently, the weighted regional operator is bounded entrywise by

    R_direct = sum_{j=0}^n M_j tau^-j P_j.

For a fixed comparison composition, its exact count probabilities may
be substituted before applying the caps. Nonnegative matrices preserve
the resulting entrywise inequality under the ordered regional product.
The existing outer count and variance partition can therefore use
R_direct^256, without the separate S^(256n) factor. At a zero-width
mean cell, this recovers the count-ratio expression, except that the
additional probability cap of one can tighten it. The intended gain
is in interval width and coverage cost, not a new encoder parameter.

The selected-interval tool enables this bound with `--direct-counts`.
The saved witness records the Boolean `regional_direct_counts`; an
absent flag retains the old count-ratio interpretation. New regional
proposals in the scalar cover use direct counts. Replay still checks
the saved rational MGF witnesses and does not rerun their optimization.
Exact small-instance tests cover binomial maxima, Poisson-binomial
count masses, and the full ordered regional product. They include
deterministic types and noncommuting local matrices.

For three updates and distance 9.9%, fresh 256-bit checks give:

| Reference activity | Mean-cell center | Radius | Log2 upper bound |
|---|---|---|---:|
| 0.175 | 7/183 | 1/10000 | -70.9996193604 |
| 0.20 | 3/67 | 1/10000 | -99.0941787733 |

These checks use 16 variance parts, tilted atom bounds, the finer MGF
tilt grid, and finite-type tilted variance. The radius is 100 times
that of the earlier selected intervals. Independent 384-bit evaluation
agrees to the displayed digits. These bounds cover only the stated mean intervals and
q>=97, not the full dense domain or the sparse prefix.

### Implementation and Current Scope

`regional_count.py` implements checked count ratios, direct count masses, and ordered
polynomial matrices. Its proposed witness records every rational MGF
dual alongside the complete variance partition. Replay checks matching
intervals and counting duals, recomputes every MGF majorant, and does not
rerun the LP when a saved witness is supplied. The scalar cover's
`--regional-count` option requires birth classes and variance partitions.
Ordinary cells retain the existing cheaper bound.

Optional integer witness fields `regional_lazy_density_through` and
`regional_joint_return_through` sharpen the local matrices before their
ordered polynomial product. The first retains density through selected
lazy steps; the second uses a freshly cross-checked joint return census.
Their state-domination argument and selected 10% checks are recorded in
[the density note](BIRTH_DENSITY.md#retaining-density-through-lazy-updates).
The verifier applies the return refinement before the density change.
Omitted fields retain the original matrices, so earlier records remain
replayable without changing their interpretation.

Tests compare polynomial powers against exact ordered placement and
enumerate small comparison compositions, including deterministic packet
types. The full weighted regional sum is checked against those exact
compositions. Other tests reject mismatched partitions, invalid dual
signs, and an attempt to replay the option under a different model.
Additional tests enumerate tilted Bernoulli count laws and check the
new atom factor, including deterministic types and both tilt signs.

At distance 9.25% and q>=97, three mean intervals of radius 10^-6 give
the following log2 upper bounds:

| Reference activity | Mean-cell center | Log2 upper bound |
|---|---|---:|
| 0.15 | 9/281 | -467.0660567919 |
| 0.175 | 7/183 | -421.5893617014 |
| 0.20 | 3/67 | -630.0255413580 |

Independent 256- and 384-bit evaluations agree to the displayed digits.
These are complete bounds on the stated intervals and all their variance
parts, not coverage of the entire mean domain. No whole-code 9.25%
claim follows yet. At 9.5%, floating diagnostics still leave positive
bounds at two difficult activities even with 32 variance parts. That
failure is not evidence of a low-weight codeword.

The tilted-atom refinement improves those 9.5% bounds. With eight
variance parts, fresh 256-bit checks give -17.3739 at reference activity
0.15 and +91.1114 at 0.20. With 32 parts and the finer MGF-tilt grid,
the bounds are -107.7590 and -9.0172. These use the same mean-cell centers
and radius as the table. The second interval still misses 40 bits, and
neither this pair nor the earlier three intervals covers the full
composition domain. The encoder remains unchanged.
The 64-part follow-up also checks reference activity 0.175, which the
two-point trial omitted. Its interval bound is still positive,
+50.5177. Thus even the selected 9.5% frontier is not closed; the
negative endpoints must not be interpolated across this gap.

An output-tilt scan at that interval finds only a modest improvement:
scaling the original tilt by 39/40 gives floating log2 proposal +41.7663.
With the finite-type tilted-variance minorant as well, a fresh 256-bit
evaluation gives +37.3798. This still does not meet the target. It is
not evidence that the encoder has a codeword below 9.5% distance.
