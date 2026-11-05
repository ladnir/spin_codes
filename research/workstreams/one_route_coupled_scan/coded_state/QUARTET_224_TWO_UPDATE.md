# Closing the 224/32 family without changing the construction

2026-10-26. Bounded follow-up to
[the quartet trajectory diagnostic](QUARTET_NATIVE_TRAJECTORIES.md).

## Result and scope

The specified family closes at 46% distance with a directed exponent below
-104.10. No outer, route, inner map, acceptance rule, or implementation changed.
The gap came from two counting relaxations, not a demonstrated code defect.

More precisely, use the retained construction over J=GF(2^32): K=2^20,
N=2^21, 8192 fixed RS[256,128] outer rows, 256 independent uniform regional
permutations, and 512 state coordinates over 4096 updates. Each inner update
has an independent uniform ordered quartet partition and independently
accepted field-multiplication maps. Acceptance requires A and A+delta I
to be superregular, not the stronger mixed-diagonal rule. Setup is sampled
once and fixed for all messages; Q_0=0 and the terminal state is free.

Let E be the event that some nonzero J-message satisfies all of:

- exactly 1974 outer rows are inactive;
- wt(Q_t)=224 for every t=1,...,4096;
- exactly 32 coordinates leave the state support at each t>=2;
- wt(Y)<=floor(0.46 N)=964689.

The frozen replay proves

    Pr_setup[E] <= 2^(-4578692708060433/43980465111040)
                 < 2^-104.10.

The union includes all row sets, messages, and labeled support histories
in this family. It is not a fixed-message estimate. It is not a whole-code
certificate, a 49% result, or coverage of other inactive-row counts.
The eventual whole-code proof over J transfers minimum distance to the
extension-field encoder. An individual extension-field message need not have
the same support history as one of its J-components, so this restricted
family statement is explicitly over J.

## What the two-update analysis establishes

Retain the local constraints and the conditional first-moment argument from
[the joint quartet count](RANDOM_PAIRING_SUFFICIENCY.md). For fixed positive
markers, average one transition over the fresh quartet partition and the
independent-label proof experiment. Its weighted row sum depends on the
incoming support only through its size. All incoming supports have size 224.

Consequently, if K_1 is this exact restricted row sum, composing two updates
and summing over their shared support gives K_1^2. Keeping that support's
labels, without additional information, cannot improve the bound. This is
an identity for the relaxed weighted count, not a statement that realized
state vectors or consecutive inputs are independent.

Likewise, two successive batches inside one region contain distinct outer
rows. The route places each row once in each region. Two batches across a
regional boundary can expose two evaluations of a row; two evaluations
alone do not impose an RS[256,128] parity relation. A stronger rowwise parity
argument would need information retained across many regions.

However, the previous transition estimate was not K_1: it replaced three
coefficient extractions by positive evaluations. The resulting loss repeats
4095 times. Removing this loss is the first useful refinement.

## Exact incoming, birth, and death counts

For one quartet, let U be its incoming support, V its next support, and T
its inactive-row input positions. Define a=|U|, b=|V minus U|, d=|U minus V|,
and j=|T|. Let P_(U,V,T)(ell) be the preceding polynomial counting permitted
selected output zeros. It uses only the retained A and A+delta I constraints.
Put c=j for a live quartet; for an empty quartet put c=|T intersect V|.

For a canonical U of size a, the local weight is

    binom(4,a) p^j (1-p)^(4-j)
      * (P_(U,V,T)(ell) / q^c)^theta.

Sum over V,T with the specified (a,b,d,j). Round the fractional power upward
to a 48-bit dyadic; this gives rational coefficients B_(a,b,d,j) dominating
the original local weights. There are 35 possible triples (a,b,d).

Introduce positive markers H,X,Y for a,b,d, respectively. Define

    C_(a,b,d,j) = B_(a,b,d,j) H^a X^b Y^d,
    Z = sum C_(a,b,d,j),
    rho_(a,b,d,j) = C_(a,b,d,j)/Z.

The artificial law rho is a coefficient-evaluation device, not an
independent-partner approximation to actual quartet occupancies. For 128
independent draws from rho, let mu be the probability that their a,b,d totals
are 224,32,32. Uniform-partition counting then bounds the transition kernel by

    K_1 <= Z^128 mu / (binom(512,224) H^224 X^32 Y^32).

The denominator binom(512,224) gives the actual incoming-support normalization.
The old evaluation bound replaced mu by one.

The frozen theta is 7/40. Only proof markers change; theta is not a code
parameter. A positive integer dynamic program proves

    mu <= 287823 / 2^31.

Each rational coefficient rounds upward on the 2^-31 grid. Each stage uses
integer multiplication and addition, followed by upward rounding. Truncating
degrees above (224,32,32) is exact because degrees are nonnegative. A scalar
maximum bound is checked before every stage to exclude uint64 overflow.
The DP uses three arrays of at most 225*33*33 uint64 entries, below 6 MiB.

Compared with the evaluation bound at the same frozen markers, this saves
52,682.90 bits over the 4095 transitions. It leaves a directed upper exponent
+1057.39. Exact local counts alone do not close the family.

## Retaining the regional inactive-row quotas

Exactly r=1974 rows are inactive, and each contributes one position to every
region. Thus each region contains exactly r inactive labels. The preceding
proof generated independent Bernoulli-p labels and conditioned on these
256 quotas, paying pi^-256, where

    pi=binom(8192,r)p^r(1-p)^(8192-r).

It then discarded the quota indicators in the numerator. We retain them.
This is row-membership information already present in the construction;
it introduces no new randomness or independence assumption.

Condition on every quartet's (a,b,d) in the positive polynomial product.
Its inactive-label count has the rational law

    nu_(a,b,d)(j) = C_(a,b,d,j) / sum_k C_(a,b,d,k).

These label-count laws factor across quartets in that product. Let

    s_(a,b,d) = sum_(j=0)^3 nu(j) nu(j+1).

A rational dual witness verifies all 35 inequalities

    s_(a,b,d) >= lambda_0 + lambda_a a + lambda_b b + lambda_d d.

Summing over a step with the exact quotas gives

    sum_groups s >= 128 lambda_0 +224 lambda_a +32 lambda_b +32 lambda_d
                  = s_* > 6.11513678717.

The helper stores the four rational dual coefficients and verifies their
inequalities exactly. Optimizer output is not trusted in the frozen replay.

### A uniform bound on a regional count

For an integer-valued variable with law nu and characteristic function phi,

    |phi(t)|^2
      =1-2 sum_(i<j) nu(i)nu(j)(1-cos((j-i)t))
      <=1-4s sin^2(t/2).

Keeping just adjacent pairs in the sum yields
|phi(t)|<=exp(-2s sin^2(t/2)). For independent variables the exponents add.
Fourier inversion therefore bounds every atom of their sum by

    (1/(2 pi)) integral_(-pi)^pi exp(-2S sin^2(t/2)) dt.

For |t|<=1, sin(|t|/2)>=(23/24)|t|/2. Extend the resulting Gaussian integral
to the real line and use pi>3. For the remaining interval, use
sin(|t|/2)>=23/48. The resulting elementary upper bound is

    cap(S) = min(1, 24/(23 sqrt(6S)) + 2^-floor(529S/1152)).

The helper rounds the square root outward with integer comparisons. Ordinary
regions have 16 transitions and hence S>=16s_*. The first region has only 15:
its first batch activates the zero state. The initial batch's remaining label
count is an independent convolution factor, which cannot increase the maximum
atom. Thus the numerator quota indicators contribute at most

    cap(15s_*) cap(16s_*)^255.

This holds uniformly over the quartet-type sequences satisfying each step's
exact quotas. It therefore remains valid after summing those sequences.
The exact-quota DP probability and this cap can be multiplied. The original
pi^-256 conditioning cost is still paid; it has not been removed or paid twice.

The frozen caps are

    cap_first <= 195622534269 / 4398046511104,
    cap_ordinary <= 757642817361 / 17592186044416.

Their directed log2 contribution is at most -1161.4953565386. Combining the
caps with the exact transition counts gives the claimed exponent -104.1074.

## Row-by-row shortened supports

The stronger condition that each active outer row has at least 129 covered
positions remains valid but is not needed for this closure. To retain it
exactly, condition on the inactive placements and the covered active positions
in every region. Active row identities are uniform bijections to the remaining
positions, independently across regions. For one fixed active row, the
coverage generating function is

    product_(regions j) ((1-v_j)+v_j z),

where v_j is the covered fraction among the 6218 active positions in region j.
Different rows are not independent within a region. Raising a one-row tail
probability to the 6218th power without justification would be invalid.

The relaxed tilted local distribution has mean covered-active count about
2.5964 per quartet, corresponding to roughly 219 positions per active row.
This is diagnostic only, before exact-count conditioning; it is not a
uniform statement about all geometries. It suggests the individual 129-floor
is not the first missing constraint here. The useful retained information
was the exact number of inactive rows per region.

## Directed result and validation

| Refinement at the frozen theta/markers | Log2 upper bound |
|---|---:|
| Positive evaluation of all transition-count markers | +53,740.2904 |
| Exact transition-count upper coefficient | +1,057.38794 |
| Also retain regional inactive-row quotas | -104.10742 |

The global accepted-map density allowance remains 2^theta, as in the earlier
family proof. All support, row-set, and selected-zero multiplicities remain
present. Logs use rational atanh-series intervals after dyadic normalization;
no floating-point log enters the certified exponent.

Run from the repository root:

    python -B research/workstreams/one_route_coupled_scan/coded_state/quartet_two_update.py --frozen

The replay performs no optimization and asserts the exact stored final
exponent. Checks include 26 integer-DP fixtures, six two-update row-sum
identities, 1050 enumerated partition/transition cases for the three-variable
coefficient identity, four rational characteristic-function inequalities,
six binomial atom bounds, 26 logarithm checks, every algebraic root interval,
all 35 dual inequalities, and integer overflow guards. The preceding frozen
trajectory replays also pass. Runs take about four seconds, are serial and
file-free, and have a 60-second watchdog and bounded arrays. At least 20 GiB
free space was retained. No benchmark or construction change was made.

Next test neighboring activity, turnover, and inactive-row counts using
these same refinements. The roughly 104-bit reserve covers this family only;
it cannot be assigned to all missing profiles without an explicit union.
The successful fixed-construction repair favors this bounded coverage pass
before reconsidering the inner.
