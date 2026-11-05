# Regional rank and partial-state resets

For the fixed three-shear model, every history with at most five defective
incoming batches now has a restricted 46% certificate below 2^-78.
The route exceptions are included. Six-defect coverage remains incomplete,
but its first sparse residual closes after charging the short resets jointly.
Neither result is a whole-code distance certificate.

Definitions and the construction are unchanged from
[UNIFORM_PAIR_DEFECTS.md](UNIFORM_PAIR_DEFECTS.md). In particular, q=2^32,
N=2^21, K=2^20, b=512, L=4096, and M=8192. The RS outer has 256 coordinates
per row and dimension 128. Each region independently permutes its M rows.
Maps and routes are sampled once for all messages.

## Why the fifth defect exposed another counting loss

The old no-full-reset bound fails near h=5,a=1,c=6,T=2058. It permits
48 proper incoming states, leaving at least 2010 full batches. Its
deterministic parity-check rank bound is 968584. That yields a positive
log2 contribution of 427855 at 46%.

Let G consist of these full batches. Only the regions intersecting interval
endpoints have fractional occupancy. There are at most twelve such regions.
Retaining those occupancies bounds the shortened outer dimension on G by
3330 outside an explicit route exception. The parity-check rank is then
at least 2010*512-3330=1025790. The same class's good-route contribution
is at most 2^-1402737. A simultaneous version of this argument, rather than
a selected-point calculation, supplies the certificate below.

There was also a geometric incompatibility. The previous bound separately
maximized the number of proper states and the number of regions reached
by short early intervals. Increasing an early interval's length reduces
the remaining terminal interval's length. The new count retains this
dependence and the defects required by each early interval.

## An actual-route bound for shortened dimension

Fix a set G of complete batches. Suppose it is a union of at most m
intervals and contains A batches. Write mu=A/16 for the mean number of
included coordinates per RS row. At most j_max=2m regions are partial.
In an actual occupancy profile, let C regions be full and j regions partial,
with inclusion fractions p_1,...,p_j in {1/16,...,15/16}. Then

    sum_l p_l = mu-C.

For row i, let f_i be its number of included coordinates. The shortened
dimension on G is d(G)=sum_i(f_i-128)_+. Equivalently,

    d(G) = sum_(ell=129)^256 #{i: f_i>=ell}.

Fix a threshold ell and put u=ell-C. If 1<=u<=j, choose a moment order
k in {1,...,u}. The indicator that a row hits at least u partial regions
is bounded by binom(S,k)/binom(u,k), where S is its number of partial hits.
Expand binom(S,k) as a positive sum over subsets of regions.

For s specified rows, multiply these polynomials. In each region, a monomial
requiring r distinct rows has probability (|G_l|)_r/(M)_r<=p_l^r.
Independence across regions and summation of the positive expansion give

    Pr[all s rows have f_i>=ell]
      <= [e_k(p)/binom(u,k)]^s
      <= [binom(j,k)/binom(u,k) * ((mu-C)/j)^k]^s.

The second inequality follows by pairwise averaging in the elementary
symmetric polynomial, as in the four-defect argument. Choose k before
taking the worst actual (C,j) profile. The resulting rational bound is
called beta_ell. The cases u<=0 and u>j give bounds one and zero.

There is no independence assumption between row events. Consequently,

    Pr[#{i:f_i>=ell} >= s] <= binom(M,s)*beta_ell^s.

The verifier chooses a cap for each threshold whose failure is at most
2^-B. Summing those caps gives d(G). All fractions and comparisons are
exact. Logarithms are bounded upward by raising fractions to the 256th
power and comparing integer bit lengths; every chosen cap satisfies the
resulting directed inequality. A cap M is deterministic and costs no event.

The complement provides a second rank bound. With H_G the outer parity
check restricted to G,

    rank(H_G) = |G|-d(G) = K-d(complement G).

The complement has the same partial-region count. Apply the threshold
argument to both sides and take the stronger rank lower bound. This is
useful on both sides of the half-length crossover.

## Make the route statement simultaneous

The full incoming batches form at most h+1 intervals. Retain an initial
subset if only a lower bound on their number is known. Truncating the last
retained interval does not increase the interval count. For memoization,
the verifier rounds that number downward to a multiple of eight batches.
This can only weaken the rank lower bound.

The possible input/output support I consists of live intervals expanded by
their activation batches, and possibly the terminal batch. It has at most
h+2 intervals. To bound d(I), the verifier pads its batch count upward to
a multiple of eight. Extend existing intervals into gaps, merging them
when necessary. This does not increase their number or decrease dimension.

For at most m intervals there are fewer than 2m*(L+1)^(2m) endpoint choices,
including every length. For each set, union over at most 128 thresholds
and the three uses of the dimension bound. Nine additional bits suffice.
The per-threshold budget is selected to leave 80 bits after this entire
union. Thus the event is uniform over setup-adaptive interval choices.
It is not a bound for one support fixed independently of the codeword.

On the complementary route event, use the previously proved inverse bound
with the new rank. The forward bound uses d(I) rather than its worst-case
deterministic value. The density inflation for all pair maps is still below
two. No further randomness is inserted in the construction.

## Couple the early lengths to their defect budget

Retain the branch with no full incoming state followed by a zero incoming
state before the terminal state. If any full state occurs, the final live
interval continues through Q_(L-1). All earlier intervals are proper.

Of h defects, c-a activate the nonzero intervals from zero. Put e=h-c+a.
For an early interval of length l, the clean-run argument gives

    l <= 8+9e_i,

where e_i counts its internal defective batches. Therefore that length
requires at least ceil((l-8)_+/9) extra defects. Let S be the combined
early length, and let u be the sum of these minimum requirements.
The final interval has at most 8+9(e-u) proper states. Hence

    P <= S+8+9(e-u).

An integer dynamic program retains (S,u) jointly with the maximum region
cover of the early intervals. Their first start is aligned when a=1.
The final interval has an aligned end and reaches at most
ceil((T-S+1)/16) regions after including its activation batch.

The actual interval placement determines its early lengths and u. Thus
the endpoint union already covers these classes; no additional independent
union over S or u is needed. The verifier checks all possible states to
obtain a uniform bound for each placement. It does not discard histories
without a full state: the previous bound handles them, and the helper
explicitly retains any unresolved such branch.

## Complete five-defect result

The no-full-reset replay checks 45020 count cells, excludes 21792 by the
outer distance, and checks 29991 joint refinements. There is no residual.
The full-reset replay checks the same 45020 cells, excludes 21466, and
refines 1085. On the good route, it takes the best of the inverse count,
the previous direct bound, and the full-reset probability with the new d(I).
There is no residual in this branch either.

| Branch | Bad-word contribution on good routes | Route exception |
| --- | ---: | ---: |
| No full reset | <=2^-1253 | <=2^-80 |
| At least one full reset | <=2^-1599 | <=2^-80 |

Each branch totals less than 2^-79. Their union is less than 2^-78 at
W=floor(.46N). This covers every J-valued message whose reconstructed
history has at most five defective incoming batches. Higher defect counts
and whole-code closure over J or its extension field remain open.

## Six defects: a discriminating failure and its first repair

The same broad replay does not close h=6. The no-full-reset branch leaves
137 coarse cells, with log2 upper 6491. The full-reset branch leaves 25,
with upper 9147. These are failed upper bounds, not bad codewords.

The first no-full-reset residual is h<=6,a=0,c=6,T=1893. Its allowed input
dimension is bounded by 21, but the existing forward bound charges no
short resets when its random-coordinate count is below the weight limit.
The geometry determines what that count omitted.

All six defects activate an interval from zero. There are no one-sided
pairs or internal births. Each of the five early intervals has length
at most eight. If their combined length is S, they reach at most ten
regions; the final interval reaches at most ceil((1894-S)/16). Nonzero
RS input requires 129 regions. Since S>=5, this forces S=5: all five
early intervals have length one and straddle distinct regional boundaries.

The resulting I reaches exactly 129 regions. It has at most eleven partial
regions and total mass (1893+6)/16=118+11/16 in row-coordinate units.
Each partial region has at least 1/16 occupancy. These facts force exactly
eleven partial regions, each with occupancy 1/16; the others are full.

Fix s active RS rows. Each eligible row has a one-dimensional shortened
code and is nonzero in every one of the 129 allowed regions. Eligibility
of s specified rows has probability at most 2^(-44s), by the actual
regional falling-factorial bound. Each early activation batch therefore
contains exactly s nonzero inputs. A one-row message would activate a
one-sided pair, which this class forbids. Thus s>=2.

Each early episode has one incoming state followed immediately by zero.
For a fixed message and route, its s activation coordinates are known
before exposing that update. The pair-image law charges q^-s, with the
same total density factor below two. The five resets jointly charge
q^(-5s). There are at most q^s messages on the chosen rows. Consequently
the per-placement contribution is bounded by

    2 * sum_(s=2)^M binom(M,s) * (2^-44*q^-4)^s.

With y=M*2^-44*q^-4=2^-159, this is at most y^2/(1-y)<2^-317.
Even the broad h<=6 endpoint union of 186 bits leaves a restricted bound
below 2^-131 for this entire profile. An exact check of all 32768 early
length vectors and the rational sum bound reproduces it.

This does not close the remaining six-defect cells. It identifies a
transferable next interface: count the initial active coordinates of
every closed partial episode, jointly with the outer row choices and
their routing. Charging only a fully active reset misses this resource.

More generally, for a fixed message and fixed episode placement, an
activation support of weight w is known before the future pair maps.
Resetting it within l updates requires an active-to-zero transition on
each of those w coordinates. Choosing their times and exposing maps in
order gives at most 2*l^w*q^-w. For multiple disjoint episodes, multiply
the time-choice factors and q powers, but use the total density factor
two only once. Turning this local statement into a useful all-message
count requires the joint source/route accounting illustrated above.

## Reproduction and next decision

Run serially from coded_state:

    python -B diagnose_pair_region_rank.py --test
    python -B diagnose_pair_region_rank.py
    python -B diagnose_pair_region_rank.py --joint-cover
    python -B diagnose_pair_region_rank.py --full-reset-cover
    python -B diagnose_pair_region_rank.py --joint-cover --defects 6
    python -B diagnose_pair_region_rank.py --full-reset-cover --defects 6
    python -B diagnose_pair_region_rank.py --sparse-reset

Tests include 1500 exact joint factorial-moment checks, 18 exact tail-cap
checks, 21 exhaustive joint geometry tables, and the 32768 sparse-reset
length vectors. Independent five-defect replays agree; they take about
6.4--8.2 seconds each. Six-defect diagnostics take about 12--15 seconds.
Coverage has a 45-second cap; selected checks have 20/30-second caps.
All arrays and caches are bounded. Only source and compact notes are kept.
C: had more than 23 GiB free. No raw output files or Python caches were made.

The previous defect helper remains unchanged, with SHA256
7238fec9e5553e7b4d1ce3edb3089cec06abcc34cb090c029d0f3d6b54bbaed3.
Both generic frozen proof-helper hashes also match their checkpoints.
No benchmark, encoder/paper edit, commit, push, or promotion occurred.

Continue the proof, not a redesign yet. Next extend the source-conditioned
partial-reset count to neighboring six-defect cells and assess whether it
can join the pin-transfer argument for widely distributed defects. Do not
equate another small-defect certificate with whole-code progress unless
the new interface also explains how to address the remaining families.
