# Regional episode bound through 64 active outer blocks

For the unscaled post-diagonal inner, the complete p=32,...,64 source class
has distance greater than49.1% except with probability at most2^-47.
Combining this with the [p<=31 bound](SPARSE_PROOF.md) gives failure less
than2^-46 for every nonzero message using at most64 outer blocks.
The probability now includes both the actual regional permutations and the
independent diagonal setup. Messages may be selected after setup.
Messages using65 or more outer blocks remain outside these certificates.

The construction is unchanged: Y=X+D A Q and Q'=Q+B Y, with D uniform
diagonal over GF(2^128), b=128, s=64, and fixed RS[256,128] outer blocks.
The field size is q=2^128. Write N=2^21, M=8192, and n=N/b=16384.

## Keep the terminal case separate

An episode consists of an activation batch followed by live batches until
a reset or the end of the scan. For p<=64, every nonzero input batch
activates an empty state. Thus every nonzero input belongs to an episode.
Let U be a specified union of c such intervals, with v batches in total.

There are two representations. Put delta=1 when the last interval reaches
the end of the scan, and delta=0 when it stops earlier. The number of
prescribed reset endpoints is a=c-delta. All closed intervals have length
at least2; the last terminal interval may have length1. The last state is
unrestricted in the terminal representation, including a possible reset
at the very last batch.

The number of representations at (v,c,delta) is at most

    R_delta(v,c) = binom(v-a-1,c-1) * binom(n-v+c-1,c-delta).

For delta=0, v>=2c and v<=n-1. All intervals have minimum length2, and
the final outside gap has minimum length1. For delta=1, v>=2c-1 and v<=n;
the terminal interval has minimum length1 and the final gap is zero.
Counting lengths and outside gaps gives the two binomial factors.
Adjacent episodes are allowed in this upper bound.

At most r=min(2c-delta,256) regions are partially occupied by U. The terminal
case saves one partial boundary; the nonterminal case pays an additional
reset. Ignoring both distinctions was the largest loss in the first probe.

## Jointly charge resets and defects

Fix the route, a candidate input supported in U, and all a reset endpoints.
The reset probability is at most q^(-64a). Use E and Z from SPARSE_PROOF.md:
E counts holes in nonzero refreshed expanded states, while Z counts zero
outputs at coordinates with nonzero feedback inputs.

Let Z_clean exclude output coordinates in the a prescribed reset batches.
At most128a coordinates were removed, so Z<=Z_clean+128a. An exact episode
trajectory has v-c live batches and at most c activation-created states.
Consequently output weight at most W implies

    E+Z_clean >=128v-191c-128a-W.

Define

    m = ceil(max(0,128v-191c-128a-W)/2),
    eta =1-log_q N=107/128,
    T =64a+eta*m.

At least one of E or Z_clean is at least m. A prescribed hole witness for E
cannot occur at a prescribed reset: its state is required nonzero. Thus
the chronological conditional-probability argument charges q^-64 at each
reset and q^-ell for ell holes at each other selected refresh step.
This gives q^(-64a)*(N/q)^m after counting witnesses.

For Z_clean, the selected diagonal entries occur in batches different from
all prescribed reset batches. Expose each batch chronologically and its
entries individually. The same joint upper bound follows. These two defect
events may be dependent; taking their union gives

    Pr_inner[specified episodes and wt(Y)<=W] <=2 q^-T.     (1)

For m=0 this is simply a harmless factor2 over the reset bound. No uniformity
of a state conditioned on future output events is assumed.

## Count shortened source words before averaging the route

Fix a set P of p outer blocks. Let t_i be the number of coordinates from
block i that the route places inside U. For every i in P to be nonzero,
t_i must be at least129. On this event the shortened source dimension is

    D_P(U)=sum_(i in P)(t_i-128).

There are fewer than2q^(D_P(U)-1) projective input representatives.
The construction and all events being counted are invariant under nonzero
scalar multiplication of the input. Combining this count with (1), the
conditional probability of a bad word with this P,U is at most

    1[all t_i>=129] * min(1,4 q^(D_P(U)-1-T)).              (2)

Fix 1<tau<=q and theta=ln(tau). Applying min(1,z)<=z^(theta/ln(q)) to (2)
leaves the prefactor tau^(1/64-1-T) and the route moment

    E_route product_(i in P) [1[t_i>=129] tau^(t_i-128)].

This fractional moment is applied after the all-message union for fixed
P,U. It does not replace GF(2^128) by a smaller field.

## Positive route moments for selected source blocks

In each region, positions in U select a uniform subset of outer-block
labels of a fixed cardinality. For any nonnegative increasing function g,
condition on the other regions and write g(t_i)=A_i+B_i Z_i, with A_i,B_i>=0.
Expanding the product over selected labels and using
(u)_d/(M)_d <=(u/M)^d shows that replacing this region by independent
Bernoulli indicators can only increase the expectation. Repeat by region.

Apply this to g(t)=1[t>=129] tau^(t-128). The block occupancies then become
independent. For one block, let mu=|U|/M=v/64. Its positive moment is bounded
using the convex function tau^((t-128)_+):

    F <= min(1+E, E*tau/(tau-1)),

where E bounds E[tau^((t-128)_+)-1]. Pad the partial-region count to r and
average its Bernoulli probabilities, as in the retained regional proof:

    E_r(mu,tau) = max_f E[tau^((f+Bin(r,(mu-f)/r)-128)_+)-1].

The integer f ranges over 0<=f<=256-r and 0<=mu-f<=r. This is a finite
sum of nonnegative terms. The padding/equalization applies to the convex
excess function, not to the discontinuous positive-indicator function.
Put F_r=min(1+E_r,E_r*tau/(tau-1)). The route moment in (2) is at most F_r^p.

If ceil(mu)+r-1<=128, no block has enough selected coordinates and the
contribution is zero. Otherwise a complete per-cell upper bound is

    binom(M,p) R_delta(v,c) tau^(1/64-1-T) F_r(mu,tau)^p.   (3)

One may minimize (3) over fixed markers. The minimization is over upper
bounds on the same event and changes neither the setup nor the code.

## Complete finite replay

`route_episode_certificate.py` covers all p=32,...,64 and every feasible
v,c in both terminal classes. A floating pass proposes disjoint rectangular
boxes and a rational marker per box. Every box then receives an independent
outward-rounded replay of (3), at80 and120 decimal digits.

For a box, binomial counts use upper entropy bounds at safe endpoint values.
The moment uses the largest v and partial-region count. The coefficient
-1-T uses the smallest reset cost and defect witness count in the box.
If F_r>=1, its power uses the largest p; otherwise it uses the smallest p.
The excess E_r and F_r are exact rational numbers. All logarithms and signed
multiplications are enclosed with the retained directed-arithmetic routines.

Result at W=floor(.491N)=1,029,701:

- 15,709 accepted boxes, from31,416 box evaluations.
- 4,429,185,024 feasible discrete cells, including both terminal cases.
- Every cell contributes at most2^-80 by its box bound.
- Largest directed log2 upper bound: -80.173269845194262... at both precisions.
- Fewer than2^33 cells, giving failure at most2^-47 for p=32,...,64.

The p<=31 failure bound is less than2^-50. Their sum is less than2^-46.
This is a combined sparse-source certificate, not a whole-code certificate.

## Independent checks and remaining obstruction

`test_route_episodes.py` checks128 small interval-count cases against explicit
enumeration. It checks2500 exact route-moment comparisons by enumerating
all selected-label subsets in four regions, without assuming label
independence. It also verifies pairwise disjointness of the complete box
cover on all66 source-count/terminal slices and matches an independent
triangular-domain count of4,429,185,024 cells.

Both precisions and the checks passed serial bounded runs. The initial
large-marker floating probe overflowed in a deterministic exponential;
a separate stable evaluator fixed that diagnostic without editing the
retained helper. No generated reports or coefficient dumps are needed.

Beyond64 active blocks, a nonzero input batch may lie in ker(B), producing
a nonzero output without activating the state. Such a batch has weight
at least65. The next proof must include these quiet, heavy batches rather
than assume that episodes cover all nonzero inputs. The exact refresh
and inverse lemmas remain valid. The older one-/two-tap designs and the
certified window construction remain unchanged.

Reproduction, run serially with Python -B:

    route_episode_certificate.py
    test_route_episodes.py

The first command needs NumPy for the floating proposal. Both commands
write no output files. The sources include bounded caches and runtime caps.
