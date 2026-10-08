# Mixed intervals reject the regional two-state candidate

The proposed two-state additive inner does not reach 40% minimum distance
with a 40-bit setup-failure guarantee. At K=2^20 and rate one half, the
construction below fails 40% distance with probability greater than 499/500.
The probability is over the regional permutations. The conclusion holds
for every valid choice of field coefficients, including route-dependent choices.

The previous cut used only safe intervals. This cut permits unsafe holes
inside an interval, charging one rank dimension for each hole. Joining
safe runs can save two endpoint dimensions while introducing fewer unsafe
dimensions. An exact interval-selection rule exploits this tradeoff.

## Construction and the shared check space

Work over F=GF(2^16). Let b=64, k=32, R=32768, K=kR, and N=bR.
Apply an injective [b,k] constituent independently to each of R message rows.
Place coordinate j of row i at position jR+i. Let B denote this outer map.
Its full-rank parity check H0 has R groups of k rows; column jR+i is
supported only in check group i. The proposed GRS constituent has this form.

Let Q be the composition of two input-weighted prefix scans, with nonzero
weights and zero initial states. Both scans carry state across region
boundaries. Their inverse Q^-1 has lower bandwidth two: column t is
supported on t,t+1,t+2, where those positions exist.

Let P independently permute each of the b regions of length R. Each
regional permutation is uniform. After routing, apply J, initialized with
s1=s2=0 and defined at each position by

    y = s1 + f*s2 + c*x,
    s1' = s1 + d1*x,
    s2' = s2 + d2*x.

The proposed setup requires c,d1,d2,f nonzero and c+d1+f*d2 nonzero.
All N outputs are retained, and the final states are discarded. States
continue across region boundaries. The final code has parity check

    H = H0 Q^-1 P^-1 J^-1 = M J^-1,
    M := H0 Q^-1 P^-1.

Choose a=8192 source-row groups, specifically the first a groups. Let U
be their check-coordinate space, of dimension h=ka=262144. Before routing,
each region has A=a-2=8190 designated safe positions: positions 0 through
a-3. The inverse precode sends each safe coordinate only to the first a
source-row groups. Its column in H0 Q^-1 therefore lies in U.

After routing, the safe positions in each region form a uniform A-subset.
These subsets are independent across regions. Every safe column of M
lies in U, regardless of all scan coefficients. The source-row groups
remain fixed; only the intervals selected below depend on the realized route.

## The mixed-interval rank bound

The inverse postscan has two states. Writing s=(s1,s2), d=(d1,d2),
and l=(1,f), its recurrence is

    x  = c^-1 * (y - l*s),
    s' = (I - d*c^-1*l)*s + d*c^-1*y.

Subtraction is addition in F. In particular, J is invertible because its
input/output matrix is triangular with nonzero diagonal c.

Fix an interval I of output coordinates. An input to J^-1 supported in I
has zero response before I. Its response after I depends only on the two
states at the end of I. The inverse response outside I thus has rank at
most two. This claim holds for all coefficients, including rank-deficient
state transitions; no generic-rank assumption is used.

Suppose I contains m safe and u unsafe positions. Multiplication by M
sends the within-I safe part into U. The within-I unsafe part lies in the
span of u columns of M. The remaining tail contributes at most two dimensions.
Consequently, for r disjoint intervals whose union is S,

    rank(H_S) <= h + unsafe(S) + 2r.

This is a global rank bound. The intervals need not reset the encoder
states, and their tails need not be independent. Overlapping tail spaces
only reduce the rank. It is unnecessary to solve independent local
cancellation problems or to assume that endpoint constraints are independent.

Define the gain and cost of S by

    G := safe(S) - 2r,
    C := unsafe(S) + 2r.

Then |S|=G+C and rank(H_S)<=h+C. Set e=ceil(2N/5)-1=838860.
Whenever

    G > h,                 C <= e-h-1,                      (1)

the columns H_S are dependent. Some subset of at most h+C+1 columns is
dependent, so a nonzero vector in ker(H) has support at most h+C+1<=e.
Because H is a full-rank check for J P Q B, this vector is a codeword of
the actual encoder. Thus (1) implies failure of the 40% distance target.

## A fixed interval policy with exact moments

Partition every region into consecutive blocks of length L=256. The
selection rule operates independently on each block and does not join
intervals across block boundaries. This restriction can only weaken the cut.

For a realized safe-bit string, choose disjoint intervals maximizing

    3G-C = 3*safe(S) - unsafe(S) - 8r.

Use the fixed tie rule implemented in `theory/hd_two_state_cut.py`.
This rule is part of the witness construction, not an assumption that
the encoder itself performs optimization.

The interval dynamic program retains an outside score O and an inside
score I. Before processing the next bit, let D=I-O. The next bit has
weight w=3 if safe and w=-1 otherwise. Then

    O' = O + max(0,D),
    D' = w + min(0,max(-8,D)).

Opening a new interval costs eight score units. Initialize O=0 and D=-8;
after the final bit, take the larger of the inside and outside scores.
The reachable difference D is an integer between -9 and 3. Ties retain
the outside solution at D=0 and start a new interval at D=-8.

The program aggregates all block words by their Hamming weight and D.
It also retains the first and second moments of the gain and cost, with
the cross-moments needed for the affine updates. No expectation is passed
through a maximum. The result is an exact integer sum over every block word,
not a Monte Carlo estimate or a stationary approximation.

To specify the without-replacement calculation, let f be either G or C
for one block. For 0<=j<=L, define

    F_j := sum of f over all length-L words with j safe bits,
    Q_j := sum of f^2 over the same words.

Write D0=binom(R,A) and m=R/L. The probability of any particular block
word with j safe bits is binom(R-L,A-j)/D0. For two disjoint specified
blocks, the probability of any specified word pair is
binom(R-2L,A-j-k)/D0. Hence

    mu = sum_j F_j * binom(R-L,A-j) / D0,
    nu = sum_j Q_j * binom(R-L,A-j) / D0,
    xi = sum_j,k F_j*F_k * binom(R-2L,A-j-k) / D0.

The total statistic within one region has mean m*mu and variance

    m*nu + m*(m-1)*xi - (m*mu)^2.                            (2)

Equation (2) includes the covariance between distinct blocks in the same
fixed-cardinality region. It does not assume independent block occupancies.
Independent regional permutations let us multiply each regional mean and
variance by b. Every calculation uses integer arithmetic and exact fractions.

## Finite K20 rejection

The exact computations give the following values, displayed approximately:

| Statistic | Mean | Variance | Slack to the adverse threshold |
|---|---:|---:|---:|
| G | 277761.3096984583 | 234772.6673994744 | E[G]-262144 = 15617.3096984583 |
| C | 518633.6326814108 | 1685404.511861882 | 576715-E[C] = 58081.3673185892 |

Here e-h-1=576715. Apply Chebyshev's inequality to G<=h and C>576715,
then take a union bound. No independence between G and C is needed:

    Pr[distance <= e]
      >= 1 - Var(G)/(E[G]-h)^2 - Var(C)/(576715-E[C])^2
      > 499/500.

The exact lower bound is approximately 0.9985378145964314. The comparison
with 499/500 uses exact fractions; displayed decimals do not establish it.
The randomness is solely the independent uniform regional permutations.
The favorable cut event implies failure for every valid coefficient choice,
even if those coefficients are selected after observing the route.

For 128-bit records comprising eight parallel F lanes, put the scalar
codeword in one lane and zero the others. Scalar support then equals record
weight, so the same failure bound applies to the requested record distance.

## Checks and consequence

`theory/hd_two_state_cut_tests.py` checks the interval policy against every
interval union on short strings. It enumerates block words to verify all
first and second moment sums. Separate exhaustive fixed-cardinality checks
verify (2), including the covariance between blocks.

GF4 tests build the forward two-state matrix from its recurrence and invert
it independently. They verify mixed-interval ranks with arbitrary safe and
unsafe check columns. A complete small striped outer/precode/route/postscan
fixture constructs a nonzero codeword from the claimed column dependence.
Its states continue across every region boundary.
All eight retained tests pass. No native build or benchmark was run.

The rejection uses only two inverse-tail dimensions, the regional safe
subsets, and the bandwidth-two inverse precode. It does not depend on a
special collapse of the varying observation coefficient. A new candidate
must change at least one of those structural properties. Measuring this
rejected encoder is unnecessary for the stated distance goal.

## More states and removal of the prescans

For an inverse with at most s tail dimensions, the same interval bound
uses G=safe(S)-s*r and C=unsafe(S)+s*r. A separate bounded numerical screen
examines stationary interval policies at safe density 0.05. It mixes two
policies on disjoint large blocks until their limiting gain equals 0.025N.
The policies maximize G-lambda*C for lambda in {1/17,1/18,1/19}.

| Inverse-tail dimensions s | Numerical limiting support threshold |
|---|---:|
| 3 | 0.359631 N |
| 4 | 0.391820 N |
| 5 | 0.412004 N |
| 8 | 0.443587 N |

These are stationary floating-point screens, not finite failure bounds.
The calculations do not include finite block-boundary losses or establish
strict finite gain slack. Dimensions three and four remain plausible cut
targets. Five is the smallest tested dimension not rejected by this bounded
screen; that statement supplies no distance guarantee or exhaustive optimum.
The helper is `theory/hd_two_state_cut_screen.py`.

Removing both prescans gives a different, stronger obstruction. Consider
the direct striped GRS[64,32] outer, followed by regional routing and any
causal s-state linear postscan. Assume zero input and zero state produce
zero output, as they do for the proposed additive-state recurrence.

In each routed region, retain its first W=13107 positions. Their combined
size is bW=838848<0.4N. For outer row i, let t_i count its coordinates
inside these prefixes. Its shortened outer-code dimension is

    d_i = max(0,t_i-32).

Let D=sum_i d_i. The b regional permutations are independent, so each
t_i has distribution Bin(64,W/R). Distinct rows need not be independent.
The exact expectation is

    E[D] = R * sum_{t=33}^{64} (t-32)*binom(64,t)
                                *(W/R)^t*(1-W/R)^(64-t)
         = approximately 2845.5563156975436.

Restrict the outer message to the resulting D-dimensional subspace, whose
routed input vanishes outside the prefixes. Requiring zero postscan state
at each prefix end adds at most bs linear equations. If D>bs, a nonzero
message remains. Starting with zero state, each region's trailing zero
input then has zero output and carries zero state into the next region.
The whole codeword has weight at most bW.

Since 0<=D<=K, the elementary bounded-variable estimate gives

    Pr[D>bs] >= max(0,(E[D]-bs)/(K-bs))

when bs<K. Exact fraction calculations give:

| States s | Failure lower bound for the direct construction |
|---|---:|
| 3 | 0.002531091962198530 |
| 4 | 0.002470196424467285 |
| 8 | 0.002226539901854795 |
| 44 | 0.000028263000781770 |

All exceed 2^-40. With s=45 this particular expectation bound becomes
zero, which is not evidence of distance. Thus dropping the prescans and
adding three, four, or eight additive states does not meet the requested
setup-failure target. The helper `direct_window_bound` also checks the
claim without assuming independent row occupancies.

The next useful task is a structural review of a topology that avoids
both the interval-hole cut and the direct prefix-window cut. Increasing
the state count modestly, or deleting the prescans alone, is not yet a
credible implementation path.

```text
python -B research/workstreams/char2_field_spin/theory/hd_two_state_cut.py
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p hd_two_state_cut_tests.py -v
python -B research/workstreams/char2_field_spin/theory/hd_two_state_cut_screen.py
```
