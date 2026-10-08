# A complete obstruction for one uniform permutation and a causal inner

The proposed full-uniform-permutation convolution cannot meet the new
high-distance target. With repeated GRS[16,8] blocks, its 30% target fails
with probability at least 1-2^-282 at K=2^16. The bound strengthens with
length. It holds for every field, every convolution width, and every tap
assignment. Uniform GF256 taps and native polynomial taps have the same
obstruction.

This is a bound on the complete code's failure probability. It is not a
local transfer calculation or a large first moment of bad messages.

## The exact shortened-space identity

Fix a field F and an [b,k] MDS constituent. Apply that constituent to L
independent message blocks. The resulting outer code C_0 has length N=bL
and dimension K=kL. Let P be one uniform permutation of the N coordinates.
Let A be any invertible lower-triangular F-linear map on F^N. A convolution
with initial state zero and diagonal one is an example, for arbitrary taps
and arbitrary width. Define the final code C=A P C_0.

For an integer e, let J be the suffix consisting of the final e coordinates.
Let F^J denote the vectors supported inside J. Both A and A^-1 preserve F^J.
Therefore

    dim(C intersect F^J) = dim(P C_0 intersect F^J).

For outer block ell, let X_ell count its coordinates whose permuted positions
belong to J. MDS shortening gives the explicit identity

    dim(C intersect F^J)
       = sum_(ell=1)^L max(0, X_ell-(b-k)).                 (1)

This holds for each fixed permutation and each fixed A. The taps may even
be chosen after the permutation is known. If any X_ell>=b-k+1, equation
(1) gives a nonzero final codeword of weight at most e.

For GRS[16,8], nine coordinates from one block in J suffice. Choose the
nonzero outer word shortened to those coordinates and set all other message
blocks to zero. Its permuted source vanishes before J. The causal inner
cannot create any earlier nonzero output.

The companion's cut theorem identifies the same obstruction. Write H_0 for
the block-diagonal outer parity check and D=A^-1. The final check is
H=H_0 P^-1 D. Its columns indexed by J satisfy

    H_J = (H_0 P^-1)_J D_(J,J).

The square matrix D_(J,J) is invertible. Hence

    rank(H_J) = sum_ell min(X_ell,b-k).

J is a terminal interval container whose erosion is all of J. There is no
width-dependent expansion at this endpoint. Counting that cut once removes
message multiplicity, but its rank defect is real and survives every finite
coefficient specialization. The stronger counting framework therefore
rejects this construction rather than certifying it.

## An exact finite failure lower bound

Put h=b-k+1, and let p be the probability that one outer block supplies at
least h suffix coordinates. Since P is uniform,

    p = sum_(j=h)^b binom(e,j) binom(N-e,b-j) / binom(N,b).   (2)

As usual, binomial coefficients outside their valid range are zero.
Let Z count the blocks satisfying X_ell>=h, and put mu=L p.

The block indicators have nonpositive pairwise covariance. To prove this,
condition on X_1=j. The second block's suffix occupancy then samples from
N-b positions containing e-j suffix labels. Its upper-tail probability
decreases with j. The indicator 1[j>=h] increases with j. The covariance
of an increasing function and a decreasing function of one variable is
nonpositive. Consequently

    E[Z^2] <= mu + L(L-1)p^2,
    Pr[d_min(C)<=e] >= Pr[Z>0]
                    >= mu/(mu+1-p).                       (3)

The last step is the second-moment inequality. Equation (3) counts a
routing event that necessarily creates a bad code; it does not count
field-valued messages. Every quantity in (2) and (3) is rational.

An independent arithmetic check also computed the exact pair probability

    p_2 = sum_(j=h)^b sum_(l=h)^b
            binom(b,j) binom(b,l)
            (e)_(j+l) (N-e)_(2b-j-l) / (N)_(2b),

where (a)_r is a falling factorial. All checked finite cases satisfy
p_2<=p^2. Using mu^2/(mu+L(L-1)p_2) slightly strengthens (3); the tables
below use the simpler bound (3), except for the stronger dense bound next.

## A strong bound when many blocks fit the suffix

Set delta=e/N. First include each coordinate independently with probability
delta. Under this temporary sampling law, block occupancies are independent
Binomial(b,delta) variables. Let

    p_B = sum_(j=h)^b binom(b,j) delta^j (1-delta)^(b-j).

Conditioning the sample size to equal e recovers the uniform e-subset law
of the permutation. Since e is a mode of Binomial(N,e/N), that conditioning
event has probability at least 1/(N+1). It follows that

    Pr[Z=0] <= min(1, (N+1)(1-p_B)^L).                    (4)

A dyadic consequence uses exact rational arithmetic. Define

    a = floor((10/7)L p_B),
    B = a-ceil(log2(N+1)).

The inequality ln(2)<7/10 gives (1-p_B)^L<=2^-a. Therefore, for B>0,

    Pr[d_min(C)>e] <= Pr[Z=0] <= 2^-B.                    (5)

The exponents reported below use (5), not floating logarithms. For example,
ln(2)<7/10 follows because the fourth-order Taylor lower bound on exp(7/10)
already exceeds two.

For every fixed b,k and every fixed delta>0, p_B is a positive constant.
Equation (4) then tends to zero as N grows. Thus this ensemble has vanishing
relative minimum distance in probability for every fixed local constituent.
Increasing the causal scan width does not change this conclusion.

## Evaluated whole-code failures

A target d_min>=ceil(delta N) forbids weights through
e=ceil(delta N)-1. For the power-of-two half-rate lengths below, this equals
floor(delta N) at both displayed distances.

At rate one half, the 30% target has the following complete failure lower
bounds. They hold over P for every realization of the taps.

| Message length K | Output length N | Forbidden e | Expected bad blocks mu | Failure probability at least |
|---:|---:|---:|---:|---:|
| 65,536 | 131,072 | 39,321 | 210.23585 | 1-2^-282 |
| 262,144 | 524,288 | 157,286 | 841.19547 | 1-2^-1181 |
| 1,048,576 | 2,097,152 | 629,145 | 3364.99922 | 1-2^-4785 |

Even the 10% target cannot have a useful setup-failure guarantee for this
uniform-permutation architecture. These failure lower bounds are rounded
downward from exact rational evaluations of (3).

| K | N | Forbidden e | Failure probability at least | Maximum possible failure-margin bits, approximately |
|---:|---:|---:|---:|---:|
| 65,536 | 131,072 | 13,107 | 0.046189459652 | 4.4363 |
| 262,144 | 524,288 | 52,428 | 0.162482389074 | 2.6217 |
| 1,048,576 | 2,097,152 | 209,715 | 0.437067300185 | 1.1941 |

Here failure-margin bits means -log2 of the setup-failure probability.
The displayed maximum is an upper bound derived from the suffix event;
other bad-code events can only make the margin smaller. This does not
contradict the binary SPIN control or the regional field constructions:
their route is not the full uniform P assumed here.

For GRS[16,12], rate is 3/4 and five suffix coordinates suffice. Natural
message lengths round the requested K upward to a multiple of twelve.
Already at the 10% target, the obstruction gives:

| Requested K | Actual K | N | Forbidden e | Failure probability at least |
|---:|---:|---:|---:|---:|
| 65,536 | 65,544 | 87,392 | 8,739 | 1-2^-115 |
| 262,144 | 262,152 | 349,536 | 34,953 | 1-2^-511 |
| 1,048,576 | 1,048,584 | 1,398,112 | 139,811 | 1-2^-2101 |

These are the natural rate-3/4 lengths, not the half-rate output lengths
from the preceding tables. The bound at 30% is stronger still.

### How large would a uniform-route outer need to be?

The same exact calculation gives a necessary-only length screen. For each
requested K, round up to complete MDS outer blocks, keep the displayed rate
exact, and increase the outer length until (3) no longer rules out a 2^-40
failure target. The first lengths passing this obstruction are:

| Requested K | Half rate, 30% distance | Rate 3/4, 10% distance |
|---:|---:|---:|
| 2^16 | 346 | 320 |
| 2^18 | 362 | 336 |
| 2^20 | 378 | 352 |

These are not sufficient lengths or certified codes. They only stop failing
one necessary test. Even MDS[256,128] at K20 has a setup-failure lower bound
above 5.99898e-8 for 30% distance, limiting its margin to less than 24 bits.
A fixed small outer plus a wider causal inner cannot evade that test.

    python -B research/workstreams/char2_field_spin/theory/uniform_suffix_gate.py --find-minimum
    python -B research/workstreams/char2_field_spin/theory/uniform_suffix_gate.py --find-minimum --rate 3/4 --delta 1/10

## Which construction change addresses the obstruction

Regional routing changes the necessary event. If every outer coordinate
has its own equal-sized region, a short suffix cannot contain arbitrarily
many coordinates from one outer block. For GRS[16,8], a suffix of length
30% of N intersects at most five of the sixteen regions. For GRS[16,12],
a suffix of length 10% of N intersects at most two. These are below the
respective shortening thresholds nine and five.

This observation removes the single-suffix witness for those regional
geometries. It is not a distance theorem: unions of windows and other
shortened subspaces still require complete analysis. The GF256 regional
direction construction is therefore a relevant next candidate; replacing
its coordinate spread with a full uniform permutation would destroy the
property just identified.

Simply closing the scan into a cycle is also insufficient when its state
dimension is small. If an outer subspace of dimension D is confined to an
interval, imposing zero final state adds at most w equations for a w-symbol
state. When D>w, a nonzero word survives with state zero on both sides.
It remains confined even in an invertible cyclic realization. A cyclic
candidate must therefore control shortened dimensions, not only repair
the one-word endpoint example.

The next proof should retain regional coordinate spread and test every
relevant confined subspace before optimizing local coefficient transfer.
The polynomial moment in [HIGH_DISTANCE_GATE.md](HIGH_DISTANCE_GATE.md)
remains valid, but the full-uniform-permutation prototype there is withdrawn.

## Reproduction and audit

[uniform_suffix_gate.py](theory/uniform_suffix_gate.py) computes the exact
hypergeometric lower bound, checks proposed failure targets with rational
arithmetic, and returns the dyadic exponent from (4). Its companion
[tests](theory/test_uniform_suffix_gate.py) exhaust small subset laws,
pair covariances, MDS support witnesses, every small tap assignment, and
natural-length rounding. Six tests pass.

An independent audit recomputed both single-block and two-block probabilities
using the falling-factorial expressions above. It also recomputed the
Bernoulli-conditioning exponents and the natural rate-3/4 lengths. These
calculations agree with the saved evaluator. No native code, remote work,
or encoding benchmark was used.

    python -B research/workstreams/char2_field_spin/theory/uniform_suffix_gate.py --blocks 16 --delta 3/10
    python -B research/workstreams/char2_field_spin/theory/uniform_suffix_gate.py --blocks 16 --delta 1/10
    python -B research/workstreams/char2_field_spin/theory/uniform_suffix_gate.py --blocks 16 --rate 3/4 --delta 1/10
    python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_uniform_suffix_gate.py -v
