# One-route BAA: two scalar accumulators cannot reach 40%

For the tested local outer lengths 16, 32, 64, 128, and 256, two weighted
scalar accumulators have a deterministic distance ceiling below 40% at K20
and rate 1/2. This holds for every field, coefficient assignment, and
intervening permutation. Changing the first outer layout from contiguous
blocks to coordinate-major order does not remove this ceiling.

The minimal next topology considered here uses two consecutive weighted
accumulators before the route and one after it. Its inverse has group degree
at most six, which permits 40% but excludes 45%. This necessary condition
is not a lower-distance guarantee or a timing claim.

## Model and the finite obstruction

Let K=2^20 and N=2K. Repeat a full-rank linear [b,b/2] outer independently
across g=N/b rows. The outer may be MDS, but the following obstruction does
not require that property. Fix an arbitrary order of its N coordinates.
Let B denote this ordered outer map.

A weighted scalar accumulator maps u to y with zero initial state and
recurrence

    y_i = a_i y_(i-1) + d_i u_i,

where every a_i and d_i is nonzero. Its final state is discarded, but all N
outputs are retained. The accumulator is invertible; its inverse D has at
most two nonzeros in each column. Input-only random weighting, used for fast
prefix implementation, is the special case a_i=1.

For two such accumulators and a coordinate permutation P, the encoder is

    E = A_2 P A_1 B.

The source dimension is K. A 128-bit application record may contain sixteen
GF256 lanes or eight GF65536 lanes sharing this scalar map. Distance counts
nonzero records/field symbols, not bits. The obstruction holds already for
one nonzero lane.

Choose a parity-check matrix H_0 of B. It has g disjoint row groups, each
containing rho=b/2 checks; every column touches only one group. The arbitrary
outer ordering changes column locations but not this property. A parity
check of E is

    H = H_0 D_1 P^(-1) D_2.

Each final column reaches at most two middle coordinates through D_2 and
at most two outer coordinates per middle coordinate through D_1. Hence it
touches at most Delta=4 check groups. Nonzero diagonal scalings do not enlarge
these supports. No independence or randomness assumption is used.

The exact group-count lemma is in
[width-coefficient-tradeoff.tex](../../../../gen-BAA/field-universality/paper/width-coefficient-tradeoff.tex),
section `sec:group-ceiling`, lemma `lem:check-group-ceiling`. In this checkout
the source is `C:/Users/peter/repo/gen-BAA/field-universality/paper/width-coefficient-tradeoff.tex`.
Its finite statement is:

> Suppose N=2*rho*g and each column of a parity-check matrix touches at most
> Delta groups of rho rows. If Delta<=m<=g and
> `2*C(m-1,Delta-1) > C(g-1,Delta-1)`, then its kernel has distance at most
> `rho*m+1`.

To prove it, select m check groups uniformly. A column touching j<=Delta
groups is confined to the selected checks with probability
`C(m,j)/C(g,j) >= C(m,Delta)/C(g,Delta)`. The displayed inequality makes
the expected number of confined columns exceed rho*m. Some selection
therefore confines at least rho*m+1 columns to only rho*m rows. Their linear
dependence gives the claimed nonzero codeword. The proof works over any field.

Binary search on the exact integer inequality gives:

| Outer | g | Minimal m | Deterministic distance upper bound |
|---|---:|---:|---:|
| [16,8] | 131072 | 104033 | 832265 |
| [32,16] | 65536 | 52017 | 832273 |
| [64,32] | 32768 | 26009 | 832289 |
| [128,64] | 16384 | 13005 | 832321 |
| [256,128] | 8192 | 6503 | 832385 |

The 40% target requires distance at least 838861. Thus every listed
two-scalar-accumulator encoder fails, not merely a fraction of setups.
As g grows, the ceiling approaches `2^(-1/3)/2 = 0.396850263...`.
Larger coefficient alphabets, a stronger local MDS implementation, and a
different single permutation cannot change the inverse group degree.

## A second obstruction that explains the layout dependence

The following calculation preceded the deterministic degree check. It remains
useful for distinguishing topologies and for future wider scans.

Put the outer in coordinate-major order. Each coordinate region contains
the g source rows in their common row order. Select a adjacent source rows
and set all other source rows to zero. Their source dimension is ak. Force
the first scalar accumulator to zero at the exit of each of the b resulting
bands. This costs at most b homogeneous equations. The surviving space has
dimension at least ak-b.

Each band exit is itself a known zero output. The intermediate subspace is
therefore confined to a fixed container of size C=b(a-1), not ba. If a final
suffix of e coordinates receives T of these container positions under a
uniform global P, shortening costs at most C-T further dimensions. A
nonzero middle word confined to that suffix exists whenever

    T > (b-k)a.

The final causal accumulator preserves a suffix, so this is a bad-code
event for every choice of accumulator coefficients. Disjoint adjacent row
groups give disjoint containers. If p is their common hypergeometric event
probability and M=floor(g/a), their indicators have nonpositive pairwise
covariance under uniform subset sampling. The second-moment bound gives

    Pr[at least one bad container] >= M*p / (1+(M-1)*p).

The helper uses eight exact positive terms from the hypergeometric tail as
a lower bound on p. It selects a through a bounded floating search, then
evaluates the reported lower bound as a rational number.

| Outer | a | C | Required hits | Setup-failure lower bound | Reliability at most |
|---|---:|---:|---:|---:|---:|
| [16,8] | 5 | 64 | 41 | 0.694118 | 0.527 bits |
| [32,16] | 6 | 160 | 97 | 0.001232287 | 9.664 bits |
| [64,32] | 6 | 320 | 193 | 1.03382e-9 | 29.849 bits |
| [128,64] | 6 | 640 | 385 | 1.96498e-21 | 68.786 bits |

The last row merely survives this particular 40-bit obstruction; the group
ceiling above still rejects its 40% distance. Distant repetition `[x;x]`
also has short containers: three adjacent source symbols leave at least
one dimension supported on four middle positions after the two band exits
are zeroed. A distant repeat does not by itself supply an effective first
random interleaver.

For a genuine order-r recurrence whose state consists of its last r outputs,
zeroing every band exit removes its final r output coordinates. The analogous
container has size b(a-r), dimension at least ak-br, and the same hit
threshold `(b-k)a+1`. Independent unmixed scalar channels do not receive this
benefit: one can select source rows from a single channel and repeat the
scalar attack.

These hypergeometric probabilities require one uniform permutation of all
N coordinates. Independent permutations within coordinate regions have a
different law. The numerical table must not be reused for that routing law.

## Exact q-ary accumulator enumeration

Write q for the field cardinality, so the nonzero multiplier alphabet has
q-1 elements. For a fixed nonzero source word, independent input multipliers
make each live increment uniform on the nonzero field elements. The same
support transition law holds with fresh nonzero feedback multipliers.

If an incoming state and an input are both nonzero, the new state is zero
with probability p=1/(q-1); otherwise it remains nonzero with probability
r=(q-2)/(q-1). Zero plus zero remains zero, and exactly one nonzero summand
gives a nonzero state. Thus the output-weight-marked matrices are

    T_0(y) = [[1, 0], [0, y]],
    T_1(y) = [[0, y], [p, r*y]].

Let P_N,q(w,h) be the output-weight probability for a uniform w-subset of
input positions, with independent nonzero live increments. The initial
state is zero and the terminal state is free. For w,h>0,

    P_N,q(w,h) = 1/C(N,w) * sum_c C(N-h,c) * p^c *
      [ C(h-1,c-1) C(h-c,w-2c) r^(w-2c)
       +C(h-1,c) C(h-c-1,w-2c-1) r^(w-2c-1) ].

Invalid binomial coefficients are zero; omit terms with negative exponents.
The first term has c closed nonzero
runs; the second has c closed runs and one open final run. The number of
q-ary input words in this input/output weight class is
`C(N,w)*(q-1)^w*P_N,q(w,h)`. The zero-input case has output weight zero
with probability one. These identities are exact, including q=2.

For a fixed structured first scan, let G_h be the expected number of source
messages producing middle weight h. A uniform intervening permutation and
fresh final coefficients give the valid whole-code first-moment inequality

    Pr[d<=e] <= 1/(q-1) * sum_(h=1..min(N,2e)) G_h
                                      * sum_(j=1..e) P_N,q(h,j).

The divisor counts projective message lines. Both accumulators and the outer
are injective, so every bad line has exactly q-1 nonzero representatives.
This formula has only one permutation normalization. Inserting an additional
`1/C(N,w)` before the first scan would incorrectly assume another interleaver.

For `[x;x]`, an exact four-state transfer retains the common source support
between the two halves. It is

    M(y) = T_0(y) tensor T_0(y) + (q-1) T_1(y) tensor T_1(y).

The moment is the sum of entries `(M(y)^K)[(0,a),(a,b)]` over a,b in {0,1}.
This boundary gluing preserves the first-half terminal state as the second
half's initial state.

For a coordinate-major MDS outer, the analogous full transfer has 2^b states,
but permutation symmetry gives a positive transfer on b+1 Hamming weights.
Let A_w count outer codewords on one prescribed support of size w. For an
incoming support of weight u, choose a births, c returns, and d live-to-live
active inputs. The outgoing weight is v=u-c+a. The row-transfer entry is

    M_y(u,v) = y^v * sum_(a,c,d: v=u-c+a)
      C(b-u,a) C(u,c) C(u-c,d) A_(a+c+d) p^c r^d.

The MDS counts A_w are exact. If a positive vector z satisfies M_y z<=lambda*z,
the complete region-boundary gluing is bounded by
`2^b * max(z)/min(z) * lambda^g`. This gives a finite positive upper bound
without replacing the structured first scan by a random interleaver.

The new floating screen applies this bound and the last accumulator's
generating function. For [128,64], GF65536, intermediate densities .4 and .6
give margins approximately -4.781 million and -3.321 million bits. These are
failed upper bounds, not further counterexamples. The first-scan mark lies
near q^(-1/2); localized source subspaces can make counting messages much
looser than counting support constraints and retained dimensions.

## Minimal next topology and its precise open questions

The intended next candidate is

    GRS[64,32], coordinate-major
      -> two independent weighted accumulators in the same order
      -> independent within-coordinate-region permutations
      -> one weighted accumulator.

The two consecutive inverse bidiagonal maps have bandwidth two and at most
three nonzeros per column. This is a genuine coupled order-two map, not two
independent scalar channels. With the final inverse bidiagonal map, each
output check column touches at most Delta=3*2=6 outer groups. The group lemma
still applies to regional permutations, because it needs no routing law.
At K20 with [64,32], its distance upper bound is 912865, approximately
43.5288%. It permits 40% but excludes 45% for every setup.

Three consecutive scans before the route and one after it give Delta=8.
Their exact [64,32] ceiling is 949761, approximately 45.2882%, leaving only
a narrow necessary-condition margin for a 45% target. Two scans on each
side give Delta=9 and ceiling 961569, approximately 45.8512%. These ceilings
do not establish achievable distance.

The double-prefix support law already has an exact four-state local transfer.
Its state records whether each prefix value is nonzero. First update the
first prefix with the current input, then update the second prefix with that
new value, using independent fresh multipliers. Mark only the second output.
The new helper implements this transfer and verifies that it is coupled.

Before native adoption, the next proof work is:

- derive regional, rather than global, support placement after this coupled
  order-two pre-accumulator;
- repeat the dimension-window test with two pre-state constraints per band
  and the final accumulator's exit constraints;
- compose the four-state transfer with the outer support/rank count in one
  common quotient, or justify a stronger first moment over all source classes;
- cover sparse classes, the middle range, and the terminal boundary with a
  finite union bound before directed replay.

The current exact scalar IOWE cannot be used as the first-stage law of the
double prefix. Likewise, the global hypergeometric band probabilities above
cannot be relabeled as regional probabilities. These are the immediate
proof obligations for the next candidate; no long search was started here.

## Artifacts and verification

Only new `hd_baa*` theory files and this note were created. Nine tests cover
small-field exhaustive IOWEs, large-field normalization, distant-repeat
boundary gluing, the coupled double-prefix law, exact tiny MDS row transfers,
finite group ceilings, band lower bounds, and a small exact second-moment
check. The finite group inequalities use integers, and band lower bounds use
rational arithmetic. Moment searches are floating and cover selected classes
only. No build, benchmark, remote command, or complete certificate was run.

```text
python -B research/workstreams/char2_field_spin/theory/hd_baa_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_baa_screen.py --ceilings
python -B research/workstreams/char2_field_spin/theory/hd_baa_screen.py
python -B research/workstreams/char2_field_spin/theory/hd_baa_moment.py --density .4 .6
```
