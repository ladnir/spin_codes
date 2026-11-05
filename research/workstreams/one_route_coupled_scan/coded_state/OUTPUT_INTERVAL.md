# Count outputs before bounding the inverse randomness

The dense mixed example in [MIXED_RANK.md](MIXED_RANK.md) defeats the forward
message count, even with perfect quiet rank. An output-side count closes
that example and an entire surrounding class without changing the encoder.
The whole-code near-half certificate remains open.

## Fix the trajectory pattern, not the diagonal setup

Use the column-scaled post-diagonal inner over F=GF(2^128):

    Y_t=X_t+D_t A Q_(t-1),     Q_t=Q_(t-1)+B R_t Y_t.

The maps, outer, route, and independent setup distributions are as in
[QUIET_BRANCH.md](QUIET_BRANCH.md). Fix the route and every nonzero R entry.
Only the independent uniform full-field entries of D remain random.
The proof below holds for every such fixed route and R, not just typical ones.

Specify c disjoint episodes occupying v batches, including each activation
batch. Let V denote this batch union. Let J be a disjoint set of j nonzero
quiet batches. All other output batches are zero and have zero incoming
state. Each episode starts with zero incoming state and then has nonzero
incoming state until its reset or the end. Let A_live be the coordinate
union obtained by removing the first batch from every episode. It contains
128(v-c) coordinates. The final state is unrestricted.

For a fixed output y realizing this pattern, its state sequence is determined
by y and R, independently of D. Let S(y) contain the nonzero coordinates of
A Q_(t-1) in each batch. The MDS property of A gives

    S(y) subset A_live,           |S(y)| >=65(v-c).

Over D, the inverse input is uniform on y+F^S(y). Write d_C(T) for the
dimension of the routed outer shortened to coordinate set T. Its intersection
with that affine space has either zero or q^d_C(S(y)) elements. Consequently

    Pr_D[inverse(y) in C]
       <=q^(d_C(S(y))-|S(y)|)
       <=q^(d_C(A_live)-65(v-c)).                         (1)

This upper bound removes the need to enumerate internal feedback-hole masks.
It uses containment of their coordinate set, not independent state bits.

## Quiet outputs contribute rank to the output count

For t in J, the equation B R_t y_t=0 holds, and y_t is nonzero. Every
64 columns of B R_t are independent. Thus wt(y_t)>=65; on any fixed
possible support of size w_t, the kernel has dimension w_t-64.

Put u=128(v+j) and let W be the total output-weight threshold. On a fixed
global output support of weight w<=W, the j quiet equations constrain
disjoint variables and have total rank 64j. There are at most q^(w-64j)
eligible output assignments. Other trajectory constraints can only reduce
that count. Summing over at most 2^u supports gives

    #eligible outputs <=2^u q^(W-64j).

Combining this count with (1) gives the conditional all-message bound

    min(1, 2^u q^(W-64j+d_C(A_live)-65(v-c))).              (2)

This is a union over fixed output candidates, so it covers messages chosen
after setup. No uniformity conditioned on future output events is assumed.
For 1<tau<=q, fractional averaging of (2) over the route gives

    tau^(u/128+W-64j-65(v-c)) E_route[tau^d_C(A_live)].      (3)

The existing interval moment bounds apply to A_live, a union of at most c
intervals. Quiet placements contribute a factor at most binom(n-v,j),
and episode placements use the retained terminal/nonterminal interval counts.
Equation (3) is a general interface, not yet a complete sum over all patterns.

## Complete dense one-episode class

At K=2^20, n=16384 and each region contains 64 batches. Consider every
trajectory with exactly one episode, of length 1<=v<=8192, and at most
32 zero-output batches outside that episode. All other outside batches are
nonzero quiet batches. Write z for the number of outside zero batches;
then j=n-v-z. Retain only feasible pairs satisfying 65j<=W.

The live interval has L=v-1 batches. Regardless of its starting offset,
it touches at most

    h(L)=floor((L+126)/64),     h(0)=0,

regions. Each outer block contributes at most one coordinate per region.
For every regional route realization,

    d_C(A_live) <=8192 max(0,h(L)-128).

The log2 bound from (2), for each specified interval and outside-zero set,
is therefore at most

    128(v+j)+128[W-64j+8192 max(0,h(v-1)-128)-65(v-1)].    (4)

At W=floor(.491N)=1029701, exact integer enumeration covers all 252978
feasible (v,z) pairs. The largest bound is -44672, at
(v,j,z)=(8131,8221,32). There are at most n^2 choices of interval,
33 choices of z, and n^32 choices of its outside-zero positions. Their
combined count is less than 2^482. Thus this entire class has failure
probability at most 2^-44190 at >49.1% distance, uniformly in route and R.
The large margin belongs only to this class, not to the whole code.

For the former obstruction, take the episode in the final half and all
8192 earlier batches quiet and nonzero. Its live coordinates occupy at most
128 regions, so d_C(A_live)=0. Equation (2) gives log2 failure at most
-1359104, rather than the positive exponent from the forward argument.

## Rank from prescribed state holes

Longer or more fragmented live trajectories are not covered above. To avoid
replacing every live batch by its worst 65-coordinate feedback support,
retain constraints that its holes impose on the output variables.

Fix output support sets W_t. To prescribe h holes of the nonzero state
expansion A Q_t, let H_t be those h rows of A, where 0<=h<=63. The new
variables y_t[W_t] have coefficient matrix

    A[H_t,:] B[:,W_t] diag(R_t[W_t]).

Its rank is at least max(0,h+min(64,|W_t|)-64), by the full row rank of
A[H_t,:], the MDS columns of B, and the rank inequality for a matrix product.
To prescribe Q_t=0 instead, its rank on these new variables is
min(64,|W_t|). Earlier state constraints contain no y_t variables. A
block-triangular minor therefore adds these lower bounds across steps.
If their sum is rho, at most q^(w-rho) output assignments realize the
specified supports and state constraints. Do not constrain Q_n: it is free.

Zero-output steps satisfy Q_t=Q_(t-1) exactly. Their expanded-state hole
sets are identical; re-counting them independently would be unnecessary
overcounting. Any refined enumeration must preserve that dependency.

The product-rank inequality is sometimes strict and sometimes attained.
In the GF16, b=8, s=4 analogue, all h=w=3 submatrices of AB have ranks
2 or 3: 552 of rank 2 and 2584 of rank 3. Thus replacing every minor by
the generic full rank is false, but treating all of them as worst-rank
also loses structure. The next step is to retain the geometry and count
of these exceptional minors, then combine it with (3) or a finer inverse
support bound. No new randomizer is selected by this observation.

## Checks and reproduction

[output_interval.py](output_interval.py) performs the full integer class
check and the explicit old-obstruction calculation. No floating-point
arithmetic is needed. [test_output_interval.py](test_output_interval.py)
also checks:

- 512 exact interval/region-count comparisons.
- All 65535 nonzero length-eight GF4 outputs, grouped into 129 exact
  trajectory/weight classes, against the quiet-constrained output count.
- 23808 exact submatrix ranks in the GF16 map; the product-rank inequality
  is strict in 7604 cases.
- 512 joint state-constraint rank fixtures with a free terminal state,
  including zero-output state persistence.

The tests validate the stated local interfaces; the integer computation
covers precisely the dense class above. Run serially with no output files:

```text
python -B research/workstreams/one_route_coupled_scan/coded_state/test_output_interval.py
```

The earlier sparse-source and entirely quiet certificates remain unchanged.
Other mixed trajectories still require a complete all-message proof.
