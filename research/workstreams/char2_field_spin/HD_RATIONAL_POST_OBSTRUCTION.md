# A run cut rejects the regional rational scalar postscan

The proposed rational scalar postscan does not repair the regional
construction. At K=2^20 and rate one half, the construction below fails
40% minimum distance with probability greater than 19/20 over its regional
permutations. The conclusion holds for every valid coefficient assignment.
It therefore rules out a 40-bit setup-failure guarantee over GF(2^16).

This is an actual code obstruction, not a failed first-moment upper bound.
The postscan has a dense inverse, but its inverse tails have rank one.
Long safe runs allow that one-dimensional tail cost to be amortized.

## Construction and probability space

Fix a field F with at least three elements. Let b=64, k=32, R=32768,
K=kR, and N=bR. Apply an injective [b,k] constituent independently to
each of R message rows. Store coordinate j of row i at position jR+i.
The argument applies to the proposed GRS constituent, but does not require MDS.

Let B be this striped outer map. Its parity-check matrix H0 has R row
groups of k checks each. Column jR+i is supported only in check group i.

Apply two input-weighted prefix scans in the same order. Write their
composition as Q. Each scan has zero initial state and nonzero input
weights. Thus Q is invertible, and Q^-1 has lower bandwidth two:
column t is supported on positions t,t+1,t+2, where those positions exist.

Independently in each of the b regions of length R, sample a uniform
permutation. Let P denote the resulting block-diagonal permutation map.
After P, apply the rational postscan J:

    state[-1] = 0,
    state[i]  = state[i-1] + d[i]*v[i],
    y[i]      = state[i-1] + c[i]*v[i].

Require c[i],d[i] to be nonzero and distinct. All scans carry state across
region boundaries; they do not restart at those boundaries. Output all N
entries of y. The final code is the image of J P Q B, with parity check

    H = H0 Q^-1 P^-1 J^-1 = M J^-1,
    M := H0 Q^-1 P^-1.

The probability below is over the independent uniform regional permutations
only. The bound holds for arbitrary valid coefficients, even if they are
chosen after observing P. In particular it holds under the proposed
independent uniform nonzero coefficient sampling.

## Safe columns and one tail direction per interval

Choose an integer a and let U be the check-coordinate space of the first
a outer row groups. Its dimension is h=ka. Before routing, each region
has A=a-2 designated safe columns: its positions 0 through a-3.
The bandwidth-two inverse precode maps each such column only to outer row
indices below a. Therefore its image under H0 lies in U.

After routing, the safe positions in each region form a uniform A-subset.
The subsets in different regions are independent. Every safe column of M
lies in U, for every coefficient assignment.

The inverse rational scan uses one state. Solving its equations gives

    v[i]     = c[i]^-1 * (y[i] - state[i-1]),
    state[i] = (1-d[i]/c[i])*state[i-1] + (d[i]/c[i])*y[i].

Over GF(2^16), subtraction is addition, giving the proposed XOR recurrence.

Consider input to J^-1 supported on an interval I=[l,r]. Before l its
response is zero. After r the response depends only on the single scalar
state at r. Consequently, the part of (J^-1)_I outside I has rank at
most one. This remains true even when the interval crosses a region boundary.

If every position in I is safe, multiplication by M maps the within-I
part into U. Hence the columns H_I lie in U plus one additional direction.
For a union S of rho disjoint safe intervals,

    rank(H_S) <= h+rho.                                      (1)

The same U is shared by all intervals, so its dimension is charged once.
There is no independent-tail assumption. The tail directions may overlap,
which can only reduce the rank.

## Select safe runs of length at least two

Within each region, retain all maximal safe runs of length at least two.
Treat region boundaries as run boundaries. Let S be their union, let rho
be their number, and define

    T := sum over retained runs I of (|I|-1).

Thus |S|=T+rho. Also T is exactly the number of adjacent safe-safe pairs
within regions. A retained run starts either with a boundary pair 11 or
with the pattern 011. Singleton safe runs contribute to neither statistic.

Let e=ceil(2N/5)-1=838860 be the largest forbidden output weight. If

    T > h,                 rho <= e-h-1,                    (2)

then |S|>h+rho. Equation (1) implies a dependent set of at most h+rho+1
columns of H_S. Their dependence gives a nonzero codeword with weight at
most h+rho+1<=e. Conditions (2) therefore imply failure of the 40% target.

## Exact means and a conservative variance bound

For a uniform A-subset of a region of length R, every adjacent pair is
safe with probability A(A-1)/(R(R-1)). There are R-1 pairs. Counting the
boundary 11 event and the R-2 interior 011 events gives

    E[T]   = b*A*(A-1)/R,
    E[rho] = b*(R-A+1)*A*(A-1)/(R*(R-1)).                   (3)

Both statistics change by at most two when a safe and an unsafe position
are exchanged within a region. For T, deleting a safe bit removes at most
two adjacent pairs, and inserting it adds at most two. The net change
therefore lies between -2 and 2.

For rho, changing one bit changes the number of non-singleton runs by at
most one. An inserted bit can extend a run, promote a singleton, join two
runs, or create a singleton. Each case changes rho by at most one.
Deletion is the reverse operation, so an exchange changes rho by at most two.

Reveal the safe bits successively, region by region. Given a revealed
prefix, couple the two possible next-bit conditional laws by exchanging
that bit with one appropriate unrevealed bit. To couple a next zero with
a next one, choose uniformly one of the remaining safe positions and
exchange it with the next position. The resulting remainder is uniform
with one fewer safe position, as required.

Under this coupling, conditional expectations of either statistic differ
by at most two. Its Doob martingale therefore has increments of absolute
value at most two. Orthogonality of martingale increments, over at most N
reveals, gives the conservative bounds

    Var(T) <= 4N,             Var(rho) <= 4N.                (4)

This argument samples without replacement. It does not replace a regional
permutation by independent Bernoulli bits. Deterministic reveals contribute
zero variance, but counting them in N preserves the bound.

## K20 rejection with rational arithmetic

Set a=floor(13R/25)=17039 and A=17037. Then

    h = 545248,
    e-h-1 = 293611,
    E[T]   = 18575509248/32768,
    E[rho] = 292229911489536/1073709056.

The two positive expectation gaps are approximately 21631.55 and 21442.38.
Chebyshev's inequality and a union bound applied to the complements of
(2) give

    Pr[failure] >= 1 - 4N/(E[T]-h)^2
                      - 4N/(e-h-1-E[rho])^2
                > 19/20.                                  (5)

The program evaluates the displayed lower bound as an exact Fraction.
Its approximate value is 0.9638277172631955. The comparison with 19/20 is
exact; the decimal is not used to establish the claim.

For 128-bit records processed as eight independent GF(2^16) lanes, put the
nonzero scalar codeword in one lane and set the others to zero. Its record
weight equals its scalar support, so the same obstruction applies to the
requested record-distance target.

The same argument shows why more fixed-order local prescans do not solve
this problem. With inverse precode bandwidth L, replace A=a-2 by a-L.
For fixed L and a/R tending to alpha,

    E[T]/N -> alpha^2,
    E[rho]/N -> alpha^2*(1-alpha),
    h/N -> alpha/2.

For alpha just above 1/2, the first inequality in (2) has positive linear
slack, and the resulting distance upper bound approaches 3/8. The O(N)
variance bound makes failure tend to one for every fixed target above
3/8. This is an asymptotic probabilistic ceiling, not a deterministic
finite bound of exactly 3N/8.

## Verification and next decision

`theory/hd_rational_cut.py` contains the exact means and finite bound.
`theory/hd_rational_cut_tests.py` independently checks:

- The run identities and exchange Lipschitz bounds on every bit string through length nine.
- Exact means by enumerating every fixed-size subset at small lengths.
- Rank-one inverse tails and selected-run ranks over GF4, including all length-three coefficient setups and support strings.
- Separated runs with arbitrary safe/unsafe check columns, plus an actual two-prescan matrix and regional permutation.
- The K20 strict integer thresholds and exact comparison in (5).

These finite tests support the written argument; they do not replace its
proof for arbitrary lengths and fields. No native implementation or
benchmark is needed to reject this topology.

All seven tests pass with:

```text
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p hd_rational_cut_tests.py -v
```

The next design must change the rank-one interval structure, the regional
safe-prefix geometry, or the bounded-band inverse precode. A wider
full-Cauchy state is a distinct candidate; a denser scalar inverse alone
is insufficient. Its proof and cost remain separate work.
