# Rank-three bound: retaining the all-one column constraint

2026-09-24. This note concerns the uniform-shared-coordinate g=4,c=4
distribution in RANDOM_BLOCKS_GFNI.md, at K=2^20, N=2^21, and output weight
at most 209715. The construction, production code, and paper are unchanged.
The initial screens use binary64. Subsequent 192-bit and 384-bit outward
replays verify a 40.643-bit bound for all rank-three messages confined to one
group, including every possible active-group location. The implementation
and coverage of this replay are described below. Ranks one and two retain
their previously verified bounds.

## Which relaxations matter

The preceding worst-shape envelope gives only 30.22 bits for rank-three
messages in one group. Two inexpensive refinements barely change it:

- `basis_lattice.py` uses attainable basis weights rather than all integers,
  and rounds cheap-basis counts to multiples of h!. It passes 73 exhaustive
  small-code inequalities. The improvement at the dominant BCH supports
  is generally below one bit.
- `random_group_backward.py` maximizes each complete shape's action on the
  continuation vector instead of maximizing its matrix entries separately.
  The rank-three diagnostic improves only to 30.235 bits.

These tests point away from matrix-entry splicing as the main problem.
The envelope permits arbitrarily many all-one columns. Higher-rank BCH
tuples cannot have that pattern everywhere on their union support.

## Penalizing all-one columns

Fix a rank-h four-row outer tuple, before the independent lane permutations.
Let u be its union-support size and m its
number of columns equal to (1,1,1,1). If this pattern lies in the column
space, the linear combinations annihilating it form an (h-1)-dimensional
subcode. That subcode has support size u-m: every other nonzero column
is detected by some annihilating combination. If the pattern is outside
the column space, then m=0.
Lane permutations preserve column weights, so they preserve this count.

Consequently m <= u-d_(h-1), where d_(h-1) is a lower generalized weight
of the BCH code. The valid values used here are d_1=38, d_2=57, d_3=67.
`random_group_penalty.py` checks this inequality by exhaustive enumeration
on small codes.

For 0<rho<=1, weight each local shape transfer by rho to its number of
all-one columns before taking the worst-shape envelope. The resulting
product bounds the output moment multiplied by rho^m. Since m<=M implies
rho^m>=rho^M, multiplying by rho^(-M) gives a valid upper output moment.
Support size and first occupied macroregion remain conditioned exactly.

A coarse grid improves the rank-three screen to 35.39 bits. Refining the
tilts and rho values improves it to 38.06 bits. The remaining gap suggests
coupling the column constraint to the BCH multiplicity bound, rather than
using a worst-case M for every tuple with the same union support.

## Count the annihilating subspace as well

Let U be the three-dimensional span of the four row words. When the
all-one pattern belongs to its column space, let H be the two-dimensional
subcode annihilating that pattern. Write u=|supp(U)| and v=|supp(H)|.
Then the all-one pattern occurs exactly u-v times.

For any binary r-dimensional subcode of support size w, the sum of its
nonzero word weights is 2^(r-1)w. Thus the four words in U outside H have
total weight 4u-2v. At least one has weight at most floor(u-v/2).

Let H_v be an upper count of two-dimensional BCH subcodes with support
exactly v. The weight-triple method in RANK_TWO.md supplies such shell
caps by summing its per-type subspace caps. Let A_<=a bound the number of
nonzero BCH words of weight at most a. The number of pairs (U,H) with
|supp(U)|<=u and |supp(H)|=v is at most

    H_v A_<=floor(u-v/2).

Indeed, each pair has a representative word outside H below this weight
threshold; H and any such word determine U. Counting all candidate words
can only overcount the pairs.

For fixed U,H, exactly 168 spanning four-row assignments map the nonzero
functional annihilating H to the all-one pattern: after fixing that image,
the remaining two independent images have 14 and 12 choices. Therefore

    F_v(u) = min(T_3(u), 168 H_v A_<=floor(u-v/2))

is an upper CDF for the corresponding tuples. Set F_v(u)=0 for u<max(67,v).
Here T_3 is the existing upper CDF for all rank-three tuples.

For each v, the penalized transfer calculation now uses the exact count
m=u-v. Taking the best bound over the fixed finite tilt/penalty grid gives
a failure upper bound p_v(u). Over a chosen sparse interval u<=L, replace
it by a decreasing majorant in u. Summation by parts combines it with
F_v using only nonnegative coefficients, then sums over v. This avoids
interpreting CDF differences as exact shell counts.

There is a second class. Eight of the fifteen three-dimensional subspaces
of F_2^4 exclude the all-one vector. Their spanning-row assignments form
exactly 8/15 of every support-size class. For them, discard all shapes
containing a weight-four column and use floor(8 T_3(u)/15) as the CDF cap.
This class is bounded over the full range u<=256, not only the sparse part.

For tuples in the other seven image subspaces with u>L, use the earlier
unrestricted-shape moment bound and floor(7 T_3(u)/15). The restricted
class has CDF zero through L; beyond L its CDF is at most this cap.
These two tail classes and the sparse flags cover all rank-three tuples.
Finally union-bound over the 2048 possible active groups.

`rank_three_flags.py` tests the coset-weight identity, the 168 factor, and
the 8/15 partition against exhaustive small codes. All outer multiplicity
calculations use integers. The moment calculations still use binary64 and
an explicit tiny-value floor, so their displayed margins are diagnostics.

## Numerical screen

With L=120, sixteen fixed tilts and eight positive penalty values, plus
the no-all-one case, the final binary64 screen gives:

| Contribution, union over all active-group locations | Negative log2 upper-bound diagnostic |
|---|---:|
| Sparse flags and the full no-all-one class | 43.6795 |
| Dense remainder in the seven other image subspaces | 40.8305 |
| All rank-three tuples in one group | **40.6430** |

The no-all-one class alone gives 55.63 bits. The main sparse contributions
come from annihilating-subspace support sizes around 57 through 66.
The dense remainder dominates the combined bound. This screen identified
a parameter grid for the outward replay; it is not itself a certificate.

## Outward verification and the remaining failure budget

`rank_three_verify.py` reconstructs the local orbit census and uses Arb
throughout the moment calculation. The tilt grid and the split L=120 are
fixed in the verifier. Penalty values are exact dyadic rationals. No
binary64 transfers or tiny-value floor enter this replay.

The verifier first minimizes the penalized conditional moment over tilts.
For each annihilating-subspace support v, it multiplies by rho^(-(u-v)),
minimizes over positive penalties, and caps the probability at one. All
reused bounds are rounded upward. Decreasing majorants are exact dyadic
values, so their nonnegative differences can safely weight the outer CDF
caps. The no-all-one class uses rho=0; the dense remainder uses rho=1.

Independent runs at 192-bit and 384-bit precision agree and verify

    rank-three failure upper < 5.824288e-13 < 2^-40,
    margin approximately 40.6429836351 bits.

The sparse flags contribute less than 7.096705e-14. The full no-all-one
class contributes less than 1.790025e-17, and the dense remainder contributes
less than 5.114439e-13. These contributions are union-bound terms, not
measured failure probabilities.

The run authenticates 163 BCH dependency files, verifies 67 exact rational
shortening witnesses, and rechecks the rank-two shell-count factors. It
also checks the attainable-basis-weight refinement, small-code hyperplane
counts, conditional transfer identities, and sparse/tail partition.
The final comparison is against 2^-40 directly, not a rounded margin.
Authentication does not replay every historical BCH proof.

Combining conservative rounded upper bounds from the three completed
one-group rank classes gives

    rank 1: 2.410310e-14,
    rank 2: 8.519283e-17,
    rank 3: 5.824288e-13,
    sum:    6.0661709283e-13.

The sum has margin greater than 40.584 bits. It leaves approximately
3.028776e-13 below the overall 2^-40 budget. Rank four and messages spanning
multiple groups must fit within this allowance together, equivalent to
about 41.586 bits for their combined contribution. Giving each missing
class a separate 40-bit bound would not be sufficient.

The full-distance certificate and the 2x performance target remain open.
In addition to the remaining proof classes, the next engineering step is
to revisit the kernel's measured cost. The subsequent code-equivalent
implementation in VECTOR_ROUTE.md reaches about 1.76--1.78x, not 2x.
Production defaults and the paper remain unchanged.

Reproduce the screens from the repository root:

    python -B research/workstreams/permutation_locality/basis_lattice.py
    python -B research/workstreams/permutation_locality/random_group_backward.py
    python -B research/workstreams/permutation_locality/random_group_penalty.py --fine
    python -B research/workstreams/permutation_locality/rank_three_flags.py --maximum 120
    python -B research/workstreams/permutation_locality/rank_three_verify.py --precision 192
    python -B research/workstreams/permutation_locality/rank_three_verify.py --precision 384

The scripts authenticate the BCH inputs and reconstruct the same local maps
as the preceding calculations. They write no experiment data files.
