# Nonlinear input permutations across RS columns

Status: exact algebraic proposal for a new construction. No native implementation,
regional provider run, or whole-code certificate is claimed. Existing receipts
remain attached to their original construction.

The static H* mixer admits quiet four-row affine planes. A deterministic input
permutation that changes with the RS column prevents one such plane from remaining
quiet throughout the word. The exact checks below establish at least five
nonquiet columns for every four-row group and at least two for every five-row group.

## Construction and quiet columns

Identify the sixteen mixer coordinates with GF(16), using polynomial
`x^4+x+1` (`0x13`) in the binary label basis. This label field is distinct from
the GF(256) field of symbol amplitudes. Define `inv(0)=0` and

\[
P_j(i)=\operatorname{inv}(i+j),\qquad 0\leq i,j<16.
\]

Here label addition is XOR. Column j sends original RS row i to mixer
coordinate P_j(i). The sixteen source rows remain the same RS[16,8] codewords
throughout the word. The scheduled column map is

\[
y_j=D_{2,j}H_*D_{1,j}P_jx_j.
\]

The diagonal entries are independent uniform nonzero GF(256) values sampled once
at setup, as in the ideal model for the existing H* construction. The schedule
itself has no randomness. The inverse table is

```
0,1,9,14,13,11,7,6,15,2,12,5,10,4,3,8
```

Write V for the sixteen-by-five binary affine-evaluation matrix whose row at
label i is `(1,i_0,i_1,i_2,i_3)`. The existing identity
`H*=I+VBV^T`, with invertible B and full-column-rank V, gives

\[
H_*u=u\quad\Longleftrightarrow\quad V^Tu=0.
\]

A column is **quiet** when its pre-H* vector u satisfies this equation. Empty
columns are quiet. Diagonal scaling preserves support, so the support restrictions
below hold for every realization of the nonzero diagonals.

## Four active rows: at most four quiet-capable columns

Fix a four-element original row set A. A nonempty quiet vector supported on
P_j(A) must use all four coordinates: the affine-check code has minimum weight
four. It exists exactly when the four mapped labels sum to zero.

The following polynomial argument bounds the number of such columns. Define

\[
p_A(X)=\prod_{a\in A}(X+a)
      =X^4+e_1X^3+e_2X^2+e_3X+e_4.
\]

For j outside A, the mapped-label sum is

\[
\sum_{a\in A}\frac1{a+j}
=\frac{p_A'(j)}{p_A(j)}
=\frac{e_1j^2+e_3}{p_A(j)}.
\]

For j in A, the zero convention omits that summand. Factoring the simple root
j from p_A gives the quiet condition

\[
e_1j+e_2=0.
\]

The denominator in this last identity is p_A'(j), which is nonzero because A
has distinct elements. If e_1 is nonzero, the derivative has one root in GF(16),
outside A. There is at most one additional quiet column in A. If e_1 is zero,
squarefreeness implies e_3 is nonzero, so no column outside A is quiet-capable.
In this case either all four columns in A qualify (e_2=0), or none do.
Thus every A has at most four quiet-capable columns.

The complete exact enumeration agrees with this argument:

| Quiet-capable columns | Labelled four-row sets |
|---:|---:|
| 0 | 120 |
| 1 | 1,200 |
| 2 | 480 |
| 4 | 20 |

Each active RS row has at least nine nonzero columns. At most four of those
columns can be quiet. Hence every group with exactly four active RS rows has
at least five nonquiet columns. This is a deterministic support statement,
not a probability estimate or a distance bound for the concatenated code.

## Five through seven rows: exact incidence test

Fix an original row set A of size a. For i in A, let b_i count the columns j
for which coordinate P_j(i) can be nonzero in a quiet vector supported on P_j(A).
The exact rank criterion is

\[
\operatorname{rank}V_{P_j(A)}
=\operatorname{rank}V_{P_j(A\setminus\{i\})}.
\]

Deleting a column of V^T preserves its rank precisely when that coordinate is
not a coloop. Binary rank equals rank after extension to GF(256). For each
column, all such coordinates can be nonzero simultaneously: fewer than 256
proper coordinate hyperplanes cannot cover the quiet kernel. This per-column
existence statement does not impose the RS equations across columns.

An all-quiet RS word requires b_i at least nine for every active row. More
generally, any such group needs at least

\[
\max\{0,9-\min_{i\in A}b_i\}
\]

nonquiet columns. Complete enumeration over all A gives:

| a | Distribution of min_i b_i (value: count) | Sets checked | Universal nonquiet-column floor |
|---:|---|---:|---:|
| 5 | 3:960, 4:2208, 5:960, 7:240 | 4,368 | 2 |
| 6 | 9:1920, 10:4560, 11:1200, 12:328 | 8,008 | 0 |
| 7 | 12:720, 13:3520, 14:6960, 15:240 | 11,440 | 0 |

Every five-row set therefore excludes an all-quiet RS word. The incidence test
excludes no six- or seven-row set. Passing this necessary test does not prove
that an all-quiet RS word exists.

## Quiet dimensions across the complete word

For a fixed original row set A of size a, define

\[
d(A)=\sum_{j=0}^{15}\bigl(a-\operatorname{rank}V_{P_j(A)}\bigr).
\]

This is the dimension of the product of the sixteen quiet-column kernels before
imposing the RS constraints. A separate exact enumeration checked all 49,946
labelled row sets of arities four through nine:

| a | Distribution of d(A) (value: count) | Maximum d(A) |
|---:|---|---:|
| 4 | 0:120, 1:1200, 2:480, 4:20 | 4 |
| 5 | 5:1488, 6:1440, 7:960, 8:240, 9:240 | 9 |
| 6 | 16:480, 17:2640, 18:4080, 19:640, 20:120, 21:48 | 21 |
| 7 | 32:7600, 33:3840 | 33 |
| 8 | 48:12390, 49:480 | 49 |
| 9 | 64:11440 | 64 |

In particular, an eight-row set maps to an affine hyperplane in at most one
column. The other columns have affine rank five. Its maximum quiet dimension
therefore falls from the static worst case 64 to 49. For nine rows, every
column has rank five, independently of the schedule.

These counts can support a separate setup-expurgation lemma, but the dimension
helper does not compute such a probability. That lemma must account for the RS
counting measure, setup diagonals, source positions, and the row-set union.
Neither the table nor the incidence test is a whole-code distance certificate.

## Global setup event for too few nonquiet columns

The dimension counts also give a stronger necessary activity condition for
arities five through eight. This condition holds outside one small setup event,
simultaneously for every source group and every nonzero message in that group.

Fix A of size a and let r_j be its quiet-check rank in column j. The whole
column space has q^a elements, and its quiet subspace has q^(a-r_j) elements.
For q=256, the exact number of enlarged arrays with at most h nonquiet columns,
excluding only the global zero array, is

\[
C_{A,h}=q^{d(A)}[Y^{\leq h}]
       \prod_{j=0}^{15}\left(1+(q^{r_j}-1)Y\right)-1.
\]

Here `[Y^{<=h}]` means the sum of coefficients through degree h. The factors
partition each column into quiet and nonquiet vectors. In particular, this
formula does not union overlapping choices of columns allowed to be nonquiet.

For a valid nonzero RS row with prescribed support T, shortened-code counting
gives at most 255^(|T|-8) full-support codewords. The independent D1 entries
assign any prescribed scaled amplitude vector probability 255^(-|T|).
Consequently its counting measure is at most 255^-8 per scaled row array.
The product bound for a fixed A is 255^(-8a). Adding arrays that violate RS
support floors only enlarges the positive count. Global nonzero scalar multiples
have the same event, so one projective divisor 255 applies.

The probability that any of the 8192 source groups has such a message is thus
at most

\[
\frac{8192}{255^{8a+1}}\sum_{|A|=a}C_{A,h}.
\]

For a=4 the ranks are three or four. For a=5 through nine they are four or
five. The saved dimension histogram determines their multiset, which is enough
for this unweighted count. No additional rank enumeration or regional provider
is used.

| Arity a | Excluded nonquiet count h or less | Directed log2 probability upper bound |
|---:|---:|---:|
| 4 | 4 | -77.529371232863 (raw count; actually impossible) |
| 5 | 4 | -69.670579896229 |
| 6 | 3 | -77.743510457400 |
| 7 | 2 | -80.043520941525 |
| 8 | 1 | -61.671478482257 |
| 9 | 0 | -45.106001459030 |

After removing the deterministically impossible four-row event, the union
through arity eight is at most

\[
2^{-30832912669877/500000000000}
=2^{-61.665825339754}.
\]

Its complement enforces nonquiet-column floors
`(9,9,9,5,5,4,3,2)` for arities one through eight. Arity nine can remain
unrestricted. If the zero-nonquiet event for arity nine is also excluded, the
combined log2 ceiling is -45.105986525291. The exact rational probabilities
and matching Arb256/384/512 ceilings are saved separately.

Increasing the floor by one is not free. The same raw method bounds the
arity-six event with at most four nonquiet columns by 2^-36.736410118299,
and the arity-seven event with at most three by 2^-37.926582228798. These
majorants exceed a 2^-45 setup budget; they do not prove the underlying events
are that likely.

The useful proof decomposition is unconditional: add the setup-event charge,
then bound bad codewords satisfying the new per-group floors under the original
setup distribution. Do not condition the operators on successful expurgation
and assume their previous independence properties survive. The setup lemma
applies to arbitrary numbers of active groups, but does not itself prove the
required output distance.

## Proof integration and implementation boundary

For a fixed original A, a column bound Q_A under the old mixer can be used at
column j with label set P_j(A). The whole-word bound must retain the ordered
product

\[
Q_{P_0(A)}Q_{P_1(A)}\cdots Q_{P_{15}(A)}.
\]

It cannot be replaced by a power of an average column matrix. The original A
remains fixed until its whole-word contribution is combined with other row sets.
The existing source-group routing and inner operator do not change.

The independent cached control `hd_inverse_column_control.py` uses the same
schedule exactly. Its saved FLOAT comparison covers r=18, F<=1, and mixed
arities seven through nine only. It neither proves the sparse-row benefit nor
provides a directed receipt for this new construction.

A native implementation could load source row `inv(r) XOR j` when forming mixer
coordinate r in column j. This changes producer indices without an output scatter.
No runtime saving is asserted: compilation, register pressure, and column-loop
unrolling still require correctness and performance checks.

The next bounded proof test should retain a nonquiet-column flag with the floors
above, charging the arity-at-most-eight setup event once for the whole code.
This directly targets the old all-quiet obstruction. No further
uniform-information-cover checker is needed for that test.

## Reproduction and pins

From the theory directory:

```
C:/Python314/python.exe -m unittest hd_bit_syndrome_column_permutations_tests -v
```

Six exact tests passed in 1.249 seconds. The complete five-through-seven
incidence calculation took 0.913 seconds on its first run. The saved receipt
contains a fresh 0.932-second replay, all histograms, witness row sets, the
complete permutation table, and the independent root-schedule equality check.

| File in theory/ | SHA-256 |
|---|---|
| hd_bit_syndrome_column_permutations.py | fecd7c914a221901aad426a7d93f93b603805642377288431a26c38a9cc8163f |
| hd_bit_syndrome_column_permutations_tests.py | a457bbe8b4575e62c3c590ae9f76d1266898ae923dd8ad73161d20f1b15bc676 |
| hd_bit_syndrome_column_permutations_exact.json | 1d7b0e5b86c378e32f28ce48730bee51444ece43b5792ce495a5939ba512d6bd |
| hd_bit_syndrome_column_dimensions.py | 765e27e8da66a6c2946eafb888653815f9f50b70e6042e1ce99de091981bf3fc |
| hd_bit_syndrome_column_dimensions_tests.py | e1c1a604b9de38024acabd1ad03327b849e96eee807587a140339ca13b747164 |
| hd_bit_syndrome_column_dimensions_exact.json | c78f3c8caa755e6d9fcbbfdd7cd01df56341a4191294373ca75e0ea907e16e1c |
| hd_bit_syndrome_nonquiet_expurgation.py | 787b19dbd87665fa9a23d1fd172f956f6348dacfd7e6c39a953042cb59f783b0 |
| hd_bit_syndrome_nonquiet_expurgation_tests.py | dae1996ecba1979da2bfc0739af6835cac85f9ae965955285f571f5b54fa925e |
| hd_bit_syndrome_nonquiet_expurgation_exact.json | e8c8b21914d57817989e9980327965410de34da330b7ff88fe1df16ff89568d2 |

The receipt pins the affine-rank dependency and the independent cached-control
source. Hashes bind these artifacts; they do not replace the mathematical
argument or exact tests.

The separate dimension helper's four tests passed in 0.213 seconds; its full
enumeration took 0.876 seconds. Reproduce those tests with
`C:/Python314/python.exe -m unittest hd_bit_syndrome_column_dimensions_tests -v`.
The original schedule source, incidence tests, and incidence receipt were not
modified when adding this dimension check.

The nonquiet-event helper's six tests passed in 0.058 seconds. They include
complete GF(4) toy-array enumeration, recovery of the earlier zero-nonquiet
count, and comparison with the padded-subset union. Reproduce them with
`C:/Python314/python.exe -m unittest hd_bit_syndrome_nonquiet_expurgation_tests -v`.
