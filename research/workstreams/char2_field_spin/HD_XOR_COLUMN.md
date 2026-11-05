# Cheap XOR column mixers: exact support interface

The eight-row XOR mixer is worth a cost screen and a bounded proof screen.
It has a useful expected-measure comparison that retains the number of active
constituent rows. Its two-row cancellation cases are substantial codeword
families. They cannot be discarded as rare setup failures.

The comparison below does not yet prove 35% whole-code distance. No native
benchmark was run for this note. The existing GF256 Toeplitz construction
and its 72.3975-bit certificate remain the certified control.

## Construction and arithmetic

Let F=GF256, q=256, g=8, n=32, and k=16. Encode eight message rows with the
same RS[n,k] map. At each column, independently sample nonzero diagonal
entries in D1 and, for the second candidate, D2. The columns and groups
have independent setup randomness. Define

    H = I + 11^T,
    Hx = x + (sum_i x_i) 1.

The characteristic is two and g is even, so H^2=I. The candidate columns
are H D1 x and D2 H D1 x. Both maps are invertible and F-linear. The routing,
256 regions, and certified inner stay as in HD_GROUPED_SUBFIELD.md.

The scalar shared by both diagonals can be removed: divide D1 by its first
entry and multiply D2 by that entry. Sampling the normalized first diagonal
and an independent full second diagonal preserves the matrix distribution.
The proof may therefore analyze independent unnormalized diagonals.

| Column family | Variable products per eight outputs | Binary XORs |
| --- | ---: | ---: |
| Two current 4-by-4 Toeplitz columns | 18 | 30 |
| H D1 | 8 | 15 |
| Normalized D2 H D1 | 15 | 15 |

The new rows also use RS[32,16] instead of RS[64,32]. These operation counts
are not timing predictions. Normalized D2 H D1 is the stronger proof target:
its output amplitudes are independent nonzero field elements conditional
on the output supports. H D1 retains correlations between amplitudes in
different regions.

A nearby candidate replaces H by H_c=I+c11^T, for c outside {0,1}.
It is also an involution and adds one fixed multiplication per column.
Every matrix entry is then nonzero, so a single active coordinate spreads
to all eight outputs. Its two-coordinate equal-value cancellation remains.
The counts below concern c=1; they must not be reused for H_c without a
separate calculation.

## Two exact failures of stronger shortcuts

A fixed common H without per-column random diagonals can be absorbed into
the constituent message basis, because every row uses the same RS map.
Invertibility alone therefore does not provide the current atom bound.

For two active input coordinates x_1,x_2, the H output has weight two when
x_1=x_2 and weight eight otherwise. Under independent nonzero input
diagonals, those probabilities are 1/(q-1) and (q-2)/(q-1). An output
diagonal does not change the weights. In particular, D2 H D1 cannot satisfy
the current dimension-eight pointwise column cap.

There is also a deterministic outer-code example. Fix any diagonal setup,
two constituent rows, and any 31 of the 32 columns. The conditions that
their scaled values agree in these columns are 31 homogeneous equations
in 32 message coefficients. A nonzero solution exists. Neither row can
vanish identically: a nonzero RS[32,16] row cannot have 31 zero coordinates.
At the selected columns, the H output has at most two nonzero coordinates.
The last column has at most eight. Thus every setup has a nonzero outer
word of weight at most

    31*2 + 8 = 70.

This is at most 27.344% of the outer block. It is not a counterexample to
35% distance after the recursive inner. It does rule out a proof that
charges all two-row cancellations to rare bad setup.

## Exact expected membership and support counts

Write B_w for the number of RS[n,k] words on one prescribed support of
size w. Put d=n-k+1. The MDS support formula is

    B_0 = 1,
    B_w = 0                                     for 0<w<d,
    B_w = (q-1) sum_(j=0)^(w-d) (-1)^j binom(w-1,j) q^(w-d-j)
                                                for d<=w<=n.

For H D1, fix a nonzero complete output y and compute z_j=H y_j columnwise.
Let w_i be the weight of row i in z. Independence of the nonzero diagonals
gives the exact expected membership

    E[number of messages producing y] = product_i B_(w_i)/(q-1)^(w_i).

Each encoder is injective, so the expectation is also its membership
probability. This formula includes codewords chosen after setup.

For D2 H D1, define K(S,T) as the number of vectors z in F^g with support
exactly S and with supp(Hz)=T. Define P(S,T)=K(S,T)/(q-1)^|S|.
These are finite support kernels; P(empty,empty)=1.

Here is a direct calculation of K. Put s=|S|. The number of nonzero s-tuples
whose sum is zero is

    Z_s = ((q-1)^s + (q-1)(-1)^s)/q.

These tuples contribute Z_s to K(S,S). For a nonzero sum t, let B be the
subset of S whose entries equal t. Then the output support is [g]\B.
After dividing the entries by t, the remaining r=s-|B| entries avoid 0
and 1, and sum to 1+(|B| mod 2). Define

    L_0(0)=1,  L_0(1)=0,
    L_r(0)=L_r(1)=((q-2)^r+(q/2-1)(-2)^r)/q       for r>0.

For every B contained in S, add

    (q-1) L_(s-|B|)(1+(|B| mod 2))

to K(S,[g]\B). Addition inside the L argument is in characteristic two.
The zero-sum and nonzero-sum contributions are disjoint. Some can have the
same output support, so their counts must be added.

Fix a nonempty active constituent-row set A. For a prescribed sequence of
output supports T_1,...,T_n, its exact expected message count is

    sum_(S_1,...,S_n subset A)
        product_(i in A) B_(w_i) * product_j P(S_j,T_j),
    w_i = number of j with i in S_j,

where each w_i is required to be positive. For a prescribed output vector
with these supports, divide the expression by (q-1)^(sum_j |T_j|).
The output diagonal justifies this division even when the original row
messages are correlated.

This gives a finite generating function. Introduce variables x_i for i in A
and marks t_T for output masks, and define

    F_A(x,t) = sum_(S subset A) (product_(i in S) x_i)
                              sum_T P(S,T) t_T.

Raise F_A to n, extract each x-degree vector w, multiply by product_i B_wi,
and retain only positive row weights. Substituting t_T=z^|T| gives the
exact expected outer weight enumerator. For an ordered prescribed mask
sequence, use the product of the corresponding mask-specific polynomials
instead. A whole-code analysis must retain sufficient regional information;
an ordinary outer weight enumerator alone does not do that.

## A smaller comparison for a first certificate screen

The full eight-variable enumerator is not necessary for the initial screen.
Fix a nonempty active-row set A and put a=|A|. Let U_A be uniform on the
coordinate space F^A, embedded in F^g with zeros outside A. Let nu_A be the
law of D2 H U_A, using an independent nonzero output diagonal.

For one active RS row, independent nonzero coordinate scalars give expected
membership at most (q-1)^(-(n-k)) at every prescribed row vector. To prove
this, condition on an information set in its nonzero support. The prescribed
values then determine at most one shortened-code message. Each of the n-k
remaining coordinates independently costs at most 1/(q-1). An impossible
support or the zero target contributes zero.

Consequently, the expected measure for that active row is dominated by
q^n/(q-1)^(n-k) times the uniform measure on F^n. Multiply these comparisons
over the rows in A, fix zero rows outside A, and push the resulting measure
through H and the independent output diagonals. For every complete group
output y, the resulting comparison is

    mu_A(y) <= beta_a * nu_A^tensor(n)(y),
    beta_a = q^(na)/(q-1)^(a(n-k)).

Here mu_A counts messages with active constituent-row set exactly A.
The comparison includes zero coordinates and zero rows artificially; it
does not assert that actual active rows can be zero. At q=256,n=32,k=16,
log2(beta_a) is approximately 128.09034*a, rather than about 1024 per group
regardless of a.

The support law of nu_A is particularly small:

    Pr_nu_A[supp(Y)=T] = N_A(T)/q^a,
    N_A(T) = sum_(S subset A) K(S,T).

Conditional on T, every output amplitude is independent and uniform in F*.
The n columns are independent in this comparison. For a=8, U_A is uniform
on all of F^8, and nu_A is exactly uniform. Thus beta_8 and its comparison
are the old grouped-outer interface, without an extra density penalty.

For A={i}, the masks are empty with probability 1/q and [8]\{i} with
probability (q-1)/q. For A={i,j}, the complete mask distribution is

| Output mask | Probability |
| --- | ---: |
| empty | 1/q^2 |
| [8]\{i} | (q-1)/q^2 |
| [8]\{j} | (q-1)/q^2 |
| {i,j} | (q-1)/q^2 |
| [8] | (q-1)(q-2)/q^2 |

This retains the pair cancellation as an ordinary comparison branch.
The lower-dimensional class also pays a much smaller message count.

To use an existing local inner bound formulated for uniform values including
zero, a further positive comparison is available. Uniform nonzero values
on a mask T are dominated by (q/(q-1))^|T| times uniform F^T. Use the mixture
weights

    lambda_A(T) = N_A(T)/q^a * (q/(q-1))^|T|.

The sum of these weights need not equal one; its excess must be retained
in the calculation. For a=8, use the exact full-space comparison instead
of this looser decomposition.

The remaining difficulty is support correlation, not amplitude correlation.
One column produces a mask across eight regions. After independent regional
permutations, the occupied positions are uniform conditional on each regional
count, but the eight counts remain correlated. Retain the column-mask
generating function when averaging these counts. Replacing their joint law
by independent region occupancies is unjustified. The source-row-set union
must also include mixtures of different A across groups.

## Initial numerical screen and validation

For two fixed active rows, the exact expected outer weight polynomial has
a compact form. For row weights w_1,w_2 and overlap u, the number of ordered
support pairs is

    binom(n,w_1) binom(w_1,u) binom(n-w_1,w_2-u).

Multiply this by B_w1 B_w2 and by

    z^(7*(w_1+w_2-2u)) * ((z^2+(q-2)z^8)/(q-1))^u,

then sum over w_1,w_2>=17 and feasible u. This polynomial applies to either
XOR candidate because output diagonals preserve support. Dividing its
coefficients by q-1 counts projective scalar orbits in expectation.

The helper evaluates the integer formula with floating probabilities for
q=256. The following values concern one fixed pair in one group.

| Outer weight at most | log2 expected projective count, FLOAT |
| --- | ---: |
| 34 | -98.8341 |
| 40 | -76.3471 |
| 52 | -36.3149 |
| 64 | -6.3726 |
| 70 | 5.9270 |
| 89 | 39.7073 |

These are expectations, not probabilities or minimum-distance receipts.
The total expected message count agrees with (q^16-1)^2.

The exact helpers and finite-field tests are
[theory/hd_xor_column_probe.py](theory/hd_xor_column_probe.py) and
[theory/hd_xor_column_probe_tests.py](theory/hd_xor_column_probe_tests.py).
Run from the theory directory:

    python -B -m unittest hd_xor_column_probe_tests -v
    python -B hd_xor_column_probe.py

The tests enumerate small-field support laws, verify H^2=I, compare direct
setup membership with the formula, check the prescribed MDS support counts,
and verify the uniform-comparison law including the full-eight-row case.

The next proof screen should start with a=1 and a=2 using nu_A and beta_a.
Then test a=7 as the class closest to the old full-space envelope. Retain
the exact two-row enumerator as a refinement if its simpler comparison
does not close. Only after these screens succeed should the calculation
expand to all row-set mixtures and a directed whole-code union. Native
cost comparisons should include both H D1 and normalized D2 H D1; the latter
is the recommended first target for certificate work.

## Scalar-route screening with correlated masks

The next screen uses the unchanged scalar route and inner. It finds useful
negative floating bounds for one active group with one, two, or seven active
constituent rows. The one- and two-row cases retain exact RS support counts.
No complete coverage or directed certificate is claimed. These results do
not transfer to the separate alternating-pair route.

### A product comparison that keeps its density cost

For 0<p<=1, define f_p on F: sample a uniform field element with probability
p, and output zero otherwise. Thus f_p(0)=1-p+p/q and f_p(x)=p/q for x!=0.
The correlated column law nu_A is dominated by c_a(p) times f_p^tensor(8),
where a=|A| and

    c_a(p) = max_T [N_A(T)/(q^a (q-1)^|T|)]
                   / [(p/q)^|T| (1-p+p/q)^(8-|T|)].

This pointwise maximum replaces correlation by a proved domination. It does
not declare the original regions independent. The group density becomes
beta_a c_a(p)^32. Summing over the binom(8,a) possible row sets gives a common
product comparison for that homogeneous row-count class. The artificial
potential positions are independent Bernoulli thinnings of the originally
uniformly routed positions.

The existing packet operator averages this thinning before its fractional
power. Each potential block occupancy l uses a Binomial(l,p) average of
its local quantities. Independent thinning across the four emission blocks
preserves their positive product recurrences. A separate positive recurrence
isolates empty thinning and uses the inverse-refresh cap on that event.

The product comparison works well for seven active rows. It is too loose
in the tested one- and two-row classes. Selected FLOAT log2 bounds follow;
each includes the group-set union and all row-set choices in its class.

| Active groups | a=1 | a=2 | a=7 |
| --- | ---: | ---: | ---: |
| 1 | 28.36 | 34.89 | -15.45 |
| 16 | 322.04 | 429.54 | -2126.56 |
| 128 | 1755.65 | 2521.84 | -1697.40 |
| 512 | not screened | not screened | -2633.16 |
| 2048 | not screened | not screened | -21514.87 |
| 4096 | not screened | not screened | -90943.59 |
| 8192 | not screened | not screened | -1236657.28 |

The first three rows use regional coefficient extraction. Later rows use
the global marker. These are proposed witnesses, not exhaustive optimization
or statements about a failing construction.

A common product comparison for all row-set mixtures is also valid: use
the group density sum_a binom(8,a) beta_a c_a(p)^32. It does not close the
tested counts. Its FLOAT logarithms at 1,128,512,2048 groups are about
59.09,5321.46,19095.66,41600.12. Thus this simple aggregation does not remove
the need to retain the smaller-row classes more carefully.

### Exact one-group operators

For one active group, each region receives either no nonzero input or exactly
one. The route places the latter uniformly among 8192 positions. This special
case admits a small exact support calculation.

Put d=q^16-1 and b=d/(q-1)^16. Let P_0 be the existing empty-packet operator,
including its inverse-refresh live entry. Let Jbar_1 and V_1 be the established
return and live caps for one uniform potential coordinate. A forced nonzero
input X has C X!=0, because every column of C is nonzero. Starting at zero,
the packet therefore cannot remain at zero and its emitted weight is one.
A valid forced-nonzero operator is

    P_* = [[0,                b z],
           [kappa Jbar_1/d,   kappa b V_1]],
    kappa=q/(q-1).

The flat-state entries follow by conditioning the uniform input on being
nonzero. Their caps retain kappa. The first row uses the actual nonzero-input
behavior. A uniform-input comparison would reintroduce impossible quiet
inputs and loses much of the sparse bound.

Average over the occupied packet's location in a region:

    R_0 = P_0^128,
    R_* = (1/128) sum_(j=0)^127 P_0^j P_* P_0^(127-j).

The single-input packet operator already averages its 64 input positions.
For an output mask T within an eight-region column block, define

    Q_T = product_(i=0)^7 R_(1[i in T]),

where R_1 means R_* and the product is chronological. State carries across
regions. Independent regional permutations and fresh inner coins justify
this expectation for the comparison inputs. Actual outer coordinates are
not assumed independent before their support and output-diagonal randomness
have been accounted for.

For an input mask S, define M_S=sum_T P(S,T) Q_T. Its conditional amplitudes
are already nonzero, so no additional kappa factor belongs in M_S. The
factors needed for the local return/live caps are inside P_*.

For a fixed one-row set A={i}, define G_A(x)=M_empty+x M_{i}. Its expected
weighted message count is at most

    e_0 [sum_(w=17)^32 B_w [x^w] G_A(x)^32] 1.

For A={i,j}, use

    G_A(x,y)=M_empty+x M_{i}+y M_{j}+xy M_{i,j},

and sum its coefficient matrices against B_u B_v for 17<=u,v<=32. In
particular,

    M_{i,j}=(Q_{i,j}+(q-2)Q_[8])/(q-1).

The pair-cancellation branch remains in this polynomial. The matrix product
retains chronological order and regional positions. Only a degree-32
polynomial in one or two variables is required. Its coefficients already
count ordered support sequences; no extra binomial coefficient belongs here.

Sum over the row sets A of the declared size. Multiply by 8192 for the
active-group choice, divide once by q-1 for global projective scalar orbits,
and multiply by z^(-734003) for forbidden output weights. This uses the
ordinary first moment, with alpha=1. The route average precedes Markov.

The resulting FLOAT screens are:

| One active group | Weight mark z | log2 bound |
| --- | ---: | ---: |
| One active constituent row | 0.9952831508888066 | -1043.6402 |
| Two active constituent rows | 0.9982002694097376 | -300.5626 |

Replacing the forced-nonzero operator by uniform-amplitude domination gives
approximately -133.54 and -23.24. Retaining its exact first row is therefore
material. Neither calculation drops codewords chosen after setup.

Seven active rows can use the smaller nu_A comparison directly. Average its
exact mask mixture over each eight-region block before taking the entrywise
fractional power. Jensen bounds the route average at that block boundary.
The one-group FLOAT result, including all row-set and group choices, is
approximately -95.34. Without exact MDS support floors, the corresponding
one- and two-row results remain about 5.16 and -8.18. Their support restriction
is essential in these screens.

### Validation and remaining coverage

The new code is [theory/hd_xor_column_screen.py](theory/hd_xor_column_screen.py),
with [its tests](theory/hd_xor_column_screen_tests.py). It reuses the existing
scalar-route coefficient calculation and inner bounds.

    python -B -m unittest hd_xor_column_probe_tests hd_xor_column_screen_tests -v
    python -B hd_xor_column_screen.py --active 1 2 --one-group-exact --iterations 30
    python -B hd_xor_column_screen.py --active 7 --one-group-correlated --iterations 55
    python -B hd_xor_column_screen.py --active 7 --rows 512 2048 4096 8192 --iterations 65

All fifteen tests pass. They check direct thinning averages, agreement with
the unthinned control, noncommuting matrix-polynomial coefficients, and
exclusion of quiet behavior for a single
forced nonzero input. An independent audit found no blocker in the one-group
normalization, conditional operators, ordered polynomial, projective divisor,
product-density comparison, or thinning recurrence. All numerical bounds
above still require directed replay before acceptance as receipts.

An optional input-support mark retains sum_i w_i>=17a in the product comparison.
Replace N_A(T) by sum_(S subset A) u^|S| K(S,T), and multiply the group density
by u^(-17a), for u>=1. The `--input-mark` screen implements this refinement.
A bounded search improved the one-row class at 128 groups from about 1756
to 191 but did not close it. It does not replace the exact sparse calculation.

The next step is to extend the support and nonzero-amplitude treatment to
several active groups, then cover all constituent-row mixtures and integer
group counts. A region can then contain multiple inputs, whose feedback can
cancel. P_* cannot be reused as if every occupied packet had one input.
Retain group or row-support marks and derive small-occupancy operators for
nonzero amplitudes. Follow successful floating screens with directed sparse
receipts and a complete union. These results justify continuing the scalar
XOR2 proof; they do not yet certify its encoder timing or any new route.
