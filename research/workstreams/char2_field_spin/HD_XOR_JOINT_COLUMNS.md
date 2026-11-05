# Joint column weights and occupancy envelopes for XOR2

For exactly five active groups, each with two nonzero RS[32,16]
constituents, a directed first-moment bound is below

    2^(-129.9217876158025605).

The bound covers every active-group choice and all 28^5 constituent-pair
choices in the unchanged scalar-route D2 H D1 construction. A polynomial-time
occupancy envelope also closes every class with six through sixteen active
groups. The directed log2 bound on that eleven-class union is below
-63.04480338044788. All classes in this note have exactly two active
constituents per active group.

The proof retains
the joint column-category contractions from
[HD_XOR_COLUMN_POTENTIAL.md](HD_XOR_COLUMN_POTENTIAL.md). It replaces the
separable marked count T(f_1,f_2)^5 by a positive polynomial tail.

At the same z and state potential, the earlier separable FLOAT bound was
about -13.6184 bits. The improvement is in the proof's support counting,
not in the construction. The unrestricted two-constituent classes above
sixteen groups, other constituent counts, and mixed classes remain open.
A scoped seventeen-group control below covers only groups sharing one
constituent pair.

## Keep the column categories before counting supports

Fix r active groups and a constituent pair in each group. At one outer
column, each group has one of four source-support patterns: empty, first
constituent only, second constituent only, or both constituents. Let s count
singleton groups and d count double groups in that column.

Fix a weight mark 0<z<1 and a state potential v=(1,t)^T with 0<t<=1.
For the column's source-support tuple S, let M(S) be its chronological
eight-region operator, averaged over the exact XOR cancellation mixture.
This is the operator defined in HD_XOR_COLUMN_POTENTIAL.md. Choose positive
category caps c_(s,d) such that

    M(S) v <= c_(s,d) v                                   (1)

for every source-support tuple in that category. The inequality is
componentwise. In particular, choose c_(0,0)=1 exactly.

The directed checker verifies (1) over all group multisets. It preserves
the eight labelled regions and all independent cancellation coins, including
repeated constituent pairs. At r=5 there are 749,398 universal checks.
These checks cover source-support membership; their count is not a factor
in the message union.

For a fixed message, iterate (1) across its 32 columns. The initial state
is zero. Bounding the terminal all-ones vector by v/t costs 1/t once:

    E[z^W] <= (1/t) product_(j=1)^32 c_(s_j,d_j).           (2)

The expectation is over the construction's setup coins: D1, D2,
independent regional permutations, and the inner randomness. No message
or RS support is sampled from a surrogate distribution. The category
weights can be used for counting only after (1) is verified universally.

## Bound prescribed RS supports, then relax their floors

For a prescribed support S of size w>=17, let B_w count RS[32,16] words
whose support is exactly S. Shortening the RS code to S gives dimension
w-16. An information set of that shortened code determines its word
injectively. Every coordinate in that information set is nonzero for a
word counted by B_w. Therefore

    B_w <= 255^(w-16),    17<=w<=32.                       (3)

This bound is used only for valid nonzero constituent words.

A fixed r-group message has 2r constituent supports with weights w_i>=17.
Write J=sum_i w_i. Their prescribed-support count is at most

    product_i B_(w_i) <= 255^(J-32r),
    J >= 34r.                                             (4)

First apply (4) to valid support arrays. Then enlarge the nonnegative sum
to every binary 2r-by-32 support array whose total weight is at least 34r.
The enlarged set may include a row of weight less than 17, including zero.
Equation (3) is not asserted for those added rows. Their artificial positive
weight 255^(J-32r) only enlarges the sum after the valid-message bound.

This relaxation drops the individual support floors but retains their sum.
It keeps the joint category weight in every column, which the separable
T(f_1,f_2)^r bound had replaced by a product of per-group marks.

## A univariate positive polynomial counts the enlarged set

For the fixed labelled constituent pairs, there are

    binom(r,s) binom(r-s,d) 2^s

source-support patterns in category (s,d). Each has total column support
s+2d. Define

    P(X) = sum_(s,d>=0; s+d<=r)
             binom(r,s) binom(r-s,d) 2^s c_(s,d) X^(s+2d).

The coefficient of X^J in P(X)^32 sums the products of category caps over
all support arrays of total weight J. This is a combinatorial identity;
it requires no independence hypothesis about RS supports.

By (2)--(4), the marked message count for the fixed constituent pairs is
at most the following polynomial tail:

    C_r = 255^(-32r) sum_(J=34r)^(64r)
                        [X^J] P(X)^32 255^J.              (5)

The implementation evaluates P(255X)^32 and sums its coefficients from
34r through 64r, then divides by 255^(32r). All coefficients are positive.
For r=5, the polynomial degree is 320 and the retained tail starts at 170.
No ten-variable coefficient table or Chernoff relaxation is needed.

There are binom(8192,r) active-group choices and 28^r constituent-pair
choices. All nonzero global scalar multiples have the same output weight.
Applying Markov to projective message classes gives

    Pr[there is a nonzero message in this class with W<=734003]
       <= binom(8192,r) 28^r C_r z^(-734003) / (255t).       (6)

The denominator has one projective factor of 255 for the whole message.
The counting relaxation in (5) is an upper bound on the scalar-invariant
sum over actual messages; it does not introduce another normalization.

## Category caps without source-multiset enumeration

Fix a category (s,d), and let K count double groups whose D1 values cancel.
For every fixed source-support tuple in that category,

    K has distribution Binomial(d,1/255).

Condition on K=k. Each uncancelled double occupies every region. Each
cancelled double occupies two regions, and each singleton occupies seven.
Every resulting labelled occupancy sequence therefore satisfies

    d-k <= m_i <= s+d,   0<=i<8,
    sum_i m_i = 7s+8d-6k.                                 (7)

This box-and-sum condition can admit sequences that no given source tuple
can produce. Enlarging to those sequences is a valid upper bound.

For a lower occupancy l and upper occupancy h, define vector bounds F_j(u)
for products of j regional matrices. The index u counts occupancy above
the lower baseline jl. Initialize F_0(0)=v and omit infeasible indices.
Use the recurrence

    F_j(u) = componentwise max_(0<=a<=h-l; u-a feasible)
                 R_(l+a) F_(j-1)(u-a).                    (8)

For every labelled sequence with occupancies in [l,h] and shifted total u,
its product applied to v is componentwise at most F_j(u). To prove this,
remove the first regional matrix and apply the induction hypothesis to the
suffix. Multiplication by the nonnegative first matrix preserves the bound.
The two coordinates of F_j may be maximized by different suffixes; that
only enlarges the vector bound. The recurrence does not assert that any
single sequence realizes both maxima. It keeps chronological products and
does not assume that regional matrices commute.

In (7), take l=d-k and h=s+d. The shifted total is

    sum_i m_i - 8l = 7s+2k.

The following vector bounds the column's complete cancellation mixture:

    U_(s,d) = sum_(k=0)^d binom(d,k) 254^(d-k)/255^d
                  F_8(7s+2k; l=d-k, h=s+d).

Choose c_(s,d)=max(U_(s,d)[0], U_(s,d)[1]/t). This choice proves (1)
for every source-support tuple, without enumerating those tuples. The exact
cancellation-count probabilities are averaged before the final potential
ratio. The empty category is checked separately and assigned c_(0,0)=1.

Directed evaluation takes the upper endpoint of each candidate in (8),
then compares those exact endpoints. It does not order overlapping Arb
balls. These endpoint maxima remain valid upper bounds under subsequent
nonnegative matrix multiplication. The final category maximum uses the
same endpoint rule.

All categories reuse tables indexed by (l,h). Through r groups there are
binom(r+2,2) tables. Their final layers contain

    sum_(w=0)^r (r-w+1)(8w+1)

two-state vectors. Work is O(r^4) and retained table storage is O(r^3),
with the eight-region length fixed. At r=16 there are 153 tables and
6,681 vectors. The positive support-count polynomial still has degree 64r.
The evaluator now permits r<=32 after checking the underlying formulas:
packet profiles cover all degrees through 64, every block occupancy stays
at most 16, and regional normalization is defined through 8192. This guard
extension does not itself certify any additional class.

## Directed receipts and validation

The saved witness uses

    r=5, z=99829/100000, log2(t)=-499/8.

The 21 literal rational exponents log2(c_(s,d)) are saved in
[theory/hd_xor_joint_columns.py](theory/hd_xor_joint_columns.py) under
`SAVED_WITNESSES`. Replay does not run a floating optimizer. It checks
every inequality (1) with Arb, then evaluates (5) and (6) with Arb.

The 256-bit result is

    [-129.9217876158025605116594248224667007488248212332812505603630050535137816
      +/- 2.65e-71].

The fourteen tests in
[theory/hd_xor_joint_columns_tests.py](theory/hd_xor_joint_columns_tests.py)
also replay the five-group witness at 384 bits. That replay agrees with the displayed
256-bit interval and encloses the same log2 bound with radius 3.41e-110.
They independently check:

- every binary one-column pattern for one through five groups;
- the information-set envelope at every relevant RS support weight;
- all 4,096 binary support arrays in a two-group, three-column instance;
- the weighted message sum from independently generated GF4 RS[3,2] words;
- rejection of an invalid empty cap and a false nonempty cap.

The additional envelope tests exhaust noncommuting toy products, check
(7) against every two-group cancellation branch, and compare envelope caps
with exact three-group caps. They also test overlapping-ball endpoint
maxima, the six-group streaming extension, and the directed resource-family
union. An independent read-only audit found no gap in the counting or
resource-envelope arguments and their implementations.

The small support-array test explicitly includes newly admitted arrays
whose individual row floors fail. It checks the intended aggregate-floor
relaxation, rather than silently treating that relaxation as equality with
the original RS message set.

The regional replay uses the existing four-plus-four factorization. It
stores 1,296 prefix matrices and suffix vectors, processes source cases
in bounded batches, and caches at most 16,384 occupancy actions. This
avoids a full 6^8 table and a complete five-group source-case table.

A six-group control used the full literal checker before the occupancy
envelope was introduced. It checked 5,245,786 source multisets and
157,773,595 cancellation branches, giving a directed log2 bound of
-127.840307156752122 at 256 bits. The FLOAT screen took about 31 seconds;
the directed replay took about 189 seconds. These are proof-evaluator costs,
not native encoder benchmarks. The control and the five-group literals
remain saved; the exponential six-group replay is not repeated by default.

The resource-envelope witnesses below are literal rational z and log2(t)
values. Their category caps are derived by directed recurrence (8), not
proposed by a floating optimizer during replay.

| r | z | log2(t) | Directed log2 bound |
| ---: | --- | --- | ---: |
| 6 | 49907/50000 | -3997/64 | -127.8454999577 |
| 7 | 997990/1000000 | -4002/64 | -125.2305057312 |
| 8 | 997840/1000000 | -4007/64 | -122.4090081647 |
| 9 | 997740/1000000 | -4004/64 | -120.3725052217 |
| 10 | 997590/1000000 | -4001/64 | -122.8830366773 |
| 11 | 997440/1000000 | -3998/64 | -126.1400254569 |
| 12 | 997290/1000000 | -3995/64 | -130.0777906853 |
| 13 | 997140/1000000 | -3992/64 | -134.5117327562 |
| 14 | 996990/1000000 | -3997/64 | -132.8230974180 |
| 15 | 996840/1000000 | -4002/64 | -110.3640309424 |
| 16 | 99669/100000 | -4063/64 | -63.0448033804 |

All eleven receipts agree at 256 and 384 bits. Their directed union at
384 bits is

    [-63.0448033804478875669277837338750520427251584669072903770428166501237209226109794790277578478069134229122930553
      +/- 4.81e-110].

Combining this union with the five-group receipt and the earlier r=1..4
column-potential union gives a directed log2 bound below
-44.76797103805864 for all homogeneous two-constituent classes with
1<=r<=16. This combines distinct active-group counts; it does not include
messages with mixed constituent counts across groups.

The full envelope family replay took about three seconds at each precision.
Witness selection used five local FLOAT points per r=7,...,16. The initial
r=16 point reached only -39.65 bits; five further potential values at fixed
z found the saved -63.04-bit witness. No optimization over an unbounded
parameter domain or extension beyond r=16 was attempted.

Run from the theory directory:

```text
python -B -m unittest hd_xor_joint_columns_tests -v
python -B hd_xor_joint_columns.py --groups 5 --replay --precision 384
python -B hd_xor_joint_columns.py --groups 5
python -B hd_xor_joint_columns.py --groups 6 --envelope --replay --precision 384
python -B hd_xor_joint_columns.py --replay-envelope-family --precision 384
```

## Next proof step

The aggregate-floor count and resource envelope remove the exponential
support enumeration for this homogeneous two-constituent family. The first
observed difficulty was the smaller r=16 margin under inherited parameters;
a bounded adjustment of t resolved it. This is not evidence of a
construction defect or an impossibility at larger r.

The bounded extension below identifies a geometry-counting issue before
a suffix-frontier refinement becomes useful. The existing r=1..4
column-potential receipts and r=5 literal receipt remain unchanged.

## Bounded extension and the first unresolved class

The next round capped numerical work at eight minutes. It used one
inherited point and at most forty local optimization evaluations for each
of r=17,20,24,32. Five additional z values checked the upper boundary of
the r=32 search box. Total local-search evaluation time was about one minute;
no exponential source enumeration was resumed. The best recorded FLOAT
values in these bounded searches were:

| r | z | log2(t) | FLOAT log2 bound |
| ---: | ---: | ---: | ---: |
| 17 | .9965482169227191 | -63.5 | 9.37406970 |
| 20 | .9961345155811148 | -63.5 | 229.42040667 |
| 24 | .9955894533856386 | -63.5069298254 | 516.60732365 |
| 32 | .9945399873746938 | -63.6140820122 | 1043.70432014 |

These are failed bounded screens, not optimized lower bounds or evidence
of encoder failure. No unrestricted new class received a directed
certificate. The previously saved receipts remain valid.

At the inherited r=17 point z=.99654 and log2(t)=-63.504375, the marked
support tail was dense: its largest term had total
support 1082 out of a possible 1088. The aggregate support floor was 578.
Thus the individual-floor relaxation is not the dominant loss at this
point.

The focused probe
[theory/hd_xor_geometry_probe.py](theory/hd_xor_geometry_probe.py) evaluates
the 28 geometries in which every double group uses the same pair. For
each fixed pair, it averages every cancellation-count branch before taking
the potential ratio. At z=.996548217 and log2(t)=-63.5, it obtains

    maximum repeated-pair row-0 cap: -259.94396536054614,
    maximum repeated-pair row-1 cap: -259.94400817669117,
    resource-envelope c_(0,17):     -259.94396536054614.

These are FLOAT log2 values. The maximizing pairs are the last two and
first two labelled regions, respectively. The dominant cancellation counts
in the envelope's two row sums are 17 and 0: rare complete cancellation
competes with ordinary dense live-state evolution.

The dominant category already matches a genuine repeated-pair geometry
under the existing regional operators. A Pareto frontier of suffix vectors
cannot reduce that required category cap at these parameters. It would
target coordinatewise-max slack elsewhere, not the observed bottleneck.

A counterfactual diagnostic deleted selected entries of the regional cap
matrices. Deleting positive-occupancy zero-to-zero entries had no visible
effect. Deleting live-to-zero entries improved a small potential scan only
to about -15.66 bits. These modified matrices are not valid operators for
the construction and certify nothing. The diagnostic does not identify
an actual encoder failure or prove that a return event attains its cap.

## Directed common-pair control at seventeen groups

Restrict the message class so that all r active groups use one shared
constituent pair. There are 28 choices for that pair, instead of 28^r
independent pair choices. Every category inequality and the support
polynomial above still upper-bound this subclass. Replace only the pair
factor in (6):

    common-pair failure probability
       <= binom(8192,r) 28 C_r z^(-734003)/(255t).           (9)

At r=17, z=996548217/1000000000, and log2(t)=-127/2, the directed result is

    [-67.5436393308537773319903487954042307917343736601654380272229851945722548911300847181289327990063023030283122942
      +/- 5.18e-110].

The 256-bit replay agrees. The improvement is exactly 16 log2(28), about
76.92 bits, relative to charging these caps at 28^17 choices. The receipt
sets `constituent_pair_geometry` to `one shared pair` and
`full_homogeneous_group_class` to false. A test confirms that the same
marks with the unrestricted 28^17 factor fail the 40-bit target.

```text
python -B hd_xor_joint_columns.py --groups 17 --envelope --common-pair --replay --precision 384
python -B -m unittest hd_xor_geometry_probe_tests -v
```

The four focused geometry tests validate exact multiset multiplicities,
match repeated-pair operators against the ordered two-group checker,
check a small Hölder screen, and reproduce the dominant r=17 category.

## Geometry moments as the next structural step

The common-pair control shows that charging the worst geometry at every
pair choice can lose useful mass. It does not bound the remaining
geometries. A partition by collision multiplicity or a geometry moment
must account for those choices explicitly.

For a fixed source-support array, let sigma encode the constituent-pair
geometry shared by all 32 columns. If c_j(sigma) is its scalar column cap,
Hölder gives

    E_sigma product_(j=1)^32 c_j(sigma)
       <= product_(j=1)^32 (E_sigma c_j(sigma)^32)^(1/32).   (10)

Here sigma is a normalized index of the message union, not additional
setup randomness. To make singleton categories symmetric, represent each
unordered pair with a uniformly ordered labelling. Swapping that labelling
also swaps its two constituent words, so the auxiliary 2^r representations
cancel. The original union factor remains 28^r. A fixed relative singleton
index is then uniform over the eight rows; a double source uses a uniform
unordered pair.

For a category with s singleton and d double sources, its multiset moment
weights are therefore

    (s!/product_i n_i!) / 8^s
       times (d!/product_e m_e!) / 28^d,

where n_i and m_e are singleton-row and double-pair multiplicities.
Uniform weighting of the multiset cases would be incorrect.

The focused FLOAT implementation evaluates these exact assignment counts
only for r<=3. At fixed marks, replacing universal maxima by the order-32
geometry moments changes the joint bound from -137.5741 to -140.1543 for
r=2, and from -134.6070 to -137.2204 for r=3. These are diagnostic screens,
not directed receipts or evidence that the improvement scales to r=17.
Large-r moments need a collision stratification or a rigorously bounded
moment sum; no new exponential enumeration is proposed.

One cheap counting ingredient is available. Let mu be the largest number
of active groups using the same pair. The number of ordered pair assignments
with maximum multiplicity at most mu is

    r! [X^r] (sum_(j=0)^mu X^j/j!)^28.

Subtract the expression with mu-1 to count the exact stratum. This count
is small to compute. The missing proof ingredient is a useful category or
geometry-moment cap for each stratum. For example, among k cancelled
double groups, the sum of occupancies at any two regions is at most
k+mu. With s singleton groups and d-k uncancelled doubles, the corresponding
bound is 2s+2d-k+min(k,mu). Incorporating this constraint must retain the
labelled region order and must not silently replace a fixed source geometry
by independently resampled geometries in each column.

## Mixed constituent counts: a possible interface, not a certificate

The whole XOR2 outer also contains groups with active-constituent counts
a other than two, and messages mixing different a values. Extending only
the homogeneous two-constituent family does not cover those messages.

A weight-only local law can provide a common interface. Suppose one column
has b nonzero source coordinates, with independent D1-scaled values x_i
in F*. Put tau=sum_i x_i. If tau=0, the output support is the source
support itself, of size b. Its probability is

    Pr[tau=0] = (255^b + 255(-1)^b)/(256*255^b),  b>=1.

If tau is nonzero, normalize by tau. The output vanishes exactly at source
coordinates whose normalized value equals one. If there are j such
coordinates, the output weight is 8-j. For u=b-j>0, the number of ways
to choose the remaining values from F\{0,1} with the required sum, either
zero or one, is

    E_u = (254^u + 127(-2)^u)/256.

This follows by additive-character summation: nontrivial additive characters
chi satisfying chi(1)=1 contribute (-2)^u; those with chi(1)=-1 contribute zero.
For u=0, the required sum is satisfied exactly when j is odd. Consequently
the joint probability of tau nonzero and exactly j normalized values equal
to one is

    255 binom(b,j) E_(b-j) / 255^b,

with the stated u=0 parity convention. The empty source has output weight
zero with probability one. These branches define probabilities p_(b,w)
for the local output weight w; branches giving the same w are added.

For a group with a active constituent rows, the polynomial

    A_a(X,Y,Z) = sum_(b=0)^a binom(a,b) X^b
                         sum_w p_(b,w) Y^w Z^(1_(w=8))

records source support size, output support size, and whether every region
is occupied. Convolving these polynomials for a fixed mixed profile gives
coefficients indexed by source-column weight J_col, output-column weight L,
and the number F of full-output groups. Every resulting occupancy sequence
satisfies F<=m_i<=r and sum_i m_i=L. Thus the same directed resource
recurrence can bound its action. This bound is uniform over output-mask
locations; it does not resample fixed constituent geometry across columns.

If A is the total number of active constituent words, the corresponding
information-set bound is 255^(J-16A) and the total RS support floor is
J>=17A. The row-set factor is product_g binom(8,a_g), with the appropriate
multinomial factor for active-group type choices. These replace the
two-constituent factors, rather than being multiplied into them.

This is a proposed counting interface, not an implemented mixed-class
bound. A dense three-index convolution can be memory-heavy, and the
weight-only location envelope may be too loose. Summing over mixed global
profiles is a separate coverage problem. Geometry stratification or
Hölder moments should be designed together with this interface before
investing in a large mixed-profile evaluator. No mixed-message certificate
or construction change is claimed here.

No encoder change, native benchmark, remote operation, or commit was made.
