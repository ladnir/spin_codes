# Routing complete outer groups in packets

The proposed route keeps four adjacent outer groups in one source packet.
This changes the sampled code: four potential coordinates now move together
in each region. The existing 35% / 72.3975-bit scalar-route certificate does
not certify this distribution. A fixed permutation of final output coordinates
is different: it preserves each sampled code's distance exactly.

The reduction below gives valid conservative bounds for packet routes.
Its implementation is a FLOAT screen, not a directed certificate. Initial
four-group screens fail on medium families of fully occupied source packets.
Two-group packets with alternating cyclic matchings admit stronger h4
FLOAT screens. No complete packet-route claim is made.

## Construction and comparison measure

Work over F=GF256. There are 8192 independent outer groups, each encoding
128 inputs into 256 coordinates. The certified outer membership bound is
255^(-128) at every fixed nonzero target group. It applies to the independently
sampled invertible Toeplitz column mixers described in HD_GROUPED_OUTER.md.

First consider four-group packets. Packet p contains groups 4p,...,4p+3,
so each of the 256 regions contains 2048 packets. Each region samples an
independent uniform permutation pi of these packets. Sixteen consecutive
destination packets supply one inner step, in the order

    X[16*lane+j] = region[4*pi(packet_index+j)+lane].

The inner has 16 state coordinates and 64 outputs per step. Its state starts
at zero and its final state is unrestricted. For local rank h=4, the emission
map A is four RS[16,4] blocks, C=A^T, and CA=0. Each step computes

    Y = X + A Q,
    Q_next = D2 H16 D1 (Q + C X).

The two diagonals have independent nonzero entries. All refreshes, regional
permutations, and outer group setups are independent.

Fix an active set U of r outer groups. The expected outer counting measure
for messages whose active group set is exactly U is dominated by

    beta^r * Uniform(F^(256r)),  beta = 256^256 / 255^128.

Coordinates outside U stay zero. This is a comparison of measures, not an
assertion that encoded coordinates are independent. Under the comparison
measure, however, every potential coordinate is independently uniform in F.
Routing complete packets preserves this independence of values. It correlates
the locations of the potential coordinates.

## A local operator for correlated locations

For one inner step, let l_b count potential inputs in emission block b.
The existing common-refresh operator M_z(l_0,...,l_3) is valid for each fixed
potential set with these counts. Its entries act on zero-state mass and
a density bound over nonzero states. In the notation of the existing proof,

    M_z(l) = [[K, b(T-K)], [Jbar/d, bV]],
    d = 256^16-1, b = d/255^16.

For an empty potential set, replace its live-to-live entry by
min(bV,L0(z)), using the established inverse-refresh bound. All other
entries are unchanged. This local argument does not require independent
potential positions. It requires the independent uniform values supplied
by the comparison measure and independent refresh randomness.

Fix 0<z<1 and 0<alpha<=1. Let M_z(l)^[alpha] denote entrywise powers.
The fractional-moment argument expands the product over the two state types
and uses (sum_i a_i)^alpha <= sum_i a_i^alpha for nonnegative a_i.
Consequently a fixed sequence of potential sets contributes at most the
corresponding product of these entrywise powers. Composition of two step
operators before taking powers gives a stronger bound by the same inequality.

## Coefficients indexed by packet occupancy

Each source packet has a mask S contained in {0,1,2,3}: its active groups.
Keeping all sixteen mask counts gives an exact multivariate route enumeration.
The following four-count envelope is much smaller.

Let n_s be the number of source packets containing s active groups, for
1<=s<=4, and put n_0=2048-sum_s n_s. Then

    r = sum_s s*n_s,
    H(n) = multinomial(2048; n_0,...,n_4) * product_s binomial(4,s)^n_s

counts all active-group sets in this class. There is no fractional power on
this support count: the proof takes a union over these sets.

For a single step, let c_s count packets of size s among its sixteen slots,
including c_0. To bound all possible masks, sort the four lane counts l in
decreasing order. A binary incidence matrix with column sizes described by
c and row sums l exists exactly when

    sum_b l_b = sum_s s*c_s,
    sum_(b=1)^k l_b <= sum_s min(k,s)*c_s  for k=1,...,4.

These are the Gale--Ryser degree constraints. They retain the correlations
between lanes without treating individual scalar positions as permuted.
The common-h4 operator is invariant under permuting its four emission blocks.
Define E_z(c) as the entrywise maximum of M_z(l) over these feasible counts.
Different entries may attain their maxima at different l; this only enlarges
the positive operator. It does not posit an actual route with those maxima.

Set x_0=1 and define the matrix polynomial

    F(x_1,...,x_4) = sum_(c_0+...+c_4=16)
        multinomial(16; c_0,...,c_4) E_z(c)^[alpha]
        * product_(s=1)^4 x_s^c_s.

Write A(n)=multinomial(2048;n_0,...,n_4). For one region, a valid averaged
operator is

    R(n) = [x_1^n_1 ... x_4^n_4] F(x)^128 / A(n).

To justify the denominator, fix U. Its sequence of packet sizes after a
uniform permutation is uniform among the A(n) size sequences. Conditional
on a size sequence, each local E_z bounds every possible remaining mask
assignment. No independence of the masks within or across steps is assumed.
The independent regional permutations then give R(n)^256.

Let W=734003 and let B_n denote the event that some nonzero message in class
n has output weight at most W. Over the ideal independent setup randomness,

    Pr[B_n] <= H(n) * (beta^r / (255*z^W))^alpha
                    * (1,0) R(n)^256 (1,1)^T.

The factor 255 accounts for projective scaling of nonzero codewords.
The bound includes the actual refresh density loss inside each local
operator. The comparison measure may assign positive mass to all-zero
potential values; deleting that mass without a new argument is invalid.

The global marker bound substitutes positive x_s into F and divides by
product_s x_s^(256*n_s) and A(n)^256. The screen also extracts regional
coefficients exactly as polynomials in a single marker for homogeneous
classes. The arithmetic used for this extraction is still floating point.

## Four-group screens and their limitation

The implementation is theory/hd_grouped_packet_screen.py. It includes every
local profile; it does not truncate occupancy tails. Representative h4
results at 35% are:

| Source counts | Active groups | Composition | Regional FLOAT log2 bound |
|---|---:|---:|---:|
| n1=1 | 1 | One step | -114.08 |
| n1=4 | 4 | One step | -483.09 |
| n4=1 | 4 | One step | -521.81 |
| n1=128 | 128 | One step | -1977.68 |
| n4=16 | 64 | One step | +652.62 |
| n4=16 | 64 | Two steps | +632.99 |
| n4=128 | 512 | One step | +6195.04 |
| n4=128 | 512 | Two steps | +4851.09 |
| n1=512 | 512 | One step | +459.35 |

For homogeneous n4 classes, every active packet has all four lanes active.
A step containing j active packets therefore has l=(j,j,j,j). The envelope
over masks is exact in these rows. The positive values cannot be attributed
to the four-count relaxation. Fifteen dispersed optimization starts reproduce
the one-step n4=16 and n4=128 values. This is evidence about the available
bound; it neither proves global optimality of the numerical search nor
provides a low-distance codeword.

Replacing h4 by two RS[32,8] blocks gives +654.51/+6261.66 at n4=16/128.
A global RS[64,16] block gives +655.57/+6306.57. These screens use the
conservative principal-nullity bound <=h in their return terms. The h8/h16
n4=16 results also survive fifteen dispersed starts. A cyclic permutation
of the four lanes cannot change the fully occupied mask.

## Two-group packets

The second route partitions the 8192 groups into 4096 pairs. Each region
independently permutes these pairs. One inner step receives 32 pairs.
Even destination pair slots feed emission blocks 0 and 1; odd slots feed
blocks 2 and 3. If every active pair is full, and j0/j1 count active even/odd
slots, the h4 occupancies are

    (j0,j0,j1,j1),  with count binomial(16,j0)*binomial(16,j1).

The marker degree is j0+j1. Regional normalization uses binomial(4096,n2).
This enumeration is exact for full pairs. For h8 the occupancies are
(2*j0,2*j1); for h16 the sole occupancy is 2*(j0+j1). The proposed maps retain
width 64 and state 16, but h8/h16 are different inner constructions.

| Active full pairs | Active groups | Inner | Composition | Regional FLOAT log2 bound |
|---:|---:|---|---:|---:|
| 32 | 64 | h4 | One step | -219.50 |
| 32 | 64 | h4 | Two steps | -248.22 |
| 256 | 512 | h4 | One step | +967.35 |
| 256 | 512 | h4 | Two steps | +83.07 |
| 256 | 512 | h8 | Two steps | +124.29 |
| 256 | 512 | h16 | Two steps | -159.09 |

The last row is a useful design lead. It is only one active-pair class;
singleton pairs, mixed occupancies, and complete coverage remain open.
The paired h16 screen gives -143.28, -271.54, and -38.44 at n2=1,32,128,
but its n2=512 search gives +430.13, including fifteen dispersed starts.
These values are not a whole-code
union bound. At h16, full-pair profiles collapse to a single potential count
2j with multiplicity binomial(32,j); the screen uses this exact simplification.
Mixed singleton/full-pair classes can likewise use two markers and potential
count c1+2c2, because the global-block bound does not depend on lane orientation.
Reproduce it from the theory directory with

    python -B hd_grouped_packet_screen.py --histogram 0 256 0 0 --regional --route-width 2 --rank 16 --steps 2 --iterations 100

Exact four-step composition resolves the sampled h16 n2=512 gap. The screen
enumerates all 33^4 sequences of step counts using two arrays of 1089 paired
operators. It multiplies the operators before taking entrywise powers and
groups only afterward by total marker degree. Per-entry log normalization
keeps the FLOAT summation scaled. At the unchanged two-step witness, the
bound becomes -843.06. A short search improves it to -985.33, at
z=0.4815867884177727 and alpha=0.0754332185491407. Reproduce the search with

    python -B hd_grouped_packet_screen.py --histogram 0 512 0 0 --regional --route-width 2 --rank 16 --steps 4 --exact-four --initial-witness 1.05402699 0.083416613 --iterations 55

This improvement identifies substantial slack from applying fractional
powers before intermediate state-path sums. A different, cheaper four-step
bound first averages all internal profiles conditional on their total degree.
Jensen makes that bound valid, but it loses the route fluctuations needed
here. Its 4/8/16-step searches remain positive at approximately
+10971/+7378/+5605. The CLI omits `--exact-four` for these averaged screens;
their metadata explicitly marks the internal profile average.

For h4, retaining all detailed half profiles across four steps is larger.
The implemented bounded alternative first averages M(j0,j1) conditional on
each step's degree j=j0+j1, then retains all four degrees until composition.
Jensen justifies the conditional average. It preserves more route information
than averaging an entire four-step block. This variant gives -200.67 at
n2=256 but +2083.07 at n2=512. Its metadata sets
`conditional_half_profile_average=true`; the h16 four-step result needs no
such average. These are distinct bounds even though they share the
`--exact-four` composition option.

Mixed h16 classes also admit an exact two-marker global screen. Its local
profile is (c1,c2), with potential count c1+2*c2 and multinomial multiplicity
for 32 source pairs. Two-step composition retains both profiles until after
the matrix product. Initial global results are:

| n1 | n2 | Active groups | FLOAT log2 bound |
|---:|---:|---:|---:|
|64|224|512|+2406.50|
|256|128|512|+818.85|
|128|448|1024|+2893.01|
|512|256|1024|-865.02|

These global bounds pay both marker-conditioning losses. The screen does
not implement regional extraction for mixed classes. A positive value is
therefore not an obstruction to a sharper regional argument.

## A possible cyclic pairing route

Another proposal chooses the parity of adjacent pairs independently in
each region. Parity zero pairs (0,1),(2,3),..., while parity one pairs
(1,2),(3,4),...,(8191,0). It then independently permutes the 4096 pairs.
This is another setup distribution and requires another bound.

For a fixed active set U in the cycle of 8192 groups, let r=|U|. Let e0 and
e1 count edges with two active endpoints in the two matchings. Matching
sigma has n2=e_sigma and n1=r-2*e_sigma. If R_sigma is a valid regional
operator for that matching and active set, the independent random parity
gives regional operator (R_0+R_1)/2. The full route bound uses its 256th power.

For 0<r<8192, let k count maximal runs of active groups around the cycle.
Then e0+e1=r-k, because each active run of length l contains l-1 active
adjacencies. The number of labeled supports with parameters r,k is

    8192/k * binomial(r-1,k-1) * binomial(8192-r-1,k-1).

One can count them by fixing a starting active run, choosing positive active
and inactive run lengths, then dividing by the k possible active starts.
For an exact joint count in r,e0,e1, define a 2-by-2 matrix indexed by bits
a,b, with T_u(a,b)=t^b*u^(a*b). The required count is

    [t^r x^e0 y^e1] trace((T_x T_y)^4096).

Both matchings can have almost r/2 full pairs only when k is small. Such
supports have much less multiplicity than arbitrary collections of full
pairs. At h4, the bound must also retain or conservatively handle singleton
orientations; at h16 only their total counts matter.

Deterministic alternation of the two matchings is another construction.
It samples independent packet permutations but fixes the matching parity
of each region. Starting with parity zero, its route operator is
(R_0 R_1)^128. This avoids the fluctuations in the number of regions using
each matching. A source support that is full under the odd matching instead
requires the reversed product. The following screens cover only supports
that are full under the even matching.

## Alternating h4 screens with exact orientation counts

Suppose U is a union of a even-matching pairs. Contract those pairs to a
cycle of length 4096, and let c count active runs in that cycle. There are
exactly

    H(a,c) = 4096/c * binomial(a-1,c-1) * binomial(4096-a-1,c-1)

such labeled supports. In the odd matching, each active run creates one
singleton of each orientation. Its packet-type counts are therefore

    (n00,n10,n01,n11) = (4096-a-c, c, c, a-c).

The orientation counts are exact, not an independence approximation.
The implementation in theory/hd_shifted_pairs_screen.py retains them.
For a half-step of sixteen pairs, its potential lane-count generating
polynomial is (1+x*s+x*t+u*s*t)^16. The two halves provide four lane counts,
to which the established h4 operator applies. Equal markers for the two
singleton orientations are a permitted witness choice; normalization still
uses both counts c separately.

At a=512, initial one-step global bounds for alternating matchings give
+3032.47, +2160.56, and -196.97 for c=128,256,448. Exact regional extraction
and two-step composition for the full even regions improve these to
+1273.39, +520.00, and -1675.43. The odd regions still use a one-step global
coefficient bound. Random matching parity is substantially weaker in this
screen: its initial global bounds are +3124.38, +3644.18, and +4480.78.

## Positive alias sums for the odd-region coefficients

The remaining global bound discards the odd region's three count constraints.
A root-of-unity filter can retain them without storing the complete degree
4096 polynomial. Let P(x,y,u) be the 2-by-2 polynomial for one odd-region
step, with coefficients indexed by its two singleton counts and full-pair
count. Its coefficients include the entrywise alpha power of the local h4
operator. The regional polynomial is P^128.

Fix positive markers m=(m10,m01,m11). For numerical scaling, let lambda and
v>0 be the Perron eigenvalue and right eigenvector of P(m), and D=diag(v).
Define

    Q(x,y,u) = lambda^(-1) D^(-1) P(m10*x,m01*y,m11*u) D.

Then Q(1,1,1) has row sums one. Every coefficient of Q is nonnegative.
For a power-of-two grid size g, evaluating Q at all triples of gth roots
of unity, taking its 128th matrix power, and applying the inverse DFT gives

    Alias_g(t) = sum_(j>=0, j congruent t mod g) [x^j] Q(x)^128.

Each entry majorizes the desired coefficient because all summed entries
are nonnegative. Undoing D and lambda, then dividing by the positive marker
monomial and the exact multinomial normalization, gives a regional operator.
This is a mathematical upper bound even when other degrees alias to t.

The current implementation performs the FFT and matrix powers in binary64
complex arithmetic. It adds a normalized 1e-12 cushion after extraction,
but that cushion is not a proved rounding-error enclosure. Small imaginary
residuals and agreement across grids are diagnostics, not certification.
It rejects target entries below 100 times the larger of 1e-12 and the
selected cushion. A retained regression fixture exhibits relative error
3.47e-5 in a tiny coefficient despite absolute error only 1.11e-19; the
diagnostic now rejects that fixture. This numerical rejection threshold
does not replace the missing rigorous error enclosure.

At the refined global witnesses above, the alias screens give:

| a | c | Grid | Alternating h4 FLOAT log2 bound |
|---:|---:|---:|---:|
|512|128|64^3|-756.59696|
|512|128|128^3|-818.24728|
|512|256|64^3|-1647.40758|
|512|448|64^3|-3793.41070|

The larger c128 grid removes additional positive aliases. Its minimum
normalized extracted entry is approximately 4.39e-6; the target imaginary
residual is below 2.48e-20. Across these cases, Perron row-sum errors are at
most 2.22e-16. None of these observations supplies a rigorous error bound.

For reproducibility, the c128 witness is

    z = 0.48180002666992056
    alpha = 0.09468173632710614
    marker_log_odds = [-2.807354922057604, -4.337800371916814, -2.4207471537600247]

Call `alias_screen(512,128,z,alpha,marker_log_odds,grid=128)` from
hd_shifted_pairs_screen.py. The even-region coefficient uses exact polynomial
extraction in FLOAT and two-step composition. The odd-region filter uses
the three source-type coordinates (n10,n01,n11).

The necessary next step is to cover support profiles beyond sets that are
full under one matching. Retain at least total weight r and the two matching
counts (a,b), with whatever singleton-orientation information their operator
requires. Independently, establish a rigorous complex-rounding enclosure
or another directed implementation of the positive alias bound. Only then
can a coverage driver attempt a new whole-code certificate.

## Route-independent dense tail and validation

The existing dense-tail argument remains valid for either packet width.
Fix the route and inner setup. Each candidate output has a unique inverse
input because the systematic inner map and the route are invertible. Its
active outer-group set is determined, so no subset union is necessary.
Outer membership charges 255^(-128r), followed by projective scaling.
This argument does not require a uniform scalar permutation.

The unchanged directed dense receipt covers r>=7665, with log2 bound at
most approximately -655.53. It is a terminal class bound, not a packet-route
certificate for the omitted group counts.

The dedicated tests enumerate all 65,536 four-position mask sequences and
compare their attainable lane counts to the Gale--Ryser envelope. They also
check exact rational h4 entries, the inverse bound for empty steps, route
normalization, cancellation of the regional marker, and equality of one-step
and two-step composition at alpha=1. The one-active-group comparison agrees
with the established scalar-route screen up to its explicit dyadic rounding.

The 26 focused tests pass with

    python -B -m unittest hd_grouped_packet_tests hd_shifted_pairs_screen_tests hd_shifted_pairs_tests -v

Next, retain the scalar-route implementation as the certified reference.
For structural work, broaden the alternating h4 route analysis to (r,a,b)
classes beyond supports full under one matching. Establish a rigorous error
enclosure for the complex arithmetic before directed replay. A new coverage
check and whole-code union bound are required before promoting this route.
