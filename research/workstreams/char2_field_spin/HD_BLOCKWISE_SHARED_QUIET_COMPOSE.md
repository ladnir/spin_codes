# Quiet composition through regional value types

The shared beta64 construction admits an all-density bound for quiet words
whose inverse images have sufficiently informative regional value types.
This bound uses actual outer-code membership, not just outer support sizes.
It leaves a precise low-information complement and does not certify the
whole code.

At the current K20 geometry, the following class has failure probability
below `2^-80.8712`: the sum of the largest 128 regional type-information
values is at least **7,693,904 bits**. The probability is over the 256
independent regional permutations. No randomness of message amplitudes or
state mixers is needed for this quiet-class statement.

## Model and the base-field reduction

Fix outer GRS[256,128], R=8192 rows, N=256R=2^21 symbols, packet width 64,
and state dimension 4. Each coordinate region receives an independent uniform
permutation of its R rows. The shared feedthrough uses beta64 and eight
FFT8 tiles per packet. Let F denote its packet feedthrough and C its
four-row feedback map. These fixed maps are over GF256, embedded in GF65536.
The facts used below are

    F^2=I,       CF=TC for an invertible T,
    ker C is an MDS[64,60,5] code over GF256.

There are P=N/64=32768 packets. A word is quiet if every emitted state is
zero. Initial state is zero, and the terminal state is free. Thus quiet inputs
satisfy C X_i=0 for i<P-1, with no constraint on the final packet. The quiet
output is Y_i=F X_i. Equivalently, C Y_i=0 for i<P-1.

For every fixed route, the minimum symbol distance of this quiet code is
unchanged by extension from GF256 to GF65536. To see this, decompose an
extension-field word in a GF256 basis. Both components satisfy the outer and
quiet equations, and each component's output support lies within the original
output support. At least one component is nonzero. The converse follows by
embedding base-field words. This argument, established in
[HD_BLOCKWISE_SHARED_CUTS.md](HD_BLOCKWISE_SHARED_CUTS.md), lets us partition
**base-field witnesses** when bounding quiet failure.

The type class below is a class of those GF256 witnesses. Descent need not
preserve an extension-field word's type class. This causes no issue when
the high-information and low-information base-field classes are eventually
bounded together. Nonquiet words are not covered by this reduction.

## A finite quiet-output ball bound

Write q=256 and e=838860. Then weight at most e is exactly symbol distance
strictly below 40% at N=2^21. For one constrained packet, let

    A(z)=sum_(w=0)^64 binom(64,w) M_w z^w,

where M_w counts codewords with a prescribed exact support of size w in
MDS[64,60,5]. Explicitly,

    M_0=1,
    M_w=0                         for 1<=w<=4,
    M_w=(q-1) sum_(j=0)^(w-5) (-1)^j binom(w-1,j) q^(w-5-j)
                                  for 5<=w<=64.

The unconstrained final packet has weight polynomial B(z)=(1+255z)^64.
For 0<z<1, the number of nonzero quiet-output projective lines of weight at
most e is bounded by

    L_e <= z^(-e) A(z)^(P-1) B(z) / 255.

Each nonzero scalar orbit has 255 elements and preserves symbol weight.
Counting the zero vector in the polynomial only enlarges this bound.
The same orbit reduction applies to regional type information, because
nonzero scaling only relabels the 256 value multiplicities.

Directed 160-bit arithmetic at the fixed rational mark z=2/765 gives

    log2 L_e <= 7693823.1287572013310355600912018628354165293620
                with numerical enclosure radius below 8.40e-41.

This is an upper bound on a superset of actual quiet outputs: outer membership
has not yet been imposed. The final-packet correction is included explicitly.

## Exact regional type information

Fix a nonzero GF256 output Y in this quiet ball before sampling the route.
Its inverse input X=FY is then fixed. For region j and value a in GF256,
let m_j(a) be the number of coordinates of X_j equal to a. Define

    H_j(X)=log2(R!) - sum_(a in GF256) log2(m_j(a)!).

A uniformly permuted X_j is uniform on its multiset orbit. Its probability
of equalling any particular ordered vector in that orbit is exactly
`2^-H_j(X)`. It has probability zero of equalling a vector of another type.

Let T(X) be the sum of the largest 128 values among H_1(X),...,H_256(X).
Ties can be broken by region index. The fixed-word membership bound is

    Pr_route[X belongs to the routed outer code] <= 2^-T(X).

To prove the bound, select the 128 largest-information regions as parity
regions. Expose the permutations of the other 128 regions first. Any 128
coordinates of each GRS[256,128] row form an information set, so this exposure
uniquely determines all required parity values. The remaining regional
permutations are independent. Each parity match costs either zero probability
or its exact inverse orbit size. Their product is at most `2^-T(X)`.

The parity selection depends only on fixed X, not on sampled permutations.
This proof does not condition amplitudes to be independent and does not
replace regional routing by a global permutation.

For any threshold h, a union bound over projective quiet outputs therefore
gives

    Pr_route[there exists a nonzero quiet Y with wt(Y)<=e and T(FY)>=h]
        <= L_e 2^-h.

Two directly verified integer choices are:

| Required T(FY), in bits | Failure upper bound |
| --- | --- |
| 7,693,864 | less than 2^-40.8712427986 |
| 7,693,904 | less than 2^-80.8712427986 |

The second threshold averages `480869/65536`, about 7.33748 bits per position,
across the selected 128 regions. This is the actual class condition. In
particular, the remaining class has **low top-half type information**, not
merely one unusually frequent field value.

## A simple sufficient type condition

Suppose at least 128 regions have every GF256 value multiplicity at most 47.
The least informative type under this cap has 174 values of multiplicity 47
and one value of multiplicity 14. Convexity of log factorial gives this
extremal type: moving one count from a smaller nonempty bin to a larger
nonfull bin cannot decrease the denominator of the orbit size.

Consequently each such region contributes at least

    log2(8192! / (47!^174 14!))
        = 60307.482679861565881285036205756707863423013163 bits.

For this sufficient class, the directed failure bound is below

    2^-25534.654265079101768924543134995771101616323.

The multiplicity cap is only a convenient corollary. Replacing 47 by 48 in
this worst-type argument produces an uninformative upper bound above one;
it does not establish failure or negate the actual T(FY) criterion. A column
with maximum multiplicity 48 can still have much greater type information
than this worst-type envelope.

The finite alphabet matters here. R=8192 exceeds q=256, so an all-distinct
column is impossible. The maximum-information type has all 256 counts equal
to 32. No bound in this note uses the impossible all-distinct value R! as
the size of an actual regional orbit.

## Why exact outer values remain necessary

There is a checked local obstruction to an all-density proof that retains
only outer supports and packet feedback ranks. Choose one output coordinate
in each FFT8 tile:

    S={0,8,16,24,32,40,48,56}.

The restriction of C to S has rank 4. Its kernel is MDS[8,4,5], and its
full-support vectors number

    256^4 - 8*256^3 + 28*256^2 - 56*256 + 35 = 4162570275.

The all-ones vector on S is one explicit quiet output. For beta64 its inverse
input repeats the following eight nonzero values across the packet:

    (16,105,170,110,165,217,144,80).

Thus output weight is 8, input weight is 64, and both feedback syndromes are
zero. The helper also checks this phenomenon for beta8.

Repeating the packet witness yields output weight N/8 and full input support.
After any regional route every source row has all 256 coordinates present.
All support-only outer minimum-weight tests therefore pass. The actual outer
**value equations have not been imposed**, so this is not a codeword or a
construction counterexample. It shows why local ranks and support counts
alone do not settle the remaining dense quiet class.

## The remaining interface

The open quiet class consists of low-weight inverse outputs with
T(FY)<7,693,904. Sparse-row bounds from the cuts note already cover a portion
of that class. The uncapped quiet first-moment obstruction from that note
also remains: its rare route contributes more than 2^4110 expected
projective quiet lines. The new high-information theorem does not transfer
that expectation into a failure probability.

A type-sensitive continuation can preserve more than the information-set
upper bound. For a fixed tuple of regional multiplicities m, let A(m) count
ordered R-row outer-code matrices with exactly those column types. The exact
fixed-input membership probability is

    A(m) / product_j (R! / product_a m_j(a)!).

The numerator is the corresponding coefficient of the R-th power of the
outer code's complete weight enumerator. This identity exposes the needed
joint value consistency; it does not yet give a scalable computation.

The next proof step should bound the low-T inverse-output classes jointly
with this outer consistency, and cap rare-route failure events before
counting all words in a large kernel. Uniform outer messages do give
128-wise independent uniform coordinate columns by the MDS property. That
fact cannot be reused as conditional amplitude independence after fixing
quiet constraints, low output weight, or inverse types.

No bound on this low-information complement or on nonquiet histories is
claimed here. Therefore these results are not a full 40%/40-bit certificate.

## Reproduction

The new helper is
[theory/hd_blockwise_shared_quiet_compose.py](theory/hd_blockwise_shared_quiet_compose.py).
Its six tests check the exact small-field feedback enumerator, the extremal
capped type, and all 4096 inverse words for a two-row GF4 GRS[3,2] code under
all eight regional routes. They also check both local witnesses and the
directed positive and uninformative type receipts.

From the theory directory:

    python -B -m unittest hd_blockwise_shared_quiet_compose_tests -v
    python -B hd_blockwise_shared_quiet_compose.py

No native implementation, benchmark, remote run, or existing proof file was
changed for these results.
