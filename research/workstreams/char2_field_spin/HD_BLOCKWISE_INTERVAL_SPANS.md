# Complete-row interval traps for outer256, width16, state4

At the new GRS[256,128] geometry, the dimension-forced complete-row connected-window family has
route probability below 2^-528.39821182. This bound uses exact order statistics
instead of a union over window endpoints. It is uniform over the inner mixers.

The same geometry still has an actual-code first moment above 2^740.75317482
for low-weight projective lines. These statements are consistent: a rare route
can support a large bad subspace. The capped route bound counts that route once.
Neither result is a whole-code distance certificate.

## Construction and event

There are R=8192 independent GRS[256,128] rows over GF(65536), with K=2^20 and
N=2^21. Each of the b=256 coordinate regions has an independent uniform row
permutation. The inner has packet width t=16 and state dimension s=4:

    Y_i = X_i + A_i Q_i,
    Q_(i+1) = M_i Q_i + C_i X_i,     Q_0 = 0.

The terminal state is free. The route and inner setup are independent. All
claims below hold for every fixed realization of the causal linear inner.
The identity feedthrough makes the emitted map invertible: each packet input
is recovered from its output and the state determined by earlier packets.

Define B_full as the following route event. There exist a nonempty row set J
of cardinality a and one nonempty packet-aligned interval in every region such
that all b*a coordinates of J lie in these intervals. Their total length is at
most T=838848, the largest packet multiple strictly below 0.4N. Let c be the
number of connected components in their union along the global stream. Set
nu=1 if the last component reaches the terminal boundary, and nu=0 otherwise.
Require

    128a > 4(c-nu).

On B_full, use all message dimensions of the rows in J and impose a zero-state
exit after each nonterminal component. At most 4(c-nu) equations leave a nonzero
message. Its emitted support stays inside the intervals, so it violates 40%
distance. This is a sufficient failure mechanism, not a characterization.

The containment condition includes every coordinate position of J, even when
a particular outer word is zero there. A low-weight word may instead exploit
its zero entries or partially shortened source rows. Those cases are not
covered. Neither are multiple intervals per region or failures caused by an
exit system whose rank is lower than this dimension-forcing threshold.

## Why an uncapped first moment still fails

Fix ten source rows. Choose 52 regional windows of length 3264 and 204 of
length 3280. Use prefixes in the first 255 regions and a suffix in the last.
The total is T=838848. For this fixed row set and these fixed windows, exact
sampling without replacement gives the containment event E_J probability

    Pr[E_J] = product_j C(W_j,10)/C(8192,10),
    log2 Pr[E_J] = -3387.246825176145718298... .

On E_J, zero-state exits cost at most 255*4 equations. Thus a message subspace
of dimension at least d=128*10-255*4=260 emits only within these windows.
Let Z count scalar projective code lines of weight at most 838860. Then

    E[Z] >= Pr[E_J] * (65536^260-1)/65535
         > Pr[E_J] * 65536^258
         > 2^740.753174823854281701... .

No union over row sets is used, and all outer parity checks are retained.
An upper certificate based on the uncapped first moment of actual code lines
therefore cannot close. This does not show a substantial failure probability:
the displayed route event itself has probability about 2^-3387.

## Count component patterns, not window endpoints

For a fixed connected pattern, let v count its singleton components. There are
c-v components spanning at least two regions. All b regions are used. For
v<c, the number of ordered component-size patterns with specified singleton
slots is

    H_0(c,v) = C(b-c-1,c-v-1).

This follows by giving each long component two regions, then distributing the
remaining regions among those components. The all-singleton case is possible
only when c=v=b, and has H_0=1.

Every interior region of a long component must be full. Its first region is
a suffix and its last is a prefix. A singleton not fixed to the terminal
boundary can use an arbitrary interval. The following table gives the number
H of component patterns, the number f of forced-full regions, the number h of
anchored variable regions, and the number z of free-span regions.

| Terminal treatment | H/H_0 | f | h | z | nu |
| --- | --- | --- | --- | --- | ---: |
| Charged exit | C(c,v) | b-2c+v | 2(c-v) | v | 0 |
| Free final singleton | C(c-1,v-1) | b-2c+v | 2(c-v)+1 | v-1 | 1 |
| Free final long component | C(c-1,v) | b-2c+v+1 | 2(c-v)-1 | v | 1 |

Invalid binomial coefficients or negative counts contribute nothing. Patterns
whose forced-full regions and at least one packet per variable region exceed
T are also omitted. The charged case can overcount some terminal-reaching
windows; retaining that overcount is safe for the upper bound.

For each pattern, it suffices to consider exactly

    a_0 = floor(4(c-nu)/128)+1

rows. Every larger J witnessing B_full contains a subset of this cardinality
that still satisfies the dimensional inequality and the same containment
condition. Thus the source-set union costs C(R,a_0), with no sum over larger
cardinalities. The feasible patterns here use only a_0=3,...,9.

## Exact regional length laws

Fix a row set of size a. In each region its positions form a uniform a-subset
of {0,...,R-1}, independently across regions. Let L_A be the smallest prefix
length containing those positions. Reflection gives the same law for a suffix.
Let L_S be the length of the smallest arbitrary inclusive interval containing
them. For a>=2,

    Pr[L_A=l] = C(l-1,a-1)/C(R,a),
    Pr[L_S=l] = (R-l+1) C(l-2,a-2)/C(R,a),     a<=l<=R.

The first formula fixes the last selected position. The second fixes both
extreme positions and chooses the remaining a-2 positions between them. When
a=1, L_S=1 deterministically. These minimum lengths need not be packet-aligned;
using them relaxes the containment event and therefore gives an upper bound.

For a mark u in (0,1), define the exact polynomial moments

    A_a(u) = E[u^L_A],     S_a(u) = E[u^L_S].

If windows of the specified connected pattern exist with total length at most
T, then their regional minimum lengths satisfy

    fR + sum_(h anchored regions) L_A + sum_(z free regions) L_S <= T.

Markov's inequality and regional independence give, for a fixed J and pattern,

    Pr[containment by such windows]
      <= u^(fR-T) A_a(u)^h S_a(u)^z.

No separate window-start or endpoint factor appears. The minimum spans already
account for every interval that could contain the selected positions.

Combining the component-pattern and source-set unions proves

    Pr[B_full]
      <= sum_(feasible c,v,terminal)
           H * C(R,a_0) * u^(fR-T) * A_(a_0)(u)^h * S_(a_0)(u)^z.

## Directed receipt and scope

At the fixed rational mark u=1023/1024, 160-bit Arb interval arithmetic sums
31,469 pattern groups and obtains

    log2(upper bound)
      = [-528.39821182502750854061177814156802479209160
         +/- 5.13e-42].

Consequently Pr[B_full] < 2^-528.39821182. The largest group has c=160,
v=89, a free final singleton, f=25, h=143, z=88, and a_0=5.

This receipt is specific to outer[256,128], t=16, s=4. It does not carry over
the old width-64/state-16 number. Indeed, directly rerunning the older endpoint
union at the new geometry gives an uninformative log2 upper of +2005.6111.

For the full distance event Bad, a valid next decomposition is

    Pr[Bad] <= Pr[B_full] + Pr[Bad and not B_full].

The second term remains open, including partial-row shortening and arbitrary
output histories. Conditioning on not B_full leaves independent mixer setup
unchanged, but changes the regional routing law. A later proof must retain
that conditioning or use a valid dominating unconditioned quantity.

## Reproduction

The new helper is [theory/hd_blockwise_interval_spans.py](theory/hd_blockwise_interval_spans.py),
with [three tests](theory/hd_blockwise_interval_spans_tests.py). They exhaustively
verify both order-statistic distributions for R<=8, replay the actual-code
first moment, and check the directed capped receipt. The helper reuses the
previously tested component-pattern enumeration from hd_blockwise_global.py;
its exhaustive small-window count test was also rerun.

From the theory directory:

    python -B -m unittest hd_blockwise_interval_spans_tests
    python -B hd_blockwise_interval_spans.py

The combined run of these three tests, the eleven quiet-subcode tests, and the
small component-pattern enumeration test passes all fifteen tests. No native
files, remote setup, or performance benchmarks were touched.
