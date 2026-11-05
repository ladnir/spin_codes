# Confined states: count the outside projection of the RS outer

The joint pin transfer is correct, but its first full-size implementation
loses too much when conditioning on the number of active outer rows.
At R=192, state weight 256 and 4096 active rows, that FLOAT bound is
approximately +6041139.9 bits. A positive bound does not exhibit a bad word.

This continuation counts all messages together through an RS weight
enumerator. It gives a restricted certificate for the random-pairing
construction, without changing the construction or assuming uniform states.

## Family covered

Retain GF(2^32), RS[256,128], 8192 rows, 256 independent regional routes,
512 state coordinates, and independent accepted SL2 maps on uniform
matchings sampled at setup. The target output weight is
W=floor(0.46 * 2^21)=964689.

Choose R trailing regions. The state is zero before this suffix. At its
first update it acquires a support H, which is then exactly the support
of every remaining state, including the free terminal state. Thus this
family has one activation, no support changes, and no final reset.
Write w=|H| and L=16R for the suffix length in updates.

The new directed result covers R=192,w=256 and R=256,w=256. It includes
every H of that size and every message with such a trajectory. Their
combined bad-word probability is at most 2^-95 over setup. It does not
cover changing supports, repeated episodes, or the whole code. It also
does not establish an extension-field certificate over GF(2^128).

## A pairing must reach the precise nonzero inputs outside H

Fix a message and route, and expose previous inner maps. Suppose the
incoming state has support H and the next state must have support H.
Outside H, the update requires (A Q)_i=X_i. A crossing pair has precisely
one nonzero input, and superregularity makes both its image coordinates
nonzero. An empty pair has zero image. Therefore the inactive crossing
endpoints must equal the support of X outside H, not merely have its size.

Let P_w(k) be the crossing-count law in [RANDOM_PAIRING.md](RANDOM_PAIRING.md).
For any specified k-subset of the inactive coordinates, its probability
of being exactly the set of crossing endpoints is

    P_w(k) / binom(512-w,k).

Let m_w be the maximum of this expression over feasible k. Conditional
on a matching, each specified outside value costs a field equation.
Requiring selected outputs inside H to be zero adds equations on disjoint
coordinates of live pairs. The accepted-map density loss is paid once,
globally, as a factor less than 2.

For w=256, m_w is approximately 2^-251.17355. Its maximum occurs when
all 256 inactive coordinates cross, not at the zero-crossing term that
dominated the earlier marginal inverse moment. For w=2, m_w=1/511.

Let B(x) count nonzero inputs outside H in the suffix. The initial input
must itself have support H; discard messages that violate this condition.
For all other messages, its outside contribution at that step is zero.
There are (L-1)w positions where an output can vanish after activation.
Low total output weight requires at least z=max(0,Lw-W) such zeros.
A union over z selected positions gives the bound

    2 * m_w^(L-1) * binom((L-1)w,z) * q^(-z-B(x)).

It is valid to sum q^-B(x) over all shortened outer words, including
words excluded at activation and the zero word: that only enlarges it.
This is a fixed-message argument followed by a union over messages, so
messages chosen after setup are included.

## Sum the outer words before maximizing their supports

The zero prefix shortens each row to an RS[R,k] code, where k=R-128.
No nonzero row remains for R<=128. For R>128, let b_i count positions
of row i routed into H in these regions, and put m_i=R-b_i.
The row's projection outside H has length m_i.

When m_i>=k, this projection is an RS[m_i,k] code. Its nonzero weight-j
coefficient is at most binom(m_i,j) q^(k-m_i+j), and its minimum weight
is m_i-k+1. This follows by fixing the zero positions and using MDS rank.
Consequently its enumerator at q^-1 is at most 1+2^m_i q^(k-m_i).
When m_i<k, projection is surjective with kernel size q^(k-m_i), and
the exact enumerator is q^(k-m_i)(2-1/q)^m_i. The same upper bound holds.

Thus, for the fixed route,

    sum_x q^-B(x) <= product_i U_R(b_i),
    U_R(b) = 1 + 2^(R-4096+31b).

This includes all active-row counts. It keeps the RS cost of zero
evaluations instead of replacing each row by its minimum distance alone.

## A simultaneous good-route event

For fixed H, each region sends a uniform set of 16w rows into H.
Regions are independent. Hence each b_i has distribution Bin(R,w/512),
but the b_i are not independent across rows.

Uniform subsets are negatively associated. Independent unions preserve
that property, as do increasing functions on disjoint coordinate sets.
These standard properties justify the product upper bound on the moment
of row functions; see [Joag-Dev and Proschan (1983)](https://doi.org/10.1214/aos/1176346079).
The marginal distributions alone would not justify the product bound.

For the directed replay, fix an integer c and rational r>1. Put
D_c=sum_i (b_i-c)_+. With B distributed as Bin(R,w/512),

    Pr[D_c>d] <= r^-d * (E[r^((B-c)_+)])^8192.

Also, log2 U_R(b) is increasing with slope at most 31. Therefore

    sum_i log2 U_R(b_i)
      <= 8192 log2 U_R(c) + 31 D_c.

Choose d so the displayed probability is at most
2^[-80-17-ceil(log2 binom(512,w))]. Union over all masks H leaves at
most 2^-97 for this R,w. The 17-bit reserve can cover later unions over
suffix lengths and weights; it is not evidence that those profiles pass.

On this good-route event, the log2 bad-word contribution is at most

    1 + log2 binom(512,w) + (L-1) log2 m_w
      + log2 binom((L-1)w,z) - 32z
      + 8192 log2 U_R(c) + 31d.

The route event is uniform over H before any message is chosen. All
random-map probabilities are then bounded conditional on a good route.

## Directed results and checks

`certify_random_pairing_confined.py` uses integer powers for directed
logarithms and exact rational binomial moments. No FLOAT proposal is
trusted by this replay.

| R | w | c | r | Excess cap d | Good-route log2 upper |
| --- | --- | --- | --- | --- | --- |
| 192 | 256 | 120 | 9/8 | 4282 | -638009 |
| 256 | 256 | 124 | 33/32 | 91040 | -434936 |

Each row has route failure at most 2^-97 after the mask union, and total
restricted failure at most 2^-96. Union of the two rows gives 2^-95.
The replay takes about 0.49 seconds and writes no files.

The separate FLOAT screen passes all 512 support weights for R=192 and
R=256. Its worst good-route bounds are about -27384 and -36369 bits,
both at w=2. Those 1024 profiles have not received directed replay.
A 30-profile screen across five suffix lengths also passes. No screen
is counted as certified coverage.

`test_random_pairing_projection.py` checks 64 punctured RS enumerators
over GF(8), 12 exact regional moment inequalities, and 540 one-step
constraints using the actual accepted SL2 distribution. These are
independent finite checks of the counting interfaces, not substitutes
for the arguments above.

## Extension to changing supports

Fix incoming support H and next support H'. Let b=|H' minus H| count
births, d=|H minus H'| count deaths, and U=H union H'. Only coordinates
outside U have both state values zero. On these coordinates, the crossing
endpoints must equal the nonzero input positions. Crossing choices at the
b birth coordinates are unrestricted by this necessary condition.

For a fixed required set of a crossing endpoints outside U, its probability is

    sum_j binom(b,j) P_w(a+j) / binom(512-w,a+j)
      <= min(1, 2^b m_w),

with infeasible terms omitted. Death coordinates impose d additional
field equations, disjoint from the outside-U equations and selected
zero-output equations on H intersect H'. Thus the local forward bound is

    density * min(1,2^b m_w)
      * q^[-wt(X outside U)-d-number_of_selected_zero_outputs].

The same test checks 3510 changing-support constraints over GF(8), including
complete resets. All pass. This is a local bound, not a composed certificate.

Next compose it with the RS outside-projection enumerator, paying for the
choices of changing supports jointly with births and deaths. The required
route event now depends on unions H_t union H_(t+1), rather than one mask.
Repeated activations and reset boundaries must also be included. The
constant-support certificate is preserved as a control. Native performance
work remains deferred until this broader proof model is credible.
