# Cap the conditional count before averaging pairings

The conserved route bound still failed at 4096 swaps. Its first-moment
count can be large when many outer messages share the same rare pairing
event. Multiplying that event probability by the number of messages loses
the fact that its conditional failure probability cannot exceed one.

Keeping this cap gives a different proof of the same construction. It
certifies every half-full-state history with at most 8192 swaps in either
of the two retained suffixes, at 46% distance. Their joint contribution is
at most 2^-119. This remains restricted coverage, not a whole-code or
GF(2^128) extension certificate. No implementation was changed.

## Condition on the pairing geometry

Use the setup and notation from [the turnover proof](RANDOM_PAIRING_TURNOVER.md).
The suffix contains R regions and L=16R updates. After its initial
activation the state has weight w=256 at every step. There are B births
and B deaths across its T=L-1 subsequent transitions. The terminal state
is free, subject to this family's weight condition.

Fix a complete support history before sampling setup. Let U_t be the union
of its incoming and outgoing supports. Enlarge U_0 to the whole initial
batch, discarding some constraints on the initial input. Write

    A = sum_t |U_t| = Lw+B+(512-w).

Fix the route and matchings. Let J_t contain inactive incoming coordinates
outside U_t whose partner is active. Superregularity implies that X is
nonzero at these positions if the history occurs, and zero outside U_t
union J_t. Let J be their union across time.

Let a_i count the positions of RS row i in U union J. Shortening the row
to those positions leaves dimension (a_i-128)_+. Thus the candidate message
space has dimension d=sum_i(a_i-128)_+. Ignoring the required nonzero
inputs on J only enlarges this space.

Fix z selected zero-output positions in the intersections of successive
supports. Each fixed message incurs |J|+B+z field equations: the outside
crossing positions, deaths, and selected output zeros. These coordinate
sets are disjoint within each update. The existing accepted-map point
bound applies causally, with its single global density factor less than 2.

For this route and geometry, the probability of a bad message with the
specified history is therefore at most

    min(1, 2 C q^[d-|J|-B-z]),

where C bounds the number of selected zero-output sets. The outer count
is a union over fixed messages; the cap is applied after that union.
It is not a claim that a first moment equals a failure probability.

Define the parity deficit D=sum_i(128-a_i)_+. Since sum_i a_i=A+|J|,

    d-|J| = A-128*8192+D.

This identity retains the cost of outside crossings rather than replacing
all messages by the dimension of the outer shortened to U alone.

## Compare crossing endpoints with independent coins

For an incoming support of weight w, the probability of any specified
external crossing set of size k is P_w(k)/binom(512-w,k). Define

    m_w = max_k P_w(k)/binom(512-w,k),
    rho_w = 2^(512-w) m_w.

Its full crossing-set distribution is pointwise at most rho_w times
independent fair coins on the inactive coordinates. Marginalizing the
birth coordinates preserves this density bound on J_t. At w=256,
log2 rho_w is about 4.82645 bits per update, not one field factor per
coordinate. Independent matchings give a total factor rho_w^T.

This comparison concerns only support geometry. It does not replace the
accepted SL2 arithmetic by independent bits. The arithmetic density factor
remains inside the capped conditional bound.

## Average routing without a simultaneous exceptional-route event

Under the coin comparison, each position outside U joins J independently
with probability one half. For a fixed history, region r contains a fixed
number a_r of U positions. Represent their routed row labels by independent
Bernoulli indicators of parameter p_r=a_r/8192, then condition on their
exact column total. This is the actual uniform regional routing law.
The conditioning loss is at most (8193)^R, as in the previous proof.

Before conditioning, different rows are independent. A row is covered
in region r with probability (1+p_r)/2. Consequently its covered count
is a sum of independent Bernoulli variables. For r>1, the function
r^((128-c)_+) is convex and decreasing. Bernoulli smoothing bounds its
expectation by the equal-mean binomial law. Monotonicity also allows a
lower coverage parameter to give a uniform upper bound.

For 0<theta<=1, min(1,x)<=x^theta. Set r=q^theta and apply this inequality
to the conditional count before averaging. Let H count the support
histories under consideration, and let A_min bound their area below.
If e bounds A-128*8192-B-z above, the resulting bound is

    H rho_w^T (8193)^R
      * r^[e+(1+log2 C)/32]
      * (E[r^((128-Z)_+)])^8192,

where Z is Bin(R,(1+A_min/(8192R))/2). Here q=2^32.
There is no union over a separately defined bad-route event. In particular,
the history count is not paid again inside a route-tail budget.

All histories are fixed when applying the distributional comparison, then
unioned. This covers messages and histories chosen after seeing setup.

## Uniform bound for B at most 8192

Set B_max=8192. As before,

    H <= binom(512,256) S_256(x)^T / x^B_max,
    S_256(x)=sum_b binom(256,b)^2 x^b,
    x=1/8192.

Keep only z_0=max(0,Lw-W) output zeros, with W=964689. Every history in
the family needs at least this many zeros, and its intersection has at
most Tw coordinates. Hence C=binom(Tw,z_0) is a valid uniform choice.

Take A_min=Lw+256. Crucially, the extra union area B cancels the B death
equations before averaging. Thus e=Lw+256-128*8192-z_0, independent of B.
Larger B only improves the coverage probability used in the deficit moment.

`certify_random_pairing_capped.py` verifies the following using exact
rational moments and integer-power logarithm bounds:

| R | History-count log2 upper | r | Directed log2 contribution |
| --- | --- | --- | --- |
| 192 | 124214 | 3 | -182846 |
| 256 | 129942 | 8 | -57995 |

The conservative quoted bound is 2^-120 per family, hence 2^-119 jointly.
The replay includes all masks, placements of up to 8192 swaps, and messages
in these two families. It supersedes their 256-swap coverage but does not
invalidate or remove the prior proof.

## Numerical diagnosis beyond the certificate

The exact-B FLOAT screen at B=4096 gives -244656.9 and -137794.4 bits for
R=192 and 256. Previously these bounds were positive. At B=16384 the new
screen still passes, at about -103937 and -14454 bits. At B=32768 it fails,
at about +41208 and +106989. Intermediate larger B values also remain open;
very large turnover again benefits from its forced nonzero outputs.

These selected screens are not directed coverage. In particular, the
proof does not yet cover changing state weights or repeated resets.

## Checks, execution, and next step

`test_random_pairing_capped.py` passes 652 exact endpoint-domination checks,
280 capped-average checks, 295 deficit identities, and 30 exact comparisons
of regional coverage moments. The arithmetic point bounds and local
birth/death tests remain the preserved controls.

Local free space fell from23.97GiB to4.15GiB, so the local guard prevented
the numerical run. A read-only check found a44.35GiB allocated pagefile,
but did not establish the cause of the space change. No files were deleted
and no system settings changed.

The configured Peach connection was available with63GiB free in /tmp and
no numerical job under the account. It lacked SciPy, so the new helpers
use only the Python standard library. Source was streamed through stdin;
no remote files or packages were created. Each job had a30s timeout,
256MiB address-space limit, Python -B, and its own20GiB free-space check.
Runs were serial: screen0.10s, directed replay0.28s, exact checks0.10s.

Next examine the active-outer-row count jointly with this cap. The present
dimension count includes many words in the same support subspace. For an
exact active row set S, every J position must belong to S and every active
row must have at least129 covered positions. Retaining these constraints
may improve the intermediate-turnover range. Test that interface before
expanding to varying weights and reset boundaries. No native timing yet.
