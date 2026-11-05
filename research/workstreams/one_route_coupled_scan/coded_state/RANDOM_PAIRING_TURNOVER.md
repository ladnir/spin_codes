# Composing support changes at fixed weight

This continuation keeps the construction fixed and extends the
[outside-projection argument](RANDOM_PAIRING_PROJECTION.md). The support
may change at every step, but its size remains w after activation.
This family does not include changing state weights or complete resets.

The initial run was prevented by the 20 GiB disk guard. Space subsequently
returned to 25.33 GiB, and bounded serial runs resumed. The conserved-route
continuation below now certifies all histories with at most 256 swaps,
for w=256 and R=192 or 256. Their combined bad-word contribution at 46%
distance is at most 2^-118. Moderate-turnover and whole-code closure remain open.

## Count changes and output constraints together

Fix a region-aligned suffix of R>128 regions and L=16R updates. The state
is zero before this suffix; every state Q_1,...,Q_L has weight w. Put
T=L-1. At transition t between these nonzero states, let b_t count births.
There are also b_t deaths, since the weight stays fixed. Let B=sum_t b_t.

The exact number of support histories in this class is

    H(R,w,B) = binom(512,w) [x^B] S_w(x)^T,
    S_w(x) = sum_b binom(w,b) binom(512-w,b) x^b.

Indeed, from any support of size w, choose b coordinates to remove and b
inactive coordinates to add. This count is independent of its labels.
For any x>0, S_w(x)^T/x^B gives an upper bound on the coefficient.
The FLOAT helper optimizes this positive marker; directed replay would
need rational markers and an independently verified upper bound.

At each update, outputs vanish outside the union of its two state
supports. The total size of these unions, including initial activation,
is A=Lw+B. The total intersection size after activation is I=Tw-B.
Birth and death coordinates have nonzero output. Hence an output of
weight at most W requires at least z=max(0,A-W) zeros in the intersections.
If z>I, the class is impossible. Otherwise choose z of the I positions.

The earlier local lemma gives at most 2^b_t m_w for the matching constraint
and q^-b_t for the deaths. The selected output zeros impose additional
equations on disjoint live coordinates. For a fixed history and route,
the sum over messages is therefore at most

    2 H_free * m_w^T * 2^B * q^(-B-z) * binom(I,z),

where H_free denotes the outside-projection enumerator of the outer for
this history. The factor 2 is the global accepted-map density loss.
There is no factor H(R,w,B) yet: this formula concerns one history.

## Account for time-varying union masks in the route

Fix a history before sampling the route. In region r, its union masks
select a deterministic number a_r of physical positions. Those positions
receive an independent uniform a_r-subset of the 8192 source rows.
For row i, let c_i count selected positions across regions. The RS bound is

    H_free <= product_i U_R(c_i),
    U_R(c) = 1+2^(R-4096+31c).

The regional probabilities p_r=a_r/8192 need not be equal, but their sum
is A/8192. Each c_i is a sum of independent Bernoulli variables with those
probabilities. The different row counts remain negatively associated.

For lambda>0, f(c)=U_R(c)^lambda is convex and increasing. Replacing two
Bernoulli probabilities by their average increases E[f(sum)] at fixed sum
of probabilities: the difference is a nonnegative multiple of the second
finite difference of f. Repeated averaging gives the binomial envelope
with p=A/(8192R). Negative association then gives

    E[2^(lambda sum_i log2 U_R(c_i))]
       <= (E[U_R(C)^lambda])^8192,    C distributed as Bin(R,p).

Let h be an upper bound on log2 H(R,w,B). A cap V satisfying

    -lambda V + 8192 log2 E[U_R(C)^lambda] <= -(120+h)

fails for any history in the class with probability at most 2^-120 after
the history union. This event is simultaneous over histories; none is
assumed independent of setup after a bad message has been chosen.

On the good-route event, the log2 bad-word contribution is bounded by

    1+h+T log2 m_w+B-32(B+z)+log2 binom(I,z)+V.

This includes both the history choices and the route cost of making the
bound simultaneous. The 40-bit reserve above 80 is for possible later
unions over R,w,B, not a claim that those classes have been checked.

## Preserve the total occupancy in the route moment

The first screen passes small B but fails at moderate turnover. For
R=256,w=256,B=4096, its log2 upper bound is about +4700129. The simultaneous
route cap contributes about 8165263 bits. This is not an exhibited codeword.

One lost constraint is sum_i c_i=A, which holds exactly for every route.
The negative-association product bound permits all rows to gain occupancy
together. A scalar tilt can retain the constraint, but applying a decreasing
factor inside negative association would not be justified.

Instead, start with independent Bernoulli membership entries, with
probability p_r=a_r/8192 in region r. Condition each column to contain
exactly a_r selected rows. This is the actual routing distribution.
The conditioning probability is at least (8193)^(-R): the mass at the
integer mean of Bin(8192,a_r/8192) is a mode, hence at least 1/8193.

Before conditioning, rows are independent. For any lambda>0 and tau,
insert 2^[-tau(sum_i c_i-A)], which equals one on the conditioning event.
The row function U_R(c)^lambda 2^(-tau c) is convex, even though it need
not be increasing. Thus Bernoulli smoothing still applies, giving

    E_actual[2^(lambda sum_i log2 U_R(c_i))]
      <= (8193)^R 2^(tau A)
         (E[U_R(C)^lambda 2^(-tau C)])^8192,
    C distributed as Bin(R,A/(8192R)).

This argument uses independence before conditioning, not a false
negative-association claim for nonmonotone functions.

For R=256,w=256,B=4096 the resulting FLOAT cap falls to about 4397327 bits;
the final bound improves to +932193. It still fails. At B=450000 the old
bound fails at +3897867, whereas the conserved bound passes at -5728231.
For R=192,B=4096 the bound improves from +381995 to +236648 and remains open.
Both versions remain available in the diagnostic helper.

## Directed certificate for at most 256 swaps

Let B_max=256. For a rational 0<x<=1, the number of histories with
B<=B_max is at most binom(512,w) S_w(x)^T / x^B_max. Use x=2^-20.
To give every history the same occupancy total, pad its union masks to
A_max=Lw+B_max positions, using a deterministic rule independent of routing.
Padding can only enlarge the outside-projection enumerator bound.

Use the directed RS envelope

    sum_i log2 U_R(c_i) <= 8192 log2 U_R(c_0) + 31 sum_i(c_i-c_0)_+.

For rational r>1 and t>0, the conserved moment for the final sum is

    Pr[sum_i(c_i-c_0)_+ > d]
      <= (8193)^R r^-d t^(A_max-8192 a)
         (E[r^((C-c_0)_+) t^(a-C)])^8192,

where C is Bin(R,A_max/(8192R)) and a=floor(A_max/8192). Centering the
tilt at a avoids a large cancellation in the directed logarithms.
The tilted row function is convex, so the preceding proof applies.

The replay pays 120 bits plus the entire history-count bound in this tail.
After the history union, route failure is at most 2^-120. For the inner
bound, retain only z_0=max(0,Lw-W) output zeros, discard all death penalties,
and use the worst birth factor 2^B_max. This is uniform over B<=B_max.

| R | History-count log2 upper | c_0 | r | t | Excess cap | Good-route log2 upper |
| --- | --- | --- | --- | --- | --- | --- |
| 192 | 5905 | 120 | 3/2 | 1 | 14798 | -306404 |
| 256 | 5997 | 124 | 3/2 | 4/3 | 78520 | -817359 |

`certify_random_pairing_turnover.py` verifies these bounds with exact
rationals and integer powers. Each family has total failure at most
2^-119; their union has failure at most 2^-118. The initial r=9/8 proposal
for R=192 failed its directed check and was replaced by r=3/2. This changes
a proof weight, not the code. The successful replay takes about 0.72 seconds.

This covers every mask sequence, swap placement, and message in these
two fixed-weight families, including sequences chosen after setup. It does
not cover varying weights, complete resets, or GF(2^128) extension claims.

## Is the route obstruction visible in sampled instances?

A bounded diagnostic samples one full regional route, with fixed seed
20261026. It alternates between choosing half-support masks on time blocks
and choosing the rows favored by those masks. This searches for large
shortened outer subcodes; it does not sample inner maps or find codewords.

For block lengths 1024,256,128,32,16 updates, the selected histories use
378,1919,3944,16256,32552 swaps. Their exact shortened RS dimensions are
32216,45455,55269,86289,117894. Even using their actual route enumerator
envelopes, the selected-history FLOAT bounds remain below -1000000 bits.
For example, 3944 swaps give an enumerator envelope of about 2253537 bits,
well below the simultaneous cap of 4397327 for the nearby B=4096 class.

This is evidence of slack in the simultaneous cap, not a bound on the
worst route or history. The search does not establish that the construction
has no obstruction. Only sources and compact results are retained.

## Checks and next discriminating test

`diagnose_random_pairing_turnover.py` proposes a bounded screen at w=256,
R=192 and 256, and B ranging from zero to 450000. It retains the exact
history polynomial and uses the binomial envelope above. It writes no
data files. Add `--conserved` for the fixed-total-occupancy bound.

`test_random_pairing_turnover.py` exhausts 12044 small support histories,
checks 12 history-coefficient identities, 405 convex-smoothing inequalities,
and 48 conserved-conditioning inequalities. All pass. The original
constant-support certificates and local accepted-SL2 checks remain intact.

Next improve the joint count of routing and changing histories, rather
than treating every history as an unrelated possible bad route. Start
with B=4096 at R=192 and 256. Retain their negative sampled-instance
controls and positive universal bounds. Then extend to changing weights
and reset boundaries. Native implementation work remains deferred.
