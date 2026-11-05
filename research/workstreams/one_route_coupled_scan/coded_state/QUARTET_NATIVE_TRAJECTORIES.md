# Field quartets: native cost and varying-state trajectories

2026-10-26. Follow-up to [the mechanism comparison](INNER_MECHANISMS.md).
The outer and global route remain fixed. No library or paper changes.

Follow-up: [the bounded 224/32 analysis](QUARTET_224_TWO_UPDATE.md) closes the
gap reported below at 46% without changing the construction. This note
preserves the preceding diagnostic and timing checkpoint.

## Outcome

The nine-product quartet circuit costs 5.54--5.78 ms for the inner alone
with partition-ordered state buffers. The pair control costs 4.20--4.50 ms.
These are not end-to-end encoding times.

The proof now covers selected varying-weight histories. However, a newly
tested constant-weight-224 family does not close at 46% distance; it closes
at 45%. Neither result covers all histories or gives a whole-code certificate.

## Measured computation

[quartet_native.cpp](../native/quartet_native.cpp) processes N=2^21 records
of 128 bits, with 512 state records and 4096 updates, corresponding to
K=2^20 at rate one half. Each record contains four GF32 coordinates. The
quartet multiplier acts across four records, applying the same scalar map
to each payload lane; it is not a uniform GF128 randomizer of each record.

Pairs use three shears; quartets use the explicit nine-product nested
quadratic circuit. Recurrence timings include the input addition and
Y_t=Q_t+delta Q_(t-1), with the existing AES-subfield delta=2.

The host is Peach's Ryzen 7950X, affinity CPU 2, boost disabled, GCC 15.2.0,
`-O3 -std=c++20 -march=znver4`. After warmup, ten repetitions alternate the
order of compared kernels. No concurrent benchmark was launched. The table
gives final-run medians, not a best-case selection.

| Timed work | Pairs | Field quartets |
|---|---:|---:|
| Arithmetic on contiguous SIMD-grouped records | 3.314 ms | 4.286 ms |
| Recurrence with state/input gathers and state/output scatters | 9.498 ms | 12.583 ms |
| Recurrence with partition-ordered buffers | 4.498 ms | 5.782 ms |

The preceding run measured 4.196 and 5.540 ms for partition-ordered buffers.
Treat 4.2--4.5 and 5.5--5.8 ms as observed ranges, not improvements between
runs. The quartet costs about 29--32% more than pairs in this layout.

In partition order, each step writes its next state sequentially to a second
buffer. The following step gathers only that state; input and output accesses
are sequential in the current partition order. This requires composing the
fixed per-step coordinate relabeling into the outer route during setup.
Output relabeling preserves weight. The benchmark does not implement or time
that composition or the resulting outer-plus-route locality. The underlying
sampled local maps need not change.

### Scope and validation

Setup, the outer, global routing, hashing, and allocation are excluded.
The harness samples sixteen accepted keys per family and reuses them in
expanded coefficient arrays to bound setup work. It does not sample a full
independently keyed code from the proved ensemble. This is an arithmetic and
memory cost workload, not a sampled-code distance experiment. Hot kernels
are branch-free in coefficient values. Expanded coefficient streams occupy
48 MiB for pairs and 72 MiB for quartets: the pool does not replace them with
cache-resident lookups. Each update has its own precomputed shuffled partition.
The PRNG is for reproducible timing only.

Every pool key passes acceptance. Scalar replay checks all payload lanes of
the first two entire scattered updates and the first three partition-ordered
updates. The scalar quartet oracle uses direct field multiplication, not the
prepared nine-product circuit. Final warning-enabled build/checks pass.

The native GF32 tower uses alpha=0x20000000, checked to have absolute trace
one, and u^2=u+alpha, v^2=v+alpha*u. This is a separate basis instantiation;
we do not assume that basis changes preserve minors. Its checked witness is

    (0xedb82667, 0xb336d1ca, 0xc32dfddc, 0xbd3652ce).

All square minors of L_c+delta D are nonzero for every binary diagonal D.
In particular both matrices required by the retained acceptance rule are
superregular, justifying the same 280/q per-map rejection bound here.

## Varying-state bound

Fix weights w_1,...,w_L and deaths d_t=|supp(Q_(t-1)) minus supp(Q_t)|,
t>=2, before sampling the code. Set

    A=sum_t w_t, D=sum_(t>=2) d_t,
    z=max(0,A+D-W), W=floor(0.46 N).

Successive support unions contain A+D positions. Cancellations occur only
on intersections, so at least z selected zeros are required. If

    z > sum_(t>=2) (w_(t-1)-d_t),

the schedule cannot produce output weight at most W, regardless of setup.

Otherwise use the local polynomials L_a from RANDOM_PAIRING_SUFFICIENCY.md,
allowing next-weight and death markers to vary by transition type. For
incoming weight u, next weight v, and d deaths, the transition bound becomes

    K(u,v,d) = [h^u](sum_(a=0)^4 binom(4,a)L_a h^a)^128
               / (binom(512,u) next^v death^d).

The family bound is

    binom(8192,r) binom(512,w_1) (1-p)^w_1 pi^-256
      * 2^theta * 2^[32 theta (A-K+128r-z)] * zero^(-theta z)
      * product_(t>=2) K(w_(t-1),w_t,d_t),

where pi=binom(8192,r)p^r(1-p)^(8192-r). Incoming-support symmetry permits
composition across changing weights. This unions labeled supports, compatible
zero selections, messages, and inactive-row sets consistent with the fixed
schedule. It does not union all weight schedules or assume independence
between messages.

[quartet_trajectories.py](quartet_trajectories.py) bounds the incoming
coefficient by positive evaluation at h, optimizes markers with theta=3/16,
then replays rational approximations. Powers P^(3/16) round upward through
integer sixteenth-power comparisons. Rational logarithms use outward integer
bit-length bounds, treating powers of two exactly. Combining the regional
conditioning probability and each transition kernel before bounding their
logarithms avoids unnecessary loss from separately cancelling terms.
The two alternating families and the hard family's 45% result have frozen
rational markers in the helper; `--frozen` replays them without optimization
and checks their exact reported exponents.

## Restricted results

All rows use exactly 1974 inactive outer rows and 46% distance. Directed
log2 upper bounds below are rounded upward for display.

| Fixed weight schedule | Deaths per transition | Log2 upper bound |
|---|---|---:|
| 128 throughout | 16 | -1,723,772 |
| 224 throughout | 32 | +57,303: inconclusive |
| 256 throughout | 32 | -74,946 |
| 288 throughout | 32 | -242,007 |
| 384 throughout | 32 | -912,842 |
| Alternating 240,272 | 16 then 48 | -117,352 |
| Alternating 224,288 | 16 then 80 | -313,045 |
| 256 throughout | 0 | -359,171 |

Weight 256 with 128 deaths every transition, and alternating 0,256 with full
resets, are impossible at this threshold by support counting. This does not
rule out every reset-containing history.

For the 224/32 family, changing only distance to 45% gives the directed bound
-1830657/32, below -57,208. At 46%, optimizing theta improves the FLOAT bound
from +56,020 to +53,740, still inconclusive. Free-theta FLOAT bounds at 47%
and 49% are about +153,583 and +321,948. These positive upper bounds are not
lower bounds on bad-code probability or observed low-weight codewords.

## Diagnosing the new gap

The 224/32 family requires only about 20.48 zero outputs per update at 46%,
versus about 52.49 in the previous half-full family. Partially active and
empty quartets still matter. Closing the half-full family did not establish
a uniform margin across activity levels.

As an ideal control, use independent uniform raw 512-by-512 inner matrices
over GF32 with the same outer and route. For these nonempty histories, every
subsequent fixed-state transition has probability q^-512. A deliberately
loose count of messages, state values, supports, and zero selections yields

    2^[8192+N+A-w_1]
      * q^[K-128r+A-w_1-z-512(L-1)].

Its exponent is -11,931,328 at the hard 46% family. Thus the outer and route
do not force this particular obstruction. This expensive ideal control is
not a certificate for the actual quartet inner.

Two narrow countermeasures did not resolve the gap:

- Requiring every L_c+delta D to be superregular strengthens setup without
  changing arithmetic. The witness above makes every minor polynomial
  nonzero. A union bound gives rejection at most 2240/q per map and total
  density at most 128/93<2. Joint zero-state/output constraints then exclude
  a local choice when |T intersect ((group minus V) union Z)| >= |U| for
  nonempty U. The optimized hard-family FLOAT bound improves by only 1.21
  bits. This stronger acceptance rule is not adopted.
- The necessary aggregate active-row coverage condition is
  C_active >=129(8192-r). An additional positive marker at least one bounds
  this event in the joint count. Its optimizer chooses one and gives no
  improvement. This does not test the stronger row-by-row occupancy condition.

These checks argue against a long pass confined to local minor tests or
scalar-marker tuning. They do not distinguish an actual quartet weakness
from loss of multi-update or rowwise information.

## Next discriminating test and retained checks

Before a broad certificate sweep, test the 224/32 family with a two-update
count retaining surviving coordinates and rowwise outer constraints. Compare
its loss with the dense-inner control. If that bounded test does not recover
the gap, test modest cross-quartet mixing on this family before implementing
another design. Preserve this candidate and all earlier certified controls.

Native builds/runs were serial, with 60-second timeouts and 2 GiB virtual
memory limits. Python runs are file-free, single-threaded numerics with a
60-second watchdog; enumerations are bounded by quartet masks and at most
two transition types per tested schedule. Disk checks retain at least 20 GiB
free. The previous 943 local checks, restricted replay, and 10740 pairing
checks pass. No raw timing files or local build tree are retained. Frozen
helpers and earlier proofs are unchanged.
