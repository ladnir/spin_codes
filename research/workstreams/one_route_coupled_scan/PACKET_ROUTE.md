# Packet routing for the persistent-state construction

Status: exact route interface, complete finite certificates for three
profiles, and checked native implementations. Serial timings
remain in the scalar implementation's performance range. This changes the routing
distribution in [Persistent MDS state](PERSISTENT_MDS.md). Its scalar-route
margins do not automatically transfer.

The proposed route moves p adjacent outer coordinates together. For p=4,
four 128-bit records form one 64-byte packet. The intended benefit is
contiguous memory access. The cost is additional correlation between outer
rows, which can weaken the finite distance bound even when the inner is
unchanged.

## Outer rows and the sampled route

Let F have q elements. Fix an even B, put k=B/2, and use a fixed [B,k]
MDS outer code over F. A subfield RS code extended to F is sufficient.
Fix N divisible by B and write

    M = N/B,       G = M/p = N/(pB).

Require p to divide M. Partition the M independent outer message rows into
G fixed groups of p rows. The partition is the same in every region.
Write c[g,a,j] for coordinate j of row a in group g, with
0<=g<G, 0<=a<p, and 0<=j<B. Each row has its own k message symbols;
the code remains a direct sum of M constituent codes.

There are B chronological regions, each containing M scalar coordinates.
For every region j, sample an independent uniform permutation pi_j of
{0,...,G-1}. At destination packet h, place

    X[jM+ph+a] = c[pi_j(h),a,j],       0<=a<p.              (1)

Thus one permutation moves whole groups within each region. The internal
order a is fixed. All permutations are sampled once and fixed for every
message, independently of the inner setup. Concatenating the regions is
one global permutation of the N outer coordinates, not a second routing
pass.

The p rows in (1) are independent code coordinates. They are not the four
GF32 payload lanes within one 128-bit record. Those payload lanes can
still share one sampled GF32 code, as in the persistent-state construction.

## Exact shortening on complete packets

Use inner batches of b scalar coordinates. Require b to divide N and p
to divide b. Region starts and batch boundaries are then packet-aligned.
In particular, each union U of complete batches selects whole packets.
Deleting the first batch of each interval from U also leaves whole packets;
call that set A.

Fix a route realization and a coordinate set S that selects whole packets.
Let T_g(S) count the regions in which S contains the packet from group g.
Every row in that group has exactly the same selected constituent positions.
Write C_out for the routed outer code and

    D(S) = dim_F(C_out intersect F^S).

An MDS [B,k] row shortened to T positions has dimension (T-k)_+.
The rows have independent message spaces, so

    D(S) = p sum_(g=0)^(G-1) (T_g(S)-k)_+.                 (2)

This identity holds for every route, not only on average. No randomness
of the constituent generator is needed.

The complement S^c also selects whole packets, and
T_g(S^c)=B-T_g(S). Since B=2k and p sum_g T_g(S)=|S|,

    D(S) = |S|-N/2+D(S^c).                                (3)

The scalar-route complement identity therefore survives unchanged.
Equation (2) is not asserted for sets that cut through packets.

## Route moment with the actual packet distribution

Fix S independently of the route, and suppose S is a union of at most c
intervals. Set s=|S| and rho=min(2c,B). At most rho regions are partial.
In a partial region j, let m_j be the number of selected destination
packets. Its selected source groups form a uniform m_j-subset of G groups.
Thus the selection probability for each group is m_j/G, and

    mu = sum_j m_j/G = s/(pG) = sB/N,                     (4)

where full regions contribute one to the sum.

For t>=1 define the same padded-binomial envelope used in the
[regional comparison](STRATIFIED_FRACTIONAL_BOUND.md):

    H_rho(mu,t)
      = max_f E[t^((f+Bin(rho,(mu-f)/rho)-k)_+)],           (5)

where f is integral, 0<=f<=B-rho, and 0<=mu-f<=rho.

The fixed-cardinality comparison now operates on G group indicators,
not M independently routed row indicators. Replacing each regional
fixed-size selection by independent Bernoulli indicators increases the
product moment in (2). The factor for each group is nonnegative and
nondecreasing in every regional indicator. Convex equalization then
bounds its Poisson-binomial sum by (5). Padding with deterministic zero
or one indicators covers fewer than rho actual partial regions. Hence

    E_route[tau^D(S)] <= H_rho(sB/N,tau^p)^G,   tau>=1.    (6)

Both changes matter: the base is tau^p and the exponent is N/(pB).
Keeping the old base tau would discard p-fold source dimensions.

Complement before applying (6), as in the scalar-route proof. Define

    Moment_(p,c)(s,tau)
      = tau^((s-N/2)_+)
        H_rho(min(s,N-s)B/N,tau^p)^(N/(pB)).               (7)

Equations (3) and (6) give E_route[tau^D(S)]<=Moment_(p,c)(s,tau).
The complement has the same partial regions, so rho remains min(2c,B).

In the fractional-moment argument, tau=q^nu with 0<nu<=1. Therefore
1<tau<=q is still the relevant restriction. The base tau^p in (5) is
only a route-moment parameter and need not be at most q.

## What remains unchanged in the whole-code bound

Retain an F-linear persistent inner satisfying the inverse-uniformity and
MDS reset lemmas in [Persistent MDS state](PERSISTENT_MDS.md). The inner
setup remains independent of the packet route. This section changes
neither its state update nor its terminal convention.

Put n=N/b. An output with c excursions occupying v batches has

    u=|U|=bv,       a=|A|=b(v-c),       r=c-1.

The same feasible classes and representation count apply:

    2c-1<=v<=n,       (b+1)(c-1)+1<=e,
    J(v,c)=binom(v-c,c-1) binom(n-v+c,c).                  (8)

Here e is the bad-output weight threshold. Output supports need not be
unions of packets. Only U and A require that property, and batch alignment
provides it. The reset ranks and output-volume count still refer to
individual F coordinates.

Let W be the representation-wise upper bound on eligible outputs from
the persistent-state proof. For example, its two reset bounds give

    W <= min(q^(u-br),
             q^(-br) sum_(w=0)^min(u,e) binom(u,w)q^w).    (9)

For each 1<tau<=q and nu=log_q(tau), the inverse and forward envelopes
become, respectively,

    U_vol = W^nu tau^(-a) Moment_(p,c)(a,tau),

    U_reset = tau^(-br)
              (Moment_(p,c)(u,tau)-1)/(tau-1).             (10)

Their derivations are unchanged. The inverse argument uses an affine
fiber with dimension D(A). The forward argument counts projective input
lines in the shortened space of dimension D(U). In particular, D(A)=0
does not justify deleting an inverse term.

The whole-code failure probability is bounded by the sum of J(v,c) times
the smaller complete envelope in (10), also capped at one, over (8).
Different valid markers may be used for the two envelopes and each class.
The sufficient deterministic exclusion

    uB/N+2c<=k  implies  D(U)=0                           (11)

also remains valid. No such exclusion follows merely from D(A)=0.

## Finite checks and limits

A finite check must replace the route moment by (7) throughout its complete
cover, including complemented sets and the projective subtraction. Keep
the original q, batch count, interval count, and scalar output threshold.
Do not rescale those quantities as if one packet were a field symbol.

The helper `theory/packet_mds_certificate.py` preserves the frozen scalar
helper. All three completed profiles use K=2^20, rate1/2, GF(2^32), and
threshold 734003. The first two use term target 75 bits. Every proposed
rectangle is replayed using exact rational route moments and directed
arithmetic at 80 and 120 digits.

| Packet records p | Outer B | Batch=state b | Distance | Margin at least | Rectangles |
|---|---:|---:|---:|---:|---:|
|2|128|64|>35%|46 bits|1448|
|4|128|128|>35%|48 bits|720|
|4|64|128|>35%|59 bits|1652|

For the first two rows, the replay's worst log2 term bounds were at most
-103.3894958 and -116.5384571, respectively. The reported margins use the
common75-bit term target and whole-code union factors below2^29 and2^27.
No larger margin is claimed from the displayed rounded values.

The third row keeps the cheaper RS[64,32] outer. It uses term target67
and the finer bounded rational grid selected by `--fine-markers`; its
worst directed log2 term is at most -67.375900734. The union factor is
below 2^27, giving a coarse 40-bit margin. Summing the individual covered
rectangle contributions with `--sum-boxes` proves the stronger 59-bit margin.
The directed 80/120-digit replay gives the exact dyadic upper bound

    117849570452783704458019672788624293933296618182474848723848 / 2^256
      < 2^-59.

The complete cover still has 1652 rectangles. The summation retains each
rectangle's verified bound and multiplicity instead of assigning every
class the common 67-bit bound. This refines only proof accounting, not
the route, inner, or encoder work. An independent analytic review checked
the optional summation interface.

Reproduce serially from the repository root:

```text
python -B research/workstreams/one_route_coupled_scan/theory/test_packet_route.py
python -B research/workstreams/one_route_coupled_scan/theory/packet_mds_certificate.py --packet 2 --B 128 --batch 64 --field-bits 32 --term-bits 75 --seconds 90
python -B research/workstreams/one_route_coupled_scan/theory/packet_mds_certificate.py --packet 4 --B 128 --batch 128 --field-bits 32 --term-bits 75 --seconds 90
python -B research/workstreams/one_route_coupled_scan/theory/packet_mds_certificate.py --packet 4 --B 64 --batch 128 --field-bits 32 --term-bits 67 --fine-markers --sum-boxes --seconds 30
```

The small-field test passes4096 actual shortened-dimension checks,4096
complement checks, and255 exact rational moment comparisons across all16
routes. Negative controls detect partial-packet supports and omission of
the packet exponent. This validates interfaces, not the K20 certificate.
An independent analytic review checked the finite helper and its full union.

Failed floating proposals are preserved as failed upper bounds: p2/B64/b64
first rejects (v,c)=(11557,69), at log2 bound -73.3966101; p4/B128/b64
first rejects (11562,70), at -74.8278626. Neither is a low-distance witness.
The p1 control reproduces the original B64/b64 certificate with2846 rectangles.

The coarse-marker p4/B64/b128 proposal first failed at (5791,46), with
bound -66.3716641. The finer grid above closes that profile completely.
For p2/B64/b64, the fine grid at term69 first fails at (11544,57), with
bound -64.7727485; a denominator64 marker grid improves that point to
-68.4207553 but does not establish the required complete cover.
The fine-grid term40 probe also fails, at (11541,54), with bound
-28.7149277. That smaller-state profile remains uncertified and is not
instantiated by the packet benchmark.

## Native checks and serial measurements

`native/persistent_native.cpp` selects packet width at compile time with
`SPIN_OR_PACKET=1`, `2`, or `4`; one is the default. Packet variants use
the additive-coset inner. They store N/p labels and preserve each packet's
lane order. Width two combines two aligned 32-byte loads per four records;
width four uses one aligned 64-byte load. Encoding fuses those reads into
the inner, without materializing a routed array.

All three certified profiles pass the native `--check` gate. It checks
regional group bijections and actual packet sources with an independent
scalar index formula. It also compares fused and materialized routing,
checks the RS outer and full stream directly, and verifies the conceptual
global Vandermonde state after each batch. Shift checks cover basis and
dense inputs through field bit21 against direct Cauchy matrices.
K256 receives the direct stream oracle. K1024 has nontrivial packet
permutations, including the B128/p4 profile whose K256 regions contain
only one packet. A further K=B*b fixture exercises two batches per region.
The larger fixtures avoid repeating the quadratic field oracles.

The measurements use Peach CPU2, GCC15.2, and
`-O3 -std=c++20 -march=znver4`, with three warmups and medians of 31 samples.
Runs are serial.
They encode K=2^20 input records into N=2^21 output records. Setup and
allocation are excluded. Component times are not additive.

| Profile (p,B,b) | Route only | Inner, contiguous | Inner, fused route | Full |
|---|---:|---:|---:|---:|
| (2,128,64) |4.062343 ms|10.474301 ms|16.962221 ms|23.489125 ms|
| (4,128,128) |3.155800 ms|14.321332 ms|18.382182 ms|25.127033 ms|
| (4,64,128) |3.154939 ms|13.710633 ms|18.484083 ms|23.184977 ms|

Packet routing reduces isolated route cost, but the complete certified
profiles also change batch or outer size. Their full times do not establish
a clear improvement over scalar B64/b64. Repeat full medians, in table
order, are 23.447340, 25.826732, and 24.051388 ms. A same-session scalar
B64/b64 control takes 23.325273 ms. These results remain in the same
performance range or slower; do not select a winner from the lowest median.
All three packet profiles also pass UBSan with recovery disabled. The
[native checkpoint](NATIVE_PERSISTENT.md) records outer times, setup sizes,
checksums, source receipts, and reproduction commands.

Next use the component costs to choose a structural change rather than
promoting packet routing as a speedup. The historical 3.798 ms reference
remains far faster than these measured profiles.
