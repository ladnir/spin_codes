# Cheap inner mechanisms with the outer and route fixed

2026-10-26. This is a construction/proof comparison, not a production change.
The surviving candidate has a proved restricted trajectory bound. It does
not yet have a whole-code certificate or a native timing.

Follow-up: [native timings and varying trajectories](QUARTET_NATIVE_TRAJECTORIES.md)
now records the circuit's inner-only timing and broader restricted bounds.
The statements below describe the preceding mechanism-comparison checkpoint.

## Fixed problem and outcome

Keep J=GF(q), q=2^32, K=2^20, N=2^21, 8192 RS[256,128] outer rows,
and the 256 independent uniform regional permutations. The state has 512
coordinates, and L=4096 updates use

    Q_t = A_t Q_(t-1) + X_t,
    Y_t = Q_t + delta Q_(t-1),
    Q_0 = 0.

The terminal state is free. Delta is the retained nonzero, nonunit GF256
element. Every inner choice is sampled independently at setup, then fixed
for all messages. The new designs change only the inner. Local quartet
partitions do not add a global routing stage.

| Candidate | Mixer work, excluding common input/output operations | Finding |
|---|---|---|
| Accepted random pairs, control | 6 variable J-scalings per four state coordinates | Existing hard-family upper bound remains positive |
| Forward/backward weighted four-coordinate sweep | 6 variable J-scalings | Structural sparse inverse and a fixed-output-coordinate exception; cannot inherit the quartet proof |
| Global rank-one update I+u v^T | 2 variable J-scalings per state coordinate | Exact lazy/uniform law, but a whole-code obstruction prevents a near-47% guarantee at 40 bits |
| Accepted field-multiplication quartets | 9 variable J-scalings and 15 XORs per quartet | Strong quartet proof interface established; restricted hard-family certificate transfers |

The product counts are not timings. On 128-bit records, a J-scaling acts
on four J lanes. Counting it as a scalar operation does not hide that work.
The quartet circuit adds 50% more variable scalings than the pair mixer,
but uses 43.75% fewer than a direct sixteen-product four-by-four matrix.

## Recursive sweeps: spreading is not the missing proof by itself

Within four coordinates, apply

    q1 += a1*q0; q2 += a2*q1; q3 += a3*q2;
    q2 += b2*q3; q1 += b1*q2; q0 += b0*q1.

The first sweep is a weighted accumulator. Each sweep's inverse is a
bidiagonal difference map. The product of the inverses is tridiagonal,
so an interior unit output has a preimage of weight at most three.
A boundary unit output has a preimage of weight at most two. This holds
for every coefficient choice, not just unfavorable samples. Thus these
six operations never give the superregular four-by-four map used by the
strong control: that map would require input weight plus output weight
at least five.

There is also a distributional obstruction. Starting at (0,0,0,1), the
last output coordinate stays one for every choice of all six coefficients.
The first three coefficients do nothing. For independent nonzero backward
coefficients, the other three outputs range bijectively over nonzero triples.
Each such point has probability (q-1)^-3, not approximately q^-4. Randomly
ordering the four coordinates preserves this exception with probability
1/4 for each fixed singleton input; it does not supply the missing factor q.

These are local obstructions to a proposed proof interface, not low-weight
outer-compatible codewords. Longer or differently wired sweeps are not
excluded. A successful sweep proof would have to count the exceptional
directions, rather than assuming uniform quartet refresh. An unweighted
XOR butterfly also maps an all-equal input to one coordinate; it is not a
replacement for the stronger distribution without additional work.

## Concentrating randomness: useful local law, genuine global limitation

First consider fixed feedback directions: A_t=I+U_t R_t V with rank(V)=r.
Every vector in ker(V) is fixed by every update. This common kernel has
dimension at least 512-r and contains a nonzero vector of weight at most
r+1. Such a state stays localized through empty-input stretches. A changing
deterministic base map might remove this common invariant space; the claim
here is only about the identity-base design.

A more promising candidate resamples both directions at each update:

    A_t = I + u_t v_t^T,   u_t,v_t sampled from J^512.

It costs one dot product and one scaled-vector addition. For any fixed
a!=0 before exposing this update, v_t^T a is uniform in J. Conditional on
its being nonzero, u_t(v_t^T a) is uniform over J^512. Therefore the raw law is

    Law(A_t a) = q^-1 point_mass(a) + (1-q^-1) Uniform(J^512).

This is much stronger spreading than pairs for a typical fixed state.
Rejecting v_t^T u_t in {1,1+delta} makes A_t and A_t+delta I invertible;
its rejection probability is at most 2/q. However, this conditioning cannot
repair the following obstruction, which holds for every choice of keys.

### A message can neutralize every low-rank refresh

Fix a set S of w state coordinates, and allow input only in those same
positions in every batch. Let D_pi be the dimension of the routed outer
code shortened to these positions. If every rank-one dot product is zero,
the state recurrence reduces to Q_t=Q_(t-1)+X_t and stays inside S.
Imposing those dot-product equations on that prefix-sum trajectory uses
at most L-1=4095 homogeneous equations on the message. The first entering
state is already zero. If D_pi>4095, a nonzero message survives and all its
outputs are confined to S. Its output weight is at most wL.

More generally, A_t=I+U_t V_t with rank at most r can be neutralized with
at most r(L-1) equations. The same is true of r sequential rank-one updates:
require each factor to fix the common incoming state. Invertibility,
superregularity of other components, or a different key distribution cannot
avoid this dimension argument for these identity-base designs.

The actual fixed route permits an exact expectation. One outer row has
H allowed positions, with H distributed as Binomial(256,w/512), since
regions are independent. RS shortening gives (H-128)_+ dimensions. Hence

    E[D_pi] = 8192 E[(H-128)_+],    D_pi <= wL/2.

The upper bound follows from (h-128)_+ <= h/2 for h<=256, summed over rows.
Rows need not be independent. With c=r(L-1),

    Pr[D_pi>c] >= max(0, (E[D_pi]-c)/(wL/2-c)).

This is a lower bound on bad-code setup probability whenever the desired
distance exceeds w/512. The helper evaluates it with exact rational arithmetic.

| w | Output fraction at most | E[D_pi], rounded | Rank r | Bad-code probability lower bound, conservative |
|---|---:|---:|---:|---:|
| 235 | 45.8984375% | 2875.953 | 1 | This test is inconclusive |
| 238 | 46.484375% | 4248.497 | 1 | greater than 0.00031 |
| 240 | 46.875% | 5429.211 | 1 | greater than 0.0027 |
| 245 | 47.8515625% | 9533.891 | 2 | greater than 0.0027 |
| 250 | 48.828125% | 15647.031 | 3 | greater than 0.0067 |
| 251 | 49.0234375% | 17144.731 | 4 | greater than 0.0015 |

Thus even fresh full-field randomness in a rank-one correction does not
support the desired near-half ladder. At 46% it remains open, not ruled out.
The important failure is too few independent constraints per update, not
insufficient entropy in the sampled keys. A nontrivial deterministic base
map is a possible countermeasure, but needs its own trajectory analysis.

## Surviving candidate: a transitive quartet with a nine-product circuit

Partition the state uniformly into quartets at each update, as in the
[quartet control](RANDOM_PAIRING_SUFFICIENCY.md). Identify each quartet of
J-coordinates with a degree-four field extension E of J. Sample c from E
and use its J-linear multiplication matrix L_c. Reject and resample unless
both L_c and L_c+delta I are superregular. This is not a uniform matrix
distribution; the proof only requires its fixed-input law and the minors.

### Exact fixed-input law and acceptance bound

For every fixed a in E minus zero, c -> ca is a bijection on E. Before
rejection, L_c a is exactly uniform in J^4. Multiplication matrices depend
linearly on the four J-coordinates of c. Each size-j minor of L_c or
L_c+delta I is therefore a polynomial of degree at most j in those coordinates.

The helper verifies one actual-J key where all 138 nonempty square minors
of the two matrices are nonzero. Consequently none of these polynomials is
identically zero. The elementary polynomial root bound gives rejection
probability at most

    2 sum_(j=1)^4 j binom(4,j)^2 / q = 280/q.

There are N/4 matrices. Conditioning independent raw keys on acceptance
everywhere has density at most

    1/(1-(N/4)*280/q) = 1024/989 < 2.

This is exactly the density allowance used by the quartet control. Every
accepted map supplies both of its deterministic local constraints, and the
raw fixed-input point probability is q^-4. The joint support/matching/outer
count therefore transfers unchanged. No independence between different
messages is needed; the proof counts them before applying the first moment.

Checking acceptance is required setup work. Omitting rejection and merely
noting that bad minors are rare does NOT give a 40-bit setup guarantee.
Neither setup cost nor native encoding cost has been measured here.

### Explicit field and implementation circuit

For exact checking, represent J with modulus 0x10000008d. The helper verifies
irreducibility using Frobenius and a polynomial gcd. The element
alpha=0x02000000 has absolute trace one. Set

    u^2+u+alpha=0,
    v^2+v+alpha*u=0.

The first quadratic is irreducible by its trace. The absolute trace of
alpha*u in J(u) is also one, since its trace to J is alpha. The second
quadratic is therefore irreducible, giving E=J(u,v). Use basis (1,u,v,uv).

The proof representation embeds the retained AES-field element delta=0x02
as 0xa2e260ca, a checked root of X^8+X^4+X^3+X+1. This preserves the small
coefficient subfield rather than replacing delta by a degree-32 element.
A future native implementation should use the existing tower representation
and transport setup constants through a fixed field isomorphism.

The verified simultaneous-minor witness has coordinates

    (0xc423f3e3, 0xb01ba9e8, 0x5aa87c13, 0x915276c9).

For a quadratic extension z^2=z+mu, multiplying x0+x1*z by a fixed key
c0+c1*z uses precomputed coefficients c0, mu*c1, c0+c1:

    p0 = c0*x0;
    p1 = (mu*c1)*x1;
    p2 = (c0+c1)*(x0+x1);
    output = (p0+p1, p0+p2).

Recursing through the two quadratic levels costs exactly nine variable
J-products and fifteen J-coordinate XORs. Constant-key reductions are
absorbed into setup coefficients. This is an explicit circuit, not a claim
that arbitrary dense matrices have nine-product multiplication algorithms.

The proof is over J. For 128-bit payloads, extend the J-coefficient generator
to GF(2^128). A nonzero message has a nonzero J-basis component, whose output
support is contained in the full output support. Conversely J-messages embed.
Thus minimum distance is unchanged for every sampled code. We do not claim
uniform refresh over four GF(2^128) state symbols.

## What is proved against the hard trajectories

The directed replay proves the same restricted bound as the strong control:
at W=floor(0.46N), all histories with state weight 256 after every update,
exactly 32 swaps per subsequent transition, and exactly 1974 inactive outer
rows contribute at most 2^(-2842457/64), below 2^-44413. The bound includes
all messages, support histories, selected zero-output positions, and row
sets in that family. It is not a fixed-message estimate.

The twelve neighboring profiles, swaps in {16,32,64} and inactive-row counts
in {1024,1974,2560,4096}, again pass the FLOAT screen. The least negative
value is approximately -99575.69. These screens are not directed certificates
and do not cover varying state weights, arbitrary turnover schedules, resets,
or the complete distance event. The pair control remains inconclusive.

## Validation and next discriminating test

Run serially from the repository root:

    python -B research/workstreams/one_route_coupled_scan/coded_state/inner_mechanisms.py --replay
    python -B research/workstreams/one_route_coupled_scan/coded_state/diagnose_pairing_sufficiency.py --group 4 --screen

The new checks cover 943 local algebra fixtures, the exact 64-point rank-one
image law over GF4^3, rational route expectations and obstruction decisions,
and the directed hard-family replay. The quartet tests include the actual-field
minor witness, a nine-product call counter, independent reference multiplication,
all 255 nonzero input image laws over GF4^4, and invalid-map negative controls.
The small-field checks validate code paths; the algebra above proves the laws.

An initial attempt to use the polynomial-basis element as a subfield generator
failed its order assertion. The helper now finds and verifies an element of
order 255 before constructing the AES-field embedding; the final replay passes.
No unsuccessful check was interpreted as a passing certificate.

All runs are serial, file-free, under 60-second watchdogs, with at least
20 GiB free. The largest working arrays are the existing length-257 coefficient
arrays and small local tables; no large data job or native benchmark was run.
Production code, old proofs, and the frozen diagnose_pair_defects.py are unchanged.

Next broaden the quartet proof to varying activity and turnover, while timing
this specific nine-product circuit against accepted pairs. Do not assume the
arithmetic count predicts end-to-end latency. Retain recursive sweeps as an
open alternative requiring a different local count; park identity-base low-rank
refresh as the near-half candidate because of the explicit obstruction.
