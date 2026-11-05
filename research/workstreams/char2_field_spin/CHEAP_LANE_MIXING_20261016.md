# Cheap lane mixing from the 3.8-ms baseline

The useful candidate changes only the state refresh of the
[certified 35% design](HD_SHIFTED_CLOSED_CERTIFICATE.md). Matched K20 tests
put it close to the unchanged baseline: approximately 3.48 ms against
3.47 ms in this harness. This is a performance result and a collection of
local proof lemmas, **not a new whole-code distance certificate**. The
35% / 55.09-bit certificate remains attached only to the original design.

The experiment is
[cheap_lane_forward.cpp](../../../spin/experiments/char2_field_spin/cheap_lane_forward.cpp).
It leaves the baseline source and production library untouched. This
work does not promote a profile, update the paper, or publish changes.

## Construction: change the refresh, not the data movement

Let E be AES GF256. Represent each 128-bit record in

    F = E[Y] / (Y^16 + Y^5 + Y^3 + 0x32).

Write theta for the residue of Y. Byte j is the coefficient of theta^j.
The existing field oracle checks the modulus by Rabin's irreducibility
criterion. This is a specified tower representation, not an identification
with an application's arbitrary binary-polynomial representation of GF128.

Retain K=2^20, N=2^21, the four RS[64,32] constituents per superblock,
the GF256 Toeplitz outer mixers, and the independently swapped paired
route. Also retain the 64-record inner steps and 16-record state.
The maps A, C=A^T, and H are unchanged. A is the direct sum of four
RS[16,4] encoders. H is the GF256 MDS interpolation map on 16 coordinates.

For each step t and state coordinate j, independently sample a[t,j] from E.
Retain the baseline independent nonzero GF256 diagonals D1[t] and D2[t].
All samples are independent of outer and routing setup. Setup is sampled
once and then fixed for every message. Starting with Q[0]=0, compute

    Z[t]   = X[t] + A Q[t],
    Q[t+1] = diag(theta + a[t,j]) D2[t] H D1[t] (Q[t] + C X[t]).

Here X[t] and Z[t] contain 64 records, and Q[t] contains 16 records.
The state persists across region boundaries. The output uses the same
fixed native ordering as the baseline.

Every operation is F-linear. In particular, this is genuine lane mixing
that remains compatible algebraically with field-linear code consumers.
It is not merely sixteen independent GF256 encodings. Representation
conversion and the missing distance certificate still prevent claiming
drop-in compatibility with a particular PCS or VOLE implementation.

The inner is an invertible causal map on N records: each output step is
the current input plus a function of earlier inputs. The outer remains
injective over F, so the rate remains 1/2.

### Cost and implementation

For a record v, the new operation is theta*v + a*v. Multiplication by a
is a bytewise GFNI product. Multiplication by theta shifts the 16 bytes and
reduces the top coefficient using theta^16=theta^5+theta^3+0x32.
The fixed companion kernel uses two byte shifts, one shuffle, one GFNI
product, and one ternary XOR. No general GF128 multiplication is used.

The work is applied to 16 state records per 64 output records, immediately
after the existing refresh. There is no extra data-buffer traversal,
global permutation, larger state, or changed outer. The additional
precomputed bytes occupy N/4 = 524288 bytes at K20. The experimental
Mersenne Twister seeds reproduce tests; they do not certify a sampled code
or replace the ideal independent setup distribution in probability claims.

## Proofs for this exact coefficient distribution

The local probability statements fix the message and its lane span before
sampling refresh randomness. Condition on the outer/routing setup and
earlier refreshes; current coefficients remain fresh unless explicitly
exposed. These statements do not directly apply to a message or span
chosen after observing current coefficients. A whole-code proof must
handle the union over messages separately.

The [goal and proof ledger](CHEAP_LANE_MIXING_GOAL.md) records independent
audits and directed bounds: all quiet words and the one-active-group
class close at 40%; the two-through-four-group classes close at 35%.
The nonquiet remainder is still open. The encoder has not changed during
these proof refinements.

### One-coordinate and joint cancellation bounds

Let b be the existing D2 entry, sampled from E minus zero, and a the new
independent uniform byte. The combined coefficient lambda=b(theta+a) is
uniform on

    S = {c + b theta : c in E, b in E minus zero},  |S| = 256*255.

The map (a,b) -> (ba,b) is bijective. Every member of S is nonzero.
For every fixed nonzero v in F and target u in F,

    Pr[lambda*v = u] is either 0 or 1/65280.

For the joint form, fix nonzero v[1],...,v[m], a matrix M over F, and a
target z. Sample independent lambda[i] from S. If M has rank r, then

    Pr[M diag(v) lambda = z] <= 65280^(-r).

To prove this, choose r pivot columns and condition on the other
coefficients. The equations determine at most one assignment to the r
pivots. Each pivot has 65280 equiprobable values. Zero v coordinates must
be removed before computing the rank. In the encoder, this bound can be
used after conditioning on the history and D1, but before revealing the
current independent a and D2 entries.

This is not uniform F randomness. Each scalar has only about 16 bits of
entropy, and its support lies in a two-dimensional E-subspace. The old
uniform-F proof cannot be imported by changing a field-size parameter.

The family has no common proper E-linear invariant subspace: it contains
multiplication by theta, whose minimal polynomial over E has degree 16.
An E-subspace invariant under that operator is invariant under E[theta]=F,
and hence is either zero or F. This fact alone does not prove distance.

### A quotient argument for every proper input-lane span

Fix a nonzero proper E-linear subspace U of F. Consider a fixed input
stream whose records all lie in U. This includes every message with
proper E-lane span, since the unchanged outer and route preserve U.
Write h for the number of state coordinates outside U before a step.

The output modulo U is A applied to the state modulo U. Any nonzero
four-coordinate block produces at least 13 nonzero quotient records.
To see this, choose an E-linear functional on F/U that detects one of its
coordinates; the resulting scalar word is a nonzero RS[16,4] word.
Consequently every step satisfies

    wt(Z[t]) >= 13 ceil(h/4).

Now suppose h>0 and condition on the current state, X[t], and D1[t]. Set

    W = H D1[t] (Q[t] + C X[t]).

Let g count the coordinates of W outside U. Since C X[t] lies in U^16,
the MDS property of H gives h+g >= 17. This property also holds in F/U:
apply a functional detecting a nonzero quotient coordinate, then use
the scalar MDS support bound.

For w outside U, at most one a in E makes (theta+a)w lie in U.
Two distinct such bytes would imply (a-a')w lies in U, contradicting
w outside U. For w in U, the result lies outside U exactly when theta*w
does; this condition does not depend on a. Nonzero b in E does not change
membership in U.

More precisely, let d count the g outside coordinates that have an
exceptional byte, and let l count inside coordinates whose theta multiple
is outside. Conditional on W, the next quotient support has the exact law

    h_next = g + l - Binomial(d, 1/256),
    17-h <= g <= 16, 0 <= d <= g, 0 <= l <= 16-g.

The binomial trials are independent because the a bytes are independent.
This is an actual transition law for the implemented randomizer, not an
independent-uniform-state assumption. A proof may retain g,d,l or bound
over their feasible values; it must not assume a transition based only on h.

For example, let D be the number of these exceptional-byte hits. Then
h+h_next >= 17-D. The two consecutive outputs satisfy

    wt(Z[t]) + wt(Z[t+1]) >= 13 ceil((17-D)/4).

With no exceptional hit, this is at least 65 of 128 records. Since D is
dominated by Binomial(16,1/256), the conditional expected lower bound is
about 64.2109 records, or 50.1647%. This expectation is conditional on
h>0. It is neither a global distance guarantee nor a concentration bound.

### Activation and a deterministic rank-one window

When h=0, the quotient argument has not established activation. W then
lies in U^16, and a coordinate escapes precisely when W[j] is outside
U intersect theta^(-1)U. This intersection has dimension at most dim(U)-1;
equality with U would contradict the invariant-subspace statement above.
Write T for that intersection. Since H D1 is invertible over E and
preserves T^16 in both directions, h_next=0 holds exactly when
Q+CX lies in T^16. Thus activation is determined before the current
refresh samples. The outer and feedback constraints must account for it.

For U=vE with v nonzero, the intersection is zero. Thus any nonzero
feedback from a zero state activates a nonzero quotient state.
There is also a deterministic initial guarantee. Divide the entire stream
and state by v, and begin at its first nonzero feedback. For the next
15 incoming states, their polynomial degrees in theta are exactly
1,2,...,15. More explicitly, let tau be the first step with C X[tau]
nonzero. Then Q[tau]=0. For k=1,...,15, define c[k] as the coefficient
of theta^k in Q[tau+k]/v. These leading vectors satisfy

    c[1]   = D2[tau] H D1[tau] C(X[tau]/v),
    c[k+1] = D2[tau+k] H D1[tau+k] c[k],  1 <= k <= 14.

The input contributes only degree zero. Before reduction at degree 16,
it cannot cancel the leading coefficient. The matrix on the right is
invertible, which proves the induction. The MDS support argument gives
65 nonzero records per consecutive pair of these protected output steps.
If all 15 steps remain before the end of encoding, pairing 14 and using
the single-step bound on the last gives at least 468 nonzero records.
This is local output weight, not a relative-distance certificate.

The quotient argument remains available after polynomial reduction;
the degree-growth argument alone does not extend beyond this window.

## What remains unresolved

There is no whole-code certificate for this modified encoder, including
at 35%. In particular:

- Delayed activation and quiet feedback require the actual outer/routing
  counting bounds. Refresh mixing does nothing when its input is zero.
- Proper lane spans require counting input spans and constrained histories,
  not multiplying one-coordinate marginal bounds.
- Full E-rank-16 inputs have U=F. The quotient argument is then empty and
  needs a separate rank/counting argument.
- A setup-failure certificate must cover every message and terminal tail.
  Neither the 65/128 local bound nor its expectation supplies that union.

There is a useful deterministic ceiling. Restrict messages to 4096 of
the 8192 superblocks. Their dimension is 4096*128=524288 over F.
Impose C X[t]=0 for the first 32767 inner steps. This gives at most
32767*16=524272 equations, so a nonzero solution exists. Its state is zero
before every emitted output, including the last one. Its output equals
the routed outer word and has weight at most 4096*256=N/2. Thus this
construction cannot exceed 50% minimum distance. This does not rule out
approaching 50%; it shows why quiet histories must stay in the analysis.
The subsequent [audit](CHEAP_LANE_MIXING_GOAL.md#an-actual-distance-ceiling)
sharpens this ceiling to 255/512, or 49.8046875%, by averaging the
nonzero GF256 words in that quiet subcode.

The recommended next proof step is to combine the exact quotient law
with the preserved outer/route counting interface, starting with rank-one
and quiet histories. Audit full-rank inputs before launching a large
certificate search. Keep the 3.8-ms implementation as the control.

## Rejected first candidate

An outer-store experiment used L=T^b P^a independently on each record,
for two fair bits a,b. P cyclically rotates the 16 bytes; T adds byte 0
to byte 1. This is E-linear but generally not F-linear.

The family has no common proper E-linear invariant subspace: conjugating
T-I by P produces the cyclic matrix units, whose products produce all
matrix units. However, its exact maximum point probability is 1/2, not
anything close to uniform F. For example, a basis vector away from the
shear coordinates has only two possible images. A common permutation
with only an optional shear is weaker still and leaves a large invariant
subspace after factoring out that common permutation.

The implemented four-map family took about 4.00--4.02 ms versus baseline
3.47--3.60 ms in matched preliminary runs. Its preexpanded masks cost
4.5 MiB at K20. The state-refresh candidate has a better cost and proof
interface, so the four-map outer experiment is not the leading direction.

## Measurements and checks

The final three serial runs used Peach CPU 2 (Ryzen 9 7950X), GCC 15.2,
and `-std=c++20 -O3 -march=znver4`. Each process measured baseline,
candidate, candidate, baseline, with 101 calls per block and five warmups.
Setup and allocations were excluded. Input was produced once, with no
producer or hash between timed calls. The shared host remains noisy.

| Run | Baseline block medians, ms | Refresh block medians, ms |
|---|---|---|
| 1 | 3.509372, 3.428331 | 3.480499, 3.465830 |
| 2 | 3.557181, 3.545359 | 3.539838, 3.501658 |
| 3 | 3.371495, 3.360564 | 3.466532, 3.483814 |

The medians of the six block medians are 3.4688515 and 3.4821565 ms:
about 0.4% difference. The quieter third run gives about 3.2% overhead.
This supports low overhead, not a speedup claim. The historical 3.798160-ms
checkpoint used a different harness; do not attribute the absolute
difference from that number to this change.

Checks passed for all 128 record basis bits and all 256 byte coefficients,
against independent polynomial multiplication. Complete K13 and K14
encodings match scalar-refresh composition, including exact in-place use.
A complete K13 encoding commutes with multiplication by theta. The
unchanged control passes its dense oracle. Four rank-one streams also
pass degree-growth and per-step emission-floor checks. The field oracle's
six tests, including irreducibility, passed remotely.

Final driver SHA256:
`28bee7d5ec10cd7491e9f9ec69070f498ef8e5e9b4952bd4ca1b7e757cbc4969`.
Final binary SHA256:
`7a86a3dbbb90cdf5f67f4fdc349daafce595069a63044788c60295ed00a3dada`.
K20 control checksum: `7701295210980077712`.
K20 refresh checksum: `1815379804882320596`.

From the preserved source root:

```sh
g++ -std=c++20 -O3 -march=znver4 -Ispin/include \
  spin/experiments/char2_field_spin/cheap_lane_forward.cpp \
  spin/experiments/char2_field_spin/Portable.cpp -o cheap_lane_forward
taskset -c 2 timeout 60 ./cheap_lane_forward --check
taskset -c 2 timeout 60 ./cheap_lane_forward --bench-refresh 20 101
```

Only one benchmark may run at a time. This is an aligned-buffer research
driver, not a new public API.

The complete source bundle, binary, checks, and final timing logs are
preserved on Peach at
`/home/prindal/spin-research/cheap-lanes/20261016-refresh-v1`.
Its `RECEIPT.json` inventories 89 files and has SHA256
`7d466d6786fe4e2a50bad4b619e479eb08ff47d448717c51388ea85dbf43aeaa`.
All listed hashes were checked after moving the bounded scratch directory
to durable storage. Reproduction can run directly from that directory;
no local archive or build tree is needed. Preliminary v1/v2 logs are
explicitly marked exploratory because their intermediate source revisions
were not retained. The final three runs are pinned to the source above.

The preserved remote directory contains 3645835 bytes, including the
receipt. The corresponding temporary directory is gone. This batch adds
approximately 30 KiB of local source, notes, and index text, with no local
build, bulk proof data, or transport archive. The local volume still has
about 50.8 GiB free. No certified evidence or unrelated files were removed.
