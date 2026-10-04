# RS8/GF32 at K18 with the exact s22 inner

The complete RS8/GF32 construction at K=262144 has
**58.245625332 bits** of setup-failure margin at 10% relative distance.
The certificate covers every active-group count q=1,...,2048.
It keeps the exact retained t64/s22 inner and four-bit packet route.
All saved sparse support endpoints are refolded with the new outer measure;
the authenticated q>=3 endpoints transfer unchanged.

## Construction and claim

There are L=2048 groups, each mapping 128 message bits to 256 code bits.
Interpret each input as eight rows of four GF(16) symbols.
Encode each row with the same fixed systematic MDS RS[8,4] generator.
Collect the eight symbols at each code position into one 32-bit column,
using a fixed row order and binary basis.

Independently in every group and column, sample a nonzero GF(2^32) coefficient
and multiply the column by it. The binary adjoints of these multiplication
maps give the same expected-measure argument. An implementation must declare
which orientation defines its forward map.
All eight columns, including systematic columns, receive independent randomizers.

Split each column into eight four-bit packets. Independently shuffle all
64 packets within each group. Route one packet per group to each region.
Independently in each of the 64 regions, permute the 2048 groups uniformly.
Concatenate the regions, with sixteen packets per physical 64-bit step.

Let A be the retained binary 64-by-22 expansion map and let C=A^T.
Their exact map identity is
`e2467b9d5de8403ccd8263f0f55c73df128a8a2dd4c4cc405bc3a8e057fc3581`.
The authenticated receipts contain all expansion rows and feedback columns;
the new receipt copies this complete map record.
For input x and incoming state a, each step emits x+A*a and sets the next
state to M*a+C*x, with M sampled independently from GL(22,2).
The state starts at zero, persists across every step and region, and is
discarded without a flush. Two chronological steps form each proof macro;
there are 64 macros per region.

Outer randomizers, group shuffles, regional permutations, and state updates
are mutually independent. Setup is sampled once and fixed across all messages.
Let E be the resulting binary linear encoder from 262144 to 524288 bits.
Over this ideal setup distribution, the certificate bounds

\[
 \Pr[\exists m\ne0:\operatorname{wt}(E(m))\le52428]<2^{-58}.
\]

Thus a successful setup gives an injective rate-1/2 code with minimum distance
at least 52429. The stored exact upper bound has about 58.245625332 bits of margin.
The claim concerns the ideal ensemble, not the distance of a specific seed.
It does not reach the 60-bit target of the retained RS16 outer.

## Exact outer measure

Write Q=2^32. Shortening the eight-row RS code to a specified set of h columns
leaves max(h-4,0) GF(16) dimensions per row. The shortened code therefore has
Q^(h-4) words for h>=4. Inclusion-exclusion gives the column-support counts

\[
 A_h={8\choose h}\sum_{j=0}^{h-5}(-1)^j{h\choose j}
          (Q^{h-4-j}-1),\qquad 5\le h\le8.
\]

Every fixed nonzero column becomes a uniform nonzero 32-bit vector under
its independent scalar map. The same statement holds for binary adjoints:
for fixed nonzero x, the map c->M_c^T*x is a binary linear bijection.
Its kernel is zero because M_c is invertible for c!=0.
Define the exact expected packet-support counts C_v by

\[
 P(z)=(1+15z)^8-1,\qquad
 \sum_{v=5}^{64}C_vz^v=\sum_{h=5}^8 A_h\left(\frac{P(z)}{Q-1}\right)^h.
\]

These counts satisfy sum_v C_v=2^128-1. They are expectations of message
counts over outer setup, not deterministic spectra or probabilities.
For each specified packet support, the complete expected measure gives
equal weight to all nonzero packet-label assignments.
After the independent group shuffle, every vector with v active packets
has expected count C_v/({64 choose v}15^v).
Independent setup between groups permits products of these measures.
No independence between different messages under one setup is needed.

## Complete sparse folds

For v=0,...,64, let W_v be the saved outward upper bound for one group
with v active packets. For 0<=u<=v<=64, let V_(v,u) be the analogous
two-group bound. Each bound averages the uniform supports, nonzero labels,
regional route, and inner updates, at output cutoff 52428.
The two-group bound retains the common route and continuous inner state.

The saved producer computes all 65 one-group endpoints and all 2145
triangular two-group endpoints before folding the old outer counts.
Its support calculation does not use those counts.
Consequently, its endpoints at supports five through eight apply to RS8
even though those supports have zero count under the old RS16 outer.
An endpoint that remained at the trivial bound one is still valid.

The new first moments are

\[
 U_1=L\sum_{v=0}^{64}C_vW_v,
\]

\[
 U_2={L\choose2}\sum_{v=0}^{64}\sum_{u=0}^v
             (2-\mathbf1_{u=v})C_vC_uV_{(v,u)}.
\]

The binomial coefficient chooses an unordered pair of distinct groups.
Off-diagonal terms count both assignments of the two supports to that pair.
Diagonal terms count one support assignment. No factor two applies there.

The receipt separately records the q1 contribution from supports five
through eight and the q2 contribution with at least one such support.
Both charges are already included in U_1 and U_2. No sparse support is
discarded, and no separate event excludes low-support outer words.
Their margins are 59.250475311 bits for q1 and 136.918019804 bits for q2.

## Transfer of every tail endpoint

Fix a nonzero 256-bit atom before the group shuffle, with h nonzero columns.
For a specified h-column support, at most (Q-1)^(h-4) codewords have exactly
that support: the shortened MDS code is determined by any h-4 coordinates,
each required nonzero. Independent column randomization divides the expected
atom count by (Q-1)^h. Atoms with fewer than five columns have count zero.
Every atom thus has expected count at most (Q-1)^-4.
The zero atom also has count zero in the nonzero-message measure.

Relative to the uniform measure U on 256-bit strings, the new group measure
is pointwise bounded by beta_new*U, where

\[
 \beta_{\rm new}=\frac{2^{256}}{(2^{32}-1)^4}
 <\frac{2^{256}}{(2^{16}-1)^8}=\beta_{\rm old}.
\]

The independent shuffle preserves this envelope. Fix any active-group set
and route. Product domination holds before applying the unchanged inner.
Every retained ordinary uniform-input bound therefore remains valid.

For the routing-conditioned bounds, condition on the same potential slot
counts in each physical step. Outer randomizers retain their independent law.
The conditional expected bad-message count still has the old upper bound A(J).
Hence the conditional existence probability remains at most min(1,A(J)).
The retained fractional-power bound on this quantity remains valid, including
its chronological matrix products, routing weights, and group-set union.
No fractional endpoint is rescaled. The q>=3 endpoints cover every integer
through 2048, with no interpolation or missing occupancy.

The certificate sums the new outward q1/q2 endpoints and all 2046 retained
tail endpoints exactly. Floating-point logarithms are used only for display.

| Contribution | Margin in bits |
|---|---:|
| q1 complete exact-shell fold | 58.245633869 |
| q2 complete exact support-pair fold | 126.092935394 |
| Retained q>=3 union | 75.612279616 |
| Complete union | **58.245625332** |

## Authentication and reproduction

[k18_rs8_certificate.py](k18_rs8_certificate.py) first runs the independent
[K18 receipt audit](../k18_margin60/verify_fractional.py).
That audit authenticates 112 files, four original components, all 2048 old/new
selections, and both exact unions. The new producer additionally authenticates
the sparse receipts by raw SHA256, checks their exact map and route metadata,
and independently refolds their original RS16 counts as a regression.
It pins the unchanged enumerator and rational arithmetic helpers.
All eight tests pass. An independent review checked the complete support
interfaces, pair multiplicities, exact union, and routing-conditioned transfer.
An independent fold using an alternate MDS formula reproduced the same margin.

The retained whole receipt is `tmp/k18-margin60-fractional-whole-p256-v2.json`,
SHA256 `44fcdb0389d5cc7be870bdfcbeb0ec5500193706b7e690313054971d052a01d0`.
The sparse q1 and q2 receipts are
`tmp/rs16-s22-k262144-whole-fresh-p256-q1.json` and
`tmp/rs16-s22-k262144-whole-fresh-p256-q2.json`, with SHA256 values
`bc55924c4c4708b293a9bb1e1baacb81344e2a7fd7a25a0546ad85761c3f880c` and
`5e411ebefb9c5c8d48c22e20fc070fa98a4cfa01213231e7cc75045aeaed37d9`.
All inherited pins and the full per-occupancy endpoints appear in the output.

```powershell
C:/Python314/python.exe -B -m unittest discover -s research/workstreams/rs8_sizes -p test_k18_rs8_certificate.py -v
C:/Python314/python.exe -B research/workstreams/rs8_sizes/k18_rs8_certificate.py --output tmp/k18-rs8-gf32-s22-fresh.json
```

The calculation needs no new local census and runs no benchmark.
Saved local operator and fractional endpoints are authenticated, not regenerated.
Existing producers and receipts remain unchanged.
Next verify that a proposed implementation preserves this exact map,
randomizer orientation, packet route, and ideal sampling distribution.
