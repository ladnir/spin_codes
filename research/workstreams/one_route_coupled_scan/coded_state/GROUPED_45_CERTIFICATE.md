# A 45% checkpoint, then cheaper feedback

At K=2^20, N=2^21, and rate1/2 over F=GF(2^128), both retained outer/route
profiles now have complete certificates for symbol distance >45% and
setup-failure probability below2^-40. This applies first to group32 with
uniform feedback. It also applies to a simpler running-sum reset and to
the two-diagonal feedback specified below. These are research constructions,
not library promotions. Subsequent full-encoder timings and correctness
checks are recorded in [SANDWICH_NATIVE.md](SANDWICH_NATIVE.md).

## Construction and setup

The two profiles are:

| Outer over its indicated subfield, extended to F | Blocks | State and batch b |
|---|---:|---:|
| RS[512,256] over GF(2^16) | 4096 | 256 |
| RS[256,128] over GF(256) | 8192 | 512 |

Keep the existing route: output coordinate j from every outer block goes
to region j, which independently permutes its block labels uniformly.
Concatenate regions, then divide them into b-symbol batches X_t. There
is one global structured routing stage. Setup is sampled once and fixed
for every message; the minimum-distance event quantifies over all messages
after setup.

The group32 control retains the reset maps in CHEAP_RESET_COORDINATES.md,
separately on every32-symbol group of the b-symbol state. Its feedback is
an independently sampled uniform b-by-b matrix U_t over F at each batch.

The simpler variant has state Q_0=0 and recurrence

    Y_t = X_t + U_t Q_t,
    Q_(t+1) = Q_t + Y_t.

Output every Y_t and leave the terminal state free. This is the d=1 case
of the grouped-reset proof: each state coordinate accumulates its outputs.
The total state dimension remains b, not1. The inverse reconstructs Q from
Y alone, then recovers X using the same keys. Thus the encoder is linear
and injective for every setup.

## Complete counting bound

Let c excursions occupy v batches, with r=c-1 required closed excursions
and a free final excursion. A closed excursion has at least two batches.
The number of lengths and placements is at most

    C(v-c,c-1) C(n-v+c,c), where n=N/b.

On every valid output support, each closed excursion imposes at least d
independent reset equations and has weight at least d+1. For group32 this
follows from the distinct-point MDS syndrome in any active group. For d=1,
a nonzero closed output has at least two symbols, and its sum-zero constraint
has rank at least1. Different excursions use disjoint output coordinates.

For the inverse bound, at most bv output coordinates contribute, the support
weight is at most W=floor(.45N)=943718, and the reset rank is at least dr.
Every noninitial batch has nonzero entering state. Its inverse input is
uniform on F^b, independently of the other feedback keys for fixed Y.
The randomized input coordinates number a=b(v-c). Consequently the bound
uses output volume divided by q^(dr), and the same affine outer-membership
bound q^(D(A)-a) as the retained proof. Here D(A) is the outer's shortened
dimension on these a coordinates.

For the forward bound, the outer input is supported in the union U of the
v batches. Count projective nonzero outer inputs in U, then charge q^-b for
each of the r prescribed live-state resets. Average each complete bound
over the actual regional route, using the retained exact fractional moments.
The sets A and U have at most2c partially selected regions. Take the smaller
complete forward or inverse bound per (v,c), then sum; no independently
averaged probabilities are multiplied.

The verifier bounds all feasible pairs, including every placement, support,
and field amplitude. It does not sample messages or discard intermediate
classes. Group activation refinement is unnecessary at45%.

## Directed verification receipts

`grouped_reset_certificate.py` obtains a FLOAT partition from the existing
coarse grouped bound. It then independently checks adjacency on every
integer c slice, including the exact feasible lower endpoint v=2c-1.
Removed, duplicated, and shifted-boundary controls are rejected.
Every rectangle is replayed with exact rational regional moments and
outward-rounded Decimal logarithms at both80 and120 digits.

| Outer / b | Reset d | Rectangles | Feasible pairs | Worst directed log2 class bound | Conservative union margin |
|---|---:|---:|---:|---:|---:|
|512 /256|32|300|16781312|-69.7935898514629|44 bits|
|256 /512|32|309|4196352|-65.7748106558776|42 bits|
|512 /256|1|814|16781312|-66.4619787772504|41 bits|
|256 /512|1|462|4196352|-65.7748106558776|42 bits|

The final column floors the negated worst bound and subtracts a conservative
cell-count budget of25 or23 bits, respectively. We claim40 bits uniformly.
The partition hashes, in table order, are:

    169be10f1d42cd3e1edaa53b84ac5d8de4e7c2a9653ff8cd5cf8ec741966f50d
    183207d2f0bccb98bdbd2c0d33ba75419e97cecb6ab9f626736defc420901bf3
    423d8d6de48600b37cf0dec86efad21ceb491eb45c9870d6d58f70b385992909
    6df48246f68dadb282a7d64313ffb6d711c69373a7c45520fdc9a51671ef6fc0

Reproduce from this workstream, serially:

    python -B coded_state/grouped_reset_certificate.py --outer 512 --group 32 --precision 80

Repeat for outer256, for precision120, and for group1. Each proposal and
replay phase is capped at60 seconds by default (maximum120). No generated
reports, raw data, or on-disk coefficient caches are written. A cap or
failed rectangle rejects the run; a partial replay is not a certificate.

## Replace dense random feedback with two diagonal scalings

Fix a b-by-b matrix D over GF(2^16) whose every square minor is nonsingular.
For example, choose2b distinct subfield points alpha_i,beta_j and set

    D[i,j] = 1/(alpha_i+beta_j).

The Cauchy determinant gives the required minors. Equivalently, use the
fixed interpolation map between two disjoint additive cosets, with its
row/column scalings. It can be applied by small-coefficient transforms.

At each batch, independently sample u_t,v_t from (F\{0})^b, independently
across coordinates, batches, and the route. Replace U_t by

    R_t = diag(u_t) D diag(v_t).

Keep either certified reset unchanged. In particular, the simpler complete
inner is

    Y_t = X_t + u_t .* D(v_t .* Q_t),
    Q_(t+1) = Q_t + Y_t.

All coins are again fixed at setup. No per-message randomness is required.
The diagonals use full-field nonzero elements, not small-field coefficients.

### Point-probability lemma

Fix z!=0 with support size s, and fix a proposed feedback output f of
weight w. Put h=b-w. The first diagonal makes the s nonzero entries of
diag(v)z independent uniform elements of F\{0}. The h zero coordinates
required of D diag(v)z impose a matrix with rank min(h,s).

If h>=s, those zero equations force the first-scaled vector to zero, which
is impossible. Otherwise, solve h selected input coordinates in terms of
the other s-h. There are at most(q-1)^(s-h) admissible choices, so the
probability of those zeros is at most(q-1)^-h. Requiring all remaining
coordinates to be nonzero can only reduce it.

For every compatible first-scaled vector, the w prescribed nonzero output
values determine w independent entries of u. Their probability is
(q-1)^-w. Therefore, for every z!=0 and every f,

    Pr[R_t z=f] <= (q-1)^-b = eta q^-b,
    eta=(1-1/q)^-b.

This is a pointwise bound, not a claim of an exactly uniform vector.
It holds after conditioning on the past because current diagonal keys are
fresh. The local proof extends from the subfield to F without assuming
that state coordinates lie in the subfield.

### Transfer to the whole-code certificate

The forward argument uses at most n prescribed live-state outputs, so
its bound increases by at most eta^n. For fixed Y, the inverse state is
deterministic and independent of all keys; its noninitial input batches
have independent point densities bounded by eta q^-b. Summing any affine
outer fiber therefore incurs at most the same eta^n factor. Fractional
powers nu<=1 do not increase that factor. Thus the minimum of the two
complete bounds and their class union increase by at most

    eta^n = (1-1/q)^-N < 2,

where Bernoulli's inequality gives (1-1/q)^N >= 1-N/q > 1/2.
The uniform-feedback receipts above have at least41 bits even for d=1.
They therefore imply failure below2^-40 with the actual two-diagonal
feedback. This transfer is analytic; no assumption about computational
pseudorandomness replaces the ideal setup distribution.

Exact small-field checks enumerate15+4095 first-scaled vectors,252 support
laws, and74 square minors. The second diagonal is summed exactly by output
support, not sampled. Removing D or allowing zeros in the first diagonal
breaks the bound in retained negative controls. The test is
`coded_state/test_sandwich_feedback.py`; the algebra above proves the
actual field and dimensions.

## Work saved, and the next implementation gate

The running-sum variant removes the entire reset FFT/translation. Feedback
uses2N full-field products, N XORs for the output, N for the state update,
and one fixed b-symbol MDS map per batch. Prepared random keys occupy64MiB
in raw form at K20, versus8/16GiB for literal dense feedback in the two
profiles. Setup tables and any transformed key format are additional.

This also removes one MDS-map application relative to the retained exact
two-map uniform-feedback construction, at the same2N full-field products.
The subsequent native implementation takes about34.5ms end to end for
RS256/state512; see SANDWICH_NATIVE.md. The old approximately2.8ms and2ms
figures measure group32/group8 reset alone, not this encoder.

The native pass measures the fixed MDS map, full-field products, complete
inner, and alternative outer/route compositions separately.
Do not infer an end-to-end win from transform counts. Do not use small-field
random diagonals without a new proof: the counting lemma specifically uses
q-1 choices for each solved coefficient.

## Distance as an outcome

Keep the45% checkpoint; treat larger distance as a measured proof outcome.
With d32, the same coarse FLOAT cover completes at47% for both profiles
(742/843 rectangles). It fails its chosen per-class threshold at48% and49%.
With d1,45% closes but47% already fails the coarse bound. These failures
are not low-weight codewords or distance ceilings. Group8/16 also pass
the45% FLOAT screen; no directed receipts were needed after d1 closed.

Smaller running-sum states64 and128 fail this coarse45% bound under both
outers. This is a proof diagnostic, not a structural impossibility.
The next cost/distance comparison should vary state and reset coupling
together rather than demand49% from every cheap candidate.
