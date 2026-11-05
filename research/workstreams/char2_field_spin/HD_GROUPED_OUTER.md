# Grouped small-RS outers for the large-alphabet code

## Construction and objective

This work replaces the independent RS[256,128] outer of the certified
[35% checkpoint](HD35_OPTIMIZATION.md) with grouped small constituents.
The initial comparison leaves its 64-symbol, 16-state inner unchanged.
The target is faster complete K20 encoding, rate1/2, at least35% GF65536
symbol distance and at least40 bits of whole-code setup margin.

Write F=GF65536 and q=65536. Choose g in {2,4,8,16}, n=256/g, and k=n/2.
Each group contains g rows of k input symbols. Encode every row with the
same fixed MDS RS[n,k] map, using evaluation points in the GF256 subfield.
At column j, collect the g row outputs into v_j in F^g. Independently for
every group and column sample nonzero diagonal entries in D1_j and D2_j.
Emit w_j=D2_j H_g D1_j v_j, where H_g interpolates from nodes0,...,g-1
to g,...,2g-1. All its square minors are nonzero.

Split the column outputs into individual F symbols. Coordinate (j,a)
enters region g*j+a. Each of256 regions independently permutes all8192
groups. There is no extra local shuffle. The local and regional layouts
are fixed before messages; the regional permutations are independent of
all outer and inner randomness. All setup is sampled once and reused.

The native representation applies this same F-linear scalar encoder to
eight SIMD lanes per128-bit record. Thus128 input records yield256 output
records per group. Column mixing couples distinct records, not those lanes.

## Outer-measure transfer

The [diagonal-Cauchy lemma](HD_DIAGONAL_CAUCHY_REFRESH.md) applies with
state dimension g: for fixed nonzero v,y in F^g,

    Pr[D2 H_g D1 v=y] <= rho,       rho=(q-1)^(-g).

Every sampled map is invertible. Fix a nonzero desired complete group
output y and let S be its h nonzero columns. Zero columns before and
after mixing coincide. If h<n-k+1, no nonzero message can produce y.

Otherwise fix any h-(n-k) columns of S. Shortening the RS code to S
leaves dimension h-(n-k) in each constituent row, and these columns
form an information set. Condition only on their mixers. Their prescribed
outputs determine their preimages and hence at most one group message
consistent with the zero columns. If no such message exists, the target
is impossible. If it exists but another column in S is zero, it is also
impossible. In the remaining case all n-k unconditioned columns have
fixed nonzero inputs. Their independent mixers each hit their required
target with probability at most rho. Therefore

    E[number of nonzero messages producing y]
       = Pr[y belongs to the sampled group code]
       <= rho^(n-k) = (q-1)^(-128).

Equality of the expectation and membership probability uses injectivity
of every sampled group encoder. It does not assume independence between
messages. The bound also holds at y=0, whose nonzero-message count is zero.
Consequently the expected measure is dominated by beta times uniform F^256,

    beta=q^256/(q-1)^128.

This is exactly the old RS256 envelope, not a density loss per column.
Conditioning on the information columns avoids multiplying a loose
per-message bound by the entire shortened-code multiplicity.

Independent group setups give product domination for every fixed active
group set and fixed route. The existing regional coefficient and global
marker witnesses depend on the outer only through beta and the256-region
geometry. The dense inverse-output tail uses the same pointwise atom cap,
injectivity, F-linearity, and independence from the inner. These properties
are all preserved. Thus the current35% bound of76.7730 bits and the40%
bound of56.6591 bits transfer with zero loss when their respective inner
distributions are retained. No new numerical witness search is required.
This is an ideal-setup theorem, not a certificate for a deterministic PRNG
seed or a claim about implementation performance.

Unlike the binary sparse-shell proof, these certificate paths do not need
a weight-symmetric exact local enumerator. The uniform-measure comparison
already includes every vector. A random local shuffle is therefore not
required. Omitting it does not retain intact columns in global routing:
their g components enter g distinct regions with independent permutations.

## Failure modes and checks

One active constituent row occupies at least n-k+1 columns. Each such
column becomes full-support, so g=4 yields at least132 outer symbols.
Several active rows can cancel; the pointwise proof above handles them.
A common fixed cross-row mixer can be absorbed into the row input basis.
One random output diagonal around that mixer does not satisfy the claimed
bound. The tests retain this variant as a negative control.

The exact GF4 toy checks enumerate all nonzero messages and all65535 target
vectors for two RS[4,2] rows. They verify the expected-count cap including
equality cases. Native tests independently evaluate RS polynomials and
Lagrange interpolation; no SIMD FFT is used by their oracle.

The main performance risk is the extra diagonal and fixed column map.
Shorter RS circuits alone do not prove a full-encoder speedup. Benchmarks
must include both randomization stages, the actual regional route, the
complete inner, and output stores. Existing controls remain unchanged.

## Progress

### Cheaper representations of the column randomizer

Write a for the first diagonal entry of D1. Replacing D1 by D1/a and
D2 by a*D2 preserves each matrix D2 H D1 exactly. Equivalently, sample
the first diagonal with first entry1 and independent nonzero remaining
entries, and sample an independent full nonzero second diagonal. The
resulting matrix distribution is unchanged: the omitted scalar has q-1
equally likely preimages for every normalized tuple. This removes one
GF65536 multiplication per column. It does not remove an entire diagonal.
Exact GF4 checks compare every atom before and after normalization.

For g=2, another proved family multiplies by a uniform nonzero scalar
in E=F[v]/(v^2+v+0x2000). In the existing tower F=GF256[u]/(u^2+u+0x20),
0x2000=0x20*u has absolute trace1, so the extension is a field. This is
an F-linear map on two records, with image uniform over F^2 excluding0.
Its group atom cap is (q^2-1)^(-64), no larger than (q-1)^(-128).
It therefore preserves the same certificate. The native arithmetic
compares a four-product matrix circuit with a three-product Karatsuba
circuit; neither requires repacking the eight SIMD lanes.

The first grouped implementation retained region-major diagonal storage.
It took6.9--8.2ms, versus a4.72ms control. Packing the coefficients in
four-group consumption order reduced the g4/g8 full time to about5.06ms,
versus4.90ms in that run. These are preliminary101-call measurements,
not a speedup. Later301-call timings coincided with an unrelated user's
benchmark on the shared host and are not used as final evidence.
Normalization and the three-product circuit remain to be measured cleanly.
The normalized first diagonal is also packed without its identity entries;
this reduces its coefficient traffic by1/g. Full K20 independent dense
checks pass for g4 and g8 diagonal mixers and the g2 field multiplier.
Normalized and unnormalized diagonal versions have identical checksums.

### Inner co-design screen

A smaller packet/state pair (32,8) with four RS[8,2] emission blocks keeps
the same common-refresh recurrence. It would lower the fixed refresh-map
depth per output. The floating single-packet profile screen does not close:
at128/512/2048 active groups the global-marker log2 bounds are approximately
+2190/+5318/+26336. Retaining regional coefficients still gives+520/+3472
for128/512 groups. These are failed upper bounds, not counterexamples or
impossibility claims. No implementation or certificate is claimed for this
smaller inner. Retain the proved (64,16) inner while optimizing the outer.

### Validation and next experiment

All seven exact grouped-outer tests pass. Native Release checks cover
g=2,4,8,16 at K=8192 and24576, unaligned input, all four output alignments,
cached/streaming output, guard preservation, and exact input/output alias.
The normalized diagonal versions also pass ASAN/UBSAN. K20 independent
checks give:

| Column mixer | Constituent | K20 output checksum |
|---|---|---:|
| Normalized diagonal-Cauchy, g4 | RS[64,32] |1717454094421689829|
| Normalized diagonal-Cauchy, g8 | RS[32,16] |5609734393776438674|
| Exact field multiplier, g2 | RS[128,64] |15439270946252110526|

An optional final-stage fusion now leaves each RS size-two butterfly to
the column consumer. It computes two columns in registers, randomizes them,
and stores them directly, avoiding the final intermediate store/reload.
All row loops are compile-time unrolled. The final-stage coefficient is
shared across the fixed row bundle; the sixteen-register prefix remains
bounded. The g16 path retains its old schedule to avoid excess register
pressure. This is an exact-map implementation change, not a new ensemble.
The g2/g4/g8 checksums are unchanged. Release and ASAN/UBSAN checks pass,
including the three-product field kernel. Full K20 dense checks also pass.

An unrelated benchmark initially prevented measurements. The user later
authorized CPU 2 measurements while that job remained pinned to CPU 0.
All our benchmarks still run serially under both benchmark locks. No faster
code has yet been established. Fusion remains disabled by default.

The 2026-10-08 comparison used 501-call full-encoder runs at K20, streaming
stores, and the original h2 inner. Control medians were 4.734--4.808 ms.
The normalized g4 candidate measured 4.923 ms, or 4.865 ms with final-stage
fusion. Normalized g2 measured 5.027 ms; g8 measured 8.263 ms. The g2
three-product field multiplier measured 5.110 ms, or 5.359 ms when fused.
These one-invocation candidate results do not establish a win. They retain
the same certificate but include possible shared-cache contention. The
[GF256 continuation](HD_GROUPED_SUBFIELD.md) records repeated controls,
the byte-field candidates, and phase timings.

Build the existing isolated experimental CMake project and target
`spin_char2_hd_grouped_forward`. Default configuration enables normalized
diagonal mixers. Define `HDGROUP_NORMALIZE=0` for the exact-map arithmetic
control. `HDGROUP_FIELD32=1` selects the four-product extension-field kernel
for g2, and value2 selects its three-product equivalent. Other g values
continue to use diagonal mixers. `HDGROUP_FUSE_LAST=1` enables the final-stage
fusion. Run `--check` and `--large-check` before
timing `--bench 501 full 20 nt 4 2`, changing4 to the selected group width.
These are research selectors, not production API options.

### Cheaper subfield coins: initial screen and subsequent closure

The [subfield continuation](HD_GROUPED_SUBFIELD.md) now resolves the sparse
case below with an inverse-refresh bound for empty packets. The complete
GF256-setup certificate gives72.3975 bits at35%, with256/384-bit agreement.
Its byte-multiply native candidate passes full K20 checks. Repeated g4
measurements are essentially tied with the old control, not a clear speedup.
The initial failures below are retained to explain the proof change.

A separate construction screen changes every random scalar to GF256,
including both outer and inner diagonals. All fixed matrices already have
GF256 entries. This would make random scalar multiplication a single GFNI
byte multiply. It is not a drop-in certificate transfer: the field size in
the measure bounds must also change from65536 to256.

If this GF256 code has symbol distance d, its scalar extension to GF65536
has distance at least d. Expand a nonzero GF65536 message in a GF256 basis.
At least one component message is nonzero; its encoded support is contained
in the complete output support. Thus the smaller field would still serve
the original eight GF65536 lanes per record, without a union loss over lanes.

The screen also avoids charging the refresh density penalty on transitions
to zero. For a local return cap Jbar and d=q^s-1, use

    [[K,       b(T_j-K)],
     [Jbar/d,  b V]],       b=d/(q-1)^s.

This follows from the local weighted operator in HD_COMMON_REFRESH.md:
the refresh preserves zero, and only nonzero next-state densities need b.
Each matrix entry is raised to its fractional power after this factor is
included. No uniform density penalty is added again at the end.
The old GF65536 screen remains the default and its controls are tested.

At35%, q=256, width64 and state16, FLOAT proposal bounds are:

| Active groups | Two-check local map | Four-check local map | Bound used |
|---:|---:|---:|---|
|1|+86.337|+86.330|Regional coefficient|
|16|-2778.63|-2789.06|Regional coefficient|
|128|-1748.11|-1794.35|Regional coefficient|
|512|-2455.77|-3611.18|Global marker|
|2048|+2331.07|-17718.26|Two-check regional; four-check global|
|4096|-53548.75|-90060.21|Two-check regional; four-check global|

These initial log2 upper bounds were not distance certificates. The local-rank
change repaired the sampled medium-density obstruction, but not the
one-active-group bound. Small classes used multiple optimizer starts, so
the failure was not just the original low-alpha initialization. The subsequent
inverse-refresh argument supplies the missing sparse-class analysis; its
complete certificate and native kernel are documented in the subfield note.
Keep this new distribution separate from the proved GF65536 setup path.

Nine grouped-outer regression tests pass, including the subfield screen's
normalization and an equality check against the old field's global bound.
The ten existing35% proof tests also pass after the optional field parameter.

- Transfer argument and exact finite-field tests added.
- Isolated native g=2,4,8,16 prototypes use the same regional layout and
  inner as the35% control; validation passes and measurement continues.
- No production profile, paper timing, or existing certificate is replaced.

## Invertible Toeplitz column mixers

The outer has two separate jobs. The small RS constituents force a nonzero
group to occupy many columns. The column mixers then prevent a prescribed
complete output from having too much probability over setup, even when the
constituent rows are correlated. Neither job requires a dense random matrix
on the whole group, nor a uniformly random image after each column mixer.
The pointwise bound below is the interface consumed by the inner proof.

This separation lets us trade RS work against mixer work. Two RS[128,64]
rows need three products per column; four RS[64,32] rows need nine products
per column but shorter RS transforms. Both retain the same conservative
whole-group envelope. Their choice can therefore be made by measuring the
complete encoder, without paying a new proof penalty for either width.

The next candidate changes only the column mixer. Keep g constituent
RS[256/g,128/g] rows, 256 independently permuted regions, and the selected
certified inner. For every group and column independently, sample a g-by-g
Toeplitz matrix T uniformly conditioned on invertibility. The prototype
supports g=2 and g=4, over either GF256 or GF65536. Matrices are not shared
between groups or columns. Setup fixes all matrices before encoding.

Sample the 2g-1 diagonals independently from the entire field, including
zero, and reject singular matrices. Invertibility preserves zero columns
in both directions and makes the complete grouped outer injective.

### The existing outer bound transfers

For an unconditioned random Toeplitz matrix and a fixed nonzero x in F^g,
Tx is uniform in F^g. To see this, take the last nonzero coordinate of x.
Condition on all diagonals except the g diagonals that multiply this
coordinate in the g outputs. Their map to the output is triangular with
nonzero diagonal. The earlier [Toeplitz outer proof](K20_RANK_CONTINUATION.md)
uses this same argument.

The fraction of invertible g-by-g Toeplitz matrices is 1-1/q; see
[Kaltofen and Lobo, ISSAC 1996](https://kaltofen.math.ncsu.edu/bibliography/96/KaLo96_issac.pdf).
Thus, for every fixed nonzero x and every y,

    Pr[T x = y | T invertible] <= q^(-g)/(1-1/q)
                                = 1/((q-1) q^(g-1)) =: rho_T.

This is an upper bound, not a uniform-image claim after conditioning.
Since rho_T <= (q-1)^(-g), the existing conditional-information-set
argument applies unchanged. Put n=256/g and k=128/g. For every fixed
nonzero complete group output y,

    Pr[y belongs to the sampled group code] <= rho_T^(n-k)
                                             <= (q-1)^(-128).

Specifically, its h nonzero columns determine at most one source after
conditioning on h-(n-k) information-column mixers. Every other nonzero
column has an independent mixer. Invertibility excludes a nonzero preimage
of any zero output column. Each of the remaining n-k columns costs rho_T.

Independent group setups retain the same product domination. We can use
the old, weaker envelope unchanged in every numerical witness. Therefore
the GF256/h4 construction retains at least 72.3975 bits at 35% distance;
the GF65536/h2 construction retains at least 76.7730 bits at 35%. The outer
also preserves the 56.6591-bit, 40% result when its original GF65536 inner
is retained. The native prototype below does not benchmark that 40% inner.
These transfers concern ideal setup randomness, not the test PRNG seed.

### Fixed small circuits

For T=[[a,b],[c,a]], compute

    p = a(x+y),  out0 = p+(a+b)y,  out1 = p+(a+c)x.

This costs three variable products and three XORs, with coefficient sums
prepared during setup. For g=4, partition T=[[A,B],[C,A]] into Toeplitz
two-by-two blocks. Compute A(x_low+x_high), (C+A)x_low, and (B+A)x_high.
Combining these products gives both output halves. The resulting circuit
costs nine variable products and fifteen binary XORs.

The normalized diagonal-Cauchy g4 circuit instead uses seven variable
products and five fixed GF256 products. Over GF256 that is twelve GFNI
products versus nine; over GF65536 the variable products require two GFNI
products each, so the comparison is nineteen versus eighteen. The Toeplitz
schedule also streams more variable coefficients. Product counts therefore
do not by themselves imply a native speedup.

The implementation selects these circuits with `HDGROUP_TOEPLITZ=1` in
the existing grouped-encoder target. `HDGROUP_FUSE_LAST=1` keeps the
established last-butterfly fusion. `HDGROUP_SUBFIELD=1` selects the GF256
ensemble; the default field remains GF65536. No hot-path rejection,
allocation, indirect call, or runtime mixer selection is introduced.

The optional `HDGROUP_BYTE_KEYS=1` stores each GF256 outer coefficient as
one byte instead of the four-byte GF65536 prepared representation. It
applies to either column mixer and leaves the sampled map unchanged.
Four adjacent coefficients load as one word and expand with one byte
shuffle. The inner coefficient representation remains unchanged.

The scalar tests independently enumerate every mixer atom for dimensions
2, 3, and 4 over GF2 and GF4. For g4/GF4, 12,288 of 16,384 matrices are
invertible and the largest conditional atom is 64/12,288=1/192, attaining
the bound. The tests also check the complete two-row RS[4,2] expected
output measure, the recursive circuit on every GF4 coefficient tuple,
and GF256 arithmetic. All 31 new and existing targeted tests pass.

### Native comparison, 2026-10-08

The compact-key comparison uses Peach CPU 2 at its existing 4.5 GHz setting.
An unrelated single-threaded job runs on CPU 0, sharing L3. All our timed
processes run serially under both benchmark locks; no builds or correctness
tests overlap them. Each invocation measures 501 full calls after three
warmups, with setup, allocation, and input generation outside the clock.
K=2^20 input records, N=2^21 output records, and each record is 128 bits.
Output stores are streaming. These are forward-encoder measurements.

We reversed candidate order and bracketed the runs with the frozen RS256,
GF65536/h2 control. The following table includes all five compact-key runs
per candidate, not just the best invocation.

| Encoder | Median of run medians, ms | Range of run medians, ms |
|---|---:|---:|
| Frozen RS256/GF65536/h2 control | 4.757 | 4.745--4.792 (eight runs) |
| GF256/h4, four RS[64,32], compact DHD keys | 4.492 | 4.443--4.584 |
| GF256/h4, two RS[128,64], compact Toeplitz keys | 4.460 | 4.443--4.597 |
| GF256/h4, four RS[64,32], compact Toeplitz keys | 4.387 | 4.369--4.503 |

The five g4 Toeplitz medians in order were 4.502577, 4.501425, 4.387462,
4.368907, and 4.376922 ms. Compact DHD medians were 4.492047, 4.572127,
4.583979, 4.446191, and 4.443005 ms. The g2 Toeplitz medians were 4.455479,
4.596933, 4.511684, 4.443136, and 4.459987 ms. Variability remains visible;
use about 4.4--4.5 ms for the leading g4 candidate rather than a precise
isolated-host speedup. It is about 5--8% below the frozen control. Its small
advantage over compact DHD is less decisive than the compact-layout gain.

The initial Toeplitz circuits with four-byte prepared coefficients did not
win: g4/GF256 took 5.024 and 5.520 ms, while g2/GF256 took 4.937 and
5.068 ms. GF65536 Toeplitz took 4.773 ms for g2 and 5.344 ms for g4 in a
single screen. Keep these negative results: fewer scalar products did not
automatically yield less native work. Inspection found no vector spills in
the g4 GF256 mixing loop, but every variable product loaded and expanded a
coefficient. Compacting the GF256 coefficients helped both mixer families.
The old g4 DHD path, rerun in the final comparison, took 4.721 and 4.956 ms.

Prepared outer coefficients occupy 14 MiB in the old g4 DHD representation,
3.5 MiB in compact g4 DHD, 4.5 MiB in compact g4 Toeplitz, and 3 MiB in
compact g2 Toeplitz. This accounts only for coefficients read by the outer,
not the route, data, inner setup, or scalar oracle data retained by the
research driver. The compact representation preserves each sampled map;
it is not a new heuristic distribution.

Separate 301-call diagnostic medians were 1.882 ms for compact DHD's outer,
1.835 ms for compact g4 Toeplitz's outer, 1.945 ms for compact g2 Toeplitz's
outer, and 2.625 ms for g2's routed h4 inner. These do not sum exactly to
the full timing because they are separate runs with different cache states.
The routed inner is now the larger phase. No further outer proof change is
required to optimize its implementation.

Both compact candidates pass Release and ASAN/UBSAN small-size checks,
covering the native supported group widths, both local ranks, alignment,
guards, cached/streaming output, and exact input/output aliasing. Independent
K20 dense-oracle checks pass for g2/g4 Toeplitz and g4/g8 DHD. Compact and
original coefficient layouts produce the same respective K20 checksums:

| Candidate | K20 full-output checksum |
|---|---:|
| GF256/h4, g2 Toeplitz | 2687453504903534300 |
| GF256/h4, g4 Toeplitz | 3966879751720036268 |
| GF256/h4, g4 DHD | 10655490065220255253 |

Reproduce the candidate with `HDGROUP_SUBFIELD=1`, `HDGROUP_FUSE_LAST=1`,
`HDGROUP_TOEPLITZ=1`, and `HDGROUP_BYTE_KEYS=1`. Build the existing
`spin_char2_hd_grouped_forward` target, run `--check` and `--large-check`,
then use `--bench 501 full 20 nt 4 4`. For the matched compact DHD control,
omit `HDGROUP_TOEPLITZ=1`; for g2 Toeplitz change the penultimate argument
from 4 to 2. No production defaults, library profiles, or paper numbers
change. Keep compact g4 Toeplitz as the leading research candidate and both
controls for future comparisons.

The rebuilt GF65536 DHD and three-product GF32 controls also pass their
Release check matrices. Measured compact binaries have SHA256:

- g2/g4 Toeplitz: `3291bbef178350bd207c5f9e78af66a5db4672a1982b1c71bdb5f887ae72b3d2`.
- DHD: `e10963b738b6eb28eb9baf3fb5f083ef7b52b6a3edafe4f92144501442aef19d`.

The measured `HdGroupedOuter.h` has SHA256
`15386b001a35ef24463900019ddc0e65bf2614e0ec2011bcc494ca0cd29fc060`;
`hd_grouped_forward.cpp` has SHA256
`2faf26025303b99427a24d827b3a308c4d747d02d40163b864fda1ba38c0c418`.
The local sources and uploaded sources matched these hashes at that checkpoint.

## Exact-map inner optimization

The next iteration keeps the complete sampled code unchanged. The outer,
regional permutations, h4 feedback/emission maps, random diagonals, and
natural output order are identical to the 35% / 72.3975-bit candidate above.
Only their representation and evaluation schedule change. In particular,
this is not another heuristic route or a certificate transfer to new maps.

The old inner computed emission for a whole 64-symbol packet, retained its
input vectors, then computed feedback. The new schedule processes one
four-point coset at a time across the four RS[16,4] blocks. It accumulates
feedback immediately, so it does not retain a sixteen-vector packet array.
Zero and identity constant products are removed at compile time.

Two further exact rewrites reduce data movement. First, prepare the same
route indices in the SIMD-striped order during setup. This changes when
the loads occur, not their destinations. Second, retain the state in that
striped representation between packets. State coordinate i resides in
vector i mod 4 and record lane floor(i/4). The refresh circuit and its
diagonal tables use the same coordinate permutation. This removes the
state and feedback transposes at every packet boundary; the output remains
in its original order.

The output can also supply feedback. The certified maps satisfy C=A^T and
CA=0. For Y=X+AQ, therefore CY=CX. Applying C to Y instead of X gives
exactly the same next state for every input and incoming state. This lets
the short-coset kernel release the input vectors after forming the output.
It does not rely on random messages or a zero incoming state.

The final optional fusion combines the adjacent size-sixteen stages of
the inverse/forward refresh FFTs. Write c for the fixed forward coefficient.
Their action on a pair is

    (a,b) -> (a+c(a+b), b+c(a+b)).

For a striped SIMD vector, one record shuffle forms the paired sum; one
fixed GFNI product and one XOR finish the update. This replaces three
record shuffles by one in each of four vectors. It leaves the refresh
matrix unchanged. The smaller instruction count yields only a small
timing change on the tested host, not another large speedup.

The research selectors are:

| Selector | Exact-map evaluation schedule |
|---|---|
| `HDGROUP_INNER=0` | Existing full-packet inner |
| `HDGROUP_INNER=1` | Short cosets, original route order |
| `HDGROUP_INNER=2` | Short cosets, striped route table |
| `HDGROUP_INNER=3` | Mode 2 with feedback computed from Y |
| `HDGROUP_INNER=4` | Mode 2 with persistent striped state |
| `HDGROUP_INNER=5` | Mode 4 with feedback computed from Y |

`HD_SUMRANK_BYTE_KEYS=1` stores the two sixteen-byte refresh diagonals in
32 bytes rather than 128; the hot K20 refresh tables shrink from 4 to
1 MiB. This did not establish a standalone speedup for the old schedule.
`HDGROUP_FUSE_REFRESH=1` enables the last fusion for modes 4/5.
`HDGROUP_PREFETCH=0` is a diagnostic control that disables next-region L2
prefetch. It is not selected: two reversed full-encoder runs took 8.500
and 8.461 ms, roughly twice the prefetch-enabled time. All selectors remain
compile-time choices, with no allocation or indirect dispatch in encoding.

Every new schedule passes the independent dense full-encoder oracle for
K=8192, K=24576, and K=2^20 with both g2/g4 outer variants. The small-size
checks also cover the rank-two fallback. They check cached/streaming output, unaligned input
and output, guards, and exact input/output aliasing. Modes 3/4/5 additionally
check every one of the 10,240 binary basis directions of a 64-record input
and sixteen-record incoming state for three setups/routes. Both the output
and next state match the original packet kernel in all 30,720 cases. The
full g4 K20 checksum remains 3966879751720036268.

### Measured result and retained checkpoint

The same CPU 2, 4.5 GHz, serial-lock protocol and CPU 0 shared-L3 caveat
apply. Full timings use 501 calls; separate phase timings use 301 calls.
The final comparison alternated the rebuilt candidate with the frozen
Toeplitz outer checkpoint and also reran the original RS256 control.

| Encoder | Final full-encoder run medians, ms |
|---|---|
| Frozen RS256/GF65536/h2 control | 4.738, 4.815 |
| Previous g4 Toeplitz, old inner schedule | 4.374, 4.480, 4.401 |
| New g4 Toeplitz, mode 5, compact refresh, fused stages | 4.018, 4.066, 3.973 |
| New g2 Toeplitz, same inner | 4.085 (one screen) |

Two earlier reversed runs of the selected g4 binary gave 3.974 and 3.980 ms.
Across all five runs its medians span 3.973--4.066 ms. Use about 4.0 ms as
the research timing: approximately 9% less time than the preceding grouped
checkpoint, and about 15--17% less than the RS256 control. Do not interpret
the last decimal as a stable isolated-host result.

The final separate phase medians are 1.869 ms for the outer and 2.040 ms
for the routed inner. The old routed inner measures 2.540 ms in the same
sequence, so this phase takes about 20% less time. Separate medians
need not sum exactly to the complete encoder's median.

The earlier screens explain the selected schedule:

| Change | Full median(s), ms | Outcome |
|---|---:|---|
| Compact refresh alone, old packet schedule | 4.496 | No standalone win established |
| Short cosets, natural routes | 4.204 | Promising first screen |
| Short cosets, striped routes | 4.226 | Similar first screen |
| Striped routes, feedback from output, natural state | 4.131; rebuilt 4.246, 4.269 | Improvement, but not the best state layout |
| Persistent striped state, feedback from input | 4.048, 4.041 | Removes more layout work |
| Persistent striped state, feedback from output | 4.043, 3.992; later 3.968, 4.006 | About 4.0 ms without the final fusion |
| Fused refresh, feedback from input | 3.996, 4.007 | Close to the selected version |
| Fused refresh, four-byte refresh coefficients | 3.908, 4.135 | Variable; not a demonstrated improvement over compact keys |

Thus the main win is the short-coset and persistent-state schedule, not
compacting refresh keys or the final butterfly fusion in isolation.

The selected mode 5 build passes ASAN/UBSAN, including all 30,720 packet
basis cases and the small-size alignment/alias/guard matrix. Its Release
build passes the same checks and independent full K20 checks. Rebuilt
default GF65536 grouped controls pass their small and K20 checks. All 31
targeted exact-proof and certificate-interface regression tests also pass.

Reproduce this checkpoint with the outer flags above plus
`HDGROUP_INNER=5`, `HD_SUMRANK_BYTE_KEYS=1`, and `HDGROUP_FUSE_REFRESH=1`.
Leave next-region prefetch enabled. Run the existing target with `--check`,
`--large-check`, and `--bench 501 full 20 nt 4 4`. No production default,
library profile, setup distribution, certificate, or paper number changes.
The original code remains selectable with `HDGROUP_INNER=0`.

The selected binary has SHA256
`6c0f42d94ef06b987dec159f092ea30ee286b0544f15e757fdf502fb43154264`.
The final clean rebuild and the measured screening binary have this same
hash. Source hashes, identical locally and on the benchmark host, are:

- `HdGroupedInner.h`: `5eb60b286dae2364df5010d971976fbcca06d7aa1f96b92d84d3754bd09a190d`.
- `HdSumrankInner.h`: `995effc7621bbea8eae0cd9ca8f8d721bee9747def7e98f1585db0c7a79a4798`.
- `hd_grouped_forward.cpp`: `216811f81b9a49895795ba425513b454343d7a755183f8ad171758961e16c5ac`.

The outer and inner now have comparable costs. The next bounded experiment
was to tune the outer's group tile size: the present 64-group tile reloads
128 KiB of inputs across its two RS cosets. Smaller tiles could improve
reuse, but may worsen the regional-store schedule. Preserve the fixed leaf
circuits, compare complete encoding, and retain the current candidate.

## Outer tile sweep: no material additional win

The tile experiment changes only the order in which independent outer
groups are evaluated. It preserves the RS leaf circuits, column mixers,
route tables, output layout, and inner. `HDGROUP_OUTER_TILE` selects a
positive multiple of four at compile time; its default remains 64. A tile
of T groups revisits 2048T bytes of input across the two RS cosets. The
driver prints the selected tile with its benchmark metadata.

The initial sweep tested 4, 8, 16, 32, 64, and 128 groups. Contrary to the
small-tile hypothesis, 128 looked strongest. A second sweep added 12,
256, 512, 1024, and 8192. The size-12 case also exercises a partial final
tile, rather than assuming that the tile divides the number of groups.
All eleven variants pass the small-size oracle/alias/alignment/guard
matrix, the 30,720 packet basis checks, and the independent g2/g4 K20
full-encoder oracle. Every g4 K20 output retains checksum
3966879751720036268. The 128-group ASAN/UBSAN build also passes the small
matrix and all packet basis checks, including the partial final tile at
K=24576. All 31 targeted proof-interface regressions pass.

Measurements use the same serial locks, CPU 2, 4.5 GHz, and shared-L3
caveat as the previous checkpoint. Each full run has 501 timed calls;
setup, allocation, and input production remain outside the timed region.
No own build, check, or second benchmark runs concurrently on Peach.

| Outer tile, groups | Full-encoder run medians, ms |
|---|---|
| Frozen 64-group checkpoint | 4.011, 4.025, 3.953, 4.135, 4.000, 4.004, 4.027 |
| 4 | 4.037 |
| 8 | 3.952, 4.034 |
| 12 | 4.024 |
| 16 | 4.113 |
| 32 | 4.084 |
| Rebuilt 64-group control | 4.048, 4.111, 4.098 |
| 128 | 3.875, 3.949, 3.861, 3.954, 3.953, 3.964, 4.118, 4.016 |
| 256 | 3.893, 4.002, 4.072 |
| 512 | 4.033 |
| 1024 | 4.083 |
| 8192 | 4.080 |

The median of the eight 128-group run medians is 3.954 ms, versus
4.011 ms for the seven frozen-control runs: about 1.4% less time. However,
the final reversed comparison gives 4.118 / 4.027 / 4.016 ms for
128 / frozen / 128, so the advantage is not consistent. The rebuilt
64-group control also varies from the frozen binary. Retain 128 as a
tested option, not as an established speedup or a new default. The
research timing remains about 4.0 ms. Neither a smaller tile nor a full
outer pass established a material improvement.

The final separate 128-group phases measure 1.845 ms outer and 1.991 ms
routed inner. A nearby frozen outer run measures 1.879 ms. These diagnostic
medians do not sum exactly to full encoding, and they do not establish a
cause for the small timing differences. The sampled code and the 35% /
72.3975-bit certificate are unchanged; no proof search, library promotion,
or paper update is involved.

Reproduce the alternative with the selected mode-5 flags above and
`-DHDGROUP_OUTER_TILE=128`, using the existing target and checks. The
measured 128-group binary has SHA256
`88751e0481800b0334932934ec0ef3a9a4929e0a997733705ce4a4984060c585`;
the rebuilt 64-group binary has SHA256
`1ab3a4b125f4f3fd5d52813d18698a7c30c4d5a8b2fddd519568f892beda9550`.
The frozen baseline retains its preceding hash. Uploaded and local source
hashes agree:

- `HdGroupedOuter.h`: `87f775f3127ece37f534c53db107e26ef82190406f48d06cd6292b5943780c14`.
- `hd_grouped_forward.cpp`: `fb4ad3469a73585461d47849ae75e9bb4260abaf29b8cf76ffac63cd94247adb`.

The next useful exact-map experiment is to fuse more of the RS transform
with the column mixer, aiming to remove intermediate scratch loads and
stores. Do not continue tile hill-climbing on differences of this size.
Preserve the frozen reference and validate any fusion against the same
independent oracle before comparing full encoding.
