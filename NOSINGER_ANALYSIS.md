# Singer-removal investigation

Status: rigorous local envelope and exact small-dimension experiments complete;
a first-active matrix ledger now closes the identity construction over a large
part of the large-`K` weight range.  Its current middle-slice relaxation still
does not close the densest range, and the complete theorem-facing ledger has
not yet been regenerated.

## What the current proof actually uses

Fresh Singer maps do not make different inputs independent.  They give the
conditional one-point law

```text
fixed v != 0  =>  A_i v is uniform over F_2^64 \ {0}, conditional on the past.
```

This collapses every live transition to the ordinary BCH spectrum.  Without
Singer, fresh split permutations still make the state exchangeable.  Given
incoming state weight `q`, input-block weight `r`, and intersection size `J`,

```text
t = wt(V) = q + r - 2J
```

and `V`, conditional on `t`, is uniform on the weight-`t` slice.  The exact
replacement kernel is therefore the input/output enumerator

```text
B[t,w] = #{v : wt(v)=t and wt(E(v))=w}.
```

## Rigorous local result

`scripts/analyze_nosinger_bch_kernel.py` authenticates the committed cyclic
BCH generator, exactly enumerates the boundary slices through input weight
eight (and their complements), and uses the rigorous caps `B[t,w] <= A[w]`
elsewhere.  A robust Bellman operator lets every capped slice choose its worst
permitted BCH distribution independently at every step.  Thus it is more
pessimistic than any fixed encoder.

At pole `z=2333/2373`, an exact rational weighted-sup-norm certificate gives

```text
log2(rho) = -0.615274494049
prefactor log2 = 0.140751026335
Chernoff crossing at N=2^21: 7525 live blocks
```

The Singer certificate has `log2(M(z))=-0.781489918058` and crossing 5923.
The difference is real in the present envelope: rare low-state trajectories,
especially state weight one, dominate the no-Singer low-output MGF.

The existing parity-block outer has only 41.31 bits of ultra-late placement
margin at 5949 blocks.  Moving the blanket late suffix to 7526 blocks weakens
that term to about 25.17 bits, so Singer cannot simply be deleted from the
current 40-bit certificate.

## True-distance evidence

`scripts/compare_riffle_identity_distance.py` constructs matched Singer and
identity instances and enumerates every codeword by Gray code.

```text
K=20, N=256, 20 trials:
  Singer mean distance   89.80
  Identity mean distance 89.05
  Identity wins/ties/losses = 8/1/11

K=24, N=256, 5 trials:
  Singer mean distance   85.00
  Identity mean distance 85.20
  Identity wins/ties/losses = 2/1/2
```

There is no observed structural distance collapse.  Singer's demonstrated
benefit is suppression of the ensemble's rare low-output trajectories.

## Promising replacement

Use an independent uniform random graph subcode of codimension `r` outside
the explicit parity-block outer:

```text
m -> (m, Rm),  R <- F_2^{r x K} uniform.
```

For every fixed nonzero ambient message, graph membership has probability at
most `2^-r`.  This is exactly the first-moment factor needed from the outer;
no independence between different codewords is required.

For `K=2^20`, `r=24`, and the extra local block required to preserve the user
dimension:

```text
N = 2,097,408
no-Singer crossing = 7526 blocks (the placement calculation uses 7527 conservatively)
ambient ultra-late term <= 2^-25.1683
graph-subcode ultra-late term <= 2^-49.1683
```

The transposed graph map can use three eight-bit subset tables.  Its exact
overhead is

```text
3 * 247 table XORs + 3*K output XORs = 3,146,469 GF(2^128) XORs,
with 3 MiB of packed fixed masks.
```

With no inner mixer, the projected construction count would become 95,432,317
XORs instead of 129,506,976, a 26.31% reduction.  The identity-only high-weight
proof does not currently close with the scalar episode envelope.  The PAP9
candidate below adds 18,581,724 XORs, for a projected total of 114,014,041:
11.96% below the current implementation.  Both candidates remove the 16.8 MiB
Singer-row setup table.

### First-active identity matrix ledger

The scalar episode decomposition is not intrinsic.  It can be replaced by a
65-state backward operator indexed by state weight.  For an occupancy pole
`x`, output pole `z`, and a terminal vector `v`, let `G_z(v)[t]` process a BCH
input of weight `t` and the random output/state split.  The exact boundary
input slices are used directly; each unknown slice uses the ordinary BCH
spectrum caps adversarially.  Let `H_x[q,t]` count input blocks, weighted by
`x^r`, whose XOR with a state of weight `q` has weight `t`.  One suffix block
is then the nonlinear positive operator

```text
F_x,z(v) = H_x G_z(v).
```

Expose the first active input block.  If its weight is `r`, the suffix has
`L` blocks and contains the remaining `H=h-r` input ones, Cauchy's bound and
the uniform interleaver give

```text
C(64,r) / C(N,h) * x^-H * e_r^T G_z F_x,z^L(1).
```

Summing this expression over `r` and every possible `L` automatically keeps
turnoff, exact cancellation, restart, and terminal-state discard.  There is
no episode cap or restart prefactor.  The current long-double diagnostic in
`scripts/probe_identity_matrix_ledger.py` gives:

```text
h range                      output z       occupancy x    ambient interval log2
2002--7858                   .9800511634       .008             -283.2755
7860--20550                  .9414884538       .025            -1769.4101
20552--30000                 .9414884538       .025            -4880.7851
30002--75000                 .8178969571       .060            -3270.3872
75002--110000                .8178969571       .100             -463.4333
110002--200000               .5382199631       .200           -47106.8614
200002--350000               .5382199631       .300           -16059.9394
```

These are complete sums over every supported (even) weight in each displayed
interval.  The `r=24` graph factor subtracts another 24 bits from their union.

There is also an exact global-consistency point that is special to the
identity construction.  At `x=1`, summing over every input block `U` makes
`S xor U` uniform on `F_2^64` for every fixed incoming state `S`.  Hence
`H_1 G_z` is computed from the ordinary BCH spectrum exactly; it must not be
formed by maximizing every input-weight slice independently.  At
`z=0.1002588437`, this exact operator gives approximately

```text
h=370000       -7813.1189
h=750000     -133110.6234
h=900000     -133110.6796
h=1048576    -133110.6649
h=1200000    -133115.8704
h=1500000    -127847.1084
```

The first-active and exact-uniform lanes overlap near `h=364000` with
thousands of bits of margin.

For `h>N/2`, write each input block as all-ones plus a defect block and run the
same 65-state operator on defect weight `g=N-h`; the XOR matrix rows are merely
reversed.  Exploratory pole rows overlap across the full complement-high
range.  Representative ambient bounds are

```text
h=1550000, (z,x)=(.3331239681,.5)       -41964.5173
h=1800000, (z,x)=(.3331239681,.5)      -100420.3947
h=1850000, (z,x)=(.5382199631,.2)      -155220.2428
h=2050000, (z,x)=(.5382199631,.2)      -112223.1136
h=2090000, (z,x)=(.5382199631,.2)        -5095.1143
h=2090000, (z,x)=(.9800511634,.005)      -6714.6466
h=N,       (z,x)=(.9800511634,.005)      -3355.3294
```

The complete even-weight interval sums for the final schedule are

```text
h range                 mode / (z,x)                         ambient interval log2
350002--364000          first-active / (.5382199631,.3)           -6106.3281
364002--1500000         exact-uniform / (.1002588437,1)            -1938.7850
1500002--1900000        complement / (.3331239681,.5)              -3173.8158
1900002--2090000        complement / (.5382199631,.2)              -5095.1109
2090002--2097408        complement / (.9800511634,.005)            -3355.3294
```

Together with the earlier table, this scans every supported outer weight from
`2002` through `N`.  The earlier positive dense result was entirely a
proof-relaxation artifact: taking the best of occupancy poles
`.7,.8,.9,.95,1` per suffix still left `+36879.4` at `h=2^20`, whereas
enforcing the exact `x=1` global consistency gives `-133110.7`.

An exact enumeration of BCH input slices `t=9` and `t=55`, followed by a
retuning of the proof poles, closes the remaining identity-only bridge:

```text
h range       output z   occupancy x   ambient interval log2
502--1000       .997        .001             -110.1072
1002--2000      .995        .002             -465.8677
```

The earlier `-58.0942` value for this row was from the PAP9 fallback and is
not used.  The identity matrix schedule now covers every supported nonzero
outer weight.  Its combined diagnostic is dominated by the refined low row:
approximately `-22.6085` ambient and `-46.6085` after `r=24`.

### `r=24` low-weight ledger probe

The finite cancellation calculation must retain the exact input slices near
both zero and 64.  With the exact rows through input weight eight and their
complements, `T=7527`, and at most 499 remaining input ones, exact rational
arithmetic gives

```text
worst finite termination log2 = -36.212938969731
valid uniform dyadic cap       = 2^-36
worst input weight             = 1
```

Each new identity episode also incurs the geometric-envelope restart factor
`prefactor/rho`, whose logarithm is 0.75603 bits.  The probe therefore uses the
conservative combined restart cap `2^-35`, not merely the one-step `2^-36`
termination cap.  With explicit episode terms through `e=32` and 250-block
refinement of the first two 4000-block placement buckets, it gives

```text
ambient low-weight total <= approximately 2^-22.6085
after random graph r=24  <= approximately 2^-46.6085
diagnostic margin beyond 40 bits = 6.6085 bits
```

This calculation deliberately upper-bounds each placement subbucket by its
largest binomial coefficient and each survival subbucket by its left-endpoint
inner bound.  It is therefore pessimistic, but it still uses floating
logarithms and is not yet the theorem-facing rational certificate.  The
All `h>=502` rows are covered by the matrix-valued identity ledger above.
The remaining tasks are arithmetic hardening and verifier packaging.

### Exact power-of-two puncturing

The `K=2^20`, `r=24` graph construction initially has length `2^21+256`:
`16385` data EBCH blocks and one componentwise parity EBCH block.  Puncture
the same one coordinate from each of 256 distinct data blocks before the
global interleaver.  The resulting outer length and inner recursion length are
both exactly

```text
N = 2^21.
```

The committed extended BCH generator is invariant under cyclic scaling of
the 127 nonzero `GF(128)` coordinates and under translation `x -> x+1`.
`scripts/punctured_ebch_outer.py` verifies both permutations directly by
row-space reduction.  Hence the code is coordinate transitive and the exact
one-coordinate-punctured spectrum is

```text
B[j] = ((128-j) A[j] + (j+1) A[j+1]) / 128.
```

It has parameters `[127,64,21]` and cardinality `2^64`.  The heterogeneous
ambient outer used by the diagnostic is therefore the exact product of 256
punctured spectra and 16130 ordinary EBCH spectra.  The parity-block floor
drops by at most one, from 44 to 43.

Regenerating every identity matrix interval at `N=2^21`, distance 188743,
gives

```text
ambient low-weight total       approximately 2^-22.537900097230
after random graph r=24        approximately 2^-46.537900097230
margin beyond 40 bits          approximately 6.537900097230 bits
```

The full diagnostic manifest is
`scripts/power2_punctured_identity_ledger_manifest.json`; its structural
verifier reports complete integer-weight coverage through `2^21`.  As with the
unpunctured identity ledger, the remaining gap to a theorem is outward/exact
arithmetic for the floating rows, not a missing weight interval.

### PAP9 global-consistency envelope (proof fallback, not the target encoder)

Nine PAP rounds admit a sharper proof than applying the BCH spectrum caps to
every input slice independently.  Let `P` be the exact 64-state accumulator
weight chain.  The boundary input slices `1..8` and `56..64` are known exactly.
On the remaining slices `9..55`, the nine-round density satisfies

```text
max_{t,s unknown} P^9[t,s] / (C(64,s)/(2^64-1)) <= approximately 2^0.003071858.
```

Combining this density bound with the exact boundary slices and the ordinary
BCH column spectrum gives the following EBCH-schedule live-MGF bounds:

```text
pole              PAP9 log2 M upper
0.980051163400       -0.922514752621
0.941488453800       -2.738487985849
0.817896957100       -8.811620857236
0.538219963100      -24.236171572132
0.320407579189      -38.333994400294
```

The same global-consistency calculation includes both state-zero turnoff and
exact support cancellation and gives `2^-62.41487575`, preserving the old
dyadic `2^-62` episode cap.  The decisive high-weight diagnostics are

```text
h=501..2000 critical total, ambient       approximately 2^-58.09
refined near prefix h=2001..350000        approximately 2^-1008.04
middle prefix endpoint                    approximately 2^-240.79
later endpoint family                     approximately 2^-442.08
ultra-late h>2000                         approximately 2^-1131.97
```

The remaining paired/complement families had tens of thousands of bits of
margin in the Singer certificate; PAP9 changes each live-block MGF exponent
by about 0.0031 bits at the scheduled poles, so these are not the diagnostic
bottleneck.  The low-weight row remains dominant at about `2^-46.61` after the
`r=24` graph factor.  These PAP9 calculations still use long-double/log-domain
arithmetic and must be converted to exact or outward arithmetic before the
paper can claim the construction.

## Reproduction

```powershell
python scripts/analyze_nosinger_bch_kernel.py --exact-radius 5 --steps 24
python scripts/probe_nosinger_r24_h500.py
python scripts/punctured_ebch_outer.py
python scripts/probe_nosinger_r24_h500.py --power2-punctured
python scripts/verify_identity_matrix_ledger.py --manifest scripts/power2_punctured_identity_ledger_manifest.json
python scripts/analyze_pap_weight_mixing.py --rounds 9
python scripts/compare_riffle_identity_distance.py --k 20 --trials 20
python scripts/analyze_nosinger_outer_tradeoff.py --codim 24
python scripts/probe_identity_matrix_ledger.py --h-min 2002 --h-max 7858 --output-pole 0.9800511634 --occupancy-pole 0.008
python scripts/probe_identity_matrix_ledger.py --h-min 7860 --h-max 20550 --output-pole 0.9414884538 --occupancy-pole 0.025
python scripts/probe_identity_matrix_ledger.py --h-values 370000,750000,1048576,1500000 --output-pole 0.1002588437228 --occupancy-pole 1
python scripts/probe_identity_matrix_ledger.py --mode complement --h-values 1850000,2050000,2090000 --output-pole 0.5382199631 --occupancy-pole 0.2
python scripts/probe_identity_matrix_ledger.py --mode complement --h-values 2090000,2097408 --output-pole 0.9800511634 --occupancy-pole 0.005
```

`scripts/enumerate_bch_input_slice.cpp` reproduces the committed exact slice
rows in `scripts/ebch128_cyclic_input_slices.csv`.

## Lightweight inner mixer probe

A permutation alone has no effect on the proof kernel: it preserves input
weight, and the existing split already makes the next state exchangeable.

The kernel script also evaluates fresh permutation--accumulate--permutation
(PAP) rounds.  One round costs only 63 field XORs in the transposed block
circuit: a reverse prefix XOR with both permutations fused into its reads and
writes.  Its exact accumulator weight enumerator composes with the rigorous
BCH slice envelope.  At the same pole:

```text
                         log2 contraction   Chernoff crossing   turnoff upper
identity                    -0.61527449             7525          2^-35.97
one PAP round               -0.64205618             7210          2^-39.81
two PAP rounds              -0.64706617             7154          2^-43.04
PAP9 (robust slice caps)       -0.64728552             7153          2^-58.47
Singer                         -0.78148992             5923          2^-62.41
```

Thus PAP is a real improvement and extremely cheap, but one or two rounds do
not recover the existing 40-bit outer-placement margin by themselves.  PAP9
is the first round count that closes the scalar high-weight diagnostic through
the global-consistency argument above; PAP8 retains about 0.060 bits of MGF
gap per live block at the first dense pole and does not close this envelope.
