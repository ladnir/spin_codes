# Rectangular bundles and shuffled layouts

Subsequent result: WIDE_STREAMING.md rules out the g=4,c=8 candidate below
using a late-activation subspace. Full zero-state rank did not detect that
obstruction. The timing remains valid, but that geometry cannot meet 10%.

2026-09-23. The 2x objective remains open. This iteration expands the distribution
search and rejects two families with explicit low-weight words. Production code
and the paper are unchanged.

## Rectangular distribution and obstruction

Use groups of g outer rows and a bundle size c dividing 256. There are 256/c
macroregions. Each contains one bundle of g*c coordinates from each row group.
Independently permute the groups in each macroregion. Each group retains its
sampled coordinate permutation, shared or independent across its rows as in
the first iteration. Within each of the c columns, independently permute the
g row labels. The mathematical model uses independent uniform draws.

The route is bijective and each row contributes c coordinates per macroregion.
For c>1 the original one-row marginal no longer applies. These experiments check
the new occupancy invariant rather than the library's original region invariant,
and compare the actual linear circuit to a materialized reference.

Suppose g*c divides t=128 and g*c<=128. A group's bundle lies in one epoch per
macroregion. Messages in that group have dimension 128*g. Requiring B*x_epoch=0
in every nonfinal epoch imposes at most (256/c)*19 linear constraints. If g*c>38,
there is a nonzero kernel message. Its state stays zero and its output weight
is at most 256*g, regardless of the transvections.

The experiment extracts a kernel message by back-substitution, independently
forms its BCH word, and runs the library's packed forward encoder. Every output
bit must equal the routed outer word. At K=2^20 and route seed 17:

| g | c | Message dimension | Rank | Nullity | Verified output weight |
|---:|---:|---:|---:|---:|---:|
| 1 | 64 | 128 | 75 | 53 | 102 |
| 4 | 16 | 512 | 303 | 209 | 422 |

Output length is 2,097,152. These are counterexamples, not loose bounds. The
dimension argument excludes these geometries for every route in these families;
the listed weights belong to the particular checked messages. Designs with
larger state or additional cross-bundle communication are not thereby excluded.

## Implementation findings

Host, workload, timing exclusions, and serial benchmark locks match RESULTS.md.
The rectangular screen uses one 31-call process per configuration, three warmups,
and seed 17. Shared g=4, c=4 or 8 gives about 8.46 ms with the existing tiled
kernel, versus a 9.64 ms original-distribution control. This roughly 12% reduction
is a screen, not a confirmed ranking of nearby parameter choices.

Keeping each group in shuffled column order permits contiguous bundle writes.
But making BCH loads look up that order is costly: an isolated BCH diagnostic
gave 6.17 ms versus 3.73 ms for the original circuit. Repacking one 16 KiB group
before the original circuit improves this to about 5.1 ms. Isolated timings are
not additive estimates of full-encoder phases.

We also replaced per-coordinate route indices by one base and packed lane choices
per bundle. The fixed-width unrolled IMT remains; there are no encode-time
allocations. Neither compact schedules nor cache-line streaming stores yielded
a winning full encoder.

### ISA correction

The first harness compiled its new inner kernels for AVX2, while production
Fast.cpp also enables AVX-512F/VL. The harness now matches these flags. Earlier
custom-kernel losses are not optimized performance ceilings. The existing tiled
library path was unchanged. Corrected single-process 51-call medians:

| Implementation | Time |
|---|---:|
| Original distribution, tiled control | 9.610 ms |
| Shared g=4, c=1, compact route and local BCH repack | 9.249 ms |
| Shared g=4, c=2, same kernel | 14.718 ms |
| Shared g=4, c=4, same kernel | 10.422 ms |
| Shared g=4, c=8, same kernel | 9.411 ms |

An identity-route cost diagnostic takes about 7.53 ms versus 9.65 ms for the
control, using medians of three 101-call processes. It is not a valid-distance
construction or a lower bound on all implementations. It shows that this kernel
organization does not approach 2x merely by making its route cache-local.
Fusion or a different execution organization could behave differently.

## Shared global-column distribution

A separate candidate shares one sampled coordinate permutation across all outer
rows. Region permutations of four-row groups and lanes remain independent. Its
one-row marginal is the original marginal; multi-row correlations are different.
No multi-row certificate transfers.

It supports a column-major intermediate matrix, where each inner region writes
only a 128 KiB column at K=2^20. A blocked transpose assembles 32-row tiles for
the unchanged BCH circuit. Correctness passes, but the implementation takes
about 12.94 ms versus a matched 9.72 ms control. A phase diagnostic gives about
7.86 ms for inner/routing and 5.16 ms for transpose/BCH. These perturbed means
are diagnostic, not replacement timings. The implementation loses; the
distribution has not been mathematically rejected.

## Validation, reproduction, and next step

Every timing process compares its candidate with the same-map reference first.
Small-size checks also use the dense transpose oracle, adjoint identity, and
in-place suffix checks. The nullspace witnesses are checked by the independent
packed forward path. A full sanitizer and length campaign remains future work.

The executable accepts an optional final `columns` argument. Serial scripts are
`run_rectangles.sh`, `run_mapped.sh`, `run_bundled.sh`, `run_ceiling.sh`, and
`run_column.sh`. Raw CSVs remain under `/tmp/spin-locality-JjY7yR` on Peach;
earlier AVX2 receipts are separately retained in `measurements-avx2` there.

Final executable SHA256:
`fe5e42d0acf8415d242b1de5ef7fa551d2c889ec9668410301f804a886acc926`.
The linked library hash remains the one recorded in RESULTS.md. All four
existing library tests pass after this iteration.

Next, compare production and experimental inner/routing instruction and memory
costs directly before another broad sweep. Keep the small rectangular families
for proof work and reject the dimension-obstructed geometries. The target needs
a substantially better execution organization, not another small improvement
to the best current tiled route.
