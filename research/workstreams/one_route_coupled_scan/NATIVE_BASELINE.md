# First native gate for the full-field baseline

Status, 2026-10-21: the certified inner has a correct native implementation,
but it is not competitive with the preserved 3.798 ms GF256 reference.
The fastest tested variant is about 43--44 ms at K=2^20. It retains
rate one half, >35% F-symbol distance, and the <2^-40 ideal-setup failure
bound through the [spectral-map proof](SPECTRAL_BATCH.md).

This is a forward encoder of 128-bit F-elements, not a transposed-kernel
measurement. No production profile, paper, or preserved implementation
was changed. All sources are confined to this workstream's `native/`.

## What changed

The first implementation computed the certified Toeplitz product by
evaluation at all 256 GF256 points. Both transform interfaces used
monomial coefficients. This required two polynomial-basis conversions
per batch. An alternate Cantor basis made these conversions XOR-only,
but did not make them cheap enough.

The spectral map changes the inner distribution, with a separate exact
proof. It treats the state directly as novel-basis coefficients and
samples independent full-F point multipliers. This removes the basis
conversions. Restricting randomness to 192 points preserves exact joint
uniformity and reduces products and setup storage. Three offset-specific
64-point transforms in each direction compute the restricted map exactly.

Finally, accepting novel-basis message coefficients removes the outer's
monomial conversion. This is an invertible change of message coordinates
within the same fixed RS[64,32] code, not a different outer spectrum.

The full-field multiplier uses the existing sparse-byte representation
F=GF256[Y]/(Y^16+Y^5+Y^3+0x32). Subfield scaling and full-F products use
the same element representation, so no field-coordinate conversion is
hidden around the pointwise products.

## Measurements

Peach, AMD Ryzen 9 7950X, CPU 2, GCC 15.2.0, `-O3 -std=c++20
-march=znver4`. Each row uses three warmups and the median of 31 timed
encodes. Jobs ran serially. The host is shared, so these are diagnostic
measurements, not a stabilized publication benchmark.

All rows have K=1048576, N=2097152, L=128, b=64, and B=64. Timings include
outer encoding, the single stratified route, and the recursive inner.
Allocation, setup, input generation, and the final checksum are excluded.
The setup is reused for all timed messages. No input copy is charged.

| Inner implementation | Outer basis | Keys MiB | Full ms |
|---|---|---:|---:|
| Toeplitz, original transform basis | monomial | 128 | 201.358067 |
| Toeplitz, Cantor transform basis | monomial | 128 | 156.983492 |
| Spectral, 256 points | monomial | 128 | 55.370134 |
| Spectral, 192 points, three-coset kernel | monomial | 96 | 45.103860 |
| Spectral, 192 points, unpruned control | novel | 96 | 48.294656 |
| Spectral, 192 points, three-coset kernel | novel | 96 | 43.407484 |
| Same three-coset variant, final-source replay | novel | 96 | 43.667288 |
| Same variant, GFNI-only full-F multiplier | novel | 96 | 46.747679 |

The final two CLMUL runs and the unpruned control have checksum
5231751927570795027. Different setup distributions and outer coordinate
bases legitimately give different checksums in other rows. Setup of the
last profile takes about 52 ms with the reproducible benchmark generator.

Component timings from the first novel-outer/three-coset run:

| Phase | ms |
|---|---:|
| Outer alone | 4.363787 |
| Route alone | 6.461352 |
| Inner, contiguous input | 29.369097 |
| Inner including routed reads | 38.681170 |
| Full encoder | 43.407484 |

These are separate measurements, not additive accounting. Key and input
traffic interact in the full kernel. Input occupies 16 MiB, two N-symbol
buffers occupy 64 MiB, route labels occupy 8 MiB, and keys occupy 96 MiB.
No key expansion or coefficient dump is written to disk.

Bounded 4-MiB arithmetic streams isolate the main costs. Normalizing to
the K20 batch count, the original monomial transform pair costs
163.195520 ms and its 256 products cost 26.415072 ms. In the final source,
an unpruned novel transform pair costs 14.339712 ms and 192 products cost
19.838472 ms. The transform microbenchmark is not the pruned full encoder;
these figures must not be summed to predict a complete latency.

The retained 3.798 ms GF256 result was not rerun in this session. It is
a different coefficient family. No matched speedup against it is claimed.

## Correctness and proof status

All tested variants passed the following independent checks before timing:

- Every monomial basis vector for transforms of sizes 64 and 256 against
  direct evaluation, followed by inverse recovery.
- Four-stream evaluation against four independent transforms.
- For Toeplitz, all 128 unit-state middle products against raw diagonals,
  plus four dense products against scalar field multiplication.
- For spectral maps, all 128 state-basis evaluations and all 256 extraction
  columns against explicit root-product polynomials and Lagrange coefficients.
- Four dense spectral products against those independently constructed maps.
- CLMUL/GFNI feedback equality, regional bijections, and the complete tiny
  outer/route/recursive encoder against the direct polynomial oracle.
- Undefined-behavior sanitizer on the final three-coset profile.

The pruned and unpruned 192-point controls agree on their full K20 checksum.
Independent analytic reviews verified the map's rank, coefficient indexing,
three-coset pruning, and certificate transfer. These implementation tests
do not replace the all-message proof.

The benchmark uses fixed `mt19937_64` seeds for reproducibility. The proof
concerns independent ideal uniform coefficients and regional permutations;
it does not certify that deterministic pseudorandom seed as a particular
good setup or make a cryptographic-PRG claim.

## Reproduction

From the repository root on a supported x86 host, build outside the checkout:

```sh
g++ -O3 -std=c++20 -march=znver4 \
  -Ispin/experiments/char2_field_spin/field128_lift \
  -DSPIN_OR_TRANSFORM_CANTOR=1 -DSPIN_OR_SPECTRAL_POINTS=192 \
  -DSPIN_OR_SPLIT192=1 -DSPIN_OR_NOVEL_OUTER=1 \
  research/workstreams/one_route_coupled_scan/native/toeplitz_native.cpp \
  -o /tmp/spin-spectral192
timeout 60s taskset -c 2 /tmp/spin-spectral192 --check
timeout 90s taskset -c 2 /tmp/spin-spectral192 --bench 20 31 all
```

Use `SPIN_OR_SPLIT192=0` for the exact unpruned control; select 256 points
with `SPIN_OR_SPECTRAL_POINTS=256` and no split. Value zero selects the
original Toeplitz distribution. `SPIN_OR_TRANSFORM_CANTOR=0` preserves
the original transform arithmetic. The `gfni` phase changes only the
full-field multiplication kernel. The benchmark accepts only powers of
two from K12 through K20; it is not the library's generic-size API.

Builds and jobs were bounded by timeouts and a 2-GiB virtual-memory limit;
more than 20 GiB free was confirmed before running. Remote source and
binaries occupied about 2 MiB, including the final baseline regression.
No local build tree or raw timing archive exists.

Final source SHA256:

```text
toeplitz_native.cpp       7becad3cb1b36a1b1e45ca3e788976f52b9cca7b99aa31c414fb744d85365702
Transform256.h           dee50875c09b213187a38057658263a8b092644b889d37030463d1043f3799fd
Transform256Polynomial.h cf6f2502bc800b9f214e9f68dd0b554ed1c716f78882ca83a67745f8ea6bee83
```

The final native binary hash is
96874e7b8fdd831afa380f42ffd2637d2e7e93741e5610f61d60c2ea270aa928.
Compiler and flags are part of this receipt; binary identity is not
expected across compiler versions.

## What this says about the next step

The proof baseline is now operational, but its performance target is not
met. Basis conversions were an avoidable cost. The remaining full-field
products and scalar routing are substantial costs in their own right.

The exact-uniform F-bilinear interface needs at least L+b-1=191 variable
field products per batch; the 192-point map nearly meets that bound.
This is not a lower bound on machine instructions or all code designs.
It does show why more transform micro-tuning is not a sufficient plan.

Next compare two explicit changes: a materially cheaper realization of
those full-field products, or an inner with fewer such products whose
inverse distribution is controlled by a sharper, nonuniform bound.
Keep the current proof as the control. Do not silently substitute small
coefficients and inherit its uniform-output probabilities. Routing will
also need packetization or fusion before the few-millisecond target is
credible. No new faster construction is certified by this diagnosis.
