# Striped outer and two scalar scans: rejected cost control

Date: 2026-10-06. This candidate is withdrawn for the 40%-distance goal.
The companion's grouped-check argument proves minimum distance at most
832321 at K=2^20, whereas 40% requires at least 838861. This bound holds
for every setup and every coefficient field. The implementation remains
a checked research control for scan and routing costs, not a new winner.
See [the proof investigation](HD_BAA_PROOF.md) and
[the cost-model review](HD_TOPOLOGY_COST.md).

## Implemented map

Each independent 64-symbol source row uses the normalized additive-basis
GRS[128,64] outer over the embedded AES GF256. Four rows share each SIMD
vector; the source indexing is the same four-row interleaving as the
earlier high-distance forward control. Output coordinate j is stored in
plane j, in increasing row order. There are 16384 rows and 128 planes.

Let z be the resulting column-major vector of N=2097152 field symbols.
The encoder computes u_i=u_(i-1)+d_i*z_i from zero state, then
y_i=y_(i-1)+e_i*u_(pi(i)) from zero state. Setup samples one uniform global
permutation pi and independent nonzero GF(2^16) scalars d_i,e_i. The
implementation uses reproducible setup fixtures and precomputed indices
and multiplication coefficients. Setup is excluded from encoding time.
Eight GF(2^16) components occupy each 128-bit record and share the same
scalar map. The distance objective counts nonzero symbols/records, not bits.

The experimental source is `spin/experiments/char2_field_spin/hd_striped_baa.cpp`.
It reuses the existing explicit 16-register FFT leaves and two-GFNI
GF(2^16) multiplier. The hot loops process fixed batches of sixteen
records, with four independent products before the XOR prefix carries.
There are no hot-loop allocations. All payload buffers belong to the caller.

## Serial forward timings

Host and policy match the RS and high-distance controls: Peach Ryzen7950X,
CPU2, userspace 4.5GHz, boost disabled, GCC15.2, Release, AVX512/GFNI,
`-mtune=znver4 -fno-lto`. The benchmark holds all five shared host locks.
Every process has three warmups. Input is produced with non-temporal stores
and a fence outside timing. Timing includes the outer, every required
movement, both diagonals, both scans, and cached output stores.
Allocation, setup, validation, and hashing are excluded. Large payloads
request huge pages; physical backing was not verified.

| Schedule | Timed calls | Median full forward time |
|---|---:|---:|
| Direct global gather, first launch |101 |12.152ms |
| Direct global gather, final launch |101 |11.987ms |
| Gather prefetch64 |51 |12.435ms |
| Gather prefetch256 |51 |12.202ms |
| Separate gather/product pass and final XOR scan |51 |13.726ms |
| Separate gather with prefetch128 |51 |14.099ms |
| Bucket route into256KiB destination groups |51 |12.075ms |
| Bucket route into1MiB destination groups |51 |13.926ms |

The 51-call controls include phase clocks; the 101-call timings do not.
These are process medians, not confidence intervals. None approaches5ms.
The small variation around12ms does not justify selecting the bucket route.

The direct-gather diagnostic measured 3.272ms outer, 1.191ms first scan,
and 7.746ms global gather plus second scan. Those medians need not sum
to the independently measured full-call median. The256KiB bucket control
measured 3.331ms outer, 3.048ms first scan plus scatter, and5.698ms local
gather plus second scan. It exchanges gather cost for scatter cost.

The bucket control preserves exactly the same sampled permutation. Setup
places each source symbol into the intermediate buffer of its destination
bucket, in source order. The first scan writes directly to these prepared
indices. The second scan gathers only within the current bucket. This
uses an extra32MiB workspace and two extra index maps; no second random
permutation is sampled. All staging is included in timing.

## Validation and receipts

Independent direct polynomial evaluation and scalar GF(2^16) arithmetic
agree with the complete forward map at K=256,512,1280. Tests cover dirty
buffers, four record offsets, unchanged input, exact input/output aliasing,
guard records, and timed/untimed paths. A K16 setup check verifies bucket
bijection, exact permutation equivalence, and destination locality across
multiple buckets. All K20 schedules produced checksum8673133842191526392.

The registered Release test passes. ASAN/UBSAN with leak checking passes
for the bucket variant, including its multi-bucket setup check. The hot
experimental code is instrumented; the shared auxiliary static library
was the Release build. No production implementation was changed.

The final direct-gather binary SHA256 is
`be781392af7c300ea811451920c985054ceaa14a98cd1e93f102c255d28a3aba`.
The complete experimental-source snapshot `tmp/hd-striped-source.tgz` has
SHA256 `89f935183d71410e8cde83935991fa3b490df4e26d70fe9c4fd5e1b6ca7eb3c0`.
Raw receipts remain ignored under `receipts/hd-fast/receipts/striped*.txt`.
The remote workspace remains `/tmp/spin-hd-fast-20261006-PTWML1`.

```text
cmake -S spin/experiments/char2_field_spin -B build/hd-fast -DCMAKE_BUILD_TYPE=Release
cmake --build build/hd-fast --target spin_char2_hd_striped_baa -j 3
ctest --test-dir build/hd-fast -R char2_hd_striped_baa --output-on-failure
build/hd-fast/spin_char2_hd_striped_baa --bench 101 full
```

Compile-time controls are `SB_PREFETCH`, `SB_SPLIT`, and `SB_BUCKET_LOG`;
all default to zero. Comparative timing requires CPU pinning and the shared
benchmark locks. Do not infer a distance claim from this implementation.

## Next construction gate

Two scalar scans cannot be repaired by optimizing coefficients or proof
constants. At least one recurrence must be stronger. A candidate for the
next round is GRS[64,32], two consecutive independently weighted prefix
scans, a regional route, and one final weighted scan. The inverse of two
consecutive scans has bandwidth two, rather than one. The grouped-check
degree becomes at most six, whose distance ceiling permits40% but not45%.
That necessary test is not a lower-distance proof. Regional routing is a
new distribution here and needs its own joint support/rank analysis.

The next step is to test this candidate's structural cuts before another
long optimization round. Retain the4K positive proof control and the frozen
high-rate winner. Neither the5ms target nor a new full certificate is met.
