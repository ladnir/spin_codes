# Fused BCH Loads and a Local Pair Mixer

2026-09-30. Neither prototype reaches the 6 ms target. Fusing independent
coordinate loads into BCH preparation gives no measurable speedup. The new
two-packet mixer takes about 7.315 ms with two IMT updates, even after making
its memory accesses local. No production code or existing proof input changed.

## Independent-Row Route: Same Code, Different BCH Evaluation

The existing independent-row GF16 construction has a complete 10%-distance
certificate with four IMT updates. This experiment changes only its evaluation.

The baseline first copies a 16 KiB four-row tile into local storage, restores
canonical BCH columns, and evaluates the GFNI BCH transpose. Mode 16 loads
the independently permuted dense coordinates directly into GFNI's preparation
network. It still copies the cold tile first, and prepares the systematic
coordinates separately. Mode 17 skips the cold-tile copy too.

The existing `gfni_bch.py` generator already emits the required
`bchTranspose4GfniIndependent` helper. The new modes select that helper; they
do not replace its fixed-width arithmetic or change the sampled code.

| Independent-row GF16, four updates | Median (ms) | Process-median range (ms) |
|---|---:|---:|
| Existing copy, repack, BCH (mode 7) | 7.490615 | 7.476093--7.516257 |
| Copy, fused BCH preparation (mode 16) | 7.547561 | 7.531537--7.565871 |
| Direct cold loads, fused preparation (mode 17) | 8.994067 | 8.969156--9.016975 |

Separate instrumented runs attribute roughly 3.7 ms to repacking plus BCH
for both the baseline and mode 16. Mode 17 raises that phase to roughly
5.2 ms. Inner plus routing remains about 3.8 ms. Removing the intermediate
dense stores and loads is insufficient; retaining a sequential cold-data
copy is important on this machine. These timings do not isolate individual
cache or instruction effects.

## Shared-Row Route with a Two-Packet Mixer

For each four-row BCH group, sample a uniform perfect matching of its 256
columns. For each matched pair of GF16 values `(x,y)`, sample independent
nonzero scalars `a,b,e` and apply

```text
u = a*x; v = b*y;
(x',y') = (u+v, e*(u+2*v)).
```

GF16 uses the polynomial basis modulo `X^4+X+1`. This replaces the original
per-packet multiplier; it does not add a second independent packet-only
multiplier layer. The remaining shared-column route and IMT maps are unchanged.
The benchmark evaluates the transpose of this construction.

The first prototype mixed the cold routed buffer in place. It took about
11.36 ms with two updates. The improved prototype copies each 16 KiB tile
sequentially into aligned local storage, mixes it there, and passes that
storage directly to BCH. There is no copy back. The fixed multiplication
by 2 uses two explicit SIMD shuffles and XORs.

| Construction | Median (ms) | Process-median range (ms) |
|---|---:|---:|
| Shared GF16, two updates | 5.822830 | 5.807276--5.834836 |
| Local pair mixer, two updates | 7.314681 | 7.307232--7.331767 |
| Shared GF16, four updates | 6.491356 | 6.482635--6.497663 |
| Local pair mixer, four updates | 8.041335 | 8.037002--8.092546 |

The local mixer plus BCH takes about 4.27 ms, versus about 2.67 ms for the
original BCH phase. The mixer removes the old GF16 work from the routing
phase, but its later mixing and memory work more than offset that saving.
A bounded one-update screen gave 6.992908 ms: dropping one update does not
bring this implementation close to 6 ms. That screen is one process with
31 calls, not a repeated headline measurement or a distance certificate.

The pair mixer is a new code distribution. Passing implementation checks
does not transfer the existing shared-route or independent-row certificates.
Its proof analysis is reported separately.

## Measurement and Correctness Scope

All measurements use Peach's Ryzen 7950X, GCC 15.2, core 15, `-O3`, and
`-mtune=znver4`. The workload is the precomputed in-place transposed encoder
at `K=2^20`, with 128-bit elements. Setup, allocations, input filling,
reference checks, and checksums are outside timed calls. Each process uses
three warmups followed by 101 timed calls, without resetting input.

Headline entries are medians of four process medians: two independent
seeds, with balanced execution order. Benchmarks run serially while holding
the three existing encoder/Hypercat benchmark locks. Instrumented phase
means use separate 31-call runs and include warmups.

Modes 7, 16, and 17 pass full encoder and suffix comparisons for updates 2
and 4, seeds 1 and 17, and three input patterns at `K=2^14`. Modes 16 and
17 also pass three-pattern full-size checks at `K=2^20`. The harness checks
the route, dense inner reference, inner adjoint, and GF16 linear map.

The pair mixer passes an independent scalar full-encoder reference, local
forward/transpose adjoint checks, and suffix preservation for updates 1,
2, and 4, seeds 1 and 17, and three input patterns at `K=2^14`. The local-copy
two-update case also passes the full-size three-pattern check. Both cold
and local-copy modes pass ASan/UBSan at the small sizes. The linked,
unchanged prebuilt BCH objects and library are not instrumented.

## Reproduction

Sources are `joint.cpp` (new modes 16/17), `pair_mixer_perf.cpp`, and
`fused_bch_run.sh`. The script reuses the existing generated BCH object and
SPIN library from its third argument. On Peach, this run used:

```text
root=/tmp/spin-fused-bch-oHiXdl
reference=/tmp/spin-joint-9n57TT
bash "$root/fused_bch_run.sh" "$root" build "$reference"
bash "$root/fused_bch_run.sh" "$root" check "$reference"
bash "$root/fused_bch_run.sh" "$root" confirm "$reference"
bash "$root/fused_bch_run.sh" "$root" profile "$reference"
bash "$root/fused_bch_run.sh" "$root" build-mixer "$reference"
bash "$root/fused_bch_run.sh" "$root" check-mixer "$reference"
bash "$root/fused_bch_run.sh" "$root" mixer "$reference"
bash "$root/fused_bch_run.sh" "$root" profile-mixer "$reference"
bash "$root/fused_bch_run.sh" "$root" mixer-r1-screen "$reference"
bash "$root/fused_bch_run.sh" "$root" build-mixer-sanitize "$reference"
bash "$root/fused_bch_run.sh" "$root" sanitize-mixer "$reference"
```

The pair mixer accepts `exponent updates seed calls [profile [local_copy=1]]`.
Use `local_copy=0` to reproduce the cold in-place control. Zero calls selects
reference checks only. The script exposes each phase separately; do not run
benchmark phases concurrently.

Raw logs are retained remotely under `measurements/` in that directory and
locally under ignored `tmp/fused-bch-results/`. Measured binary SHA256:

```text
joint: dc353e2e7e37d7aad51097513d2a6b7141d656aebd651017bf88bdd79ab3dc0e
pair mixer: 58e1c7c48aac2737286515b6cfa15b12422bb05dd7c843046c2034165924e5e8
pair mixer ASan/UBSan: 88d971f6f05850586e1e15195e5403c9cdbfa56f974d22fcfd2867b472ff40b9
```

Keep these prototypes as measured negative results. Neither is a compelling
implementation path to 6 ms / 10% without a more substantial change.
