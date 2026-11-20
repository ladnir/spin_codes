# OT-generation micro-kernels

Isolated kernels for the silent-OT sender at libOTe's SPIN configuration,
K = 2^18 OTs: 128 GGM trees x 4096 leaves (N = 2^19 noise blocks), the SPIN
transpose 2K -> K, and the correlation-robust hash (2 AES per OT). They show
what each stage can reach on its own, and whether SPIN can hide AES work,
before any of it moves into libOTe.

- `kernels.h`: VAES GGM expansion and hash. Formulas match libOTe
  (`RegularPprf.h`, `SilentOtExtSender.cpp`). Variants change batching and
  traversal only; the output stays identical, except the leaf variant marked
  UNPROVEN.
- `bench.cpp`: validates every variant against scalar references (plus the
  FIPS-197 vector), then times each one against its AES-unit floor. It also
  times the SPIN library transpose.
- `fused.cpp`, `otgen_hook.h`, `patch_spin.cmake`: builds patched copies of
  SPIN's AVX-512 route-epoch and outer-group loops. Each copy has a hook that
  runs chunks of real GGM or hash work, compiled either inlined or out of line.
  The SPIN library sources are only read. The copies live in the build tree,
  and the build fails if a patch anchor no longer matches exactly once.
- `run_peach.sh`: builds both binaries and runs them one at a time, pinned to
  one core, under the shared benchmark locks, and only on an idle host.

## Results (Peach, Ryzen 9 7950X, core 15, 4.48 GHz, 2026-11-20)

AES-unit floor: 2 VAES ymm per cycle, i.e. 4 block-rounds per cycle, which is
0.558 ns per AES-128 block. A pure register loop reaches it.

| Stage | libOTe today | Best kernel | Floor | Kernel vs floor |
|---|---:|---:|---:|---:|
| GGM expansion, libOTe PRG | 0.774 ms | 0.521 ms (W=8, breadth-first) | 0.454 ms | 87% |
| GGM, AES+round leaves (UNPROVEN) | — | 0.393 ms | 0.322 ms | 82% |
| SPIN transpose | 0.457 ms | 0.49 ms (library) | — | — |
| Hash | 0.415 ms | 0.299 ms (direct, 8 OTs per iteration) | 0.293 ms | 98% |

The libOTe column is from the earlier in-pipeline profile; the kernels are
measured in isolation with warm caches.

- **GGM:** the libOTe structure (8 trees, breadth-first) is already right.
  Wider batches (16 or 32 trees) and L1-blocked traversal do not help: the
  stage is AES-bound, not memory-bound. Blocked traversal at W >= 16 is slower.
  The 0.25 ms gap to libOTe is implementation overhead around the same AES
  work.
- **Hash:** the gain comes from building each (r, r^Delta) pair with
  `vbroadcasti128` plus xor/and, and storing the ymm straight to the output.
  libOTe instead packs through a stack batch and copies it out. B=16 is slower
  than B=8; streaming stores make no difference here.
- **Leaf PRG (parked):** using the internal-level derivation at the leaves
  saves 0.13 ms. Unlike internal nodes, leaves are output: an adversary may
  see both children of the punctured leaf's parent, one masked by Delta, and
  the children are algebraically tied (`child0 = R(child1 - x) ^ x`). Whether
  that lets Delta be recovered is open. Kept only as a speed bound.

### Overlap of SPIN with AES work (one thread)

`hidden = T_spin + T_work - T_fused`, where all three run in the same binary.

| Work injected into SPIN's loops | hidden |
|---|---|
| GGM, real buffers (8 MiB output) | -0.05 to -0.53 ms (fusion slower) |
| Hash, real buffers (4 MiB in, 8 MiB out) | -0.04 to -0.49 ms (fusion slower) |
| Same hash work on an L1-resident 512-OT buffer | -0.002 to +0.014 ms (0-3% of SPIN) |

The out-of-line hook build gives the same numbers, so the penalty is not
register pressure in SPIN's loop. It goes away when the injected work's memory
footprint goes away, so it is cache interference. With memory out of the
picture there is still no meaningful overlap: SPIN keeps the same two vector
units (GFNI, byte shuffles) as busy as AES does. A port microbenchmark
confirms it: vaesenc, gf2p8affineqb and vpshufb mixes add up, they do not
overlap. Peach runs one thread per core, so SMT co-scheduling was not
available to test.

### What this means for one core

With no overlap, the single-core limit is the sum of the stage floors:
0.454 + 0.293 + SPIN (0.49) = 1.24 ms. The best kernels together give
0.521 + 0.49 + 0.299 = 1.31 ms, about 95% of that limit; today's libOTe
(1.648 ms) is about 75%. Beyond that, only fewer AES calls (the leaf PRG),
a faster SPIN, or more cores help.
