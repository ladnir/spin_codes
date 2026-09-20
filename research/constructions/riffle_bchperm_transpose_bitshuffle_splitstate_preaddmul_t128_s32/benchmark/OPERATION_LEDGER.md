# RM2Sub s=16 selected-kernel operation ledger

The selected inner processes 2^21 transposed 128-bit blocks in 16,384
epochs of 128 blocks. There are 16,382 ordinary interior epochs. The first
epoch has zero incoming state, and the final state after the last epoch is
discarded.

## Source-level linear operations per interior epoch

| Component | Previous | Selected | Reason |
|---|---:|---:|---|
| `A` table construction | 44 | 44 | Four 4-bit tables, 11 XORs each |
| `A` form and output XORs | 512 | 455 | Exact optimal partition `[0,2,3,4]`, `[1,10,13,15]`, `[5,9,11,12]`, `[6,7,8,14]`; zero groups are omitted at compile time |
| `B` RM transpose | 448 | 346 | Backward slice to the constant, seven linear, and 21 quadratic correlations |
| `B` quadratic projection | 75 | 42 | Fifteen shared gates plus 27 output XORs |
| GF(2^16) state multiplication | 92 | 92 | Four table builds and three XORs per output bit |
| Syndrome injection | 16 | 16 | One XOR per state bit |
| **Total** | **1,187** | **995** | **192 fewer, a 16.2% algebraic reduction** |

The selected complete call performs exactly

`388 + 16,382 * 995 + 499 = 16,300,977`

128-bit block-XOR equivalents. The 388 term is the first-epoch `B`; the 499
term is the final-epoch `A`. Multiplication by alpha is correctly omitted in
the first epoch because its incoming state is zero.

## External word traffic

The selected fused loop reads every 128-bit word block exactly once. It
writes every block except the unchanged first epoch. Relative to the prior
separate `A` then `B` implementation, it removes 33,550,336 bytes of word
reloads, essentially one complete 32 MiB pass.

## Emitted Zen 4 logic instructions

GCC 15.2 with `-O3 -march=znver4 -mavx2` uses AVX-512VL ternary logic and
the extended XMM register file. One ordinary interior iteration, including
the out-of-line pruned `B`, contains:

- 274 `VPTERNLOG` instructions;
- 226 binary `VPXOR` instructions;
- 500 SIMD logic instructions total;
- 423 XMM-width, 16 YMM-width, and 61 ZMM-width logic instructions.

These are static hot-loop counts from the saved binary. A ternary instruction
replaces two binary XORs. The ZMM operations are confined to the broad RM
stages; the source retains an AVX2-only fallback.
