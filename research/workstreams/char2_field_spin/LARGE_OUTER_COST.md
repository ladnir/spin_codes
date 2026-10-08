# Larger outers and the width-32 state-eight packet

The structured GRS[256,176] outer has an upper circuit with 503 fixed GF256 products per four rows. Its straightforward coefficient-sharing schedule also has 6,440 Boolean operations, so product volume alone does not establish competitiveness. Pairing it with the new width-32 symmetric state-eight packet gives a mean estimate of `1471/704 = 2.08949` GFNI512 processed-byte equivalents per source record. This is lower than the four-state, width-16 estimate for the same outer, but neither circuit has a native schedule or timing result.

These are candidate arithmetic counts. They do not establish a whole-code distance or failure bound. The [state-capacity note](STATE_CAPACITY_REPAIR.md) retains the earlier width-16 construction; the width change is a separate construction.

## Outer circuit and checks

All evaluation points are AES GF256 byte labels. For length `b=4n`, first use GRS[4n,3n], then reconstruct the final `m=n/4` source coordinates from the first `k=3n-m`. The parity transpose of the first map has three symmetric XOR-circulant blocks of size `n`; the reconstruction transpose has eleven blocks of size `m`.

For each block, the superset zeta transform converts XOR convolution to disjoint-subset multiplication. For each transformed input separately, take a binary basis of its required coefficients together with one. A basis of dimension `d` needs `d-1` nonidentity GF256 products; identical coefficient products are shared across output blocks. The script chooses the basis in ascending byte order. The Boolean count reconstructs coefficients and sums terms using two- or three-input XOR operations. It includes the zeta transforms and the final systematic combinations.

| Outer | Base products | Shortening products | Products / four rows | Boolean upper count / four rows |
|---|---:|---:|---:|---:|
| GRS[64,44] | 82 | 19 | 101 | 647 |
| GRS[128,88] | 186 | 46 | 232 | 2,092 |
| GRS[256,176] | 400 | 103 | 503 | 6,440 |
| GRS[256,171], via [256,176] | — | 35 additional | 538 | 7,860 |
| GRS[256,170], via [256,176] | — | 42 additional | 545 | 8,210 |

Further shortening reconstructs five or six missing source values by ordinary interpolation. Each missing transpose input requires seven fixed products after binary coefficient sharing. GRS[256,171] has rate `171/256 = 0.66796875`, above two thirds. GRS[256,170] has rate `0.6640625`, below a hard two-thirds floor.

The small control reproduces the earlier 34-product GRS[32,24] circuit. Generic shortening gives eight products, rather than the separately optimized seven-product GRS[32,22] circuit. The script does not extrapolate that extra optimization to larger blocks.

For 128-bit records, four outer rows fill a ZMM register. Thus the outer volume is `P/(4k)`. The width-16 four-state packet adds `b/k`; the current width-16 symmetric eight-state packet adds `(91/64)(b/k)`.

| Outer | Full volume, state 4 / width 16 | Full mean volume, state 8 / width 16 |
|---|---:|---:|
| GRS[128,88] | `93/44 = 2.11364` | `30/11 = 2.72727` |
| GRS[256,176] | `1527/704 = 2.16903` | `1959/704 = 2.78267` |
| GRS[256,171] | `781/342 = 2.28363` | `997/342 = 2.91520` |
| GRS[256,170] | `1569/680 = 2.30735` | `2001/680 = 2.94265` |

## Width 32 with the same eight residue nodes

Keep `W=(0,1,92,93,224,225,188,189)` and residue nodes `16+W`. Combine the old phases `(0,1)` and `(2,3)`, respectively: the two packet-label sets are `H union (H+2)` and `(H+6) union (H+4)`, where `H` is the embedded GF16. Each packet is an affine 32-point space and the union of four `W` cosets. Keep the symmetric fixed maps with coefficient `(u+v)^7+(u+v)^3+(u+v)+1`.

The emission has eight inputs and four eight-coordinate output blocks. Sharing coefficients across all four blocks gives 36 and 37 scalar products in the two phases. The injection is its transpose, with the same multiplication count. Their paired processed-byte volumes are therefore 18 and 18.5 GFNI512 equivalents. Add two products for the eight eta-scaled state records, two for the fresh state diagonal, and eight for the 32 input diagonals.

| Quantity | First phase | Second phase | Equal-phase mean |
|---|---:|---:|---:|
| Products per fixed scalar map | 36 | 37 | 36.5 |
| Packet GFNI512 volume | 30 | 30.5 | 30.25 |
| Full volume with GRS[256,176] | `133/64 = 2.07813` | `1479/704 = 2.10085` | `1471/704 = 2.08949` |

The emission coefficient-sharing circuits use 180 and 167 scalar Boolean operations, counting ternary XOR as one. Expanding every ternary XOR conservatively to two binary XORs, and transposing the resulting circuit, bounds the two-map binary-XOR work by 744 and 692 operations. The latter bound includes the extra 24 additions caused by transposing an 8-input, 32-output circuit. Packet feedthrough and state accumulation are additional. This is a deliberately loose Boolean bound, not a proposed instruction schedule.

The product count benefits from sharing the eight transformed state inputs across four output blocks, and from halving state-refresh work per encoded record. Realizing that count requires packing distinct scalar products and scheduling the transposed fan-out. It must not be confused with issuing a fractional number of instructions in one packet. Outer scratch traffic, coefficient loads, transposes, and register pressure remain uncounted. The 256-point outer's Boolean count per four rows is about three times the 128-point count, or about 1.54 times per source record.

## What increasing the outer length saves in the proof

Let the outer rate be `R=k/b`, let `N=bL`, and fix encoded input density `alpha=W/N`. Suppose the active rows have mean support `rho b`. Then `A/L=alpha/rho`, and the leading entropy for selecting active row labels is

\[
 \log_2 {L\choose A}\simeq \frac{N}{b}H_2(\alpha/\rho).
\]

Doubling `b` halves this term at fixed `alpha` and `rho`. It does not halve the within-row support entropy: `A log2 C(b,rho b)` has leading term `(alpha N/rho) H2(rho)`. The shortened source dimension also has a leading term independent of `b`:

\[
 W-(b-k)A\simeq\alpha N\left(1-\frac{1-R}{\rho}\right).
\]

Thus larger geometry reduces row-label entropy and changes the regional capacity constraints; it does not remove all source entropy or repair a local return bound by itself. Comparing a fixed active-row count across different `b` is not a fixed-density comparison.

## Replay and next gate

From this directory run:

```text
python -B theory/large_outer_cost_reference.py --check
```

The check compares every base and shortening circulant entry with independent interpolation, checks each circulant circuit on all scalar coordinate basis inputs, and reproduces the existing smaller product counts. It also checks the eight-point width-32 fixed blocks on their coordinate bases. The width-32 algebra and distance screens are separate work; these arithmetic checks do not certify them.

Next, use the whole sparse and dense proof screens to choose an outer/packet geometry before implementing a register-bounded circuit. Width 32 merits that proof gate: its mean product-volume estimate is about 25% below width-16 state eight, without changing the outer rate. A complete-call benchmark would still be needed after a proof-positive choice.
