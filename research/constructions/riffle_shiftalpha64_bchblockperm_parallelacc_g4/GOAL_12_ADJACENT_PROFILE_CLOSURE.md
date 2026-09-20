# Goal 12: close the adjacent finite profiles

## Question

Goal 11 removed the finite profile ((22,22,22)). The remaining finite
envelope was dominated by

\[
(22,106,106)
\quad\text{and}\quad
(22,22,24).
\]

This goal counts both BCH relation families exactly. It then closes each
family at output threshold

\[
D=188{,}766=\lfloor0.09L\rfloor.
\]

## One low-shell enumeration

Let (B) be the extended binary BCH encoder, let

\[
L_{22}=\{x:\operatorname{wt}(B(x))=22\},
\]

and let (e) be the message for which (B(e)) is the all-one word. The
complete minimum shell has 243,840 messages and 15 free affine-coordinate
orbits of size 16,256.

The profile ((22,22,24)) is equivalent to

\[
x,y\in L_{22},
\qquad
\operatorname{wt}(B(x+y))=24.
\tag{1}
\]

Testing one first-word representative from every affine orbit against the
complete shell gives exactly

\[
27{,}765{,}248
\tag{2}
\]

ordered relations, or 13,882,624 after forgetting the order of (x,y).

The other profile uses complement symmetry. Since (B(e+x)) is the binary
complement of (B(x)), every relation has the form

\[
(x,e+y,e+z),
\qquad
x,y,z\in L_{22},
\qquad
x+y+z=0.
\tag{3}
\]

Thus ((22,106,106)) reuses the minimum-shell additive triples from Goal 11.
There are exactly

\[
1{,}365{,}504
\tag{4}
\]

ordered relations.

## A source-preserving inner transfer

Pooling the weights 22, 106, and 106 into one 384-bit slice gives a bound
larger than one. That loss is artificial: the pooled calculation forgets the
three source-block weights and then pays to condition on them again.

The replacement retains source membership through a positive-coefficient
envelope. Introduce one source variable (s_i) and one bit-weight variable
(x_i) for each BCH block. For a packet value (v\in\mathbb F_2^4), replace
the old transition coefficient by

\[
\sum_{i=1}^3 s_i x_i^{\operatorname{wt}(v)}.
\tag{5}
\]

After 96 packet transitions, extract the coefficient of

\[
\prod_{i=1}^3 s_i^{32}x_i^{w_i}.
\]

For positive (s_i,x_i), coefficient extraction is upper-bounded by
evaluation divided by this monomial. The exact normalization is

\[
\binom{96}{32,32,32}\prod_{i=1}^3\binom{128}{w_i}.
\tag{6}
\]

Equations (5) and (6) preserve the independent local block permutations and
the uniform global packet interleaving. Equal-weight source variables may be
evaluated at the same tilt; they remain distinct formal variables in the
coefficient argument.

The resulting rigorous one-word bound is

\[
K_D(22,106,106)\le 2^{-85.4188}.
\tag{7}
\]

The generic pooled calculation remains stronger for the nearly balanced
profile:

\[
K_D(22,22,24)\le 2^{-29.6521}.
\tag{8}
\]

## Closing ((22,106,106)) without an outer scan

Fix the coordinate roles and a relation ratio. Once two finite outer
coordinates are chosen, the second outer check determines the third
coordinate uniquely. With 16,385 finite coordinates, a fixed ratio therefore
occurs on at most

\[
16{,}385\cdot16{,}384
\tag{9}
\]

labeled supports. Combining (4), (7), and (9) gives

\[
\boxed{
\mathbb E Z_D^{\mathrm{finite},(22,106,106)}
\le 2^{20.3811+28.0001-85.4188}
<2^{-37.0377}.
}
\tag{10}
\]

No multiplier-ratio table is needed for this profile.

## Exact shifted-schedule scan for ((22,22,24))

For (1), record the ratio (r=y/x) between the two light coordinates. The
27,765,248 ordered relations occupy 7,375,088 distinct ratios; the largest
global ratio multiplicity is 133.

For every finite support, the scanner considers each of the three possible
locations of the weight-24 word. If its coefficient is (c), and the two
ordered light coefficients are (a,b), the second outer equation requires

\[
r=\frac{a+c}{b+c}.
\tag{11}
\]

The scan covers all 732,873,539,584 data-only supports and all 134,209,536
supports containing the first parity coordinate. It finds exactly

\[
9{,}331{,}089
\tag{12}
\]

outer words with profile ((22,22,24)). A direct 12-data-block enumeration
matches the optimized gap-normalized scan for all three weight-24 roles.

Combining (8) and (12) gives

\[
\boxed{
\mathbb E Z_D^{\mathrm{finite},(22,22,24)}
\le 2^{\log_2 9{,}331{,}089-29.6521}
<2^{-6.4985}.
}
\tag{13}
\]

The two adjacent finite profiles together contribute less than

\[
2^{-6.4984}.
\tag{14}
\]

## What the next envelope says

Removing both profiles lowers the open occupation-three envelope from

\[
2^{94.9201}
\quad\text{to}\quad
2^{92.2729}.
\tag{15}
\]

Finite supports still dominate. The leading projection terms now form a
cluster rather than one exceptional profile:

- ((106,106,106));
- ((24,24,24));
- the permutations of ((22,22,26)) and ((22,24,24));
- the complemented mixed neighbor ((22,104,106)).

This is a compression signal. Continuing with one exact ratio scan per
profile would reproduce the earlier case growth. The next step should group
the low shells and their complements, retain source-specific weight tilts in
the inner transfer, and sum their additive relations through one shared
outer-ratio interface.

## Reproduction

```powershell
python scripts/analyze_ebch128_adjacent_triples.py --counts-only
python scripts/analyze_ebch128_adjacent_triples.py --raw-ratio-222224 constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_weight22_22_24_ratios_raw.bin
python scripts/analyze_riffle_shiftalpha64_adjacent_inner.py --optimizer-maxiter 600
rustc -O -C target-cpu=native scripts/scan_riffle_shiftalpha64_weight222224_outer.rs -o scripts/scan_riffle_shiftalpha64_weight222224_outer.exe
scripts/scan_riffle_shiftalpha64_weight222224_outer.exe constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_weight22_22_24_ratios_raw.bin 16384 constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_weight22_22_24_ratio_counts.bin
scripts/scan_riffle_shiftalpha64_weight222224_outer.exe constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_weight22_22_24_ratios_raw.bin 12
python scripts/audit_riffle_shiftalpha64_weight222224_outer.py --histogram constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_weight22_22_24_ratio_counts.bin --optimized constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal12_weight222224_outer_scan_12.json
python scripts/analyze_riffle_shiftalpha64_adjacent_closure.py
python scripts/analyze_riffle_shiftalpha64_occ3_no_endpoint.py --exclude-finite-222 --exclude-p1-222 --exclude-finite-adjacent --optimizer-maxiter 100
```

Primary receipts are:

- `receipts/goal12_adjacent_triple_counts.json`;
- `receipts/goal12_adjacent_triples.json`;
- `receipts/goal12_adjacent_inner.json`;
- `receipts/goal12_adjacent_closure.json`;
- `receipts/goal12_weight222224_outer_scan.json`;
- `receipts/goal12_weight222224_outer_scan_12.json`;
- `receipts/goal12_weight222224_outer_audit.json`;
- `receipts/goal12_occ3_after_adjacent.json`.

The raw and compressed ratio artifacts are:

- `receipts/goal12_weight22_22_24_ratios_raw.bin`;
- `receipts/goal12_weight22_22_24_ratio_counts.bin`.

SHA-256 hashes for the main certificates are:

- adjacent relation receipt: `447A24C10811976212D4BF78D0CE2851354C6A702CB36B9DF4CC82B6E0B391F5`;
- source-preserving inner receipt: `69D27EF1EB9C3AB7E26617CC84314D4CF32A5D0DF19C4B6AB4207D0773BF715F`;
- combined closure receipt: `F9A5460DF7C188560437539A90D3C2C76723A9F26A4147B38D4E2FA6A29D6694`;
- raw ratio stream: `805E4EDFEC8B17F2ACC7E5F78A927E90119903747E9236120F4659C540D7A067`;
- compressed ratio table: `151335D22FC9298B0183B203671BCC8B9602014D0ECC124BB42FB7192FF4DB8F`;
- full outer scan: `9E60A1CB180778E1ACC8606DCED68384E2AADD27ECCB6150D0F746AD286DF6A9`;
- direct audit: `7136ACED9C063993A02F8A6B644987D01EFDE47A974A63D2F2E0ACFC8C72BC88`;
- post-exclusion envelope: `6CDFDDCB9828A18E51E1307C6F857984124F4FBFB84DDD2FE11A3F648A4F2A30`.

## Next goal

Build one low-shell/complement occupation-three envelope. Begin with the
shell set \(\{22,24,26\}\) and its complements \(\{106,104,102\}\). Count
additive relations by affine orbits where available, but expose them through
one ratio-table interface. Pair that outer table with the source-preserving
inner coefficient operator. The target is to close the entire clustered
family at once or identify a genuinely new dominant family.
