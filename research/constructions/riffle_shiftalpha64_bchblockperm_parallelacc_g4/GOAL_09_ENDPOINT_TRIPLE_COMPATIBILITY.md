# Goal 09: endpoint compatibility at occupation three

## Question

Goal 08 was dominated by one light BCH block and two heavy BCH blocks. The
scalar bound placed one heavy block at weight 128, the unique all-one word.
It then placed the light block at weight 22 and the other heavy block at
weight 106.

Does the actual outer code permit the profile

\[
(22,106,128)
\tag{1}
\]

on a three-symbol support?

## Reduction to one field value

Every three-symbol support of the outer MDS code contains a one-dimensional
space of nonzero words. Normalize one coordinate to the field value whose BCH
encoding is all ones. The first outer equation forces the other two BCH words
to be complements. Their weights therefore sum to 128.

The second outer equation selects their field values. For three data
positions, the selected ratio depends only on the two exponent differences.
It does not depend on the first position. This reduces the complete scan from
cubic to quadratic work.

The scanner also covers supports containing the first parity position. A
support containing the second parity position has two equal finite
components. Such a support cannot realize (1).

## Exact result

The scan covers all (733{,}141{,}975{,}040) three-position supports. It
tests every normalization that can place the all-one word on a finite
support. The result is

\[
\boxed{A_{(22,106,128)}=0.}
\tag{2}
\]

The smallest candidate weight is 32. The complementary weight is 96.
Therefore the lightest endpoint profile found by the scan is

\[
(32,96,128).
\tag{3}
\]

The exact candidate histogram contains (2{,}199{,}023{,}247{,}360)
normalized endpoint words. The count equals three times the number of finite
supports, as required.

A direct implementation independently enumerates the reduced 12-data-block
instance. Its complete candidate-weight histogram agrees with the optimized
scanner.

## Effect on the proof bound

Goal 08 bounded the dominant one-light, two-heavy pattern by

\[
2^{159.773}.
\]

Replacing its independent band maxima by the exact endpoint histogram gives

\[
2^{40.744}
\tag{4}
\]

for finite supports containing the all-one BCH word. Thus compatibility
recovers about 119 bits on the dominant endpoint calculation.

Bound (4) remains trivial. It is not the complete occupation-three shell.
It excludes supports containing the second parity position and words without
the all-one BCH block.

## Remaining loss

The current inner coefficient envelope treats all 96 packets as variable.
The all-one BCH block instead contributes exactly 32 packets with value
`1111`. Local permutation does not change those packets.

An endpoint-aware inner operator should retain the number of fixed `1111`
packets already placed. It can continue to coefficient-bound the other 64
packets. This adds one coordinate of size 33 and removes a large false family
of packet values.

## Reproduction

Compile with the Visual Studio x64 environment and run:

```powershell
cl /O2 /std:c++20 /EHsc /arch:AVX2 scripts\scan_riffle_shiftalpha64_endpoint_triples.cpp /Fe:scripts\scan_riffle_shiftalpha64_endpoint_triples.exe
scripts\scan_riffle_shiftalpha64_endpoint_triples.exe 16384
python scripts\audit_riffle_shiftalpha64_endpoint_triples.py --data-blocks 12
python scripts\analyze_riffle_shiftalpha64_endpoint_triple_transfer.py --optimizer-maxiter 80
```

Artifact hashes are:

- scanner source: `2DCF83F494546E1A00236F9056F721418F1C7FBEF449B17B41F3B7083D62DAED`;
- scanner executable: `403C50B7AAEC534890659F089CAD5EAAE50B6EB0E3F77586CB9528E0B31C322E`;
- audit script: `B3FF2785BAD772CBE0936DC84D7F9B26C68FAB7D6229341A1E14B6CD1BC4B29B`;
- transfer script: `6E64BB60C391234729D87912B5FDD9EB57918BA84761CECAC8E5F2EC474D1EBB`.

## Next goal

Build the endpoint-aware 16-state operator with exactly 32 fixed `1111`
packets. Recompute (4). If the finite endpoint shell closes, add the much
smaller family of supports containing the second parity position. Then return
to the remaining light-heavy profiles that do not contain the all-one word.

