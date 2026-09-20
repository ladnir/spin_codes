# Goal 11: the minimum triple after endpoint removal

## Question

Goal 10 closed every occupation-three word containing an all-one BCH block.
What profile dominates after removing that family?

The output threshold remains

\[
D=188{,}766=\lfloor0.09L\rfloor.
\]

## The first coarse answer is false

Removing BCH weight 128 from the three-band calculation reduces the bound
from \(2^{160.172}\) to \(2^{97.464}\). The new dominant band pattern has two
light blocks and one heavy block.

The maximizing weights are 22, 22, and 106. That assignment is impossible on
a finite support. The first outer check gives

\[
x_1+x_2+x_3=0.
\]

The BCH encoder is linear, so its three binary words also XOR to zero. Their
weights satisfy all Hamming triangle inequalities. In particular, a block of
weight at least 78 cannot equal the XOR of two weight-22 blocks.

## Exact weight envelope

The next calculation sums ordered BCH-weight triples instead of weight bands.
For each finite support, it retains only triples that satisfy the Hamming
triangle inequalities.

Fix an ordered profile \((w_1,w_2,w_3)\). The nonzero outer words on one
support form a field line. Projection onto each coordinate is injective.
Therefore the number of words with this profile is at most

\[
\min\{A_{w_1},A_{w_2},A_{w_3}\},
\tag{1}
\]

where \(A_w\) is the exact BCH spectrum coefficient. The inner calculation
continues to pool the three locally permuted blocks.

For a support containing the second parity position, the two finite field
values are equal. Its weight profile has the form \((h,h,k)\). The analogous
profile count is at most \(\min\{A_h,A_k\}\).

The exact weight envelope gives

\[
2^{96.221}.
\tag{2}
\]

Equation (2) remains trivial. Both support families are now dominated by the
minimum profile \((22,22,22)\).

## Exact additive triples in the minimum shell

Let \(L_{22}\subset\mathbb F_{2^{64}}^*\) contain the messages whose BCH
encodings have weight 22. A finite-support minimum triple requires

\[
x,y,x+y\in L_{22}.
\tag{3}
\]

The complete shell has 243,840 messages. Its affine coordinate symmetry has
15 free orbits, each of size 16,256. Testing one representative from each
orbit gives exactly

\[
1{,}365{,}504
\]

ordered additive triples, or 227,584 unordered triples. Their field ratios
\(y/x\) occupy 352,044 values. No ratio occurs more than 36 times.

The outer support fixes the ratio \(y/x\). An exact scan compares all 352,044
ratios with every finite support in the shifted outer schedule. None match.
Hence

\[
\boxed{
A^{\mathrm{finite}}_{(22,22,22)}=0.
}
\tag{4}
\]

The scan covers 732,873,539,584 data-only supports and 134,209,536 supports
containing the first parity position. A direct 12-data-block implementation
matches the optimized scan.

## Effect on the finite-support bound

Using (4) reduces the non-endpoint occupation-three bound to

\[
2^{94.920}.
\tag{5}
\]

The next finite-support profiles are \((22,106,106)\) and
\((22,22,24)\).

## Second-parity minimum profile

The second-parity problem is different from (3). It requires the pair
correlation

\[
\sum_{a<b}
\left|L_{22}\cap(a+b)^{-1}L_{22}\right|,
\tag{6}
\]

where \(a,b\) range over the 16,385 finite outer coefficients. The earlier
one-data scan proves that the 16,384 pairs containing the first parity
coefficient contribute zero.

An exact difference scan gives a reusable bound for the remaining pairs. The
16,384 data coefficients have 134,209,536 pairwise differences, and all
those differences are distinct. Thus the data coefficients form a Sidon set
in the additive group of \(\mathbb F_{2^{64}}\).

For \(d\in\mathbb F_{2^{64}}^*\), define

\[
I(d):=\left|L_{22}\cap d^{-1}L_{22}\right|.
\]

Every data-data support uses a different value of \(d\). Therefore

\[
\sum_{a<b} I(a+b)
\le \sum_{d\ne0}I(d)
=|L_{22}|^2
=59{,}457{,}945{,}600.
\tag{7}
\]

The equality counts every ordered pair in \(L_{22}^2\) by its unique ratio.
It does not assume that the shifted differences behave randomly.

The full-state positive-weight calculation gives

\[
K_D(22,22,22)\le 2^{-36.4472}.
\tag{8}
\]

Combining (7) and (8) yields

\[
\boxed{
\mathbb E Z_D^{p_1,(22,22,22)}
\le 2^{35.7912-36.4472}
=2^{-0.6561}.
}
\tag{9}
\]

Equation (9) closes the second-parity minimum profile. The margin is small,
but every input is exact or a rigorous upper bound.

For scale only, independent uniform multipliers would make the expected value
of the remaining sum about 0.433. This number is heuristic and is not used in
the proof.

After removing the minimum profile from both support families, the open
occupation-three remainder has bound

\[
2^{94.9201}.
\tag{10}
\]

Finite supports dominate (10). The second-parity remainder is
\(2^{82.5140}\).

## Consequence for the proof strategy

The weight-triangle step is scalable across all occupation-three profiles.
The affine-orbit ratio scan is scalable across the complete minimum shell.
The Sidon argument is more scalable: it converts every data-data support into
one use of a global pair-correlation budget.

Continuing with one scan for every adjacent profile would again create case
growth. The next calculation should treat \((22,106,106)\) and
\((22,22,24)\) together. Complementing the two weight-106 words converts the
first profile into an additive triple of low-weight words. This suggests one
joint low-shell enumerator rather than two unrelated cases.

## Reproduction

```powershell
python scripts/analyze_riffle_shiftalpha64_occ3_no_endpoint.py --exclude-finite-222 --exclude-p1-222 --optimizer-maxiter 100
python scripts/analyze_ebch128_weight22_triples.py --ratio-output constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal11_weight22_ratio_counts.bin
rustc -O -C target-cpu=native scripts/scan_riffle_shiftalpha64_weight22_outer_triples.rs -o scripts/scan_riffle_shiftalpha64_weight22_outer_triples.exe
python scripts/audit_riffle_shiftalpha64_weight22_outer_triples.py --data-blocks 12
scripts/scan_riffle_shiftalpha64_weight22_outer_triples.exe constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal11_weight22_ratio_counts.bin 16384
rustc -O -C target-cpu=native scripts/certify_shiftalpha64_data_difference_multiplicity.rs -o scripts/certify_shiftalpha64_data_difference_multiplicity.exe
python scripts/audit_shiftalpha64_data_difference_multiplicity.py --data-blocks 12
scripts/certify_shiftalpha64_data_difference_multiplicity.exe 16384
python scripts/analyze_riffle_shiftalpha64_lowblock_kernel.py --part-count 3 --part-weight 22 --samples 1 --positive-scaled-cost 33.7
```

Artifact hashes are:

- occupation-three analysis: `A32BDC0B2BBE03D10537177B4447E92911653DCC8D9561014E1F11340E1A61EE`;
- additive-triple analysis: `9EAAB8B0C8CC043F73052B875FAE33B59F3A1B7EF3F431D3E6927A1994AC9393`;
- outer scanner source: `BB1835D97410F183EE0D85522602E36AEAC2B3FFA050293B662BAC1E9692D739`;
- outer scanner executable: `7D17D081FCCBAA10E23FF0EF7E58B8F84F2C0D4BD461C49FF2BB245F0A303E1B`;
- direct audit: `F5E52903944F032ED7E4AEA7351450D4057093F26F6BF19F13021DD3B6C0AC28`;
- ratio table: `639CD3453122E719FFBBD4D47248D6BAAB86963232C3193882C80F107BC42716`;
- difference scanner source: `93CABE885AB1C4B66C22C8E19D9BD967CE765140482B279D849FE010C1665FDA`;
- difference scanner executable: `C837B0EBFFBFE37116A623D10AA94EFD74DA632672C6CAFBCB9AAAA7A34ED275`;
- difference audit: `DDC04A699655D050615BF16EA8CF20F04381E2EEF1C7F2755968456E845B0BEF`.

## Next goal

Build one affine-orbit enumerator for the adjacent low-shell relations behind
\((22,106,106)\) and \((22,22,24)\). First count their global additive
relations. Only construct outer-ratio tables if the global counts, combined
with the Sidon budget and the corresponding exact inner kernels, do not
already close the profiles.
