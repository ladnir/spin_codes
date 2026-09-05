# Goal 14: compress the even-complement weight-24 core

## Result

Goal 13 left four leading finite profiles. Three belonged to two shared
additive families, and the fourth was the pure weight-24 family. Goal 14
counts the relevant additive relations exactly and closes every selected
profile except `(24,24,24)`.

The main results are:

- `(22,24,24)` has zero shifted-schedule outer words;
- every selected two-complement profile closes below `2^-21.60` in total;
- the open occupation-three envelope falls from `2^90.0358` to `2^87.8525`;
- `(24,24,24)` is now the unique leading profile in the relation core.

These statements concern the finite occupation-three family. They do not
complete the end-to-end distance proof.

## Additive relation count

For an even number of complements, the low representatives satisfy

\[
x+y+z=0.
\tag{1}
\]

Let \(L_w\) be the set of messages whose extended BCH encoding has weight
\(w\). The complete shells have sizes

\[
|L_{22}|=243{,}840,
\qquad
|L_{24}|=6{,}855{,}968.
\]

The affine coordinate group preserves each shell and equation (1). It has
size `16,256`. The shell \(L_{22}\) has 15 free orbits. The shell \(L_{24}\)
has 385 orbits of size `16,256` and 147 orbits of size `4,064`.

Fix one representative \(x_O\) from each orbit \(O\). For target shells
\(L_a,L_b\), define

\[
p_O(a,b)
:=
\bigl|\{y\in L_a:x_O+y\in L_b\}\bigr|.
\]

Affine invariance makes this count constant on \(O\). Therefore the exact
ordered relation count is

\[
T(w,a,b)
=
\sum_{O\subseteq L_w}|O|p_O(a,b).
\tag{2}
\]

A vectorized scan of the 532 representatives gives

\[
T(22,24,24)=321{,}121{,}024
\tag{3}
\]

and

\[
T(24,24,24)=4{,}018{,}336{,}896.
\tag{4}
\]

Every relation in (4) contains three distinct words. Hence (4) represents

\[
669{,}722{,}816
\tag{5}
\]

unordered triples.

The computation has three independent checks. It reproduces the certified
`1,365,504` ordered \(L_{22}^3\) relations. It obtains (3) from both the
\(L_{22}\) and \(L_{24}\) orbit projections. Every \(L_{24}\) partner count
is even, and (4) is divisible by six.

## Two-complement closures

Complementing two messages preserves equation (1). The same relation counts
therefore control several actual BCH-weight profiles. The source-preserving
inner operator closes all selected high-weight copies without a schedule
scan.

| Relation family | Actual profile | Complete log2 contribution |
|---|---:|---:|
| \(L_{24}^3\) | `(24,104,104)` | `-27.2989` |
| \(L_{22}L_{24}^2\) | `(22,104,104)` | `-29.5505` |
| \(L_{22}L_{24}^2\) | `(24,104,106)` | `-23.0548` |
| \(L_{22}^2L_{26}\) | `(22,102,106)` | `-22.3250` |
| \(L_{22}^2L_{26}\) | `(26,106,106)` | `-32.8686` |

The sum of these complete contributions is below `2^-21.6093`. Each row
includes all ordered placements of the displayed weight multiset.

## Exact mixed-shell schedule zero

For a relation with one word in \(L_{22}\) and two words in \(L_{24}\), use
the ratio of the two \(L_{24}\) messages. The exact relation stream contains
`321,121,024` ratios. It compresses to `87,806,530` distinct ratios, with
maximum multiplicity `43`.

The shifted-schedule scanner tests every finite support and every placement
of the weight-22 word. None of the schedule ratios occurs in the relation
table. Both support types give zero:

\[
A^{\mathrm{finite}}_{(22,24,24)}=0.
\tag{6}
\]

A direct 12-data-block implementation independently reproduces zero for all
three placements. The generalized scanner also reproduces every prior Goal
12 count.

The product-stream generator handles the 147 short affine orbits through a
section of the group action. This section maps each orbit word exactly once.
The exact relation mass in (3) therefore contains no stabilizer duplicates.

## Updated frontier

Closing the two-complement profiles lowers the envelope from `2^90.0358` to
`2^88.3987`. Applying (6) lowers it again to

\[
2^{87.8525}.
\]

The finite-support remainder is at most `2^87.8164`. The second-parity
remainder is unchanged at `2^82.5140`.

The leading finite profile is now `(24,24,24)`. The next profiles are the
orderings of `(22,22,28)` and `(24,24,26)`, but each trails the leading
profile by about four projection bits.

## Quotient for the remaining profile

Expanding all relations in (4) would produce a 32 GB ordered-ratio stream.
That expansion is unnecessary.

Fix an unordered triple \(U=\{x,y,z\}\subseteq L_{24}\) satisfying (1), and
put \(r=y/x\). Its six labeled ratios form the multiset

\[
\mathcal R(r)
=
\left\{
r,
r^{-1},
1+r,
(1+r)^{-1},
\frac{r}{1+r},
\frac{1+r}{r}
\right\}.
\tag{7}
\]

Let \(S(q)\) be the number of finite supports whose fixed labeled
coefficient ratio is \(q\). The exact number of `(24,24,24)` outer words is

\[
\sum_{U}\ \sum_{q\in\mathcal R(r_U)} S(q),
\tag{8}
\]

where the inner sum respects multiplicity. Equation (8) permits two
simultaneous quotients:

1. enumerate the `669,722,816` unordered triples from (5);
2. canonicalize each six-element ratio orbit in (7) before schedule lookup.

The next implementation should build the canonical schedule table first.
It should then stream affine representatives of unordered \(L_{24}^3\)
triples directly into that table. This approach avoids the 32 GB raw stream
and reduces the schedule table by approximately the same sixfold symmetry.

## Reproduction

```powershell
python scripts/analyze_ebch128_weight24_relation_core.py
python scripts/analyze_riffle_shiftalpha64_goal14_inner.py --optimizer-maxiter 600
python scripts/analyze_riffle_shiftalpha64_goal14_relation_closure.py
python scripts/write_ebch128_goal14_mixed_ratio_stream.py constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_mixed_ratio_products.bin
scripts/convert_riffle_ratio_products.exe constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_mixed_ratio_products.bin constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_weight24_24_22_ratios_raw.bin 321121024
scripts/scan_riffle_shiftalpha64_equalpair_outer.exe constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_weight24_24_22_ratios_raw.bin 16384 constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal14_weight24_24_22_ratio_counts.bin 321121024 24_24_22
python scripts/analyze_riffle_shiftalpha64_occ3_no_endpoint.py --exclude-finite-222 --exclude-p1-222 --exclude-finite-adjacent --exclude-finite-goal13-initial --exclude-finite-goal13-schedule --exclude-finite-goal14-complements --exclude-finite-goal14-mixed-schedule --optimizer-maxiter 100
```

Principal certificate hashes are:

- relation core: `67F961A67406FF90A6CDB71621BFBF329643EE781177B545C81C22F3BC3BE519`;
- inner bounds: `ADBC2487B57E9856698A6AFDE2A94D16299B9F9995EBAF8CC97B7CDA93CE9D66`;
- combined closure: `A6E4920E0C55E9B4E52C0C29FAE610BD1731F8B446D2ED577F315C77606F69A7`;
- product stream: `D6B57D8EE0B132AE44B5BC008A4186E03041B068BFDBD7A4601B47E11CED3E89`;
- raw ratio stream: `31AEAA3A346FE70F2DD20EA94AB64135D649025B83DF82DB819F0C7274557CA5`;
- compressed ratio table: `3EC10760EA9518EF77C551E651D45549F51B560B0A6FB7E63383B329AB9B6A5B`;
- full mixed-shell scan: `E85441EF07B7114266E61A6DD974258900FE396902A93F35640BB84986D53233`;
- direct audit: `0E8F28B359C10B378EEFF08EBEBC5ACD5614EE39022F7135EF86658AF080B541`;
- post-exclusion envelope: `EDD6A72B1A86C792B918F2D1B8644E2C00B0121F1D26ABE17CF8D7E196AA2F88`.
