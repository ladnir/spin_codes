# Goal 13: compress the low-shell/complement cluster

## Result

This goal closes the complete finite-support profile `(22,22,26)` and every
odd-complement triple whose low representatives have weights in
`{22,24,26}`. The first closure uses exact shifted-schedule matching. The
second closure follows from one weight inequality and requires no schedule
scan.

After these closures, the rigorous open occupation-three envelope is

\[
2^{90.0358}.
\]

The finite-support remainder is at most \(2^{90.0279}\). The second-parity
remainder is unchanged at \(2^{82.5140}\). These are bounds on the remaining
occupation-three cases, not an end-to-end distance bound.

## Shared complement model

Let \(e\in\mathbb F_{2^{64}}\) be the message whose extended BCH encoding is
the all-one word. Write each selected BCH message as

\[
x_i+c_i e,
\qquad c_i\in\{0,1\},
\]

where \(x_i\) is a low-weight representative. The first outer parity check is
equivalent to

\[
x_1+x_2+x_3=(c_1+c_2+c_3)e.
\tag{1}
\]

Thus even-complement profiles use the untwisted relation
\(x_1+x_2+x_3=0\), while odd-complement profiles use the twisted relation
\(x_1+x_2+x_3=e\).

## Odd-complement lemma

Suppose that each low representative has BCH weight in `{22,24,26}` and that
the number of complements is odd. Equation (1) would make the XOR of the
three low BCH words equal to the weight-128 all-one word. Hamming weight is
subadditive under XOR, so this would require

\[
128
\leq
\operatorname{wt}(B(x_1))+
\operatorname{wt}(B(x_2))+
\operatorname{wt}(B(x_3))
\leq 26+26+26=78,
\]

which is impossible. Therefore every odd-complement triple over these three
low shells contributes zero. This single lemma subsumes the earlier exact
zero for `(106,106,106)` and, in particular, closes `(106,104,106)`.

The lemma concerns finite supports, where the first outer check forces the
BCH XOR relation. It is not applied to the distinct second-parity support
family.

## Exact `(22,22,26)` schedule closure

The affine relation enumerator finds exactly

\[
168{,}184{,}576
\]

ordered global `(22,22,26)` relations. They occupy `45,987,504` distinct
nonzero field ratios; the maximum global ratio multiplicity is `292`.

The exact shifted-schedule scan covers every finite support and every role of
the unique weight-26 word. It finds:

- `30,958,878` words on three data coordinates;
- `30,970,511` words on supports containing the first parity coordinate;
- `61,929,389` words in total.

The source-preserving inner transfer gives the pooled one-word bound

\[
K_D(22,22,26)\leq 2^{-30.6880351}.
\]

Consequently, the complete finite-support profile contributes at most

\[
61{,}929{,}389\,2^{-30.6880351}
<2^{-4.8039141}.
\tag{2}
\]

A direct 12-data-block implementation independently reproduces the optimized
scanner's `17,715` words, including its role-by-role split. A regression run
also reproduces every prior 12-block Goal 12 count for `(22,22,24)`.

## Other closures retained from the first pass

The profiles `(22,104,106)` and `(24,106,106)` reduce to the Goal 12
`(22,22,24)` relation family. Reusing that exact relation table and the
source-preserving inner bounds closes them at, respectively,

\[
2^{-25.9584}
\qquad\text{and}\qquad
2^{-34.1042}.
\]

No additional schedule scan is required for either profile.

## Remaining frontier

The envelope fell from `2^92.2729` after Goal 12, to `2^90.9275` after the
initial Goal 13 closures, and now to `2^90.0358`. Its leading finite profiles
are:

- `(24,24,24)`;
- the three orderings of `(22,24,24)`;
- the three orderings of `(24,104,104)`;
- the orderings of `(22,102,106)`.

All four displayed families have even complement parity. Their first outer
equation is therefore untwisted. The next genuinely new relation objects are
the low-shell families \(L_{24}^3\) and \(L_{22}\times L_{24}^2\), not the
twisted family that motivated this goal.

Before enumerating the complete weight-24 shell, the next goal should seek a
shared correlation, Fourier, or affine-orbit description of those two even
families. A successful compression would let one schedule calculation remove
several leading profiles at once.

## Reproduction

```powershell
python scripts/analyze_ebch128_goal13_relations.py --raw-ratio-2226 constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_weight22_22_26_ratios_raw.bin
scripts/scan_riffle_shiftalpha64_equalpair_outer.exe constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_weight22_22_26_ratios_raw.bin 16384 constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_weight22_22_26_ratio_counts.bin 168184576 22_22_26
python scripts/audit_riffle_shiftalpha64_weight222224_outer.py --histogram constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_weight22_22_26_ratio_counts.bin --optimized constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal13_weight222226_outer_scan_12.json
python scripts/analyze_riffle_shiftalpha64_goal13_schedule_closure.py
python scripts/analyze_riffle_shiftalpha64_occ3_no_endpoint.py --exclude-finite-222 --exclude-p1-222 --exclude-finite-adjacent --exclude-finite-goal13-initial --exclude-finite-goal13-schedule --optimizer-maxiter 100
```

Principal receipts are:

- `receipts/goal13_weight222226_outer_scan.json`;
- `receipts/goal13_weight222226_outer_scan_12.json`;
- `receipts/goal13_weight222226_outer_audit.json`;
- `receipts/goal12_weight222224_outer_scan_12_regression.json`;
- `receipts/goal13_schedule_closure.json`;
- `receipts/goal13_occ3_after_schedule.json`.

The raw relation stream and compressed multiplicity table are retained as
reproducible intermediate artifacts in `receipts/`.

SHA-256 hashes for the main certificates are:

- raw relation stream: `485888F81E2240A4721ABF39900E6CEB0E04A9F42E34E94D0617C2B5432D484E`;
- compressed ratio table: `A623BCD6C59B61BBE9C67117849EB16F0D4E8B3793D5B78F81CE396D6B440516`;
- full outer scan: `2DEA6679600670158DBFAB3D7A76C31EE36AC5DBFF87F772064E49F94EDDF870`;
- direct audit: `2A30B1B1CC7BF5580BB48C02DABF51DF94D522AD2F20006DBFD171F0D586860F`;
- combined schedule closure: `7C738B971CD0A7F7DF0C6A062768D99609DF27D59130B79F941428CF82C152DE`;
- post-exclusion envelope: `464A15D26F0F64E4847342D256323E6851EFCEC85226BE7FFD0F9DE7B35AC923`.
