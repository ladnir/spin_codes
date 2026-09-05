# Goal 15: close the finite `(24,24,24)` profile

## Result

Goal 15 closes the finite-support `(24,24,24)` profile exactly. The shifted
outer code contains

\[
60{,}462{,}011
\]

words with this BCH-weight profile. Combining this count with the rigorous
one-word inner bound from Goal 14 gives

\[
60{,}462{,}011\cdot 2^{-32.4042261}
<2^{-6.5547004}.
\]

Thus this profile contributes at most `0.010638` to the bad-distance union
bound. It is no longer an open case.

This result advances the finite occupation-three analysis. It does not close
the end-to-end distance proof.

## Symmetric quotient

Let (F=\mathrm{GF}(2^{64})). Consider a nonzero triple

\[
x+y+z=0.
\tag{1}
\]

Two such triples are equivalent if one is obtained from the other by a common
nonzero scalar and a permutation of the three entries. Define

\[
J(x,y,z)
=
\frac{(x^2+xy+y^2)^3}{(xyz)^2}.
\tag{2}
\]

The expression in (2) is defined because all three entries are nonzero. It is
unchanged by a common scalar. Equation (1) also makes it symmetric in
(x,y,z).

To see that (J) is a complete invariant, set (r=y/x). Permuting the three
entries generates the six transformations

\[
r,quad r^{-1},quad 1+r,quad (1+r)^{-1},quad
\frac{r}{1+r},quad \frac{1+r}{r}.
\tag{3}
\]

In these coordinates,

\[
J(r)=\frac{(r^2+r+1)^3}{r^2(r+1)^2}.
\tag{4}
\]

The generators (r\mapsto r^{-1}) and (r\mapsto 1+r) preserve (4). For a
fixed value of (J), equation (4) has degree at most six. A generic orbit in
(3) has six elements, so it is the entire fiber. The only smaller orbit
satisfies (r^2+r+1=0); its two elements form one orbit and both have (J=0).
Therefore two nonzero triples satisfying (1) have equal (J) if and only if
they are equivalent.

This invariant replaces a six-ratio table with one unordered projective key.
For (J\ne0), a matching relation triple gives one assignment to a fixed
support. For (J=0), the projective triple has three automorphisms and the
count receives a factor of three. The weight-24 relation set contains no
(J=0) triple, so the exceptional factor does not affect the final count.

## Exact schedule table

There are `16,384` data blocks. Translating a three-position support to gap
coordinates produces `134,209,536` raw gap records. The native builder maps
each coefficient triple to (2), sorts the keys, and sums the support
multiplicities. The result contains `134,201,345` distinct invariant keys.
Its masses satisfy

\[
\sum_J S_{\mathrm{data}}(J)=\binom{16384}{3}
=732{,}873{,}539{,}584
\]

and

\[
\sum_J S_{p_0}(J)=\binom{16384}{2}
=134{,}209{,}536.
\]

An independent Python implementation checks the invariant under all six
permutations on random triples. It also reproduces every record of a complete
12-block schedule table by direct enumeration.

## Exact relation enumeration

Let (L_{24}) be the complete weight-24 shell of the extended binary BCH
code. Goal 14 established

\[
|L_{24}|=6{,}855{,}968
\]

and

\[
\bigl|\{(x,y,z)\in L_{24}^3:x+y+z=0\}\bigr|
=4{,}018{,}336{,}896.
\tag{5}
\]

The entries in every relation are distinct. Dividing (5) by six gives
`669,722,816` unordered triples.

The affine coordinate group has order `16,256`. The scanner first partitions
(L_{24}) into the previously certified 532 word orbits. It then reduces the
unordered triples to 43,113 affine triple orbits. Of these, 3,682 have a
nontrivial stabilizer. The scanner expands one representative over explicit
coset representatives. The sum of the resulting orbit sizes is exactly
`669,722,816`, which independently reproduces (5).

The expanded triples contain `188,359,355` distinct (J)-keys. Only 17 keys
occur in the shifted schedule. They contain 1,846 unordered BCH relations.
The exact outer-word counts are

| Support family | Outer words |
|---|---:|
| Three data positions | 30,226,386 |
| (p_0) and two data positions | 30,235,625 |
| Total | **60,462,011** |

## Regression and audit chain

The same native path was run on the complete weight-22 shell. It reproduced
the known values

\[
1{,}365{,}504\text{ ordered relations}
\quad\text{and}\quad
227{,}584\text{ unordered triples}.
\]

It reduced these triples to 14 affine triple orbits and recovered the prior
exact outer-word count of zero. This regression exercises the BCH encoding,
affine action, word and triple stabilizers, invariant computation, sort, and
schedule merge.

The production weight-24 scan also checks the complete shell size and weight,
the certified affine word-orbit histogram, the ordered relation count from
Goal 14, the unordered orbit-size sum, sorted relation keys, and sorted unique
schedule records.

## Updated frontier

The closed `(24,24,24)` contribution is at most `2^-6.5547`. Together with the
previously closed Goal 14 profiles, the closed-profile sum is at most
`2^-6.55466`.

After removing `(24,24,24)`, the remaining occupation-three envelope is

\[
2^{87.28050}.
\]

The leading finite profiles are now the orderings of `(22,22,28)`,
`(24,24,26)`, and `(22,24,26)`.

The proof route remains tractable, but the implementation strategy must
change at this point. The complete BCH shells have sizes

\[
|L_{26}|=107{,}988{,}608,
\qquad
|L_{28}|=1{,}479{,}751{,}168.
\]

Materializing either shell as the next primary object would scale poorly.
The next case should instead keep a small shell in the projected positions.
In particular, `(22,22,28)` can be generated from pairs in (L_{22}): the
weight-28 word is their XOR and need not be enumerated as a shell. This is the
recommended next goal. The same projected-pair principle should be tested
before attempting `(24,24,26)`.

## Reproduction

```powershell
rustc -O -C target-cpu=native scripts/build_riffle_shiftalpha64_s3_schedule.rs -o constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/build_riffle_shiftalpha64_s3_schedule.exe
constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/build_riffle_shiftalpha64_s3_schedule.exe constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/goal15_s3_schedule_full.bin 16384
python scripts/verify_riffle_shiftalpha64_s3_ratio_quotient.py constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/goal15_s3_schedule_blocks12.bin 12
rustc -O -C target-cpu=native scripts/scan_riffle_shiftalpha64_equal_shell_s3.rs -o constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/scan_riffle_shiftalpha64_equal_shell_s3.exe
constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/scan_riffle_shiftalpha64_equal_shell_s3.exe constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/ebch128_weight22_messages.bin constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/goal15_s3_schedule_full.bin 22
constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/receipts/scan_riffle_shiftalpha64_equal_shell_s3.exe constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/ebch128_weight24_messages.bin constructions/riffle_shiftalpha64_bchblockperm_parallelacc_g4/goal15_s3_schedule_full.bin 24
python scripts/analyze_riffle_shiftalpha64_goal15_closure.py
```

Principal SHA-256 hashes are:

- production schedule receipt: `3E53F586B91D4FD9C11BA9AAC718DE23835B282B07F51AFD0443951E8E735304`;
- invariant and small-schedule audit: `E813583BD53AC4853300094D665194FD9CBF760DC5FB9E63F7A1666EF8939199`;
- weight-22 regression: `5A6DFEA6BC3AE3200BF8B23A13616346B87BB3E4B13EA718B21B224F77C2422F`;
- weight-24 exact scan: `52B839AF60D9B159D6F1C84147E111B269F828B7BADFBEB49BFC743EEE2887BB`;
- updated open-family envelope: `D5F334B793BA50D4238CF4A40DB1A5E66E56BD76A9C44CB26045B273CC227B87`;
- combined closure receipt: `421CEFE29EE06BA821AA15341573C1AAEA4C473778815C19754E220D17143DC2`.
