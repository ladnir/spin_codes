# One dense group with sixteen two-row groups

Fresh 256- and 384-bit replays now cover all \(28^{16}\) constituent-pair
assignments in source class \((8,0)\), cohort \((0,16,1)\). The independent
512-bit sum of rounded cell ceilings, including the exact count transfer,
gives

\[
\log_2 P\le-61000916136179/10^{12}=-61.000916136179.
\]

This is one source class, not all mixed messages or the whole code.
The original common-pair receipts below remain unchanged. The final
section records the new literals, full replays, and separate artifact pins.

## Original common-pair receipt

Fresh 256- and 384-bit replays certify the scoped class below. All 219
per-class ceilings and the raw aggregate agree across precisions. An
independent 512-bit sum of the rounded class ceilings gives

\[
\log_2 P\le-274659876277/2500000000=-109.8639505108.
\]

The three literal witnesses and disjoint assignment of all 219 dense row
sets remain unchanged from the FLOAT proposal.

This bounded experiment retains exact dense-mask mixtures and two support
floors. Its scope is one ordered source class: 17 active groups with
matching counts \((8,0)\). Sixteen groups have exactly two active RS rows;
the remaining group has a fixed row set \(A\) of size three through eight.
The sixteen constituent pairs have maximum multiplicity \(\mu=16\), so
they all use one common pair. The historical search results below are
labelled FLOAT. The directed result covers only this class, not all mixed
messages or the whole code.

The construction remains
`xor2-gf256-h4-alternating-swapped-pairs-k20-rate-half`. No encoder code,
native benchmark, or prior proof helper was changed.
Probability is over the ideal independent diagonal, route, swap, and
inner-key setup randomness. A particular pseudorandom seed is outside this claim.

## Fixed masks and the ordered resource bound

Fix a source support with matching counts \((8,0)\), its distinguished
dense group, all sixteen sparse constituent pairs, and the dense row set
\(A\). These choices remain fixed through all 32 RS columns.
The exact source count is

\[
K_{17}(8,0)=15708305994709246428578365440.
\]

There are 17 choices for the distinguished dense group and 28 choices for
the common sparse constituent pair. Each factor enters once, outside the
32-column polynomial.

For the dense group, use the counting-measure comparison proved in
[the one-group note](HD_XOR_FIXED_MASK_FLOOR.md). With \(a=|A|\), its density
factor is \(\beta_a=256^{32a}/255^{16a}\). Each comparison column samples
uniformly from \(\mathbb F_{256}^{A}\), including zero. Let \(N_A(T)\)
be the exact number of such inputs producing output mask \(T\) under the
XOR map. Thus \(\sum_TN_A(T)=256^a\) and \(N_A(\varnothing)=1\).
Independent output diagonals supply the conditional nonzero amplitudes.

Let \(\overline R_i(m)\) be the existing regional matrix bound for an
occupied subset of size \(m\), conditional on the fixed source matching
counts. The index \(i\in\{0,\ldots,7\}\) records labelled region order.
For a dense output mask \(T\), the sparse resource calculation uses

\[
R_i^T(m)=\overline R_i(m+\mathbf1[i\in T]).
\]

The same mask \(T\) is used in all eight factors of that column.
The implementation does not redraw a dense mask independently by region.

In a sparse column, let \(s\) groups have one nonzero constituent and
\(d\) groups have two. For \(k\) cancelled doubles, the necessary
occupancy constraints are

\[
d-k\le m_i\le s+d,\qquad
\sum_i m_i=7s+8d-6k.
\]

At \(\mu=16\), the previously proved pairwise upper bound is
\(m_i+m_j\le2(s+d)\), which is already implied by the coordinate bounds.
The ordered dynamic program therefore needs only the shifted sum
\(\sum_i(m_i-d+k)=7s+2k\). It returns a componentwise upper vector for
every feasible chronological matrix product applied to \(v=(1,t)^T\).
Each recurrence prepends its labelled regional matrix to the suffix.

Write \(U_{s,d,T}\) for the resulting vector after averaging
\(k\sim\mathrm{Binomial}(d,1/255)\). This upper bound holds for every
fixed sparse constituent geometry in the class. Its use across columns
does not assume independent geometry choices.

For each fixed \(A\), the two column caps are

\[
c_{A,0}(s,d)=\max_{j\in\{0,1\}}
 \frac{U_{s,d,\varnothing,j}}{256^a v_j},\qquad
c_{A,1}(s,d)=\max_{j\in\{0,1\}}
 \frac{\sum_{T\ne\varnothing}N_A(T)U_{s,d,T,j}}{256^a v_j}.
\]

Both the cancellation mixture and the dense-mask mixture precede the
potential quotient. The same potential is used for the empty and nonempty
dense branches.

## Two positive support degrees

For \(b\in\{0,1\}\), define

\[
P_{A,b}(X)=\sum_{s+d\le16}
 \binom{16}{s}\binom{16-s}{d}2^s
 c_{A,b}(s,d)(255X)^{s+2d}.
\]

The exponent of \(X\) counts nonzero entries in the 32 sparse RS rows.
Every such row has weight at least 17, so their total is at least 544.
Replacing the individual row floors by this aggregate floor is a relaxation.
The factor \(255^X\), followed by \(255^{-512}\), is the sparse shortened-code
support count. It is separate from the dense factor \(\beta_a\).

The exponent of \(Y\) in
\((P_{A,0}(X)+YP_{A,1}(X))^{32}\) counts nonempty dense output columns.
Invertibility of the outer map and a nonzero dense RS row imply at least
17 such columns. Define the positive tail

\[
T_A=\sum_{j\ge544,\,h\ge17}[X^jY^h]
 (P_{A,0}(X)+YP_{A,1}(X))^{32}.
\]

For \(0<z<1\) and \(0<t\le1\), the proposed class bound is

\[
\frac{K_{17}(8,0)\,17\,28\,z^{-734003}}{255\,t\,255^{512}}
 \sum_{A:\,3\le|A|\le8}\beta_{|A|}T_A.
\]

The projective factor \(1/255\) and the endpoint factor \(1/t\) occur once.
The sum visits all 219 fixed row sets; no additional binomial row-set factor
is applied.

The evaluator extracts the two-variable tail without a Chernoff relaxation.
Its scalar caps commute, so it uses

\[
[Y^h](P_0+YP_1)^{32}=\binom{32}{h}P_0^{32-h}P_1^h.
\]

All convolutions operate on nonnegative coefficients in FLOAT log space.
There is no FFT subtraction, coefficient floor, or interpolation.

## Initial experiment and bounded witness search

The inherited two-row witness was
\(z=996598217/10^9\), \(\log_2t=-135/2\).
Its all-219 union is +644.9849214398475 bits, a failed bound.
The dense-mask resource bank took about two seconds; the whole evaluation
took 19.45 seconds. Every row set's largest dense degree was 32.
For \(A=255\), the largest sparse degree was 1021, far above the floor.
The degree-31 dense contribution was 49.08 bits below degree 32.
Thus this failure is not dominated by empty dense columns.

A separate search used at most 36 evaluations, only for \(A=255\), with
\(z\in[.99,.9998]\) and \(\log_2t\in[-100,-20]\).
The best literal was

\[
z=8970988817957597/9007199254740992,\qquad
\log_2t=-2266954815/33554432.
\]

It gives −218.56866698083786 bits for \(A=255\), with no search boundary hit.
One all-219 evaluation at that literal gives +164.88515213189189 bits;
the obstruction shifts to arities three and four. Choosing each fixed
\(A\)'s better of these two witnesses gives +38.301069654511835 bits.
This is a deterministic choice by message class, not by sampled setup.
The largest surviving row sets were \(A=7\) at arity three and \(A=240\)
at arity four. Each received exactly 32 additional evaluations and one
all-219 evaluation at its best literal. Neither search hit a boundary.
These two searches took approximately 133 seconds in total; their all-row-set
evaluations took another 36 seconds. No other source or multiplicity class
was evaluated.

The final assignment uses the following three witnesses:

| Label | Rational \(z\) | Rational \(\log_2t\) | Assigned row sets |
|---|---|---|---:|
| `a8_fit` | `8970988817957597/9007199254740992` | `-2266954815/33554432` | 104 |
| `a3_mask7_fit` | `1121881298265941/1125899906842624` | `-66751297/1048576` | 31 |
| `a4_mask240_fit` | `8973939342378351/9007199254740992` | `-2386334719/33554432` | 84 |

The inherited witness is no longer selected. The labelled masks, not merely
their arities, determine the assignment in `hd_xor_mixed_mask_witnesses.py`.
Its `checked_classes()` interface returns the literal witnesses and checks
that the mask lists partition precisely all 219 eligible row sets.

The combined FLOAT arity unions are:

| Dense arity | Row sets | Log2 proposal |
|---:|---:|---:|
| 3 | 56 | −131.9782316258846 |
| 4 | 70 | −109.86395097346264 |
| 5 | 56 | −133.11123996704458 |
| 6 | 28 | −217.33182735941955 |
| 7 | 8 | −303.9702269364663 |
| 8 | 1 | −218.56866698083786 |

Their union is −109.86395051080125 bits. The largest individual term is
\(A=195\), at −110.66160699442929 bits. Selecting a witness by the fixed
message class does not require a witness-bank union factor. Independent
read-only review reproduced every assignment and the combined FLOAT value.

The source files are `theory/hd_xor_mixed_mask_floor.py`,
`theory/hd_xor_mixed_mask_search.py`, and `theory/hd_xor_mixed_mask_refine.py`.
Each has a corresponding `_tests.py` file. Twelve focused tests pass:
they cover labelled noncommuting offsets, cancellation weights, averaging
before the quotient, exact rational bivariate enumeration, separate floors,
normalization, fixed-class assignment, frozen literals, and search budgets. Independent
read-only review found no mathematical gap in the offset and tail construction.
The combined regression at the proposal stage passed 33 tests. These checks
alone did not turn the FLOAT proposal into a certificate.

From `theory`, the initial gate and two search records can be reproduced
into fresh output paths:

```powershell
C:/Python314/python.exe -m unittest hd_xor_mixed_mask_floor_tests hd_xor_mixed_mask_search_tests hd_xor_mixed_mask_refine_tests
C:/Python314/python.exe hd_xor_mixed_mask_floor.py --output mixed_inherited_replay.json
C:/Python314/python.exe hd_xor_mixed_mask_search.py --output mixed_a8_search_replay.json
C:/Python314/python.exe hd_xor_mixed_mask_refine.py --output mixed_a3_a4_refine_replay.json
```

The refinement command reads the two named frozen proposal artifacts in its
module. No other ordered source class or multiplicity stratum is covered.

## Fresh directed replay and independent summation

`hd_xor_mixed_mask_certificate.py` rebuilds every regional matrix in Arb
arithmetic. It retains labelled dense offsets and takes explicit upper
endpoints before each resource maximum. Suffix computations are shared only
when their fixed dense-mask suffixes agree. It discards only zero-probability
dense masks and partial sums unable to reach any requested final target.
The exact cancellation and dense mixtures still precede the common-potential
quotient. Both polynomial floors are retained without subtraction.

The initial 256-bit check replayed only the previously largest class,
\(A=195\). Its ceiling was −110.661606994429 bits. It used 27 dense output
masks and took 4.60 seconds, including 3.72 seconds for the resource bank.

After independent review of that path, full replays ran serially at 256
and 384 bits. They took 81.22 and 85.86 seconds, respectively. Both returned
the exact raw aggregate ceiling

\[
-109863950510801/10^{12}.
\]

The three banks used 256, 153, and 247 distinct dense masks. Reusing each
bank across its assigned row sets did not change the fixed-A comparison law.

`hd_xor_mixed_mask_verify.py` independently checks the exact construction
scope and the complete 219-mask assignment. It derives arity from each mask,
checks every literal witness, and requires equality of every rational class
ceiling and the raw aggregate across precisions. It then recomputes
\(\log_2\sum_A2^{u_A}\) at 512 bits, where \(u_A\) is the rounded class
ceiling. The resulting upper ceiling is −109.8639505108 bits, as stated above.
No cached aggregate or FLOAT value participates in this sum.

The verifier pins the reviewed source bytes and receipt bytes. Those hashes
identify evidence from observed fresh replays; hashes do not prove arithmetic.
Seven directed recurrence tests and four verifier tests were added. The
combined regression now passes 44 tests. Independent read-only review
confirmed the recurrence, reachability pruning, matrix order, and factors.

To replay without replacing the saved evidence, use a fresh output directory
from `theory`:

```powershell
New-Item -ItemType Directory -Path mixed-directed-replay
C:/Python314/python.exe hd_xor_mixed_mask_certificate.py --precision 256 --all --output mixed-directed-replay/hd_xor_mixed_mask_8_0_mu16_receipt256.json
C:/Python314/python.exe hd_xor_mixed_mask_certificate.py --precision 384 --all --output mixed-directed-replay/hd_xor_mixed_mask_8_0_mu16_receipt384.json
C:/Python314/python.exe -c "import json; from hd_xor_mixed_mask_verify import verify; print(json.dumps(verify('mixed-directed-replay')))"
C:/Python314/python.exe -m unittest hd_xor_mixed_mask_certificate_tests hd_xor_mixed_mask_verify_tests
```

SHA-256 pins, with paths relative to `theory`:

| File | SHA-256 |
|---|---|
| `hd_xor_mixed_mask_floor.py` | `305bd66ea8f251c3da21920e6785ad4ec16aff0da9d3e4faa4a707c86f569d25` |
| `hd_xor_mixed_mask_search.py` | `c8ff5083b4db727b385d79da182242699e47cada60ace30f5339a5b58bfc77f2` |
| `hd_xor_mixed_mask_refine.py` | `c7d3254378a0bf518ceb36b5b00bdc33c64624d97054da74309068da9d44e28f` |
| `hd_xor_mixed_mask_witnesses.py` | `6b6c4ceb63a9624ccdcbbf7d5ab78beabd19caa3224dbb193eb68509c4e09a61` |
| `hd_xor_mixed_mask_floor_8_0_mu16_float.json` | `2b3fe0af20331ca123c5370d729c4ca031e3357aab98be3ea09436d24fd3546c` |
| `hd_xor_mixed_mask_8_0_mu16_search36.json` | `88faefb4973e9a10599da7b6997158a8f596886b81b5fcb0e38ce38f0d025cc2` |
| `hd_xor_mixed_mask_8_0_mu16_refine64.json` | `0abded252784870277a1b285d93e522b9bdb717288a42448ce49469097fc1fc5` |
| `hd_xor_mixed_mask_certificate.py` | `b243ecbc00527ef39f151a276b0206323f5b4cfdd14361b380326bbf75dc8f6b` |
| `hd_xor_mixed_mask_certificate_tests.py` | `33dfc987a4c6be9ca90dc55dc8dd17ab185142cb7451cbdb80536d98f66a1de9` |
| `hd_xor_mixed_mask_verify.py` | `ff4db8203a42f3f2e6b87f4e88c29894236356c98581c342f5fbda072c0ed283` |
| `hd_xor_mixed_mask_verify_tests.py` | `1dc10b29f75a3e93f740d3e5a408d3a306568275cdd95a825cdc3c954ccf350a` |
| `hd_xor_mixed_mask_A195_receipt256.json` | `54e0b76793f9581a56ecee2d6d9e1e151372a950663f367f4b885528f05e450a` |
| `hd_xor_mixed_mask_8_0_mu16_receipt256.json` | `525871e0390b8c9d192cf041d3a049c144e3c2f89f2f42ffcd505acecab3aae9` |
| `hd_xor_mixed_mask_8_0_mu16_receipt384.json` | `99689f49085914418fdfa2c46aafe534b31ae66ef93e053d0f1275ca16b0c12a` |
| `hd_xor_mixed_mask_8_0_mu16_dual_verified.json` | `1db1308fc80eb9787ffbcd0a6eb694113b363dc17ce2d1cb596ecc18bb1d1e3b` |

## Separate transfer to all constituent-pair assignments

The preceding receipts still count exactly 28 common-pair choices. A
separate count transfer covers all \(28^{16}\) labelled pair assignments
for the same cohort \((0,16,1)\) and the same source class \((8,0)\).
Its directed log-probability ceiling is

\[
-\frac{7550725335987}{200000000000}
=-37.753626679935.
\]

This does **not** meet the 40-bit target and is not a whole-code bound.

To justify the transfer, fix the sixteen constituent pairs for all 32
columns. In one column, let \(s\) and \(d\) count singleton and double
inputs, and let \(k\) double inputs cancel outside their pair endpoints.
Write \(q_i\) for singleton inputs missing region \(i\), and \(e_i\)
for the degree of \(i\) among the cancelled pairs. The actual occupancy is

\[
m_i=s+d-k-q_i+e_i,
\qquad d-k\le m_i\le s+d,
\qquad \sum_i m_i=7s+8d-6k.
\]

These constraints hold for arbitrary pair assignments. The reviewed
directed recurrence used only these coordinate and total-sum constraints;
the multiplicity-16 pair constraint is already implied by the coordinate
upper bound. The exact binomial cancellation law also does not depend on
the pair labels. Therefore its caps hold for every fixed pair assignment,
without resampling geometry between columns. The frozen witness choice
depends only on the fixed dense row set, so the same choice remains valid.

Replacing the counted 28 choices by \(28^{16}\) multiplies the old bound
by exactly \(28^{15}\). Starting from the independently verified rounded-cell
union, the separate calculation is

\[
\log_2\!\left(28^{15}
 2^{-274659876277/2500000000}\right)
\le -7550725335987/200000000000.
\]

`hd_xor_mixed_geometry_transfer.py` performs only this positive mass
multiplication at 512 bits, cross-checking logarithmic addition. No operator
was rerun and no witness changed. Two new small tests exhaust all labelled
two-group pair assignments and column states, and independently check the
ceiling at 1024 bits. Together with the four existing verifier tests, all
six tests pass.

From `theory`, using a fresh output path:

```powershell
C:/Python314/python.exe hd_xor_mixed_geometry_transfer.py --output mixed_geometry_transfer_replay.json
C:/Python314/python.exe -m unittest hd_xor_mixed_geometry_transfer_tests hd_xor_mixed_mask_verify_tests
```

| File | SHA-256 |
|---|---|
| `hd_xor_mixed_geometry_transfer.py` | `1bdbc69fc92bb2e5e98ec2fe67dfcbc94a2942259e2a71540543bd670616f90e` |
| `hd_xor_mixed_geometry_transfer_tests.py` | `e373d4ce515e7e6151e5f7bff5e86c8e9cbd51278f213505baecb8ade02b634d` |
| `hd_xor_mixed_8_0_all_geometry_transfer512.json` | `f508b0514c658ba6fdaa7aadfc9387a6f593472e49349c24e0b79eb7911e2bb3` |

The next useful proof target is to recover at least 2.25 bits for this
same all-geometry source class, with additional margin for unioning other
classes. This transfer alone does not authorize any source-class extension.

## Bounded all-geometry witness improvement: FLOAT proposal

A new bounded search proposes a 61.000916136180436-bit bound for the same
all-geometry source class. This subsection records the **FLOAT proposal**;
the next section records its fresh directed replay. At the proposal stage,
the checked all-geometry bound was still 37.753626679935 bits. The construction, source
class \((8,0)\), cohort \((0,16,1)\), and support floors are unchanged.

`hd_xor_mixed_all_geometry_refine.py` starts with the worst fixed dense
row set, \(A=195\), at its previously assigned literal. It runs at most
32 evaluations within \(z\in[.99,.9998]\) and \(\log_2t\in[-100,-20]\).
After a material improvement, it evaluates all 219 row sets once at the
winning literal and recombines their bounds. Only if that union misses
60 bits does it repeat this procedure for the new worst row set.

The first search improved the common-pair bound for \(A=195\) from
−110.66160699442929 to −187.15732849239964. Its new full bank raised the
all-geometry union from 37.7536 to 59.3259 bits. The second search therefore
selected \(A=224\). It improved that row set from −134.51860814303598 to
−204.27369884024063 in the common-pair normalization. Neither search hit
a parameter boundary. Each used exactly 32 evaluations.

The two exact proposed literals are:

| Search row set | \(z\) | \(\log_2t\) | Assigned row sets in final schedule |
|---|---|---|---:|
| 195 | `4487137477776021/4503599627370496` | `-70552129/1048576` | 35 |
| 224 | `8975046825465429/9007199254740992` | `-2381845503/33554432` | 6 |

The remaining 178 row sets retain their previously frozen literals:
6 use `a3_mask7_fit`, 104 use `a8_fit`, and 68 use `a4_mask240_fit`.
Each assignment depends only on the fixed row set, not on setup randomness
or the fixed constituent-pair assignment. No witness-bank union penalty is
needed. The separate factor \(28^{15}\) is applied once to each old
common-pair-normalized row-set bound.

The final FLOAT union is −61.000916136180436. Arity-5 row sets dominate:
their union is −61.0009161361805, whereas arities 3 and 4 give
−108.15553350905363 and −105.84090511374325. The largest individual term is
\(A=31\), at −64.35527834210386. These are diagnostic FLOAT values,
not directed ceilings.

The complete bounded run took 169.80 seconds. Its two full-bank evaluations
took 18.23 and 18.33 seconds. The saved artifact contains all 64 trials,
both complete new banks, and the exact 219-class assignment. Four new tests
check scope, the count transfer, early stopping, and independent
recombination of that artifact. Together with the previous refinement,
transfer, and verifier tests, all 14 tests pass. Hash checks confirm that
the old helpers, witnesses, and directed receipts remain unchanged.

From `theory`, use a fresh output path:

```powershell
C:/Python314/python.exe hd_xor_mixed_all_geometry_refine.py --output mixed_all_geometry_refine_replay.json
C:/Python314/python.exe -m unittest hd_xor_mixed_all_geometry_refine_tests hd_xor_mixed_mask_refine_tests hd_xor_mixed_geometry_transfer_tests hd_xor_mixed_mask_verify_tests
```

| File | SHA-256 |
|---|---|
| `hd_xor_mixed_all_geometry_refine.py` | `6375d2ffb3983e6c5bee092ddd93d9ac051aea72ddd281935608c308ef0c86d2` |
| `hd_xor_mixed_all_geometry_refine_tests.py` | `4d32a773de77df855da02812332ad85e71b7db3c2d59a756eb18cbc08b9a5795` |
| `hd_xor_mixed_all_geometry_8_0_refine64.json` | `44fd35d8bf2e78391572a14d071cd283dc547f7bc204b3802c180132eaaab4e6` |

The proposed replay could have reused the 178 matching old receipt entries
and recomputed only 41 row sets. Instead, the following verification
recomputed all 219 row sets for simpler independent checking. No further
witness fitting or source-class extension was performed.

## Verified all-geometry replay

`hd_xor_mixed_all_geometry_certificate.py` pins the complete proposal and
calls the existing audited regional, offset-DP, mixture, and positive-tail
primitives. It neither mutates the old replay function's globals nor changes
any pinned helper. Every dense mask stays fixed through the eight labelled
regions; every dense input row set stays fixed through the 32 columns.

The recurrence's `mu=16` argument selects the universal box/sum relaxation,
not a restriction to common pairs. The preceding resource argument proves
that relaxation valid for every fixed constituent-pair assignment. Each
receipt explicitly distinguishes the base 28-choice normalization from the
derived \(28^{16}\)-choice scope. It records both cell ceilings.

Full replays ran serially at 256 and 384 bits, taking 95.60 and 102.34 seconds.
Both precisions returned the same raw all-geometry aggregate ceiling:

\[
-3050045806809/50000000000=-61.00091613618.
\]

All 219 common-pair-normalized cell ceilings and all 219 all-geometry cell
ceilings agree exactly across precisions. Both raw aggregates also agree.
For each precision, the 178 unchanged-witness cells agree exactly with
their original receipt entries.

The verifier checks the complete mask partition, derived arities, exact
literal assignment, construction, and source scope. It then sums the
rounded common-pair cell masses in a fresh 512-bit precision context and
multiplies by \(28^{15}\) exactly once. The resulting ceiling is

\[
-61000916136179/10^{12}=-61.000916136179.
\]

An independent 1024-bit sum reproduced this ceiling. No FLOAT value, cached
aggregate, or witness-bank multiplicity enters that sum. Hashes bind the
artifacts to the observed fresh replays; they do not prove arithmetic by
themselves. Five new tests cover scope, exact mask/literal coverage,
count normalization, the absence of an extra bank factor, and the observed
dual receipts. The combined focused regression passes 26 tests.
The final regression across all nine mixed-gate, mixture, search, replay,
transfer, and verification test modules passes 40 tests in 1.01 seconds.

From `theory`, replay serially into a fresh directory:

```powershell
New-Item -ItemType Directory -Path mixed-all-geometry-replay
C:/Python314/python.exe hd_xor_mixed_all_geometry_certificate.py --precision 256 --output mixed-all-geometry-replay/hd_xor_mixed_all_geometry_8_0_receipt256.json
C:/Python314/python.exe hd_xor_mixed_all_geometry_certificate.py --precision 384 --output mixed-all-geometry-replay/hd_xor_mixed_all_geometry_8_0_receipt384.json
C:/Python314/python.exe hd_xor_mixed_all_geometry_certificate.py --verify --directory mixed-all-geometry-replay --output mixed-all-geometry-replay/dual_verified.json
C:/Python314/python.exe -m unittest hd_xor_mixed_all_geometry_certificate_tests hd_xor_mixed_mask_certificate_tests hd_xor_mixed_all_geometry_refine_tests hd_xor_mixed_mask_refine_tests hd_xor_mixed_geometry_transfer_tests hd_xor_mixed_mask_verify_tests
```

| File | SHA-256 |
|---|---|
| `hd_xor_mixed_all_geometry_certificate.py` | `e39f4d8419c8361269896acf21f34b328d3b76987b2bee49690a5d73a0f59e63` |
| `hd_xor_mixed_all_geometry_certificate_tests.py` | `c2c59a8eaa9c3633d5a83da57ec5dd10352922c4eb8c82448c74afdde3ad764f` |
| `hd_xor_mixed_all_geometry_8_0_receipt256.json` | `071a50922931379a0878c78a9c0ea594ca62ccca5e1e728123fa505316ec6642` |
| `hd_xor_mixed_all_geometry_8_0_receipt384.json` | `81da27a2ddcc6807aa5c908d299630fb1c5c812bdaafb236f2228943bc071df4` |
| `hd_xor_mixed_all_geometry_8_0_dual_verified.json` | `af53e419971b68619c23489a71090c09a1d706b13288d833c1c1bb315178a095` |

The next structural question is a broader source/cohort family, rather
than additional margin for this closed class. The present receipt does not
cover another source matching count, another number of active groups, or
more than one dense group.
