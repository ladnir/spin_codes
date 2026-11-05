# Whole-code coverage plan for the paired XOR2 construction

The current XOR2 results do not yet prove whole-code distance. They close
several complete message classes, while a broad terminal argument removes
the highest constituent-row counts. The next useful calculation should
cover a family of cohorts, rather than improve one already closed cell.

This audit concerns the ideal construction
`xor2-gf256-h4-alternating-swapped-pairs-k20-rate-half`. Each of 8192 outer
groups contains eight RS[32,16] constituents over GF(256). Each column uses
independent nonzero diagonals in
\(D_2(I+\mathbf1\mathbf1^T)D_1\). The 256 regions use alternating adjacent
pair matchings, independent uniform pair permutations, and independent
fair swaps. The existing 64-input, 16-state, rank-four inner is unchanged.
Its state starts at zero and persists through every region.

Setup is sampled once and then fixed for all messages. The unchanged goal is

\[
\Pr[\exists m\ne0:\ \operatorname{wt}(Gm)\le734003]\le2^{-40},
\qquad K=1048576,\quad N=2097152.
\]

The probability concerns the ideal setup ensemble, not a selected native
seed. The native sixteen-component scalar extension needs no extra
factor-sixteen union: its record distance equals the scalar distance.
The [retained Toeplitz checkpoint](HD_SHIFTED_CLOSED_CERTIFICATE.md)
already has a complete 55.092738508223-bit certificate; it remains separate.

## What a complete XOR2 partition must cover

Fix a nonzero message. Let \(G\) be its set of active outer groups and
\(r=|G|\). For each \(g\in G\), fix the nonempty set
\(A_g\subseteq\{0,\ldots,7\}\) of nonzero constituent RS words. These
sets remain fixed across all 32 columns. Write

\[
n_1=\#\{g:|A_g|=1\},\quad n_2=\#\{g:|A_g|=2\},\quad
n_H=\#\{g:|A_g|\ge3\},\quad r=n_1+n_2+n_H.
\]

Let \(a,b\) count full pairs of \(G\) under the two ordered source
matchings. Their exact joint counts partition the indexed group supports.
For the two-row cohort, let \(\mu\) be the largest multiplicity of a
fixed constituent pair among its 28 possible pairs. Its exact stratum
counts sum to \(28^{n_2}\).

A coverage calculation must retain, or safely dominate, all of these
choices. Cohort colorings contribute \(r!/(n_1!n_2!n_H!)\); singleton row
choices contribute \(8^{n_1}\). Dense groups require their fixed row sets
or a proved joint envelope. These factors are counted once, outside the
32-column product. Neither a constituent pair nor a dense row set may be
sampled anew at each column.

The existing Toeplitz ledger indexed by \((r,a,b)\) therefore cannot be
relabelled as an XOR2 ledger. XOR2 additionally has constituent-row
structure. Count symmetry \((a,b)\leftrightarrow(b,a)\) also does not
justify reversing noncommuting regional operators from zero initial state.
Finally, passing 40 bits separately for many classes is insufficient:
the complete union must pass 40 bits.

## Directed results and their exact scope

The following are current paired-route results, not scalar-route imports.
The linked notes identify literal witnesses and replay artifacts.

| Message class | Source and constituent geometry covered | Directed log2 upper bound |
|---|---|---:|
| \(r=17,(n_1,n_2,n_H)=(0,17,0)\) | All 81 ordered \((a,b)\) classes and all fixed pair assignments | −50.891481584655 |
| \(r=1,n_H=1\) | All 8192 positions and all 219 row sets of arity 3 through 8 | −164.796015302006 |
| \(r=17,(0,16,1)\), \((a,b)=(8,0)\) | All 219 dense row sets and all \(28^{16}\) pair assignments | −61.000916136179 |

Sources are [the homogeneous paired gate](HD_XOR_PAIRED_GATE.md),
[the one-dense-group calculation](HD_XOR_FIXED_MASK_FLOOR.md), and
[the mixed exact-mask calculation](HD_XOR_MIXED_MASK_FLOOR.md).
The last row supersedes the earlier common-pair-only receipt for that
source class. It does not cover another \((a,b)\), another \(r\), or
two dense groups.

The sibling class \((n_1,n_2,n_H)=(1,16,0)\) now has a complete 324-cell
FLOAT proposal: 81 ordered source classes and four multiplicity bins.
Its union is −70.35414804665. This is not yet a directed result in this
audit. `hd_xor_sparse_mixed_certificate.py` independently reconstructs
every count and positive polynomial from the saved literals. Read-only
review found no soundness defect, and its four unit tests passed locally.

For that wrapper, \(s=\text{active singleton groups}+\text{single pair rows}\)
and \(d=\text{double pair rows}\). The multiplicity restriction still
concerns only the 16 two-row groups. The positive support polynomial has
degree at most 33 per column, aggregate floor 561, and normalization
\(255^{528}\). Source count, \(17\) cohort positions, \(8\) singleton
row choices, and the exact pair-geometry count occur once. Cancellation
branches are averaged before the potential quotient.

## A terminal bound for all cohort mixtures

The following argument avoids any source-support or cohort enumeration.
Let

\[
A_{\rm tot}:=\sum_{g\in G}|A_g|
\]

be the total number of nonzero RS[32,16] constituent words. For an integer
threshold \(A_0\), let \(E_{A_0}\) be the event that a nonzero message
with \(A_{\rm tot}\ge A_0\) has output weight at most \(L=734003\).
Then

\[
\Pr[E_{A_0}]
\le \frac{\operatorname{Ball}_{256}(N,L)}{255^{16A_0+1}},
\qquad
\operatorname{Ball}_{q}(N,L)=\sum_{j=0}^{L}\binom Nj(q-1)^j.
\]

To prove this, first use the ideal independent, unnormalized \(D_1,D_2\)
law. The native scalar gauge has the same distribution of column maps;
conditioning directly on normalized keys would not preserve the needed
independence. Condition on all routing, inner setup, and outer \(D_2\)
keys, leaving all outer \(D_1\) entries independent and unexposed.

Fix a candidate final word \(y\). Invert the inner map and routing, apply
\(D_2^{-1}\), and then apply \(H^{-1}\), where
\(H=I+\mathbf1\mathbf1^T\) and \(H^{-1}=H\). The inner inverse exists
because each emitted packet is its input plus a function of earlier
inputs. This determines a unique prescribed \(D_1\)-scaled constituent
array. Nonzero diagonal entries preserve its row supports. Thus its
zero/nonzero row classification, and \(A_{\rm tot}(y)\), are determined
before exposing \(D_1\).

For a prescribed nonzero row support of size \(w\), the shortened MDS
code has dimension \(w-16\), if such a word exists. An information set
injects its exact-support words into \((\mathbb F_{256}^{*})^{w-16}\).
Their count is therefore at most \(255^{w-16}\), not a bound with base
256. Independent nonzero diagonal entries give probability \(255^{-w}\)
per word. The row membership probability is at most \(255^{-16}\);
supports smaller than 17 have probability zero. A zero row contributes one.
Independence across rows gives \(255^{-16A_{\rm tot}(y)}\).

Sum over low-weight nonzero \(y\) with \(A_{\rm tot}(y)\ge A_0\).
This class is invariant under nonzero scalar multiplication, so divide
once by 255. There is no further support, cohort, or geometry union:
the inverse already determines those objects. Including zero in the ball
count only enlarges the upper bound. This proves the displayed inequality
even though a bad message may be chosen after setup.

Use the same tilted-ball mark as the retained terminal receipt,

\[
z=\frac{734003}{347602995},\qquad
\log_2\operatorname{Ball}_{256}(N,L)
\le N\log_2(1+255z)-L\log_2z.
\]

Fresh local Arb evaluations at both 256 and 384 bits agree exactly:

| Threshold \(A_0\) | Rational log2 ceiling | Decimal |
|---:|---:|---:|
| 61190 | −37865611218811/1000000000000 | −37.865611218811 |
| 61191 | −20721908276069/125000000000 | −165.775266208552 |
| 61192 | −146842460599147/500000000000 | −293.684921198294 |

The last row is exactly the old `terminal_receipt` ceiling: its exponent
\(128\cdot7649\) equals \(16\cdot61192\). The two preceding rows change
only this integer exponent. The calculation uses `Parameters(field=256,
state=16,limit=734003)`, `ball_log2`, and `rational_ceiling` from the existing
proof helpers. The inverse-image argument was independently challenged
and agreed with the separate rigor reviewer.

Thus \(A_{\rm tot}\ge61191\) is a broad directed tail, with no new route
assumption. It is not a bound for every message with \(r\ge7649\): mixed
groups may have far fewer active constituents. The remaining medium
\(A_{\rm tot}\) range still requires proof.

## Which older bounds transfer

Two limited transfers are valid. First, when \(r=1\), each regional support
has size zero or one. A uniform pair permutation followed by a fair swap
places its single active coordinate uniformly among all 8192 locations.
The scalar-route one- and two-row receipts therefore apply at \(r=1\).
This argument does not transfer scalar-route receipts for \(r\ge2\).

Second, suppose every active group has all eight constituents nonzero.
Its expected counting measure is dominated by

\[
\beta_8\nu,\qquad
\beta_8=256^{256}/255^{128},
\]

where \(\nu\) is uniform on 256 field coordinates. The full-row comparison
input is uniform in \(\mathbb F_{256}^{8}\) at each column. The invertible
map \(D_2HD_1\) preserves that law. This is exactly the legacy
`inherited-diagonal` outer interface of the paired Toeplitz proof. Its
regional operators, source counts, and legacy-density receipts consequently
transfer to the pure-eight-row class. There is no row-choice factor.
Discarding individual RS floors in this domination is safe; it does not
mean that the comparison law retains those floors.

The sharper `invertible-toeplitz` prefix does not transfer unchanged.
For a fixed \(r,\alpha\), conversion to the legacy prefix adds

\[
\alpha r\,96\log_2(256/255)
\approx0.5420700615496361\,\alpha r
\]

to its logarithmic upper bound. This scalar correction needs directed
evaluation, but no regional-table replay. A FLOAT audit of 3811 retained
potential rows from \(r=257\) through 4067 leaves 49 individual weights
above −40 after correction; their worst value is +65.0095 at \(r=3285\).
The union of those 3811 rows is +65.3276. Hence this is useful partial
reuse, not a completed pure-eight-row certificate. Missing weights are
1130–1135, 1200–1206, 1256–1259, 1322–1325, 1378–1383, 1441,
1646–1651, 1817–1824, 2679–2683, and 3284–3285. Other retained legacy
checkpoints, including the dense boxes, are not part of this scalar audit.

The older scalar-route homogeneous one-row certificate and scalar-route
multi-group two-row receipts do not transfer wholesale. Neither do
unswapped fixed-orientation receipts without a separate domination argument.
Construction tags must remain explicit when recording either valid transfer.

## Smallest useful next broad calculation

The first target should be all sparse cohorts \(n_H=0\) at a fixed \(r\),
then a bounded interval of \(r\). For a fixed source class, multiplicity
bin, and witness, the expensive caps \(c_{s,d}\) do not depend on the split
between \(n_1\) and \(n_2\). Reuse that cap bank for every applicable
cohort polynomial. Only combinatorial coefficients, the RS aggregate floor,
normalization, and once-only class counts change. At \(r=17\), this means
18 cohort splits rather than another isolated boundary cell. Classes with
\(n_2=0\) require no pair-multiplicity partition. Existing helper work
limits remain explicit; larger \(r\) requires a checked extension.

Two bounded controls explain why retaining the column caps matters.
Using only saved scalar-route literals, a geometry-blind single potential
envelope for the paired one-row family passes \(r=1,\ldots,7\), but gives
+338.14 bits at \(r=8\) and +1398.87 at \(r=17\). Retaining exact
\((a,b)\) source counts improves these to +64.36 and +1239.23.
These were complete FLOAT sums at each tested weight, with 13 fixed banks
and no fitting. They reject that inherited envelope, not the construction.
The whole-column contractions are therefore the next reusable proof object.

Dense cohorts with \(n_H\ge2\) remain the principal missing interface.
Powering the one-dense-group result does not preserve fixed row sets or
each group's 17-nonempty-column floor. A whole-word Hölder envelope or a
proved joint mask comparison is needed before scaling that calculation.
Such an envelope must retain each \(A_g\) through all 32 columns and pay
its choice only once. The terminal lemma removes the highest total-row
counts, but does not solve this intermediate-density problem.

Next iteration: finish the already proposed sparse mixed replay, reuse its
cap banks across the sparse cohort family, and establish the multi-dense
interface before commissioning another collection of isolated cells. A final
ledger must cover every remaining \((G,(A_g))\) and sum its directed masses;
none of the partial results above substitutes a weaker goal for that task.

## Follow-up: shared caps across all sparse cohorts at r=17

`theory/hd_xor_sparse_cohorts.py` implements the proposed reuse without
changing any pinned proof helper. It first screened source \((a,b)=(8,0)\),
then explicitly screened all 81 source classes. Both runs use only these
two fixed pairs \((z,\log_2 v_1)\), with \(v_0=1\):

\[
\left(\frac{996248217}{10^9},-\frac{127}{2}\right),\qquad
\left(\frac{4487948476156805}{4503599627370496},
      -\frac{8814329}{131072}\right).
\]

There was no fitting. Two raw regional banks were built once. For each
source and each bank, one `MultiplicityCapsFloat` instance supplied the
caps used by all 18 cohort polynomials. Multiplicity bins are
1–8, 9–14, 15, 16, and 17, truncated to \(n_2\). When \(n_2=0\), one
empty-geometry class uses the dummy cap parameter \(\mu=1\); only entries
with \(d=0\) enter its polynomial, so this parameter imposes no condition.

Each class pays its exact source count, \(\binom{17}{n_1}\),
\(8^{n_1}\), and its exact fixed-pair assignment count once. Its aggregate
RS floor is \(17(n_1+2n_2)\), and its normalization exponent is
\(16(n_1+2n_2)\). The better of the two fixed witnesses is selected per
class, before summing classes. That selection requires no bank-count factor.

For source \((8,0)\), all 18 cohorts pass a 40-bit FLOAT screen; their union
is −49.60373901560979. That run took 5.82 seconds. The complete source
screen took 396.50 seconds and retained all 2673 selected cells, including
failed controls. Its source counts sum exactly to \(\binom{8192}{17}\).
The first source artifact remains byte-for-byte unchanged, and its cohort
bounds reproduce exactly in the full run.

The following bounds sum all source classes for the indicated cohort.
Every row has \(n_H=0\) and \(n_2=17-n_1\). These are FLOAT proposals,
not directed ceilings.

| \(n_1\) | \(n_2\) | Complete-source log2 screen |
|---:|---:|---:|
| 0 | 17 | +151.7582895049 |
| 1 | 16 | +68.0292898897 |
| 2 | 15 | −35.4960930223 |
| 3 | 14 | −41.4914801693 |
| 4 | 13 | −62.3513232325 |
| 5 | 12 | −82.3673757356 |
| 6 | 11 | −135.9310724601 |
| 7 | 10 | −181.4103958771 |
| 8 | 9 | −276.7754135847 |
| 9 | 8 | −274.5940198930 |
| 10 | 7 | −335.3439008062 |
| 11 | 6 | −346.4250170210 |
| 12 | 5 | −357.7159030378 |
| 13 | 4 | −369.3851534581 |
| 14 | 3 | −381.4834510315 |
| 15 | 2 | −394.0961803915 |
| 16 | 1 | −407.3869597113 |
| 17 | 0 | −421.7682546693 |

The raw union of all 18 cohorts is +151.7582895049. The union of the 17
nonhomogeneous cohorts is +68.0292898897. Substituting the preserved
homogeneous directed bound therefore does not make this proposal pass.
These failures do not invalidate the older homogeneous receipt or the
successful \((1,16,0)\) screen: their witness schedules are not restricted
to these two pairs.

The 15 cohorts with \(n_1\ge3\) have combined screen −41.49147941122387.
For \((n_1,n_2)=(2,15)\), exactly one individual cell misses 40 bits:
\((a,b)=(0,7)\), \(\mu=15\), at −35.51892616207897. Including this entire
cohort gives −35.473654118793405 for the 16 cohorts with \(n_1\ge2\).
Thus a broad family is within reach of replay, but the one unresolved cell
and the global union budget still matter. No new fitting followed this run.

For a more useful full-family proposal, substitute the existing
homogeneous directed ceiling −50.891481584655 and the separately selected
\((1,16,0)\) FLOAT screen −70.35414804665. Keep the raw shared-bank result
for every \(n_1\ge2\). Their preliminary log2 total is
−35.473621162250815. Only the homogeneous subterm is directed at this
checkpoint; this mixed calculation is not a certificate and still fails
the 40-bit target.

The largest selected cells are listed below. Every listed cell uses the
second fixed witness above, namely
\(z=4487948476156805/4503599627370496\) and
\(\log_2v_1=-8814329/131072\). In particular, improving only the failing
\((2,15,0)\) cell would leave a narrow global probability budget.

| Cohort \((n_1,n_2,n_H)\) | Source \((a,b)\) | Multiplicity bin | Log2 screen |
|---|---|---|---:|
| (2,15,0) | (0,7) | 15 | −35.5189261621 |
| (2,15,0) | (1,7) | 15 | −41.5019554906 |
| (2,15,0) | (2,7) | 15 | −48.8658092788 |
| (2,15,0) | (0,8) | 15 | −49.9309832054 |
| (2,15,0) | (1,8) | 15 | −55.9252261478 |
| (3,14,0) | (0,7) | 9–14 | −41.5143803716 |
| (3,14,0) | (1,7) | 9–14 | −47.4973343910 |
| (3,14,0) | (2,7) | 9–14 | −54.8619776975 |
| (3,14,0) | (0,8) | 9–14 | −55.9261726313 |
| (3,14,0) | (0,6) | 9–14 | −55.9352371593 |

The first five rows are the largest cells of \((2,15,0)\); the last five
are the largest across all 15 cohorts with \(n_1\ge3\).

Eight cheap tests check counts, scope, witness selection, cap reuse, and
the no-pair boundary. Independent artifact recombination checks every
selected minimum and the full class union. Both the root review and a
separate read-only review found no defect in the original shared-cap
interface. The full-source extension also passed the root review.

There is an exact further implementation saving for a later replay. If a
final multiplicity bin is truncated to \(n_2<U\), where \(U\) is its
canonical upper endpoint, every used cancellation count satisfies
\(k\le d\le n_2\). Therefore both \(\min(k,\mu)\) and
\(\min(k,7\mu)\) agree at \(\mu=n_2\) and \(\mu=U\). All resource
bounds, hence all used cap entries, agree. Five canonical cap evaluations
can replace 17 truncated evaluations per bank. The frozen FLOAT controls
still perform all 17; the identity has a separate exact integer test.

| Artifact in `theory` | SHA-256 |
|---|---|
| `hd_xor_sparse_cohorts.py` | `aad48d68c60c2cad54504b9039f73c6ebeac5644ab3f5fa67547058e814e77ef` |
| `hd_xor_sparse_cohorts_tests.py` | `08c7e9c5e425d6a5a158de8710432b2e6d6f4e6de56bfb764f06814c8ce73c45` |
| `hd_xor_sparse_cohorts17_8_0_two_banks_float.json` | `84ff09adfddbdc619c2bddceb4dd3361292d9cb6c17eba54ea308ecada7533ec` |
| `hd_xor_sparse_cohorts17_all_sources_two_banks_float.json` | `ff786dc9a2c4a2791b759e7b244eab09206cc0213df43f7100eb71795def63ab` |

Reproduce with fresh output paths:

```powershell
C:/Python314/python.exe -B -m unittest hd_xor_sparse_cohorts_tests
C:/Python314/python.exe -B hd_xor_sparse_cohorts.py --output sparse_one_source_new.json
C:/Python314/python.exe -B hd_xor_sparse_cohorts.py --all-sources --output sparse_all_sources_new.json
```

The next bounded step is to improve the remaining \((2,15,0)\) cell and
the weakest \((3,14,0)\) cells, leaving budget for the whole-code union.
Existing literals should be checked before another fit. A successful
proposal can then share cap banks during directed replay of the complete
sparse-cohort partition. The dense-cohort obligation is unchanged.
