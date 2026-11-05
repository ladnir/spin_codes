# One dense XOR2 group: exact mask mixture and column floor

The fixed-mask comparison closes the class with one active outer group and
three through eight active constituent rows. At the inherited rational mark
`499451537/500000000`, fresh 256- and 384-bit replays give

\[
\log_2 P_{r=1,\,3\le a\le8}
\le -82398007651003/500000000000
=-164.796015302006.
\]

This is a class bound, not a whole-code certificate. The construction is
`xor2-gf256-h4-alternating-swapped-pairs-k20-rate-half`: the XOR2 outer map
is \(D_2(I+\mathbf1\mathbf1^T)D_1\), followed by the alternating paired route
with independent fair swaps and the existing rank-four inner construction.
The probability is over the ideal independent setup ensemble, not a claim
about a particular pseudorandom seed. The forbidden output weight is
\(L=734003\). There are 8192 possible positions for the active group.

## Counting measure and the comparison law

Fix the set \(A\subseteq\{0,\ldots,7\}\) of active RS rows and write
\(a=|A|\). Each row is an RS\([32,16]\) word over \(\mathbb F_{256}\).
We count messages with exactly these nonzero rows; \(A\) remains fixed
throughout all 32 columns.

For a prescribed nonzero row support \(T\) of size \(w\), let \(B(T)\)
count codewords with exactly that support. The shortened MDS code has
dimension \(w-16\) when such a nonzero word exists. An information set
injects its full-support words into \((\mathbb F_{256}^{*})^{w-16}\), so

\[
B(T)\le255^{w-16}.
\]

For a prescribed scaled row with support \(T\), independent nonzero
entries of \(D_1\) contribute probability \(255^{-w}\) per original word.
Therefore its expected counting mass is at most \(255^{-16}\).
For all \(a\) active rows, the mass of any prescribed scaled array is at
most \(255^{-16a}\). Invalid supports have actual mass zero.

Compare this measure with 32 independent uniform columns in
\(\mathbb F_{256}^{A}\), with zero coordinates outside \(A\).
Each array has comparison mass \(256^{-32a}\). Thus the actual counting
measure is pointwise dominated by

\[
\beta_a\nu_A^{\otimes32},\qquad
\beta_a=\frac{256^{32a}}{255^{16a}}.
\]

This statement concerns the sum over messages. It is not a likelihood
bound for one fixed codeword. Applying the XOR map and the independent
output diagonal \(D_2\) preserves the domination.

For an output mask \(T\subseteq\{0,\ldots,7\}\), let \(N_A(T)\) count
uniform input columns whose XOR image has support \(T\). The implementation
uses exact integer counts from `uniform_column_support_counts`:

\[
\sum_TN_A(T)=256^a,\qquad N_A(\varnothing)=1.
\]

Conditional on \(T\), \(D_2\) supplies independent nonzero amplitudes.
The existing forced-nonzero regional operator therefore applies directly.
No additional factor \(255^{-|T|}\) is needed.

## Keep the nonempty-column floor and matrix order

Each nonzero RS row has at least 17 nonzero columns. The outer map is
invertible within each column. Hence every actual active group has at
least 17 nonempty eight-coordinate output columns. This restriction
counts columns, not individual nonzero coordinates.

Let \(R_0(z),R_1(z)\) be the directed two-state regional transfer bounds
from `arb_pair_regions(z, 1)`. A single occupied source belongs to no full
source pair, for either matching. For each labelled output mask, define
the chronological product

\[
Q_T(z)=\prod_{i=0}^{7}R_{\mathbf1[i\in T]}(z).
\]

For the fixed row set \(A\), separate the empty and nonempty columns:

\[
Q_{A,0}=256^{-a}Q_\varnothing,\qquad
Q_{A,1}=256^{-a}\sum_{T\ne\varnothing}N_A(T)Q_T.
\]

All entries are nonnegative upper transfer bounds. The positive matrix
polynomial \((Q_{A,0}+YQ_{A,1})^{32}\) retains their order and the inner
state across columns. With \(e_0=(1,0)\) and \(\mathbf1=(1,1)^T\), set

\[
M_A(z)=e_0\left(\sum_{j=17}^{32}[Y^j]
 (Q_{A,0}+YQ_{A,1})^{32}\right)\mathbf1.
\]

Markov's inequality and one projective normalization give the implemented
bound

\[
P_{r=1,\,3\le a\le8}
\le\frac{8192\,z^{-L}}{255}
 \sum_{A:\,3\le|A|\le8}\beta_{|A|}M_A(z).
\]

The sum explicitly visits all 219 row sets. There is no extra binomial
factor. In particular, the evaluator does not average \(A\) before taking
the 32-column power. It uses neither a product-density comparison
\(c_a(p)^{32}\) nor an endpoint potential.

## Results and checks

Only the three existing rational marks were screened; no new mark was fit.
The FLOAT unions were +5189.1814 at `62078423/62500000`,
−164.7960 at `499451537/500000000`, and −106.1368 at `1199/1200`.
The first result is a failed bound, not a construction counterexample.

The selected mark's directed arity ceilings are:

| Active rows | Fixed row sets | Log2 upper bound |
|---:|---:|---:|
| 3 | 56 | −174.216638059926 |
| 4 | 70 | −172.760843369442 |
| 5 | 56 | −172.209760720413 |
| 6 | 28 | −172.486186664386 |
| 7 | 8 | −166.016578541808 |
| 8 | 1 | −165.646193146304 |

All 219 per-set ceilings, all six arity ceilings, and the union ceiling
agree at 256 and 384 bits. Arity eight is the largest contribution.
The largest retained marked-mass degree is 17 for every row set at this
mark; this describes the bound's terms, not an actual codeword distribution.
Each directed replay took approximately 1.7 seconds, excluding interpreter
startup. This was proof computation, not a native encoder benchmark.

Nine focused tests pass. They check the comparison law exhaustively over
\(\mathbb F_4\), shortened-RS normalization on toy supports, fixed-set
ordering, and probability normalization. Directed tests additionally
compare noncommuting rational matrix products and polynomial coefficients
with exact word enumeration, remove the zero-column contribution, and
compare all 219 saved ceilings across precisions. Independent read-only
review confirmed the counting-measure normalization and the column floor.

From the `theory` directory, reproduce into fresh output paths:

```powershell
C:/Python314/python.exe -m unittest hd_xor_fixed_mask_floor_tests hd_xor_fixed_mask_floor_certificate_tests
C:/Python314/python.exe hd_xor_fixed_mask_floor_certificate.py --precision 256 --output fixed_mask_replay256.json
C:/Python314/python.exe hd_xor_fixed_mask_floor_certificate.py --precision 384 --output fixed_mask_replay384.json
```

SHA-256 pins for the saved sources and evidence:

| File in `theory` | SHA-256 |
|---|---|
| `hd_xor_fixed_mask_floor.py` | `2dd5af28ea33015e759f43a8b7e18436a5b506a26015917c236c5e1262499044` |
| `hd_xor_fixed_mask_floor_certificate.py` | `1e06c18d0dfd2c4bc12c16723863ff697e7fcaa6b07115d3d9c2c4265ba14fdf` |
| `hd_xor_fixed_mask_floor_tests.py` | `34324454a477baf946ad01f5a49530dea5e6e3b757317acd07db47928148b110` |
| `hd_xor_fixed_mask_floor_certificate_tests.py` | `80eb64b0add9ff2daae9943816093f8858ed446679c09f4cfa30426a638007fb` |
| `hd_xor_fixed_mask_floor_r1_float.json` | `0a9d55333df1571a38099e6be02cb89165f1717402a126890dfe7425d12878b7` |
| `hd_xor_fixed_mask_floor_r1_receipt256.json` | `46d9754583dd429965501ab534f7be3c62ddf21a171561365d41f67961b4ead1` |
| `hd_xor_fixed_mask_floor_r1_receipt384.json` | `a03c455ee328186dec4152b0cb52af9c4d26f95386020f7ca619d31d52696dec` |

The remaining extension is a mixed group class, not another replay of
this one-group calculation. One dense group and sixteen two-row groups
require a separate sparse-support degree and a nonempty-dense-column
degree, while retaining each dense row set and sparse pair geometry.
