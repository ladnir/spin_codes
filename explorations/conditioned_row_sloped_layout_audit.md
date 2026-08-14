# Conditioned rows in the sloped three-band layout

## Question

The one-conditioned-row outer bound fixes one data row in every physical
tile. It then applies linear Brascamp--Lieb to the other 63 rows. The graph
formula additionally requires the punctured row to be the fixed row.

This argument needs one global set of fixed data blocks. Choosing a different
fixed row independently in each tile does not define a consistent set of
message variables. This note checks that requirement against the sloped
three-band layout.

## Frozen sloped layout

Index a data block by

\[
e=(t,l)\in \mathbb Z_{256}\times\{0,\ldots,63\}.
\]

Its three band cells occupy tiles

\[
\lambda_0(e)=t,\qquad
\lambda_1(e)=t+9l,\qquad
\lambda_2(e)=t+20l
\pmod {256}.
\tag{1}
\]

A retained set \(M\) supports the current 63-free-row formula only if

\[
|M\cap\lambda_b^{-1}(u)|=1
\quad
\text{for every }b\in\{0,1,2\},\ u\in\mathbb Z_{256}.
\tag{2}
\]

Thus each restriction \(\lambda_b|_M\) must be bijective. Condition (2)
makes \(M\) a perfect matching of the three-partite tile hypergraph.

For example, every fixed lane class

\[
M_L:=\{(t,L):t\in\mathbb Z_{256}\}
\tag{3}
\]

is a perfect matching. For fixed \(L\), each map
\(t\mapsto t+s_bL\) is a translation of \(\mathbb Z_{256}\).

## Failure of independent puncture lanes

Let \(P\) be the 128 punctured data blocks. The exact-graph conditioned-row
formula deletes one coordinate from each retained codeword. It therefore
requires

\[
P\subseteq M
\tag{4}
\]

for one matching \(M\) satisfying (2).

The independent per-tile lane sampler does not ensure (4). Consider the two
allowed punctures

\[
e_0=(0,0),\qquad e_1=(247,1).
\]

Their tile triples are

\[
(0,0,0)\quad\text{and}\quad(247,0,11).
\]

They occupy different band-zero tiles but the same band-one tile. No matching
can contain both edges. This seed has positive probability under the frozen
independent-lane sampler.

The defect is global. Exact split spectra, directed interval arithmetic, and
the unit-fugacity identity do not prove the missing matching. Consequently,
the existing exact-graph conditioned-row witnesses do not certify the sloped
construction with independent puncture lanes.

## Proof-only repair for the unchanged sampler

Fix the diagonal matching \(M_0\) and first bound the unpunctured data word.
Each final graph insertion replaces one data bit. One replacement changes one
packet class from \(j\) to some \(k\) with \(|j-k|\le 1\).

For positive packet fugacities \(t_j\), define

\[
R(t):=\max_{|j-k|\le1}\frac{t_k}{t_j}.
\]

If \(a_{\rm data}\) and \(a_{\rm final}\) are the packet profiles before and
after the 128 replacements, then

\[
t^{a_{\rm final}}\le R(t)^{128}t^{a_{\rm data}}.
\tag{5}
\]

Equation (5) gives a pointwise exact-length correction. In base-two affine
form, add

\[
128\max_{|j-k|\le1}(q_k-q_j),
\qquad q_j=\log_2t_j,
\tag{6}
\]

to the unpunctured one-row constant.

The bounded diagnostic in
`scripts/probe_packet_group_sloped_conditioned_row_repair.py` checks (6).
At the canonical `g=4` leader, the repair costs about 558.08 bits relative to
the inapplicable exact-graph row. At `g=8` leaf `h2:073`, it costs about
514.64 bits before applying the separate trivial outer cap.

## Construction repair

A smaller repair changes only the puncture sampler:

1. sample one global lane \(L^*\) uniformly from \(\{0,\ldots,63\}\);
2. choose 128 band-zero tiles uniformly without replacement; and
3. puncture the block in lane \(L^*\) of every chosen tile.

Then \(P\subseteq M_{L^*}\), so the existing one-row exact-graph formula
regroups into 256 complete BCH codewords. Uniform \(L^*\) preserves the
puncture marginal of every data block. The change introduces cross-tile lane
correlation, which the revised construction statement must expose.

This rule changes preprocessing only. It does not alter the encoder's hot
path or memory layout.

The same rule enables a second separated row. Choose a second global lane
\(L'\ne L^*\). An entrywise rank-one envelope separates the two retained
lane classes, after which both classes regroup into complete BCH codewords.
The parity-collapse and same-tile genus-two formulas remain invalid because
they do not perform this separation.

## Status

The prior `g=4` numerical certificate remains a valid outward evaluation of
its stated local inequalities. It applies to the global-lane puncture variant,
not to the independently punctured sloped layout. The `g=8` factorized receipt
has the same construction-binding limitation.

Before either result is called end to end, the construction must adopt the
global-lane rule or the verifier must replace every affected outer row with a
proof-only repair such as (6).
