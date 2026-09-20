# Sparse Expander--Accumulate Lane Freeze

Status: frozen on 2026-09-03.  The files named in
`SPARSE_EA_FREEZE_MANIFEST.json` remain in their original locations so that
their relative-path certificate bindings remain valid.

## Closed result

Let one outer constituent be sampled by choosing 512 independent uniform
weight-33 rows in \(\mathbb F_2^{256}\), applying a zero-initialized length-512
accumulator, and accepting the first full-rank constituent in at most 16
independent attempts.  The accepted constituent is reused in all 4096 outer
rows.  The routing independently samples one uniform 512-coordinate
permutation for each outer row and one uniform 4096-position permutation for
each transposed region.  RandomStepConv-M22 independently samples a uniform
binary \(23\)-by-\(23\) linear map at each of the \(2^{21}\) inner positions;
all sampled objects are then fixed as the code.

For

\[
  k=2^{20},\qquad N=2^{21},\qquad D=228{,}590,
\]

the hash-bound outward certificate proves

\[
  \Pr[d_{\min}<D\ \text{or outer setup aborts}]
  <2^{-41.3621359294}.
\]

Thus this ensemble has minimum relative distance at least
\(D/N=10.9000205994\%\), except with the displayed probability.  This is a
proved statement for the stated sparse-EA outer, uniform routing, and
RandomStepConv-M22.  It is not a theorem about the frozen Structured SPIN
interleaver, RM2Sub, or a fixed preselected sparse matrix.

The complete theorem is in
`ONE_STAGE_SPARSE_EA_RANDOMSTEPCONV_CERTIFICATE.md`.  The internal receipt
chain is bound by `ONE_STAGE_SPARSE_EA_CERTIFICATE_MANIFEST.json` and checked
by `audit_one_stage_sparse_ea_certificate.py`.

## Performance result

A matched outer-only AVX2 benchmark compares the sampled one-stage
degree-33 transform \(E^{\mathsf T}A^{\mathsf T}\) with two optimized
ExtendedBch256x128-Eq3 transposes.  Both map \(2^{21}\) packed input blocks to
\(2^{20}\) packed output blocks.  Four 31-sample processes ran sequentially
on Peach CPU 0 in sparse, BCH, BCH, sparse order.

The sparse path took 8.703 ms and the BCH path took 4.896 ms, using the mean
of the two medians for each path.  The sparse path was 1.777 times slower and
added 3.807 ms.  Its transposed circuit used 9,733 XORs and had 1,305 live
packed values, compared with 5,868 XORs and 448 live values for BCH.  This is
an exact implementation comparison for the generated circuits, but it is an
outer-only benchmark and not an end-to-end encoder measurement.

## Open depth-two route

The attempted lower-XOR outer was

\[
  C_2=A E_1 A E_0.
\]

Several degree pairs have mean spectra below every cap used by the closed
one-stage transfer.  That observation is diagnostic, not a simultaneous
high-probability spectrum theorem.  The shared first map contributes the
variance of a conditional mean.  Conditioning the square sparse mixer on
invertibility also couples its rows.  No proved conditional contraction or
signed-defect lemma currently controls that term.

Exact common-subexpression synthesis further weakens the performance case:
the sampled degree pair \((17,7)\) uses 9,556 transposed XORs, only 177 fewer
than the certified degree-33 sample.  The depth-two route is therefore both
unproved and unlikely to recover the measured BCH gap without a new idea.

## Preservation and resumption rules

`SPARSE_EA_FREEZE_MANIFEST.json` hashes the theorem, audit, cap and transfer
receipts, proof-route records, depth-two diagnostics, generated circuits, and
benchmark receipts.  `audit_sparse_ea_freeze.py` checks every binding and
then runs the semantic certificate audit.  The live attempt log and merge
summary are intentionally excluded because subsequent work updates them.

Resume this lane only for one of four explicit goals: a materially cheaper
fixed circuit, a new depth-two concentration lemma, an RM2Sub transfer, or a
proof for structured routing.  None is needed for the RM(4,9) investigation.
