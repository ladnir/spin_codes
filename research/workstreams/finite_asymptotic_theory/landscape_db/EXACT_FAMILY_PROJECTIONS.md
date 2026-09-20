# Exact-spectrum anchors and engineering projections

The finite grid stops at BCH length 128 and RM length 512. This report
provides a separate estimate of Q1 scaling beyond those lengths. It uses
only the complete exact spectra in the primary grid. Random outers are
computed directly and are not extrapolated.

For each fixed (t,s,log2(k)), let M_B be the observed Q1 margin, d_B the
outer minimum distance, A_B the number of minimum-weight words, and
L_B=k/(B/2). Form the adjusted response

    F_B = M_B + log2(A_B) + log2(L_B).

The adjustment removes the minimum-shell multiplicity and the row-position
count from the observed margin. A linear fit of F_B against d_B then
models the remaining transfer contribution. This is an empirical model;
it does not replace the full Q1 weight sum by an equality. Its input rows
must have the minimum weight as their dominant shell and an interior tilt
witness. All other rows remain in the finite grid but are excluded from
this fit.

For each prediction, fit windows containing the largest two, three, and
four available exact constituents. Restore the projected shell and row
counting costs after fitting F_B. The report keeps the median prediction
and the full window range. That range measures model sensitivity and is
not a confidence interval.

For rate-half RM(r,2r+1), the projected minimum distance and minimum-shell
count are exact:

    B = 2^(2r+1), d_B = 2^(r+1),
    A_B = 2^r [2r+1 choose r]_2.

Here the Gaussian binomial coefficient counts r-dimensional linear
subspaces of F_2^(2r+1). The displayed count enumerates affine subspaces
of codimension r. It agrees with
all four exact RM spectra in the grid and with the minimum-shell case of
`../rm511_outer/verify_rm511_low_weight_formulas.py`. The remaining spectrum
at length 2048 is still unavailable, so that code remains outside the
primary exact grid.

For BCH, both future minimum-shell quantities are estimates: fit log2(d_B)
linearly against log2(B), then log2(A_B) linearly against d_B, using the
same window. These are continuous proxies, not asserted parameters of a
specified future BCH constituent. BCH-256 bounded-spectrum evidence is
excluded from every fit.

The largest-known-size check withholds BCH-128 or RM-512 and predicts its
Q1 margin from two or three smaller exact constituents. The initial
checks have median absolute errors of about 13.7 bits for BCH and 7.7 bits
for RM, with maximum errors of about 19.0 and 14.1 bits. The median errors
are negative: the model predicts weaker margins in these checks. This
bias does not guarantee conservative predictions at a new size.

These errors are large enough that a forecast lying only a few bits above
a target cannot resolve parameter selection. In particular, Q1 forecasts
do not establish full-occupation bounds. Implementation choices should
retain that uncertainty and use full evaluated results where available.

Run `python extrapolate_exact_families.py` to reproduce the training rows,
held-out comparisons, projections, and their hashes. The outputs cover
BCH lengths 256 and 512 and RM length 2048, at the existing message lengths
and map choices where the projected geometry is native. The older
`parameter_extrapolation.json` remains historical and under review.

## Secondary BCH-256 comparison

The owning BCH worktree's `T64_S20_FULL_CLOSURE.md` reports a full outward
bound below 2^-50 at k=2^20, t=64, s=20, with decimal diagnostic margin
50.439 bits. The corresponding primary Q1 projection has median 58.411
bits and window range 44.016..72.958 bits. This is a qualitative cross-check
of the proposed scale. Quantitative calibration would require matching
the map, cutoff, and occupation coverage: the anchor uses a different
selected map, a cutoff one bit larger, and a full-occupation bound.

`secondary_bch256_comparison.json` records the source statement's hash and
these distinctions. This worktree has not independently replayed that
outward certificate. The anchor remains excluded from every fit and from
the primary exact-spectrum grid.
