# Exact-spectrum candidate for quarter-rate SPIN

Checkpoint scope: code and documentation only. The new spectrum tables,
certificate JSON, and measurement records referenced below remain local,
not committed. Data-dependent reproduction requires those inputs or their
regeneration; the reports retain the conclusions of the completed local runs.

Status, 2026-09-11: the fixed BCH-derived **[256,64,62]** generator and its
complete weight spectrum are now available. There are **512 words of weight
62**. The spectrum has 69 nonzero coefficients and total mass 2^64. It follows
from two published spectra and a verified transitive coset action, not sampling.
No outward-certified SPIN margin is claimed for this [256,64] outer.
The smaller outer below now has separate outward certificates and an optimized implementation.

Follow-up: [the fixed-inner evaluation](FIXED_INNER_EVALUATION.md) now gives
full numerical occupation coverage at K=2^16, 2^18, and 2^20 for distance
targets 0.10 and 0.19. It uses the existing optimized t128_s19 inner without
re-optimization. These numerical results are not outward certificates.

[The smaller [128,32] outer](SMALLER_OUTER.md) now has
[outward certificates](SMALLER_OUTWARD_CERTIFICATE.md) at K=2^20: 16.5%
distance with failure below 2^-40, and 19% with failure below 2^-30, using
the same t128_s19 inner. Both passed 512-bit replay; margin diagnostics are
41.08348 and 30.05251 bits.

The [quarter-rate transposed implementation](implementation/README.md) reuses the
existing optimized inner and shared routing, with a new paired AVX2 outer circuit.
See its [performance report](implementation/PERFORMANCE.md) for in-place timing,
memory, correctness checks, and reproduction commands.

## Sources

- [Fujiwara/Kasami's [256,63] table](https://isec.ec.okayama-u.ac.jp/home/kusaka/wd/EBCH/fujiwara/EBCH256_63.wd)
  is retained as `sources/EBCH256_63.wd`. It has mass 2^63, minimum weight
  64, 43,180 words of weight 64, and complement symmetry. The exact
  MacWilliams transform has nonnegative integer coefficients and mass 2^193.
  The [authors' index](https://isec.ec.okayama-u.ac.jp/home/kusaka/wd/)
  attributes the table to their earlier technical report.
- Fujiwara and Kusaka, *The Weight Distributions of the (256, k) Extended
  Binary Primitive BCH Codes with k <= 71 and k >= 187*, IEICE Transactions
  E104-A(9), 1321-1328 (2021), [DOI 10.1587/transfun.2020EAP1119](https://doi.org/10.1587/transfun.2020EAP1119).
  **Table 7, printed page 1326 (PDF page 6)** supplies the [256,71] spectrum,
  retained as `sources/EBCH256_71.wd`. Its actual minimum distance is 62,
  with 130,560 minimum-weight words. The table labels its entries A_w, A_(256-w);
  the import expands that symmetry without doubling the central coefficient.
  Its [256,63] column independently matches the retained public table.
- Earlier source: Fujiwara and Kusaka, IT2019-5 / EMM2019-5, pp. 23-28,
  May 2019. Its [official record](https://www.ieice.org/publications/ken/summary.php?contribution_id=102023&ken_id=IT&lang=en&presen_date=2019-05-23&schedule_id=6408&society_cd=ESSNLS&year=2019)
  also reports computing the [256,71] spectrum. We have not obtained its PDF.

The user supplied the final journal PDF. `sources/PAPER_TABLE_PROVENANCE.json`
records its SHA-256, page, table, and extracted coefficient-file hash. The PDF
is not redistributed: the artifact retains numerical data and provenance only.

## Concrete construction

Let P and Q be the binary primitive narrow-sense BCH codes of length 255
with designed distances 59 and 61. Their dimensions are 71 and 63. Use the
existing repository field GF(256), modulus 0x14d and primitive element 2.
Their generator polynomials p and q satisfy

    q(X) = p(X) h(X),   h(X) = X^8 + X^7 + X^2 + X + 1 (0x187).

Let C71 and C63 denote their even-parity extensions to length 256. Define
C64 as the span of C63 and the parity extension of the coefficient vector p.
The retained generator consists of the 63 extended rows q, Xq, ..., X^62 q,
followed by the extended row p. The parity coordinate is bit 255; bit i
of the remaining word is the coefficient of X^i.

`reconstruct.py` verifies rank 64, containment in C71, the all-one word,
generator divisibility, and the cyclic quotient action. This proves a
deterministic [256,64,d >= 60] constituent without using the published spectrum:
P has distance at least 59, and parity extension raises that lower bound to 60.
The reconstructed spectrum sharpens this to **d = 62**. This is a BCH-derived
intermediate code, not an assertion that a standard narrow-sense BCH code
has dimension exactly 64.

## Why two published spectra determine the intermediate spectrum

The following is our derivation, not a claim that the paper tabulates C64.
Write W_D(z) = sum over c in D of z^wt(c).

The quotient P/Q identifies with GF(2)[X]/(h): map a(X)p(X) to a(X) modulo h.
A cyclic coordinate shift acts on this quotient by multiplication by X.
The checker visits all 255 nonzero residues before returning to one. Hence
the nonzero cosets form a single orbit. Coordinate shifts preserve weight
and fix the added parity coordinate, so all 255 nonzero cosets have the same
weight enumerator, say V(z). Therefore

    W_C71 = W_C63 + 255 V,
    W_C64 = W_C63 + V = W_C63 + (W_C71 - W_C63)/255.

Table 3 of the paper also reports two coset-equivalence classes for
BCH(71)/BCH(63), agreeing with the independently checked cyclic action.

This identity is exact for the fixed generator above. It is not an ensemble
average. More generally, every intermediate subcode containing C63 and
having dimension 63+j has enumerator W_C63 + (2^j-1)V, for 0 <= j <= 8.
We make no assertion that all those higher-dimensional subcodes are
permutation-equivalent; equal enumerators suffice here.

## Artifact layout and reproduction

- `reconstruct.py`: BCH algebra, coset action, and exact spectrum checks.
- `BCH256_64.wd`: complete selected-code spectrum; omitted weights have count zero.
- `CONSTRUCTION_AUDIT.json`: generator and spectrum-independent distance bound.
- `SPECTRUM_AUDIT.json`: generator, source hashes, full spectrum, and audit results.
- `sources/`: the two published coefficient tables and PDF provenance receipt.
- `import_paper_table.py`: optional repeatable PDF extraction, requiring `pypdf`.
- `test_reconstruct.py`: fast exact regression and invalid-input checks.

Python 3.10+ and the repository's existing BCH algebra modules suffice.
These commands do not start benchmarks or enumerate exponentially many words.

```sh
python workstreams/rate_quarter_bch/reconstruct.py
python -m unittest discover -s workstreams/rate_quarter_bch -p test_reconstruct.py
```

Rebuild the complete spectrum and audit without needing the PDF:

```sh
python workstreams/rate_quarter_bch/reconstruct.py \
  --parent-spectrum workstreams/rate_quarter_bch/sources/EBCH256_71.wd \
  --parent-source 'Fujiwara-Kusaka 2021, DOI 10.1587/transfun.2020EAP1119, Table 7, p.1326' \
  --output workstreams/rate_quarter_bch/SPECTRUM_AUDIT.json \
  --spectrum-output workstreams/rate_quarter_bch/BCH256_64.wd
```

To verify extraction against a legally obtained copy of the paper:

```sh
python workstreams/rate_quarter_bch/import_paper_table.py '/path/to/paper.pdf' --check
```

The importer checks mass, parity, complement symmetry, distance support,
all MacWilliams coefficients, and nonnegative coefficientwise divisibility
of W_C71-W_C63 by 255. It then checks the reconstructed spectrum similarly.
Those checks detect transcription errors; they do not independently establish
that a supplied spectrum belongs to the stated BCH parent. Source attribution
and matching the publication's code definition remain required.

Next: implement and benchmark the complete quarter-rate encoder with the
certified t128_s19 baseline. [Inner calibration](inner_calibration/README.md)
is retained for context; t64_s16 remains numerical-only and t256 is parked.
