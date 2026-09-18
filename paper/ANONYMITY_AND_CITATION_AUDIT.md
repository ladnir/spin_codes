# Submission anonymity and bibliography audit

Checked 2026-09-17. This is an author-side audit, not material to include in
the anonymous supplement. No manuscript, bibliography entry, license, or ZIP
was changed by this audit.

## Follow-up: reviewer package prepared

The initial ZIP findings below refer to `spin-core-supplement.zip`, not the
subsequently cleaned `output/artifact/spin-core-review.zip`. Use the latter for
review. It contains 37 files (137,851 compressed bytes), with SHA-256
`ddc66d55d7ee4ef672bd114fc40feae865fb32845a56ab98f1a7bfe0b2cd9dd2`.

The reviewer package omits the source map and identifying project paths. The
rights holder explicitly authorized an anonymous copyright notice for the
forward encoder; the MIT license terms remain unchanged. Encoder source bytes
and third-party notices are preserved. The packager scans member names and
text for identifying strings and normalizes ZIP timestamps. A fresh package
reproduced the hash above. The extracted package previously passed all eight
correctness tests in a clean Linux build.

This removes the direct identifiers found in the initial archive; it does not
prevent matching byte-identical code against public sources. The PDF findings
below describe the exact earlier PDF hash, not every subsequent manuscript build.

## Initial outcome

- **Do not upload `spin-core-supplement.zip` as an anonymous supplement.** Its
  source map and license identify the Hypercat project.
- No direct author identity was found outside the bibliography in the anonymous
  submission PDF. Its author metadata is empty. The named-author draft is not
  a submission artifact.
- **All 26 bibliography entries correspond to independently traceable scholarly
  records. No fabricated reference was found.** Every entry is typeset; no
  duplicate keys or missing bibliography entries were found.
- This checks publication identity and bibliographic metadata, not whether
  every surrounding technical claim follows from its cited paper.

## Exact files checked

| File | SHA-256 |
|---|---|
| `output/artifact/spin-core-supplement.zip` | `0b6b78ef677eee4e7d4bd73267d160ce3a798765068c9c936e689c6f1bea8c40` |
| `output/pdf/spin_codes_submission.pdf` | `427db95ddae550d28b34a05101cd4a76272a864045965487e6cf10c9f808891a` |
| `paper/references.bib` | `4be19a6edc405dff4a6abf2fb2fda8a88a8241f7e4803f9cb8a6ca5dfd3b20e4` |

## Anonymity findings

### Original ZIP: not suitable for anonymous upload

All 38 archive members were inspected as text, together with member names,
archive comments, extra fields, and file types.

1. `spin-core/SOURCE_MAP.txt`, lines 6--15, names the original
   `hypercat/native/spin/` paths. This provenance makes the imported source
   project explicit.
2. `spin-core/licenses/bidirectional-MIT.txt`, line 3, identifies
   `Hypercat contributors`. The notice must not be silently stripped from an
   imported MIT-licensed source distribution. Resolve the attribution with
   the rights holders, or omit that implementation from the reviewer archive.
3. The remaining original research paths and byte-identical source files can
   assist matching against public source. Removing names alone cannot guarantee
   unlinkability to previously public code.
4. `README.md` already says the archive is not claimed to be anonymized. That
   warning accurately describes the current package; it does not make the
   package appropriate for anonymous review.

No personal name, email, home-directory path, Git remote URL, credential marker,
private-key marker, or configured remote-host name was found in the ZIP scans.
No Git metadata, binaries, build trees, logs, experiment files, or PDFs occur
in the archive. There are no ZIP comments, extra fields, or unsafe member paths.

Recommended next step: create a distinct reviewer ZIP without original-tree
provenance, resolve the identifying license notice without silently removing
required attribution, and repeat the scan on the final ZIP. Keep this full
provenance archive separate. No existing archive was deleted or overwritten.

### Anonymous submission PDF: no direct identity found

The complete text of all 78 pages was scanned for author names, affiliations,
project identifiers, personal paths, email addresses, and author repository links.
All author-name matches occur in the ordinary bibliography on pages 25--26.
Those are third-person references to real published work, not an author block.

The title page says `Anonymous submission`. `/Author` is empty; there is no
XMP metadata or embedded file. Outlines contain no tested author identifiers.
The open action is an internal page destination, not a script or external URL.
All external link annotations point to bibliography sources: DOI, IACR ePrint,
or J-STAGE. No author-identifying artifact link is present.

Metadata retains creation/modification timestamps with UTC-07:00 and standard
TeX version strings. These are environment fingerprints, not direct identity;
they could be suppressed in a final clean build. Source and phrase matching
against public work remains possible even without metadata.

The title page and both bibliography pages were also rendered and visually
inspected. The named-author `spin_codes_draft.pdf` must not be uploaded in place
of `spin_codes_submission.pdf`.

The [Eurocrypt 2027 submission instructions](https://eurocrypt.iacr.org/2027/papersubmission.php)
require anonymous submissions without author names, affiliations, or obvious
references. They permit code in a separate supplementary file. The live page
was read after removing commented-out HTML; this audit does not infer that
source-code supplements are exempt from anonymity.

## Bibliography verification

For DOI entries, publisher-deposited Crossref metadata was compared with the
local title, author list, publication venue/year, and pages where available.
DBLP, publisher, and institutional records supplied additional checks. Recent
preprints were checked directly against their IACR or arXiv landing pages.
The table below covers every entry in the typeset bibliography, in its order.

| Ref. | Local key / work | Independent record and result |
|---|---|---|
| 1 | `bms09`: The Minimum Distance of Turbo-Like Codes | [Author's Yale publication list](https://www.cs.yale.edu/homes/spielman/Research/topicList.html), [DBLP author record](https://dblp.org/pid/90/3542), and [DOI-indexed bibliographic record](https://pascal-francis.inist.fr/vibad/index.php?action=getRecordDetail&idt=20985312): Bazzi, Mahdian, Spielman; IEEE TIT 55(1), 6--15, 2009; DOI matches. |
| 2 | `benedetto98`: Serial Concatenation of Interleaved Codes | [Publisher DOI](https://doi.org/10.1109/18.669119): all four authors, title, IEEE TIT 44(3), 909--926, 1998 match. |
| 3 | `berrou93`: Near Shannon Limit Error-Correcting Coding and Decoding | [Publisher DOI](https://doi.org/10.1109/ICC.1993.397441): three authors and ICC pages 1064--1070 match. The indexed title ends `Turbo-codes. 1`; the local title omits the part suffix. This is a minor title-normalization discrepancy, not a different or invented work. |
| 4 | `fieldEA24`: Field-Agnostic SNARKs from Expand-Accumulate Codes | [Publisher DOI](https://doi.org/10.1007/978-3-031-68403-6_9), [Illinois institutional record](https://experts.illinois.edu/en/publications/field-agnostic-snarks-fromexpand-accumulate-codes/): six authors, CRYPTO 2024, LNCS 14929, 276--307 match. |
| 5 | `bcgi18`: Compressing Vector OLE | [DBLP](https://dblp.org/rec/conf/ccs/BoyleCGI18.html): four authors, CCS 2018, 896--912, and DOI `10.1145/3243734.3243868` match. |
| 6 | `ea22`: Correlated Pseudorandomness from Expand-Accumulate Codes | [Publisher DOI](https://doi.org/10.1007/978-3-031-15979-4_21), [Aarhus institutional record](https://pure.au.dk/portal/en/publications/correlatedpseudorandomnessfrom-expand-accumulate-codes/): seven authors, CRYPTO 2022, LNCS 13508, 603--633 match. |
| 7 | `bcg19tworound`: Efficient Two-Round OT Extension and Silent Non-Interactive Secure Computation | [Publisher DOI](https://doi.org/10.1145/3319535.3354255), [IACR](https://eprint.iacr.org/2019/1159): seven authors, CCS 2019, 291--308 match. |
| 8 | `bcg19pcg`: Efficient Pseudorandom Correlation Generators: Silent OT Extension and More | [Publisher DOI](https://doi.org/10.1007/978-3-030-26954-8_16), [DBLP](https://dblp.org/rec/conf/crypto/BoyleCGIKS19.html): six authors, CRYPTO 2019 Part III, 489--518 match. |
| 9 | `bcg20ring`: Efficient Pseudorandom Correlation Generators from Ring-LPN | [Publisher DOI](https://doi.org/10.1007/978-3-030-56880-1_14), [IACR full version](https://eprint.iacr.org/2022/1035): six authors and CRYPTO 2020 match. The later full-version year is not an incorrect link to another paper. |
| 10 | `blaze25`: Blaze: Fast SNARKs from Interleaved RAA Codes | [Publisher DOI](https://doi.org/10.1007/978-3-031-91134-7_5), [IACR](https://eprint.iacr.org/2024/1609): six authors, EUROCRYPT 2025, 123--152 match. |
| 11 | `flock26`: Flock: Fast Proving for Batch Boolean Computations | [IACR 2026/1329](https://eprint.iacr.org/2026/1329): Bunz, Rothblum, Wang; title, year, and report number match. |
| 12 | `CouteauRindalRaghuraman2021Silver`: Silver | [Publisher DOI](https://doi.org/10.1007/978-3-030-84252-9_17), [IACR](https://eprint.iacr.org/2021/1150): three authors, full title, CRYPTO 2021, 502--534 match. |
| 13 | `diamondPosen24`: Proximity Testing with Logarithmic Randomness | [IACR journal record](https://cic.iacr.org/p/1/1/2), [DOI](https://doi.org/10.62056/aksdkp10): Diamond and Posen, volume 1, issue 1, 2024 match. |
| 14 | `divsalar98`: Coding Theorems for Turbo-Like Codes | [JST J-GLOBAL](https://jglobal.jst.go.jp/en/detail?JGLOBAL_ID=200902100021468289): Divsalar, Jin, McEliece; 36th Allerton, 201--210, 1998 match. This supplies an independent bibliographic record despite no DBLP record being located. |
| 15 | `FujiwaraKusaka2021`: Weight Distributions of the (256,k) Extended Binary Primitive BCH Codes | [Publisher DOI](https://doi.org/10.1587/transfun.2020EAP1119): both authors, title ranges, volume E104.A, issue 9, 1321--1328, 2021 match. The local abbreviated journal name and E104-A spelling identify the same journal. |
| 16 | `brakedown23`: Brakedown | [Publisher DOI](https://doi.org/10.1007/978-3-031-38545-2_7), [JST J-GLOBAL](https://jglobal.jst.go.jp/en/detail?JGLOBAL_ID=202302255702383901): five authors, full title, CRYPTO 2023, LNCS 14082, 193--226 match. |
| 17 | `bolt26`: Bolt: Faster SNARKs from Sketched Codes | [IACR 2026/310](https://eprint.iacr.org/2026/310): Gurkan, Novakovic, Rothblum; title and report match. IACR now also identifies a CRYPTO 2026 publication; citing the actual preprint is not fabrication. |
| 18 | `kahale98`: On the Minimum Distance of Parallel and Serially Concatenated Codes | [Publisher DOI](https://doi.org/10.1109/ISIT.1998.708611): Kahale and Urbanke, ISIT 1998, page 31 match. |
| 19 | `kliewer08`: Coding Theorems for Repeat Multiple Accumulate Codes | [arXiv 0810.3422](https://arxiv.org/abs/0810.3422): Kliewer, Zigangirov, Koller, Costello Jr.; title and 2008 date match. |
| 20 | `KolesnikovEtAl2026BA`: Block-Accumulate Codes | [Publisher DOI](https://doi.org/10.1007/978-3-032-35418-1_13), [DBLP proceedings](https://dblp.dagstuhl.de/db/conf/crypto/crypto2026-8.html), [IACR](https://eprint.iacr.org/2025/1828): six authors, full title, CRYPTO 2026 Part VIII, LNCS 16807, 397--429 match. |
| 21 | `ligerito25`: Ligerito | [IACR 2025/1187](https://eprint.iacr.org/2025/1187): Novakovic and Angeris; full title, year, and report match. |
| 22 | `cryptoeprint:2026/1903`: Chosen-Block BAA Codes | [IACR 2026/1903](https://eprint.iacr.org/2026/1903): Peceny and Rindal; full title, year, and report match. The public record was approved September 10, 2026. |
| 23 | `RaghuramanRindalTanguy2023EC`: Expand-Convolute Codes | [Publisher DOI](https://doi.org/10.1007/978-3-031-38551-3_19), [IACR](https://eprint.iacr.org/2023/882): three authors, full title, CRYPTO 2023, 602--632 match. |
| 24 | `rr25`: Proving as Fast as Computing | [DBLP](https://dblp.org/rec/journals/jacm/RonZewiR25), [Publisher DOI](https://doi.org/10.1145/3721477): Ron-Zewi and Rothblum, full title, JACM 72(2), article 15, 1--54, 2025 match. DBLP/publisher omit Rothblum's middle initial used locally; the author identity is unchanged. |
| 25 | `orion22`: Orion | [Publisher DOI](https://doi.org/10.1007/978-3-031-15985-5_11), [IACR](https://eprint.iacr.org/2022/1010): Xie, Zhang, Song; full title, CRYPTO 2022, 299--328 match. |
| 26 | `ferret20`: Ferret | [Publisher DOI](https://doi.org/10.1145/3372297.3417276), [IACR](https://eprint.iacr.org/2020/924): five authors, full title, CCS 2020 match. The local entry omits pages; publisher metadata supplies 1607--1626. |

Crossref returned live metadata for 18 of the 20 DOI entries. The two rate-limited
entries were checked through the independent records listed above. Some direct
DBLP downloads returned bot challenges or timeouts; those responses were not
counted as bibliographic verification. Search-indexed DBLP records, publisher
records, and IACR/JST records supplied the stated checks instead.

The only identified bibliography cleanup is minor normalization/completeness:
the Berrou title suffix, optional author-name normalization for Rothblum, optional
Ferret pages, and an optional arXiv URL for the Kliewer entry. None is evidence
of a fabricated citation, wrong DOI target, or unrelated preprint substitution.
