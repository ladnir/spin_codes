# Spectrum sources for the rate-half comparison

## Acceptance rule

The comparison accepts a structured constituent only when its complete weight
enumerator is available and passes the following checks:

1. the coefficients are nonnegative integers;
2. their sum is (2^K) for a binary ([B,K]) code;
3. the unique weight-zero coefficient is one; and
4. the first nonzero coefficient occurs at the claimed minimum distance.

Complement symmetry is checked when the code is known to contain the all-one
word. It is not required for the shortened-XBCH constituent.

## Accepted spectra

| Constituent | Spectrum source | Authentication and checks |
|:---|:---|:---|
| extended BCH ([8,4,4]) | `spectra/ebch8_4_spectrum.csv` | The enumerator (1+14z^4+z^8) has mass (2^4). This is also the RM(1,3) enumerator. SHA-256 `ef6f4481ef2fb436a93d8a0d6f753b47db1b7d2e917bd0ac3111854414e35f35`. |
| extended BCH ([32,16,8]) | `../ebch32_16_delta8_spectrum.csv` | The exact enumerator supplied for the constituent has mass (2^{16}), minimum weight eight, and complement symmetry. SHA-256 `81aa68186b3f4ca0256893e7142914a4661236c1c52ac45a4b807c458a34a519`. It is also the RM(2,5) enumerator. |
| Philips shortened-XBCH ([64,32,12]) | `spectra/xbch64_32_philips_spectrum.csv` | All (2^{32}) words of the published systematic generator were enumerated. The output hash exactly matches `scripts/xbch64_32_philips_manifest.json`: `8b89591a98e585740626b1af493f6aa317861f506e5f8e87dbb3912f122c216d`. The generator-row hash is `68a60b201af9cb8f2918f553d540b2d491da8c9156054c113f748fd3dbd56eca`. The reconstruction receipt authenticates the BCH embedding, and weight 12 occurs 787 times. |
| extended BCH ([128,64,22]) | `scripts/EBCH128_64.wd` | The imported Okayama/Desaki table has mass (2^{64}), minimum weight 22, and complement symmetry. Its current file hash is `f633eb9a2f76c64c1d74313b066c2613adba642162e04b67466714f3c30c4b8a`. The existing table audit derives the linked primitive ([127,64,21]) table without mismatch. |
| RM(3,7) ([128,64,16]) | `spectra/rm37_128_64_spectrum.csv` | The coefficients are the published RM(3,7) distribution of Sugino, Ienaga, Tokura, and Kasami, *IEEE Transactions on Information Theory* 17 (1971), 627–628, [DOI 10.1109/TIT.1971.1054684](https://doi.org/10.1109/TIT.1971.1054684). They have mass (2^{64}), minimum weight 16, and complement symmetry. SHA-256 `89eb543fb58c193e09069cc1288b857b3336be0df1afff061cee62d347072c5a`. |
| RM(4,9) ([512,256,32]) | `scripts/rm512_256_spectrum.csv` | The repository certificate attributes the exact table to Markov–Borissov 2025, Table 6, and records a successful Gleason-form check. The spectrum has mass (2^{256}), minimum weight 32, and SHA-256 `995aab561da18f22074b5c6f5413882f492084510aa1cd19b848355f1fcd4ed7`. |

`tools/enumerate_binary64_spectrum.rs` is the retained exhaustive enumerator
for the shortened-XBCH row file. It partitions the \(2^{32}\) messages across
the available hardware threads, traverses each \(2^{20}\)-word low half in
Gray-code order, and writes only nonzero histogram entries. The generated CSV
was accepted only after its mass, minimum distance, and frozen manifest hash
all matched.

The RM(1,3) and RM(2,5) curves reuse the first two enumerators because those
RM codes have the same complete weight distributions as the corresponding
extended BCH constituents. This statement concerns the enumerators, not an
assumption that the named constructions are identical as encoded maps.

## Stopping points

The BCH-derived ladder stops after block length 128. The public extended
length-256 index lists central dimensions but does not provide the required
complete ([256,128]) table. The next primitive candidate discussed locally
has parameters ([511,259,61]), and extension gives dimension 259 rather than
256. Selecting a 256-dimensional subcode would require a specified subcode
and its complete spectrum. `scripts/bch512_spectrum_obstruction.json` records
why the available parameter-only envelopes are insufficient. No BCH-derived
curve at block length 256, 512, or 1024 is claimed.

The nontrivial rate-half RM sequence below the requested cap is

\[
  \operatorname{RM}(1,3),\ \operatorname{RM}(2,5),\
  \operatorname{RM}(3,7),\ \operatorname{RM}(4,9),
\]

with lengths (8,32,128,512). For odd (m), the middle-order code

\[
  \operatorname{RM}((m-1)/2,m)
\]

has rate one half. For even (m), no integer middle order has dimension
exactly (2^{m-1}). Hence there is no rate-half RM member at length 1024;
the next member is RM(5,11) at length 2048.

The random ladder has an analytic expected spectrum and therefore does not
stop for lack of a table. The experiment includes every doubling from block
length 8 through 1024.
