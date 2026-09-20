# SPIN code naming

The paper-facing family name is **SPIN codes**, expanded at first use as
**Single-Permutation INterleaved codes**.

Always typeset the acronym as uppercase `SPIN`. The lowercase term *spin code*
already denotes a class of quantum error-correcting codes. The paper title or
abstract should therefore include the full expansion. Do not shorten `SPIN
code` to lowercase `spin code`.

A SPIN encoder applies an outer code, one interleaver, and an inner recursive
code. The interleaver is one permutation of the complete outer-code output.
An implementation may factor that permutation into block shuffles, a
transpose, and region shuffles.

The canonical variant names are:

| Name | Outer | Interleaver | Inner | Role |
|---|---|---|---|---|
| **Accumulator SPIN** | random block code | uniform permutation | accumulator | proof warmup |
| **Random SPIN** | random block code | uniform permutation | random convolution | distance benchmark |
| **Structured SPIN** | structured block code | factored structured permutation | structured convolution | frozen main code |
| **Linear SPIN** | linear-time block code | linear-time permutation | linear-time structured convolution | future asymptotic variant |

Use one name for each variant. Do not use `Fast SPIN` as an alias for
`Structured SPIN`; speed is a measured property rather than the construction's
definition.

Exact instances use parameters after the family name when needed. The frozen
instance may be written as **Structured SPIN (B=256, t=128, s=19)**. Source
identifiers may use `StructuredSpinB256T128S19`. The construction manifest,
not the prose name, records ParityFanout-31x33 and the complete constituent
schedule.

The previous `Riffle ...` names remain as legacy artifact identifiers. Frozen
directories, receipt schemas, hashes, and historical exploration records keep
their existing paths. New paper text and cleaned implementation interfaces use
the SPIN names.
