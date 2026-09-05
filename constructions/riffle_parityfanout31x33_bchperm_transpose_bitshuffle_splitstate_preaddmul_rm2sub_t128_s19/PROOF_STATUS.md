# Proof status

## Floating-point end-to-end result

For `N=2^21` and `D=floor(0.11 N)=230686`, the current first-moment ledger
has the following components.

| Outer occupation | Margin (bits) |
|---|---:|
| 1 | 55.864642340 |
| 2 through 100 | 108.824380510 |
| 101 through 8192 | 121.019237921 |
| Combined | 55.864642340 |

The three ranges partition every nonzero message. ParityFanout-31x33 has no
exceptional all-one endpoint: its complete nonzero expected spectrum is
covered by the same Bernoulli envelope.

The envelope costs 128.006509470 bits per active outer block. The bulk
receipts are exact fixed-reference reweightings of the state-size-19 inner
transfer ledger because the Bernoulli reference probability remains 1/2.
The change in outer spectrum adds one tilt-independent likelihood constant
per active block.

## Assumptions and remaining gate

The result assumes the modeled even-floor `[256,128]` source spectrum and an
independently sampled ParityFanout-31x33 map for each outer-block instance.
It also inherits the audited fixed RM2Sub-S19 constituent assumptions.

Arithmetic in the current component receipts is nearest binary64. The final
proof gate is an independently checkable outward-rounded evaluation. Until
that gate is complete, 55.864642 bits is a strong diagnostic rather than the
final unconditional certificate.

The fine-grid distance curve is recorded in
`receipts/floating_distance_curve_d09_d11.json`. It gives 66.270087 bits at
distance 0.09 and 61.069584 bits at distance 0.10. The one-active class limits
all three distances.

## Component receipts

The active exploration receipts currently reside in the adjacent
`riffle_parityshear12_..._s19/receipts` directory:

- `parityfanout31x33_outer256_d38_expected_spectrum.json`;
- `parityfanout31x33_s19_one_active_d11.json`;
- `parityfanout31x33_s19_regular_2_100_d11.json`;
- `parityfanout31x33_s19_regular_101_8192_d11.json`.

The compact combined ledger is `receipts/floating_combined_d11.json`.
