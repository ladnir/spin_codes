# Packet-width comparison

## Question

The random-outer landscape treats `g`, `B`, and `sigma` as independent
parameters. This report locates the first observed 40-bit contour for each
packet width. The model and probability space are defined in `README.md`.

## Observed contours

| Packet width `g` | First viable outer length `B` | Margin below contour | Margin above contour |
|---:|---:|---:|---:|
| 1 | 256 | `lambda(256,1,23)=38.26` | `lambda(256,1,24)=40.16` |
| 2 | 512 | `lambda(512,2,18)=29.97` | `lambda(512,2,19)=43.12` |
| 4 | 4096 | `lambda(4096,4,19)=27.18` | `lambda(4096,4,20)=47.25` |
| 8 | none through 4096 | not near contour | `lambda(4096,8,40)=-43958.12` |

The high-memory floors explain why the preceding outer lengths fail:

| `g` | Outer length below the first viable point | High-memory floor |
|---:|---:|---:|
| 1 | 128 | 10.18 |
| 2 | 256 | 24.25 |
| 4 | 2048 | 32.62 |

Thus additional memory cannot move those cells to the 40-bit contour.

## Mechanism

Increasing `g` reduces the number of active packets produced by a fixed
number of random outer bits. A message family then pays a smaller placement
cost to move all active packets into a late suffix. The outer multiplicity is
measured in message bits and does not receive the same reduction.

At `g=8,B=4096,sigma=40`, the dominant profile has about 47,600 active
packets and 99 active outer blocks. Those packets can occupy a late suffix
whose random output contributes about 9% total weight. The inner placement
bound costs about 151,000 bits, while the outer family contributes about
195,000 bits. The resulting margin is approximately -43,958 bits.

This profile is a bulk late-start obstruction. Increasing `B` from 4096 to
8192 changes the observed floor by only about 182 bits. The interrupted
large-`B` extension therefore does not indicate a nearby `g=8` contour.

## Scope

The values are ensemble diagnostics, not certificates. The outer saddle and
support interpolation remain approximate. The matrix calculation is an upper
bound on the inner low-weight probability, and small exact tests confirm its
direction. Selected contour cells still require exact outer coefficients and
a complete support cover.
