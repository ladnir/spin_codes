# Rate-half BCH, RM, and random Q1 curves

## Experiment

Let (k=2^e), with (8\le e\le20). For every constituent, set

\[
  M(k)=e+2,
  \qquad N=2k,
  \qquad D=\lceil0.10N\rceil.
\]

A deterministic ([B,B/2]) constituent is repeated in

\[
  L=2k/B
\]

outer rows. The row-coordinate permutations, transposed-region permutations,
and RandomStepConv maps are sampled independently once and then fixed for the
code. A random control samples one uniform full-rank binary ([B,B/2]) map
once and repeats that same map in all rows. It does not resample the outer map
from row to row.

For each nonzero constituent weight shell, the calculation evaluates the
occupation-one transfer moment on a finite tilt grid. It optimizes the valid
Chernoff witness separately for each shell and then sums the shell bounds.
The displayed quantity is

\[
  -\log_2 Q_1^{\mathrm{ub}},
\]

where (Q_1^{\mathrm{ub}}) is this union-bound diagnostic. The arithmetic is
nearest binary64. Therefore the curves are neither outward-rounded bounds nor
full distance-certificate margins. Occupations (Q\ge2) remain separate
proof obligations.

## Structured results

| Family | Constituent | Largest plotted Q1 margin | Location | Matched-random maximum |
|:---|:---|---:|:---|---:|
| BCH-derived | extended BCH ([8,4,4]) | -1.282 bits | (k=2^8) | -5.566 bits |
| BCH-derived | extended BCH ([32,16,8]) | 4.107 bits | (k=2^8) | 0.866 bits |
| BCH-derived | shortened-XBCH ([64,32,12]) | 11.916 bits | (k=2^9) | 8.906 bits |
| BCH-derived | extended BCH ([128,64,22]) | 30.014 bits | (k=2^{11}) | 24.885 bits |
| RM | RM(1,3) ([8,4,4]) | -1.282 bits | (k=2^8) | -5.566 bits |
| RM | RM(2,5) ([32,16,8]) | 4.107 bits | (k=2^8) | 0.866 bits |
| RM | RM(3,7) ([128,64,16]) | 13.587 bits | (k=2^9) | 24.885 bits |
| RM | RM(4,9) ([512,256,32]) | 38.925 bits | (k=2^{11}) | 129.461 bits |

None of the authenticated structured curves reaches the 40-bit Q1 screen.
RM(4,9) is closest and misses it by 1.075 bits. The complete BCH-derived
spectra through length 128 are slightly better than their matched random
expected spectra in this transfer. The trend reverses for RM(3,7), and the
gap is large for RM(4,9). This is a statement about the weighted Q1 transfer,
not about minimum distance alone.

## Random size ladder

| Random constituent | Largest plotted Q1 margin | First plotted point at or above 40 bits |
|:---|---:|:---|
| ([8,4]) | -5.566 bits | none |
| ([16,8]) | -3.415 bits | none |
| ([32,16]) | 0.866 bits | none |
| ([64,32]) | 8.906 bits | none |
| ([128,64]) | 24.885 bits | none |
| ([256,128]) | 58.745 bits | (k=2^8), 46.514 bits |
| ([512,256]) | 129.461 bits | (k=2^9), 103.518 bits |
| ([1024,512]) | 274.879 bits | (k=2^{10}), 220.125 bits |

The first random block size that clears the Q1 screen is 256. For block sizes
256, 512, and 1024, the first closing point has exactly two outer rows. The
single-row boundary point for the 512- and 1024-bit constituents is much
weaker; the routing geometry changes discontinuously there.

The random curves use the exact ensemble expectation

\[
  \mathbb E A_w
  =(2^{B/2}-1)\frac{\binom Bw}{2^B-1}.
\]

This identity suffices for the joint occupation-one union bound over the
random outer, routing, and inner choices. It does not imply that a sampled
outer constituent simultaneously realizes the expected spectrum, and it does
not discharge the relation-type accounting required for repeated-map
occupations (Q\ge2).

## Receipts and remaining proof work

- `rate_half_family_k_margin_d100.json` is the complete machine-readable
  receipt.
- `rate_half_family_k_margin_d100.csv` is the plotting table.
- `evaluate_rate_half_family_curves.py` reproduces both receipts.
- `audit_rate_half_family_curves.py` checks every structured spectrum, the
  shortened-XBCH manifest hash, family coverage, parameter schedule, reported
  maxima, and exact reproduction of the earlier matched curves.
- `SPECTRUM_SOURCES.md` records the spectrum provenance and exact stopping
  rules.
- `spectra/xbch64_32_philips_spectrum.csv` was reconstructed by exhaustive
  enumeration; its hash matches the frozen external manifest.
- `tools/enumerate_binary64_spectrum.rs` is the retained exhaustive enumerator
  used for that reconstruction.

The next proof step should not extend the Q1 sweep. It should choose a block
size and handle occupation two with one reused outer map. If (B=256) is
acceptable, the random curve identifies it as the smallest idealized block
that clears 40 Q1 bits under the neutral schedule. To turn that observation
into a structured-code certificate, one still needs either an authenticated
complete spectrum near ([256,128]), or a realized combined-spectrum event
strong enough for the Q1 and Q2 transfer functionals.
