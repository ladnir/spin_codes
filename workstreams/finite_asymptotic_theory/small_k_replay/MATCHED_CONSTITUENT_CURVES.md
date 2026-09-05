# Matched constituent margin curves at 10% distance

## Comparison rule

For each message length \(k=2^e\), set

\[
  M(k)=e+2=\lceil\log_2 k\rceil+2.
\]

This schedule recovers the existing \(M=22\) choice at \(k=2^{20}\). It does
not depend on the outer constituent. Every curve uses output length \(N=2k\)
and cutoff

\[
  D=\lceil0.10N\rceil.
\]

Each structured curve repeats one fixed authenticated constituent in all
outer rows. Each random curve samples one uniform full-rank constituent of
the same block length and dimension, then repeats that one sampled map in all
rows. Routing and RandomStepConv use the same probability space in every
case.

The plotted quantity is the occupation-one union-bound margin
\(-\log_2 Q_1\). It is a spectrum screen, not the margin of a complete
distance certificate. Higher occupations can only add proof obligations.

## Results

| Constituent | Largest plotted Q1 margin | Location |
|:---|---:|:---|
| extended BCH \([32,16,8]\) | 4.107 bits | \(k=2^8, M=10\) |
| random full-rank \([32,16]\) | 0.866 bits | \(k=2^8, M=10\) |
| extended BCH \([128,64,22]\) | 30.014 bits | \(k=2^{11}, M=13\) |
| random full-rank \([128,64]\) | 24.885 bits | \(k=2^9, M=11\) |
| RM(4,9) \([512,256,32]\) | 38.925 bits | \(k=2^{11}, M=13\) |
| random full-rank \([512,256]\) | 129.461 bits | \(k=2^{13}, M=15\) |

No authenticated structured constituent reaches 40 Q1 bits under this
schedule. RM(4,9) comes closest, missing by about 1.075 bits at \(k=2^{11}\).
The exact BCH constituents outperform their matched random ensembles, but
their absolute margins remain below 40. The 512-bit random ensemble strongly
outperforms RM(4,9), isolating the RM low-weight spectrum as the same-size
gap; it says nothing about whether a 32-bit or 128-bit outer can close.

The curves are not monotone in \(k\). Increasing \(k\) improves the inner
length effect at first, while the number of outer rows and hence the number of
low-occupation messages also grows. Under \(M(k)=\log_2 k+2\), the latter
eventually costs about one margin bit for each doubling of \(k\) in several
curves.

## Reproducibility

- `evaluate_matched_constituent_curves.py` defines and runs the comparison.
- `matched_constituent_k_margin_d100.json` records the probability space,
  schedule, and full-precision diagnostic cases.
- `matched_constituent_k_margin_d100.csv` is the plot-ready table.

The next useful sensitivity check is a second common schedule such as
\(M(k)=\lceil\log_2 k\rceil+3\). It should be plotted separately rather than
choosing a different memory offset for each constituent.
