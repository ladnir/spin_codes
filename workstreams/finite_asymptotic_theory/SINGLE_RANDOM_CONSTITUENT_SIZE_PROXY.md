# Single random constituent size proxy

## Question

The Random SPIN proof samples an independent random local injection in every
outer block. That law gives a product of expected local enumerators. It does
not model one constituent that is sampled once and repeated in every row.

This proxy asks which rate-half constituent length first has enough sparse
margin under the repeated-code law. The target parameters are

\[
 k=2^{20},\qquad N=2^{21},\qquad
 D=\left\lceil0.109N\right\rceil=228590.
\]

## Probability space

For an even block length \(B\mid N\), set

\[
 K:=B/2,\qquad L:=N/B.
\]

Sample one uniform linear injection

\[
 G:\mathbb F_2^K\longrightarrow\mathbb F_2^B.
\]

Use this same map in all \(L\) outer rows. Independently sample a uniform
coordinate permutation for every row and a uniform row permutation for every
transposed region. Sample the RandomStepConv-M30 maps once and reuse the
complete setup for every message.

This model retains the region-permuted bit transpose. It does not use the
paper's independently resampled outer injections.

## Exact occupation-one law

Fix a nonzero local message \(u\in\mathbb F_2^K\). For a uniform injection
\(G\), the image \(G(u)\) is uniform on
\(\mathbb F_2^B\setminus\{0\}\). Therefore, the expected weight-\(w\)
multiplicity of the single sampled constituent is

\[
 \mathbb E_G[A_w(G)]
 =
 (2^K-1)\frac{\binom Bw}{2^B-1}.
 \tag{1}
\]

Linearity of expectation makes (1) exact for occupation one even though the
same \(G\) is reused in every row. The calculation sums every shell in (1),
the row and region permutations, and the RandomStepConv transfer.

The nearest-binary64 diagnostic gives:

| \(B\) | \(K\) | \(L\) | Dominant weight | Occupation-one margin |
|---:|---:|---:|---:|---:|
| 128 | 64 | 16384 | 14 | 18.7859 bits |
| 256 | 128 | 8192 | 28 | 55.3387 bits |
| 512 | 256 | 4096 | 56 | 127.9284 bits |

The shell masses in the receipt sum to \(2^K-1\), up to less than
\(6\times10^{-14}\) bits of binary64 reconstruction error.

## Interpretation

Under the current Chernoff transfer, \(B=128\) does not supply 40 sparse
bits. This failure is not a lower bound on the ensemble's true distance.

The two-block size \(B=256\) is the first tested size above 40 bits. It has
about 15.34 bits of occupation-one headroom. The four-block size \(B=512\)
has much more headroom, but it is not needed to pass this first test.

At \(B=256\), the memory scan gives:

| Memory | Occupation-one margin |
|---:|---:|
| 20 | 35.7452 bits |
| 21 | 42.1087 bits |
| 22 | 46.5928 bits |
| 30 | 55.3387 bits |

Memory 21 is the first tested value above 40 bits, but it leaves only 2.11
bits for all other occupations. Memory 22 is the safer first full-proof
target. Memory 30 remains useful when the experiment should isolate the
outer constituent from state-extinction losses.

Thus, the first proxy for the BCH--accumulator experiment should use two
128-bit BCH blocks and global 256-coordinate permutation--accumulation
layers. The target random spectrum is that of a single random
\([256,128]\) constituent, not two independently sampled \([128,64]\)
constituents.

## Occupation-two result

Equation (1) does not complete the distance proof. For two or more active
rows, their local messages can be linearly dependent. Reusing \(G\) then
correlates their outer codewords.

The independently resampled Random SPIN enumerator replaces these joint
moments by a product. That replacement is invalid here. The next calculation
must classify the active local messages by rank and retain the corresponding
joint image law under one injection.

The occupation-two calculation retains this dependence. It separates:

- equal nonzero local messages, which have rank one;
- distinct nonzero local messages, which have rank two.

At \(B=256\), memory 22, and 10.9% distance, the exact-shell binary64
diagnostic gives 64.1559 bits for occupation two. The rank-one class gives
64.1559 bits. The rank-two class gives 99.7046 bits. Thus the rank-one class
dominates, but occupation one remains the overall bottleneck.

`SINGLE_RANDOM_CONSTITUENT_Q2.md` defines the joint law and proves the
finite recurrence used by the diagnostic. Occupations three and above remain
open. Their analysis must retain the rank of the active local messages.

## Evidence and status

`evaluate_single_random_constituent_q1.py` implements (1) and the exact
occupation-one transfer. Its receipt is
`single_random_constituent_B128_256_512_q1_s30_d109.json`.
The B=256 memory receipts are
`single_random_constituent_B256_q1_s20_d109.json`,
`single_random_constituent_B256_q1_s21_d109.json`, and
`single_random_constituent_B256_q1_s22_d109.json`.

These files are diagnostics. They do not claim an outward certificate. The
size-proxy receipt covers occupation one, and
`single_random_constituent_B256_q2_s22_d109.json` covers occupation two.
