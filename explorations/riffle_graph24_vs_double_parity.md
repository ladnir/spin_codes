# Graph-24 versus Double-Parity for packet-width-eight Riffle

## Outcome

The structured outer is cheaper in the current transpose implementation.
Across two independent 31-trial campaigns, it reduced complete runtime by 3.5% to 4.9%.
It reduced the isolated outer transpose by 10.9% to 15.1%.

The engineering result does not make **Riffle Double-Parity g=8** proof-viable.
The later proof audit gives an exact lower bound

\[
\mathbb E[Z_d]\geq2^{-23.016236846906736}.
\]

Thus Double-Parity cannot meet the $2^{-40}$ first-moment target in this sequential layout.
See `explorations/riffle_double_parity_proof_audit.md`.

## Equal-length comparison

Fix $K=2^{20}$, and divide the message into

\[
B:=K/64=16{,}384
\]

blocks of 64 binary coordinates.
Both variants contain $B+2=16{,}386$ BCH blocks.
Each BCH block has length 128, so

\[
N'=(B+2)128=2{,}097{,}408=2^{21}+256.
\]

The benchmark calls the baseline **Riffle Graph-24 g=8 (unpunctured)**.
This name distinguishes it from the current power-of-two **Riffle g=8** construction.
Both benchmark variants use the same uniformly shuffled packet order at setup.
Both use the existing packet-width-eight inner transpose online.

### Riffle Graph-24 g=8 (unpunctured)

Let $m\in\mathbb F_2^K$.
Setup samples a binary matrix $R\in\mathbb F_2^{24\times K}$.
The outer map computes $r:=Rm\in\mathbb F_2^{24}$.

The construction places $r$ and 40 zero coordinates in one 64-coordinate tail block.
It also appends the XOR of all message and tail blocks.
The construction encodes every block with the existing extended BCH $[128,64,22]$ code.

The fused transpose stores three random bytes per message coordinate.
It performs three table XORs per coordinate, plus 741 XORs to build its tables.
For $K=2^{20}$, the graph step therefore uses 3 MiB and 3,146,469 block XORs.

### Riffle Double-Parity g=8

Let

\[
\mathbb F:=\mathbb F_2[x]/(x^{64}+x^4+x^3+x+1).
\]

For $\alpha\in\mathbb F$, let $M_\alpha$ denote multiplication by $\alpha$ in the polynomial basis.
Define $\alpha_i:=x^i$ for $0\leq i<B$.
The polynomial is irreducible, and $x$ is primitive.
In particular, the coefficients $\alpha_0,\ldots,\alpha_{B-1}$ are distinct and nonzero.

For message blocks $m_0,\ldots,m_{B-1}\in\mathbb F_2^{64}$, define

\[
p_0:=\bigoplus_{i=0}^{B-1}m_i,
\qquad
p_1:=\bigoplus_{i=0}^{B-1}M_{\alpha_i}m_i.
\]

The construction encodes every $m_i$, $p_0$, and $p_1$ with the same BCH code.
It stores no random outer matrix.

Let $z_0,z_1\in V^{64}$ be the adjoints obtained from the two parity BCH blocks.
Here, $V$ is the binary symbol space; the benchmark uses $V=\mathbb F_2^{128}$.
The addend for message block $i$ is

\[
z_0+M_{\alpha_i}^{\mathsf T}z_1.
\]

Multiplication by $x$ gives a sparse recurrence for the weighted state.
The implementation represents its shift as a cyclic array.
Each block uses 64 XORs to combine the two parity contributions.
Advancing the state uses three more XORs.
The outer mixing therefore uses

\[
K+3B=1{,}097{,}728
\]

block XORs.
This count is 2,048,741 below the Graph-24 count.
The specialized transpose uses 3 KiB of hot scratch and no persistent outer metadata.

Distinct coefficients give the outer block code distance three.
This property excludes nonzero codewords supported on only one or two BCH blocks.
It does not imply the random-subcode weight bound used by the existing Graph-24 argument.

## Benchmark method

The benchmark source is `libOTe_Tests/RifflePacket8Graph24DoubleParity_Bench.cpp` in the libOTe worktree.
It applies the following controls:

1. Both variants use the same $K$, $N'$, BCH code, packet order, and packet-width-eight inner transpose.
2. A master input is copied into each work buffer outside the timed interval.
3. Each trial alternates the order of the two variants.
4. Three untimed warmups precede each measurement series.
5. Each optimized outer transpose is checked against an independent full-output reference.
6. The process pins the benchmark thread to one logical processor.
7. Benchmark processes run serially.

The Double-Parity reference constructs each multiplication matrix from its columns.
The optimized reference check therefore does not reuse the cyclic-state recurrence.
The executable also rejects repeated coefficients before timing.

## Results

Each campaign used 31 paired trials.
The paired ratios below are medians of per-trial ratios.

| Campaign | Graph outer (ms) | Double-Parity outer (ms) | Paired outer ratio | Graph complete (ms) | Double-Parity complete (ms) | Paired complete ratio |
|---|---:|---:|---:|---:|---:|---:|
| A | 4.033 | 3.619 | 0.891 | 8.530 | 8.234 | 0.965 |
| B | 2.476 | 2.098 | 0.849 | 8.621 | 8.179 | 0.951 |

Absolute isolated-outer times changed between campaigns.
The paired direction and complete times remained stable.
Campaign A saved 0.296 ms in the complete path.
Campaign B saved 0.428 ms.

## Interpretation

Double-Parity exploits the 40 zero coordinates and 24 graph coordinates that already occupy a complete BCH block.
It replaces that padded block with a full structured parity block at the same unpunctured length.
The replacement removes the 3 MiB random schedule and most outer-mixing XORs.

The engineering result is favorable but the proof result is negative.
The shared inner transpose still dominates complete runtime.
The exact five-block obstruction rules out the unchanged Double-Parity candidate.

The next design step must change the positive entropy slope of contiguous eight-bit packets.
A constant number of additional parity blocks is insufficient.
