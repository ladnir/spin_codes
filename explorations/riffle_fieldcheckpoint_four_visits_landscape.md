# FieldCheckpointAccumulate with four visits per lane

The K=8 result raises a design question. Does distance require a checkpoint
every 256 bits, or does it require only four accumulator visits per state
lane before a checkpoint?

Fix implementation step (t=32). For state width (s) and checkpoint period
(K), define

\[
 m:=\frac{tK}{s}.
\]

The value (m) is the number of visits to each state lane during one epoch.
The comparison fixes (m=4) and changes only how the state is packaged.

| Candidate | Epoch bits | Mixes per 1024 bits | Middle-range margin at 9% |
|---|---:|---:|---:|
| (s=64,K=8) | 256 | (4M_{64}) | 19399.3836 bits |
| (s=128,K=16) | 512 | (2M_{128}) | 21135.6506 bits |
| (s=256,K=32) | 1024 | (M_{256}) | 21609.2409 bits |

Here (M_s) denotes the transposed circuit for multiplication by one fixed
element of \(\mathrm{GF}(2^s)^*\). Each row covers 512 through 4096 regular
active outer blocks. The calculation uses the same modeled spectrum, 9%
threshold, Chernoff grid, and exact log-domain region recurrence.

All three candidates have ample distance margin. The small differences are
not important for candidate selection. The comparison supports four visits
per lane as the controlling distance parameter. It does not support checkpoint
count by itself as the controlling parameter.

Each candidate samples approximately 256 multiplier bits per 1024 output
bits. The candidates differ in circuit packaging and state storage. Their
mixing costs are (4M_{64}), (2M_{128}), and (M_{256}). The current
transposed benchmark implements field multiplication only through width 64,
so no measured performance ordering is available for the three expressions.

A superlinear multiplication circuit favors (4M_{64}). Fewer transposes and
better batching can favor a wider map. A performance decision therefore
requires optimized width-128 and width-256 kernels, not scalar field timings.

The experiment does not cover all-one outer words or outward rounding. Those
obligations remain unchanged.
