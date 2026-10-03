# Actual 24-bit byte-native screen

The actual 24-bit candidate passes the bounded q1 gate, but its middle-q
comparison bound remains negative. The diagnostics identify temporal
concentration of potential input activity as the next proof question.
They do not establish failure of the code.

All numerical results here are floating proposals. No result has outward
endpoints, and no receipt claims a whole-code certificate. Existing
authenticated 16-bit sources were imported read-only and were not changed.

## Literal construction

Use the AES polynomial basis for the field F=GF(256), with modulus 0x11b.
The state is a=(a0,a1,a2) in F^3. Its binary encoding places a0 in the
lowest byte, followed by a1 and a2. A physical input X has eight bytes,
indexed by h=0,...,7. Define

\[
 (Aa)_h=a_0+h a_1+h^2a_2,
 \qquad
 CX=\left(\sum_hX_h,\sum_hhX_h,\sum_hh^2X_h\right).
\]

The step emits Y=X+Aa, then updates a'=Ma+CX. Each step samples a fresh
independent M uniformly from GL(3,F). The state starts at zero, persists
across every physical step and region, and is not flushed.

Both binary maps have rank 24. Every restriction of C to one, two, or three
distinct bytes has rank 8, 16, or 24, respectively. A restriction to any four
bytes has rank 24. The three-column assertion also follows from the
Vandermonde determinant at three distinct field points.

The literal composition CA is zero. Its field entries are sums of h^r
for 0<=r<=4 over the three-dimensional binary subspace {0,...,7}.
For these exponents, h^r has Boolean degree below 3, so each sum vanishes.
The code also checks every binary column of CA directly.

For each fixed nonzero state, Ma is uniform on the 2^24-1 nonzero states.
This action law, not equality of entire random-matrix ensembles, is what
the local analysis uses. A different fresh update family would transfer
only after proving the same conditional fixed-state law.

The map identity SHA256 is
`664318711a2d1a10f5add6dd28cba7accafa663c72477163cf2a173325a032ba`.
The full expansion census has 447 occurring unordered byte-weight profiles.
Its minimum nonzero binary weight is 8, attained 56 times. Its spectrum is

```text
weight:  0   8  12   16     20     24      28      32
count:   1  56 672 9820 120352 876680 3622720 7516614
```

The remaining counts are symmetric under weight w becoming 64-w.

## Exact finite formulas

Fix a weight variable 0<z<=1. For occupancy j, let X have a uniformly
chosen j-byte support and independent uniform nonzero byte labels.
Define the weighted feedback measure and fixed-state emission moment by

\[
 W_j(s)=\mathbb E[z^{\mathrm{wt}(X)}\mathbf1\{CX=s\}],
 \qquad
 M_j(a)=\mathbb E[z^{\mathrm{wt}(X+Aa)}].
\]

For each state a, the polynomial M_j(a) depends only on the unordered
profile of the eight byte Hamming weights of Aa. If n_r counts bytes
of weight r, then

\[
 M_j(a)=\frac{[u^j]}{\binom8j}
 \prod_{r=0}^8
 \left(z^r+u\frac{(1+z)^8-z^r}{255}\right)^{n_r}.
\]

The feedback census is separate. For a binary character t of the state,
let r_h be the Hamming weight of the byte character C_h^Tt. Then

\[
 \widehat W_j(t)=\frac{[u^j]}{\binom8j}
 \prod_{h=0}^7\left(1+u\frac{(1+z)^{8-r_h}(1-z)^{r_h}-1}{255}\right).
\]

The inverse Walsh transform, divided by 2^24, recovers W_j. C-character
profiles use the literal binary transposes of field multiplication;
they are not replaced by A profiles.

For each expansion profile P, accumulate

\[
 H_i(P)=\sum_{s\ne0:\,\operatorname{profile}(As)=P}W_i(s),
 \qquad \beta_i=\sum_PH_i(P).
\]

The normalized birth family is B_i(s)=W_i(s)/beta_i for s!=0, with zero
mass at s=0. Its emission moment is exactly

\[
 \sum_s B_i(s)M_j(s)
 =\frac{\sum_PH_i(P)M_j(P)}{\beta_i}.
\]

Thus profile binning preserves every birth-to-emission inner product.
It does not replace a state distribution by a weight-class maximum.

The comparison coordinates are zero, uniform nonzero, and B1,...,B8.
From a nonzero source family with emission moment T, the outgoing
weighted measure is dominated by T times uniform nonzero plus
T/(2^24-1) times the zero point mass. At occupancy 0, the zero term is
absent and the refresh is exact. Zero-source rows use W_j exactly.
Structural zeros W_j(0)=0 for 1<=j<=3 are imposed from the restriction ranks.

The only analytic relaxation here is this dominated outgoing measure.
The finite formulas themselves are exact; their implemented evaluations
are floating and therefore do not constitute certified numerical bounds.

## Bounded-memory evaluation

There are at most binomial(16,8)=12870 unordered profiles of eight byte
weights. `maps24.py` enumerates all 2^24 states in chunks and stores two
uint16 profile-index arrays, using 64 MiB total. It separately counts A
profiles and C-character profiles.

For each occupancy 2,...,8, `screen24.py` forms one float64 character
vector, applies a blocked Walsh transform, removes state zero, and bins
the result by A profile. It then discards the large vector. Occupancy 1
uses a positive direct census of the 2040 one-byte inputs.

The current vector-expression temporaries give a conservative peak
storage estimate below 0.4 GiB; resident memory was not separately measured.
No large arrays are written to disk. A complete A census took 4.71 seconds;
the A-plus-C census took about 9 seconds. Each full local tilt took roughly
13--14 seconds. These are local proof-computation costs, not encoder
benchmarks. BLAS and OpenMP thread counts are capped at one.

The transform records negative mass clipped after inversion and checks
the total against (((1+z)^8-1)/255)^j. At tilt 0.5, clipped mass was
8.85e-18 for j=2 and zero for j=3,...,8. Such checks detect gross
instability but do not replace outward rounding.

The profile-key helper is used only with at most eight packets. Its
uint32 key representation is not a supported wider-window interface.

## Outer code and ordered geometry

The outer code remains four parallel GF(16) RS[16,8] rows, independently
mixed in aligned 16-bit symbols. Each group maps 128 input bits to 256
output bits. K=65536 gives 512 groups and 131072 output bits.

Eight-bit routing has 32 regions and 512 potential byte slots per region.
Each region contains 64 physical steps with eight slots per step. The
state passes continuously through all 2048 physical steps.

For r active slots after e steps, the regional comparison matrices obey

\[
 R_{e,r}=\sum_k
 \frac{\binom8k\binom{8(e-1)}{r-k}}{\binom{8e}r}
 R_{e-1,r-k}T_k.
\]

The matrix order is chronological. Regional boundaries do not reset state.
For the middle-q screen, the pointwise outer envelope is
beta=2^256/65535^8 per active group. It compares each group with uniform
bytes, so the regional active occupancy is Bin(q,255/256). The resulting
regional matrix is raised to 32, starting from zero and retaining every
terminal coordinate. The cutoff is 13107 bits.

The q1 screen instead folds the exact outer packet-support counts through
R0 and R1. It uses the same chronological state transfer. Only one active
packet can occur in a region for q1, so the direct birth census suffices.

## Bounded numerical result

`q1_v1.json` gives **68.717321 bits** on tilts
0.00128,0.00256,0.00384,0.00512,0.00768,0.01024,0.01536.

The following middle-q values take the better value from `middle_cost_v1.json`
and `middle_grid_v1.json`. They use the actual 24-bit maps and state law,
not a changed denominator applied to 16-bit moments.

| Active groups q | Best screened margin, bits | Tilt |
|---:|---:|---:|
|64|-2818.81|0.20|
|90|-3511.19|0.35|
|119|-3810.65|0.50|
|128|-3808.77|0.55|
|160|-3512.26|0.80|
|192|-2916.21|0.80|
|256|-581.56|1.20|

The grid is bounded, not a continuous-tilt optimization or a full-q proof.
It gives no reason to launch an outward replay of the current comparison.

## What dominates the comparison

At q=119 and tilt 0.5, `trajectory_v1.json` records the full local operators
and the 32-region comparison matrix. All path statements below refer to
that positive comparison expression, normalized by its total weight.
They are not probabilities of bad paths under the actual code ensemble.

The maximum-weight regional boundary path nearly alternates zero and
uniform nonzero. Its moment is only 5.773 bits below the full moment.
Restricting every regional boundary to those two coordinates loses only
0.0618 bits. The expected number of physical nonzero-to-zero transitions
is 16.941.

Even the path family with zero state at every regional boundary gives
a negative margin of 3696.75 bits. This family allows nonzero excursions
inside regions; it is not the always-zero-state path.

`occupancy_v1.json` differentiates the log moment after multiplying each
physical occupancy operator T_j by exp(delta). This gives the expected
number of steps of each occupancy under the normalized comparison:

```text
j:       0       1     2     3      4      5       6       7       8
steps: 1499.18  1.38  3.19  7.41  20.83  63.70  144.49  193.63  114.20
```

The counts sum to 2048. Occupied steps average 6.5007 active bytes.
Expected active packets total 3567.7102, compared with 3793.125 before
weight tilting. An independent regional-marker derivative checks the
weighted packet total to within 2.2e-6.

`state_time_v1.json` separately marks zero-source rows. It gives 1517.129
zero-source steps and 530.871 nonzero-source steps. Of the zero-source
steps, 1499.177 are empty and 17.951 are occupied. Thus almost every empty
step carries zero state in the dominant comparison contribution.

At this tilt, the uniform-state emission costs 19.625 bits at occupancy 0
and approximately 20.1--20.23 bits at occupancies 2,...,8. Nevertheless,
the zero-to-zero regional moment costs only 668.31 bits. Long empty periods
and densely occupied bursts avoid paying the uniform-state cost at most
physical steps. Improving that cost alone does not address this behavior.

The next proof question is whether a route-only good event can exclude
such temporal concentration before the union over input labels. A route
event for a set of q outer groups needs a union over binomial(512,q)
sets, not over all messages in those groups. That approach requires a
separate derivation and gate; these diagnostics do not prove it succeeds.

## Reproduction and tests

Run from the repository root, choosing fresh output paths:

```powershell
C:/Python314/python.exe -B -m unittest discover -s research/workstreams/packet8_codesign/larger_state -p 'test_*.py' -v
C:/Python314/python.exe -B research/workstreams/packet8_codesign/larger_state/screen24.py --q1-only --tilts .00128 .00256 .00384 .00512 .00768 .01024 .01536 --output research/workstreams/packet8_codesign/larger_state/q1_replay.json
C:/Python314/python.exe -B research/workstreams/packet8_codesign/larger_state/screen24.py --tilts .5 --q 64 90 119 128 160 256 --output research/workstreams/packet8_codesign/larger_state/middle_replay.json
C:/Python314/python.exe -B research/workstreams/packet8_codesign/larger_state/diagnose24.py --q 119 --tilt .5 --output research/workstreams/packet8_codesign/larger_state/trajectory_replay.json
```

The five portable tests cover literal binary maps and ranks, chunked
profile identities, blocked Walsh transforms, an exhaustive GF(4)^3
birth/emission comparison, and exhaustive small boundary-path statistics.
The final test also checks the marked-matrix derivative by enumeration.
An independent read-only review agreed with the profile compression and
the candidate's literal coordinate conventions.

Every completed numerical receipt records its source hashes or an
authenticated input receipt. No production source, frozen certificate,
remote workload, encoder benchmark, or commit was changed or performed.
