# A cheaper fixed t64/s16 map

The complete fresh replay gives a setup-failure margin of
**62.8401523533463731 bits** at K=65,536 and N=131,072. In particular,

\[
\Pr[d_{\min}\le 13107] < 2^{-62.8401523533463731} < 2^{-40}.
\]

The probability is over the independent setup distribution defined below.
The resulting minimum distance is at least 13,108 outside the failure event.
The receipt is `disjoint-pairs-whole-v2-p256.json`; it covers all 512
nonzero-message occupancies with fresh 256-bit outward computations.

This candidate retains the small RS16 outer, four-bit routing packets, and
one independent state refresh per 64 input bits. It changes the fixed
expansion and feedback maps. The quadratic rows have disjoint monomial
supports: three rows use one monomial, and six use a pair. This permits
copies on expansion and six XORs on feedback for the quadratic basis
conversion. These operation counts require an implementation timing.

For z=(z_0,...,z_5) in F_2^6, order the 64 output coordinates by the integer
sum_i 2^i z_i. Define sixteen Boolean functions in this order:

\[
\begin{aligned}
f_0&=1, & f_{i+1}&=z_i\quad(0\le i<6),\\
f_7&=z_0z_1, & f_8&=z_0z_3, & f_9&=z_4z_5,\\
f_{10}&=z_0z_2+z_1z_4, &
f_{11}&=z_0z_5+z_2z_3, &
f_{12}&=z_1z_3+z_2z_4,\\
f_{13}&=z_0z_4+z_2z_5, &
f_{14}&=z_1z_5+z_3z_4, &
f_{15}&=z_1z_2+z_3z_5.
\end{aligned}
\]

Define A in F_2^(64×16) by A[z,i]=f_i(z), and define C=A^T. The rows are
independent because each quadratic monomial belongs to exactly one f_i,
and the seven affine functions are independent. Products f_i f_j have
degree at most four, below six. Their sums over F_2^6 vanish, so CA=0.
Each consecutive four-coordinate packet varies z_0 and z_1. Its restriction
has rank four because the family contains 1,z_0,z_1,z_0z_1.

Fresh enumeration of all 65,536 state images gives the following expansion
spectrum. The retained map has the same spectrum; equality of spectra alone
does not transfer its distance certificate.

| Image weight | 0 | 16 | 24 | 28 | 32 | 36 | 40 | 48 | 64 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| State count | 1 | 20 | 4,640 | 13,824 | 28,566 | 13,824 | 4,640 | 20 | 1 |

The complete code has 512 groups, each mapping 128 message bits to 256 bits
through four parallel GF(16) RS[16,8] rows. Each aligned 16-bit symbol gets
an independent invertible linear randomizer that sends every fixed nonzero
input uniformly to the nonzero vectors. Each group independently shuffles
its 64 four-bit packets. Every region independently shuffles its 512 slots.

At physical step i, the inner receives x_i in F_2^64 and state a_i in F_2^16.
It emits y_i=x_i+A a_i and sets a_(i+1)=M_i a_i+C x_i. The initial state is
zero. State continues across all 2,048 steps and is discarded after the
last output. Each M_i is independently sampled from an invertible linear
family with the same fixed-input law as uniform GL(16,2). All random maps
and routing permutations are mutually independent and sampled once at setup.
The outer is injective and the inner is triangular with identity diagonal,
so this defines a binary [131072,65536] linear code.

Uniform GL16, nonzero field16 multipliers, and their binary adjoints satisfy
the randomizer premise. Fixed invertible maps before or after these families
also satisfy it. `TRANSITIVE_FIELD16.md` gives the fixed-message law argument,
including the native outer layout's fixed coordinate permutation.

`disjoint_pair_maps.py` fixes the maps and regenerates their complete state
and packet profiles. `reproduce_disjoint_pairs.py` reads no numerical receipt.
It computes the q=1 contribution using exact outer support counts, and q=2
using all ordered support-pair moments. For q≥3 it recomputes the retained
pointwise uniform outer majorant with this candidate's local operators.
Every regional and global product retains state; none restarts from zero.
The final sum covers each occupancy q=1,...,512.

The replay uses 256-bit outward arithmetic. It sums all resulting positive
dyadic endpoints with exact integers, rounds the total upward once, and
compares that endpoint to 2^-40 with integer arithmetic. An incomplete file
keeps `whole_code_certificate=false`. A successful completed file has both
`all_occupancies_covered=true` and `whole_code_certificate=true`.

The complete replay gives q=1 margin 62.8401523697 bits and q=2 margin
118.2500739887 bits. Among q=3,...,128, the weakest individual contribution
has 89.2326258438 bits at q=5. Among q=129,...,512, the weakest has
1,386.1203134968 bits at q=129. The combined 62.8401523533-bit bound includes
the exact sum of every occupancy contribution.

Run a fresh replay from the repository root:

```text
python -B research/workstreams/k16_codesign_100us/proof/reproduce_disjoint_pairs.py --output research/workstreams/k16_codesign_100us/proof/new-whole-p256.json
```

The output filename must be unused. The claimed failure event is that some
nonzero message has output weight at most 13,107. A complete bound below
2^-40 gives minimum distance at least 13,108 except with that setup-failure
probability. It does not certify any individual deterministic setup seed.
