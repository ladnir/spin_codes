# GF256 rank history: valid transfer, unresolved six-state classes

Retaining the state-rank cap improves the tested GRS[128,86] bounds, but
does not close the six-state proof gap. At input densities .10 and .14,
the gains over the joint two-state cover are about 335 and 389 bits.
Sharper regional conditioning subsequently makes the .14 class positive.
The .10 class still has a deficit of about 17,151 bits.

The [screen](theory/gf8_rank_history_screen.py) supports the checked
width-sixteen Cauchy map restricted to four, six, or eight state coordinates.
Its field is GF256. The construction retains one regional routing stage.
The sixteen coordinates of a 128-bit record share the same scalar encoder.
The [independent audit](GF8_RANK_HISTORY_AUDIT.md) gives the local proof
contract and the cap-dependent joint return cover.

Transfer state zero denotes actual zero state. Transfer state c>0 bounds
the incoming state-form rank by c. Continuations must be nonnegative and
nondecreasing on positive caps; their value at zero is independent.
After selected inactive zero equations leave rank H, w active inputs and
i selected active zeros give successor cap min(s,H+w-i).
The transfer preserves this cap instead of restarting with the maximum
over all incoming ranks at every packet.

The return cover fixes the incoming zero flat before exposing current
diagonals. It bounds the flat rank by min(w-1,c-1,s-1-z), where z counts
selected inactive zeros. The pure and joint return covers both increase
with c, so their minimum preserves the continuation cone. The joint cover
charges each current diagonal only once.

All reported values are floating evaluations of individual integer input
weights. The first table uses the existing finite source comparison and
its regional conditioning factor (L+1)^128. No whole-code certificate or
native timing is produced. A negative margin means this upper bound fails
to establish probability decay; it does not exhibit a low-weight codeword.

For requested K=2^20, GRS[128,86] has L=12224 and N=1564672 after row
rounding. The forbidden output threshold is ceil(N/10)-1=156467.
The margin is -log2(B), where B is the evaluated upper bound.

| Input weight | Density | Two-state joint | Rank cap and joint | Rank cap with sharper returns |
|---:|---:|---:|---:|---:|
| 156467 | approximately .10 | -18314.369 | -17979.224 | -17977.829 |
| 219054 | approximately .14 | -1091.954 | -703.387 | -703.387 |

The last column uses the six-state return-output minima 15 and 13 for
input weights one and two. Those values come from the separate
[exact return-geometry enumeration](theory/gf8_return_geometry.py).
The minima strengthen both the pure output factor and the joint mask
cutoff. They improve the .10 bound by only 1.395 bits and have negligible
effect at .14.

The control comparison also evaluates the new transfer at the old
optimizer's parameters. Those margins are -17986.193 and -709.485 bits.
Thus most of the reported gain comes from the transfer itself, rather
than the optimization routine finding a different local minimum.
All four displayed new optimizations report success.

The batched aggregate sums selected active coordinates analytically before
forming the matrix. This avoids constructing seventeen packet matrices at
each objective evaluation. An exact local reference remains available for
regression tests. Thirteen tests pass, covering the old eight-state formula,
full-cap two-state masses, all continuation-cone rays, successor caps,
positive composition, source and route factors, and the batched aggregate.
The source refinements below also receive exact finite-profile checks.
Floating tail differences are checked against exact rational entries with
relative tolerance 2e-12; this numerical screen is not directed arithmetic.

```sh
python -B research/workstreams/char2_field_spin/theory/gf8_rank_history_tests.py
python -B research/workstreams/char2_field_spin/theory/gf8_rank_history_screen.py --outer 128 86 --state 6 --density .10 .14 --iterations 35 --compare
python -B research/workstreams/char2_field_spin/theory/gf8_rank_history_screen.py --outer 128 86 --state 6 --density .10 .14 --iterations 45 --return-minima 15 13 --compare
```

## Feasible source supports and modal conditioning

A nonzero GRS[128,86] row has at least 43 active coordinates. The old
source polynomial includes row weights one through 42 and counts an extra
empty row. The exact feasible row polynomial removes those terms. However,
that polynomial must still bound unequal regional support densities using
their common mean. The old proof's two-product averaging argument does
not directly establish this stronger statement.

Put b=128, p=42, and beta=tau*lambda, where lambda is the source tilt.
For a Bernoulli row with H active coordinates, the feasible source weight
is g(0)=1, g(h)=0 for 1<=h<=p, and g(h)=tau^(-p)*beta^h for h>p.
When beta>=2, this sequence has nonnegative second differences. Averaging
any two regional densities therefore increases the row expectation.
Repeated averaging bounds the heterogeneous expectation by

\[
(1-\alpha)^b+\tau^{-p}\sum_{h=p+1}^{b}
 \binom bh(\beta\alpha)^h(1-\alpha)^{b-h}.
\]

The screen evaluates this positive sum with log-sum-exp. For beta<2 it
retains the established relaxed polynomial. The restriction is necessary
for this averaging argument: when b=3, p=1, tau=3, and beta=1, densities
(1,1,0) give expectation 1/3, exceeding 23/81 at their common mean.

For even L, the existing modal-conditioning lemma replaces L+1 by
2^L/binom(L,L/2). Each regional count is compared with a Bernoulli column
at that count's own mean. This refinement reduces the loss by 827.271022
bits for the present geometry, independently of the four optimized tilts.

| Density | Refined sampled margin | Gain from feasible row supports | Gain from conditioning |
|---:|---:|---:|---:|
| approximately .10 | -17150.557106 | 0.000509 | 827.271022 |
| approximately .14 | 123.883890 | 0.000000625 | 827.271022 |

These optimizations include rank history and return minima 15 and 13.
Both report success and use beta>=2. The source gains are evaluated at
the refined optimizer's parameters. Despite the small unweighted density,
the tilted occupied-row term favors weights well above 43. Removing its
short-support terms consequently makes little difference.

```sh
python -B research/workstreams/char2_field_spin/theory/gf8_rank_history_screen.py --outer 128 86 --state 6 --density .10 .14 --iterations 45 --return-minima 15 13 --feasible-source --modal --compare
```

The GRS[128,86] next step needs to address the .10 class, whose deficit
remains about 17,151 bits. This result does not justify a certificate
campaign or native implementation of that profile. A stronger account of
state and input dependencies, or a changed packet/source geometry, remains
necessary before expanding the higher-rate search.

## Lower-rate six-state control and the operation-count gate

Reducing the outer dimension to 80 gives padded rate 5/8. The six-state
map, return minima 15 and 13, feasible source comparison, and modal
conditioning remain the same. At requested K20, this profile has 13120
rows, 1049600 padded source records, and 1679360 encoded records.
All six bounded middle-density screens are positive:

| Requested input density | Sampled margin in bits |
|---:|---:|
| .06 | 32446.604 |
| .08 | 34496.118 |
| .10 | 31161.754 |
| .12 | 30774.769 |
| .14 | 73314.553 |
| .18 | 228026.845 |

All optimizations report success. These points support further proof
screening of this profile, but supply neither interval coverage nor sparse
and dense certificates. The dimension-76 fallback was not pursued.

```sh
python -B research/workstreams/char2_field_spin/theory/gf8_rank_history_screen.py --outer 128 80 --state 6 --density .10 .14 --iterations 45 --return-minima 15 13 --feasible-source --modal --compare
python -B research/workstreams/char2_field_spin/theory/gf8_rank_history_screen.py --outer 128 80 --state 6 --density .06 .08 .12 .18 --iterations 45 --return-minima 15 13 --feasible-source --modal
```

The instruction comparison below uses the existing **systematic**
GRS[128,87], state-eight implementation, whose outer uses 516 GFNI512
instructions per four rows. It is not a comparison against the separate
coefficient-basis path with a 321-instruction outer.

A straightforward restriction of the dense packed inner costs 62 GFNI512
instructions per packet: 24 emission products, 32 injection products, two
state-diagonal products, and four input-diagonal products. Six state
records still occupy two four-record vectors. The injection and diagonal
instruction counts therefore do not fall in proportion to the state size.

For the dimension-80 outer, split the source prefix as 64+16. One
64-point translation and a 16-point interpolation to three further cosets
give a parity transpose with 434 GFNI512 instructions per four rows.
Its 978 XOR operations exclude the final systematic additions.
The new test compares all 3840 entries of this scalar circuit against
independent Lagrange interpolation. No native circuit was generated.

| Requested K20 operation count | Systematic GRS[128,87], state eight | Proposed systematic GRS[128,80], state six |
|---|---:|---:|
| Encoded records | 1548288 | 1679360 |
| Inner GFNI512 | 6773760 | 6507520 |
| Outer GFNI512 | 1560384 | 1423520 |
| Combined GFNI512 | 8334144 | 7931040 |

These full-padded-row counts give only 4.84% fewer GFNI instructions,
while the encoded record count and proportional route traffic rise 8.47%.
They exclude loads, stores, routing, coefficient handling, and spills.
They do not establish a runtime improvement. The lower rate makes the
sampled proof plausible, but this straightforward six-state packing does
not provide a compelling case for a faster complete encoder. A later
implementation decision should first identify a larger arithmetic or
packing gain that compensates for the extra routed data.
