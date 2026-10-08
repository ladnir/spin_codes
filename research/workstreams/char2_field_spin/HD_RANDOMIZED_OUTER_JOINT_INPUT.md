# Joint input/output bounds for the randomized outer

This note concerns the unchanged GF65536 randomized-outer construction.
It extends the directed quiet-word coverage to every active-row count from
1508 through 8192. Together with the earlier receipts, quiet counts 1 through
962 are also covered. Quiet counts 963 through 1507 and nonquiet counts
6 through 7547 remain open. This is not a whole-code certificate.

## The input-weight mark

Let U be the potential input positions in a packet and let X be uniform
on the coordinate space supported by U. An actual outer word with r active
RS[256,128,129] rows has at least 129r nonzero input symbols. For v>=1,
retain the factor v^wt(X) in each packet moment and multiply the global
bound by v^(-129r). This restriction applies to actual outer words; the
comparison measure still permits zero values in potential positions.

For a Cauchy tile with a potential inputs, the map [I_a;F_U] generates an
MDS[a+8,a,9] code. If M(n,k,h) counts codewords on one prescribed support
of size h, its split enumerator is

    count(input weight i, output weight w)
        = binom(a,i) binom(8,w) M(a+8,a,i+w).

This gives the exact free tile moment. The quiet moment also imposes CX=0.
Expand the output tilt positively over sets of output coordinates required
to be zero. In a tile with a potential inputs and b such zero constraints,
the residual dimension is max(a-b,0). Distinguish zero and nonzero residual
input tiles before charging the input mark. For a nonempty pattern S of
nonzero tiles, a feedback rank lower bound rho gives at most

    q^(sum_(i in S) max(a_i-b_i,0)-rho)-1

nonzero kernel vectors. This can include vectors zero in some selected
tiles, so it is an upper count. Charge v^(sum_(i in S) a_i), add the all-zero
vector exactly once, and sum the positive output-mask weights. Take the
minimum with the exact free moment and the earlier coarse marked bound.
All operations before fractional powers are exact rational arithmetic.

## A small rank refinement for the existing maps

Residual three-dimensional tile spaces have distinct feedback hyperplanes.
To see this, write the residual polynomial space as R(X) times polynomials
of degree below three. Here R is monic of degree five, with roots in the
tile's additive GF256 coset. If its next coefficient is r4, the annihilator
of the highest-four-coefficient image has first normalized coordinates
(1,r4). The value r4 is the sum of five roots and lies in that tile's coset.
Distinct tile cosets therefore give distinct hyperplanes, whose sum has
dimension four. Changing the common polynomial basis preserves this fact.

The rank dynamic program may consequently join two three-dimensional
residual images at rank four. It also retains the earlier rank-two join
for two distinct one-dimensional images. No full additivity of residual
ranks is assumed for the current maps.

## Regional coefficient caps

Condition on all eight tile occupancies before taking the fractional power
alpha. Their exact multiplicities average the fractional quiet and free
moments. At a positive Bernoulli marker p, normalize these packet
polynomials to probability polynomials U and V. The ordinary region uses
U^128; the terminal region uses U^127 V. If every coefficient of both
polynomials is at most c, coefficient extraction gains a factor c^256.
The helper checks every coefficient with Arb; c is not a normal or
saddle-point approximation.

For fixed z, alpha, v, p, and c, the final logarithmic row-class bound is
an affine function of r minus 255 log2 binom(8192,r). It is convex.
Directed endpoint bounds therefore cover every integer in an interval.
Adding log2 of the interval length bounds its union.

## Directed receipts

Each row below includes the union over its stated active-row interval.
Rounded exponents are weakened outward. All witnesses use v=1 except the
separate input-mark diagnostic described afterward.

| Active rows | z | alpha | p | Uniform coefficient cap | Union log2 upper |
| --- | --- | --- | --- | --- | ---: |
| 1508--1644 | 1/6 | 1/16 | 146211182168/577720884983 | 3/256 | -75.82 |
| 1544--1802 | 1/7 | 1/16 | 172815637333/646891988433 | 3/256 | -73.01 |
| 1658--1998 | 1/9 | 1/16 | 30758323013/106196160857 | 1 | -94.98 |
| 1933--2466 | 1/16 | 1/16 | 325686709449/951889216970 | 1 | -75.49 |
| 2395--3070 | 1/32 | 1/16 | 143989657028/347265555841 | 1 | -90.53 |
| 2779--3709 | 1/64 | 1/16 | 193629185397/402500592083 | 1 | -111.07 |
| 3609--5215 | 1/256 | 1/8 | 723745044563/1033346037698 | 1 | -180.17 |
| 4318--6585 | 1/1024 | 1/8 | 775311481896/966340120213 | 1 | -244.44 |
| 5064--7574 | 1/4096 | 1/8 | 441447453054/502077788525 | 1 | -211.38 |

The previous dense quiet receipt covers 5865--8192. Overlaps in this
table are harmless for a union bound. Directed summation of the nine
unrounded interval receipts gives a total logarithmic bound below
-72.6143614417. The helper's `--receipt` command replays this union.

The marked witness z=1/4096, alpha=1/32, v=64,
p=605679345935/1099486265534 also certifies 3979--4609 below 2^-81.39.
It is superseded by the table's stronger unmarked intervals. A bounded
large-mark search (v=2^12 through 2^64, selected z and alpha) did not close
1024, 1280, or 1536. These are failed witnesses, not counterexamples.

## Reproduction

The exact helper is [joint_input](theory/hd_randomized_outer_joint_input.py);
its tests are [joint_input_tests](theory/hd_randomized_outer_joint_input_tests.py).
Call `directed_interval(first,last,z,alpha,v,p,coefficient_cap)` with the
rational witnesses above and Arb precision at least 192 bits. A cap of one
can be omitted. The method fails rather than accepts if any proposed
coefficient cap cannot be proved.

No native implementation or timing was changed. The next question is how
to retain more state-history information in the nonquiet composition.
