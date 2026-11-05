# Pre-mixer Toeplitz pairs: comparison and next fractional gate

Replacing the input diagonal by Toeplitz2 pairs has a valid positive
comparison with the current RS counting interface. The five completed
scalar-mark controls show only small changes. They neither certify nor
reject this candidate: the unresolved proof loss concerns powered masses
of rare labelled masks, which these means do not measure.

No native code, regional proof bank, optimizer, or whole-code proof grid
was run. The next gate below is derived but not implemented here.

## Construction and comparison interface

Keep sixteen RS[16,8] rows and the current affine two-sum matrix H. At each
of sixteen columns, apply eight independently sampled, uniformly
invertible Toeplitz2 matrices on a fixed matching M, then H, then the
existing independent nonzero output diagonal. The matching may follow
the proposed XOR schedule, but its entire schedule is fixed before
sampling the keys. A constituent row set A stays fixed across all columns.

The circuit has 24 variable products for the pair maps, 15 for the
normalized output diagonal, and eight fixed products for H: 47 sites
rather than the current 39. The normalization fixes one output-diagonal
entry to one and absorbs its nonzero common scalar into all eight
Toeplitz matrices. Scalar multiplication is a bijection of the accepted
Toeplitz ensemble, so this is the same map distribution. These are
circuit counts, not a native performance prediction.

The exact pair law and its amplitude caveat are established in
[the post-mixer note](HD_TWO_SUM_TOEPLITZ_PAIRS.md). Write J_t(b) for the
safe iid-NZ majorant's mass on output support b, given input support t.
The corresponding per-vector majorant divides this mass by
255^|b|. In particular, J depends only on the input support, not its
nonzero amplitudes. Its mass can exceed one when t is a double.

This support dependence supplies the missing RS measure bridge. For
every independently nonzero diagonal D and every input x,

\[
 J(Dx,\cdot)=J(x,\cdot),\qquad
 \Pr[T x\in\cdot]\le J(x,\cdot).
\]

The second inequality is pointwise and therefore remains valid after
any positive continuation, including H and the retained output diagonal.
Insert an artificial independent input diagonal before J; the equality
shows that it changes no comparison value. For a fixed set A of a
nonzero RS rows, the established shortened-MDS argument then bounds the
expected count at every prescribed sixteen-column array by 255^(-8a).
Equivalently, a row of prescribed support size d has at most
255^(d-8) matching codewords and each matches its nonzero amplitudes with
probability 255^(-d). Independence across the rows gives the product.

Thus the actual candidate is dominated by the positive J comparison
with the same whole-word factor 255^(-8a). This does not assert that the
native candidate retains an input diagonal. The diagonal exists only in
the comparison proof. It preserves input support, so each active row's
nine-nonzero floor, or its aggregate floor 9a, may still be retained.

The helper's `exact-support` mode has a narrower meaning. Averaging a
uniform input on a coordinate space makes the exact pair output flat
within each support. This gives an exact uniform-input diagnostic, but
does not establish the diagonal-invariant bridge for actual RS words.
Use `iid-majorant` for that bridge. After H, the actual independent output
diagonal again supplies independent nonzero amplitudes conditional on
the final support.

## Completed bounded control

The control fixes A to its first a coordinates, uses XOR matching eight,
and evaluates the scalar output-weight marker z=1/2 with input marker
v=1. It reports log2(raw mass)-8a, the normalization for a uniform column
in F_256^A. This is not the whole-RS normalization or a failure bound.

| a | Existing diagonal | Post-H iid majorant | Pre-H exact uniform diagnostic | Pre-H iid majorant |
|---:|---:|---:|---:|---:|
| 1 | -7.994397 | -7.993854 | -7.994049 | -7.994049 |
| 2 | -14.955787 | -14.931756 | -14.911801 | -14.911801 |
| 4 | -15.828007 | -15.860782 | -15.905076 | -15.905076 |
| 8 | -15.905346 | -15.866053 | -15.910007 | -15.910007 |
| 16 | -15.910007 | -15.866053 | -15.910007 | -15.865709 |

For a<=8 this matching pairs every selected coordinate with an
unselected coordinate. There are no double-active input pairs, so the
iid majorant and exact diagnostic coincide. A different matching can
have double-active pairs and incur the majorant's extra mass. The table
must not be extrapolated to the entire matching schedule. Full arity's
exact uniform law is unchanged; its safe comparison pays a small tax.

The safe pre-H gains over the existing diagonal are approximately
-0.00035, -0.04399, +0.07707, +0.00466, and -0.04430 bits. A scalar
moment sums all masks before reporting its value. It cannot determine
whether the replacement substantially reduces the fractional contribution
of a rare quiet mask. In particular, this is not a negative construction
result.

## Next gate: compose each labelled mask before fractionalization

Fix A, M, an output mask T after H, and an input-support marker v>=1.
For a pair with local active-set mask a, define the raw comparison mass
on its intermediate support b by

\[
 L_a(b;v)=\sum_{t\subseteq a}255^{|t|}v^{|t|}J_t(b).
\]

The product over matching pairs is the raw mass L_A,M(B;v) of the full
intermediate support B. Conditional on B, the comparison amplitudes are
independent and uniform nonzero. Let N_NZ(k,l,w) be the exact number of
vectors with a prescribed input support B and H-output support T, where
k=|B|, l=|B intersect T|, and w=|T|. The existing closed kernel and input
Möbius inversion compute this nonnegative integer; its dependence on
these three sizes is a proved property of this H, not an assumed
permutation symmetry of H itself.

The required positive mask mass is

\[
 C_{A,M}(T;v)=
 \sum_B\frac{L_{A,M}(B;v)}{255^{|B|}}
 N_{\mathrm{NZ}}(|B|,|B\cap T|,|T|).
\]

There is no subtraction of upper bounds. Compute C using an eight-pair
polynomial. A pair e contributes

\[
 P_e(X,Y)=\sum_{b\subseteq e}
 \frac{L_{A\cap e}(b;v)}{255^{|b|}}
 X^{|b|}Y^{|b\cap T|}.
\]

If d_k,l is the coefficient of X^k Y^l in the product of these eight
polynomials, then C=Sum_k,l d_k,l N_NZ(k,l,|T|). All intermediate values
are nonnegative and the coefficient table has at most 17 by 17 entries.

For the requested first diagnostic, use alpha=1/5 and form

\[
 F_{A,M}(w)=\sum_{|T|=w}C_{A,M}(T;1)^{16\alpha}.
\]

The full composition C must precede the 16alpha power. Powering the
individual B terms first would introduce an avoidable subadditivity
loss. The target-mask multiplicity stays outside that power. Compare
these coefficients with the old fixed-A coefficients from N_A(T), not
merely their scalar first moments.

This diagnostic can avoid enumeration of all 65536 output masks. If M
contains u pairs wholly in A, c=a-2u crossing pairs, and h=8-a+u pairs
outside A, use three target-occupancy types on each internal or outside
pair, and four labelled target types on each crossing pair. There are

\[
 {u+2\choose2}{h+2\choose2}{c+3\choose3}
\]

joint types. Each has its exact multinomial multiplicity, with an
additional factor two for every singleton target in an internal or
outside pair. On a crossing pair, distinguish the selected-A leg from
the other leg. For a=8, the counts for u=0,1,2,3,4 are respectively
165, 756, 1260, 1000, and 225. These counts classify T for one fixed A;
they do not include another constituent-set factor.

Keep the same A across the sixteen matching columns. Per-column
statistics depending on u do not authorize choosing a new A, or its
worst u, independently in each column without charging that relaxation.
Likewise, F(w) is a scalar fractional diagnostic: a later inner-operator
controller must retain any required labelled-mask information or prove
its own averaging step. The whole-word normalization, RS floors, source
counts, projective factor, and class union still belong at their
established stages. If v is used to enforce total input weight>=9a,
its whole-word charge is v^(-9a), then alpha-powered with the measure.

The useful next decision is whether this fixed-A, mask-resolved
fractional comparison materially improves the quiet/dense class while
preserving a manageable comparison tax. Only then is a regional proof
controller justified. The unchanged final target is a whole-code failure
bound at most 2^-40 for distance35%; none is claimed here.

## Tests, reproduction, and pins

Six new tests plus the nine pair-law tests pass. They check exhaustive
GF4 uniform-input output atoms, pointwise domination followed by an
arbitrary positive H observable, marked normalization, a positive
coordinate-space mixture, full-arity uniformity, diagonal controls, and
scope guards. The five-column gate completed in under one second.

From `theory`:

```text
python -B -m unittest hd_two_sum_toeplitz_pairs_tests hd_two_sum_pre_toeplitz_pairs_tests
python -B hd_two_sum_pre_toeplitz_pairs.py --output NEW_OUTPUT.json
```

SHA-256 pins:

- `hd_two_sum_pre_toeplitz_pairs.py`: `58ddaa9a63c0687bb720809b509f1749d7b08a7784557f46b2d7a7c5613a60c8`
- `hd_two_sum_pre_toeplitz_pairs_tests.py`: `4f90c51194120f26dbc3e9231d94ab074c2440304f8a123f9d6b53f4a88793f3`
- `hd_two_sum_pre_toeplitz_pairs_local_gate.json`: `2bd8500bee4549d2f79706f76a69d7620fbfd86004b2046bb1c7e32732ac0e94`

The displayed logarithms are binary64 evaluations of exact rational
local calculations, not directed distance receipts. All old proof
providers, native sources, and frozen controls remain unchanged.
