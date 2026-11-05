# Cheaper resets without weakening the certificate

The [native timing pass](RESET_TIMING_AND_GROUPS.md) now measures both
profiles and starts a weaker, routing-aware reset design. Reset-only Gray
timings are about6.0/6.5ms; the normalized variants did not win. The note
below records the algebra and operation counts, not measured speedups.

The two [certified profiles](IDEAL_RESET_LADDER.md) use uniform feedback
and globally MDS reset constraints. We can retain both properties without
applying a dense Vandermonde matrix at every batch. This is an algebraic
reduction in reset work, not a native timing result. Uniform feedback
remains unchanged; its eventual efficient implementation is a separate task.

The construction extends the [additive-coset approach](../ADDITIVE_MDS.md)
to batches 256 and 512. It stores novel-basis moments rather than values
on the preceding coset. This gives a fixed transpose transform followed
by a particularly simple state translation.

## Construction and coordinates

Let F=GF(2^128), J=GF(2^32) contained in F, b=2^h, and N=2^21.
The two profiles have h=9 and h=8. Fix a Cantor chain
gamma_0,...,gamma_31 in J, with gamma_0=1 and
gamma_i^2+gamma_i=gamma_(i-1). These vectors form a binary basis of J.
Write Delta(x)=x^2+x and s_i(x)=Delta^i(x). Then s_i is monic of degree
2^i, its kernel is V_i=span(gamma_0,...,gamma_(i-1)), and s_i(gamma_i)=1.

For j<b, let beta_j be the binary expansion of j in the first h chain
vectors. Define the degree-j polynomial

    p_j(x) = product_{i: bit_i(j)=1} s_i(x),

and let C(a)[j,k]=p_j(a+beta_k). Every p_j is monic of degree j.
Consequently these moment matrices differ from ordinary Vandermonde
matrices by one fixed invertible change of rows.

Put n=N/b. Set a_(-1)=0 and let a_t be the binary expansion of
Gray(t+1) in gamma_h,...,gamma_31, for 0<=t<n. All output cosets a_t+V_h
are disjoint. The predecessor coset is reserved for coordinates, not output.
Adjacent offsets differ by

    delta_t = a_t+a_(t-1) = gamma_(h+v_2(t+1)).

There are enough quotient directions: n<2^(32-h) for both profiles.
The route, outer, message length, and emitted coordinates are unchanged.

For a vector z of b moments, define L(d) by the following h layers.
In layer i, update each coordinate j having bit_i(j)=1:

    z_j <- z_j + s_i(d) z_(j-2^i).

Coordinates having that bit zero remain unchanged. The layers commute.
Expanding each p_j(x+d) gives

    C(a+d)=L(d)C(a),
    L(d)L(e)=L(d+e),       L(d)^2=I.

The lower-triangular direction matters: this acts on moments, not on
polynomial coefficients. It is not the transpose polynomial-shift map.

Sample independent uniform b-by-b matrices U_t over F, independently of
the route, once during setup. Start Z_0=0 and encode each routed input X_t by

    Y_t = X_t + U_t Z_t,
    Z_(t+1) = L(delta_t) Z_t + C(0) Y_t.

Emit Y_t and leave the terminal state free. No additional permutation is
introduced. In particular, the output is not changed to a transform-domain
representation: distance still counts the original F-valued output symbols.

## Why the existing certificates transfer

For a fixed output, define the conceptual syndrome

    S_t = sum_{u<t} C(a_u) Y_u.

Induction gives Z_t=L(a_(t-1))S_t. Both state changes are invertible, so
they preserve zero states. Closed segments therefore satisfy the same
globally distinct-point MDS equations as the ideal reset control. A support
of w columns has rank min(b,w), and a nonzero closed segment has weight
at least b+1. This holds over F, not merely over J.

The matrices C(0) and L(delta_t) are invertible. Thus nonzero output at
zero state activates it, while zero output cannot erase a nonzero state.
For a fixed input and exposed past, a prescribed reset from a nonzero
state still has probability |F|^-b. For a fixed output, the states are
determined without the feedback keys, preserving the inverse affine-fiber
argument. Changing deterministic state coordinates transforms full uniform
feedback matrices bijectively; independence across batches is preserved.

All finite counting expressions are consequently unchanged. The existing
80-/120-digit certificates transfer without new numerical bound tuning:

| Outer | Batch and state | Distance | Setup failure |
|---|---:|---:|---:|
| RS[256,128], GF256 coefficients | 512 | >49.1% | <2^-40 |
| RS[512,256], GF65536 coefficients | 256 | >49.1% | <2^-40 |

This transfer certifies the stated uniform-feedback ensemble, not the old
small-state candidate or the running-sum reset. The exact finite tests
below check the algebra; they do not replace the whole-code proof.

## Reset arithmetic and a second iteration

Compute C(0)Y with a transpose additive FFT. At a recursion node whose
offset is a and whose level is i, its transpose butterfly is

    low = left+right,
    high = s_i(a) low+right.

This needs one coefficient product and two XORs, skipping a zero product.
The state translation needs one coefficient product and one XOR per pair.
Each map has at most bh/2 product sites, versus b^2 for a dense map.
The transpose FFT is fixed, and all its coefficients lie in GF65536
for b512 or GF256 for b256.

The Gray translation is even more favorable than a generic shifted FFT:
s_i(delta_t)=gamma_(h+v_2(t+1)-i). Most coefficients lie in small subfields.
Only 32 batches in each profile contain any GF(2^32) translation coefficient.
The occasional wider coefficient is necessary in this schedule to leave
the small subfield's limited collection of cosets.

A further deterministic diagonal change of coordinates removes all
translation products on the most frequent Gray increment. Put
n_i=gamma_(h-i), d_j=product_{i:bit_i(j)=1} n_i^-1, and D=diag(d_j).
Store W=DZ instead. Its recurrence is

    Y_t = X_t + V_t W_t,
    W_(t+1) = D L(delta_t) D^-1 W_t + D C(0)Y_t,

where the V_t remain independent full uniform matrices. Layer i now uses
s_i(delta_t)/n_i. On half of the batches, delta_t=gamma_h and every layer
uses coefficient one: only XORs remain. Charge the added D scaling of the
reset transform; it is not free. This change also preserves both certificates.

Here are exact scalar coefficient-product counts for the complete N-symbol
reset stream. Products by zero and one are omitted. Each site acts on an
F-valued element with a fixed subfield coefficient. These counts exclude
feedback, the outer, routing, XORs, and memory costs; they are not CPU counts.

| Reset implementation | Original outer / state512 | Larger outer / state256 |
|---|---:|---:|
| Literal dense reset | 1,073,741,824 | 536,870,912 |
| Direct shifted transpose FFT | 9,437,184 | 8,388,608 |
| Fixed FFT + Gray translation | 16,781,312 | 14,688,256 |
| Same, normalized state | 14,147,584 | 12,566,528 |
| GF32 products, direct shifted FFT | 4,228,608 | 4,227,840 |
| GF32 products, either Gray version | 16,128 | 8,064 |

The Gray approach uses more product sites than the direct FFT, but replaces
almost all wide coefficients with smaller ones. Full normalization reduces
the Gray product count by about 15.7% and 14.4%, respectively. It also
*increases* GF65536 product sites by 524288 in each profile. This tradeoff
needs a native comparison, not a presumed speedup. Partial normalization
of only the GF256-valued levels retains the unnormalized GF65536 count;
it uses 15159296 and 13066240 total products, respectively.

The larger-outer profile has less reset arithmetic, but its outer uses a
larger transform and wider coefficients. These counts do not select an
end-to-end winner. The prior additive-coset implementation already used
structured transforms at smaller parameters; the dense baseline here is
the new ideal proof control, not an optimized historical timing.

## Structural shortcuts screened before bound tuning

Dropping C(0) from the recurrence is not a harmless coordinate change.
The highest-index unit vector e_(b-1) satisfies L(d)e_(b-1)=e_(b-1) for
every d. Emit this vector in two consecutive batches and the state returns
to zero. That gives a weight-two closed inner trajectory. It does not
establish an outer-compatible codeword, but it invalidates the MDS reset
interface used by both certificates. Uniform feedback does not fix that
deterministic output-side identity.

More generally, for an update z'=Az+y with identity injection, choose
consecutive outputs v and Av. They form a closed segment. Retaining the
MDS interface requires wt(v)+wt(Av)>=b+1 for every nonzero v. In particular,
every column of A must have full weight. Sparse state transitions or
coordinate permutations fail this test before any routing-bound tuning.
The fixed C(0) injection is what allows a triangular translation in our
coordinates without introducing those low-weight output resets.

Making the state basis especially simple for one translation can likewise
move work into the injection. A basis adapted to that translation is not
automatically the degree-graded basis of the batch's evaluation subspace.
We retain and count the actual injection rather than assuming it stays free.

## Checks and next test

Run `python -B research/workstreams/one_route_coupled_scan/coded_state/cheap_reset_maps.py`.
The source has a 45-second cap, uses small in-memory vectors, and writes no
files. It verifies the degree32 modulus, constructs a full Cantor chain,
and checks 14 transform/translation fixtures, including dimensions256/512.
It also checks58 stream-state identities,1883 independent small-field rank
fixtures,32 coupled encoder/inverse streams,32 shortcut counterexamples,
and4 normalized-coordinate identities. Small cases use independent dense
polynomial evaluation; actual-size cases check selected columns directly.

No frozen certificate helper, implementation, paper, or prior receipt was
changed. The next discriminating test is a native reset-only comparison
of direct FFT, Gray, and normalized Gray on both profiles, followed by an
outer-plus-route comparison. Keep uniform feedback as the proof control.
Do not weaken reset rank or tune another numerical bound until these
proof-preserving arithmetic options have been measured.
