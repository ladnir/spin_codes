# Composing the local refinements

This note records two local refinements of the randomized-outer construction.
They do not change the encoder. The composed floating screen still fails to
certify the selected middle occupancies. The successful sparse class bounds
are recorded separately in the sparse-ladder note.

Use the field, maps, and three-component state envelope from
[HD_RANDOMIZED_OUTER_STATE_PROFILE.md](HD_RANDOMIZED_OUTER_STATE_PROFILE.md).
Here q=65536, D=q^4-1, epsilon=1/q, and c=1-epsilon. Packet occupancy j counts
potential input positions under the artificial uniform measure, not necessarily
nonzero inputs. The old state is independent of the current packet's positions
and amplitudes, conditional on the full packet-count path.

## Weighted returns of a residual state

Fix a nonzero old state h and a potential input subset U of size j. Write
r=min(j,4). The feedback C_U has rank r. A return to zero requires

    C_U X_U = M h.

For any target in the image of C_U, this event has probability q^-r. Conditional
on that event, X_U is uniform on an affine coset of ker(C_U). The emitted word
is consequently uniform on an affine coset of F_U ker(C_U). For 0<z<=1,
the weight moment of any such coset is at most that of the zero coset.
This follows from the nonnegative Fourier coefficients of z^wt.

Let K_U(z)=E_X[1[C_U X_U=0] z^wt(F_U X_U)]. The refresh kernel sends h to a
mixture of its scalar line, with total mass epsilon, and the uniform nonzero
state, with mass c. Thus the weighted return is at most

    K_U(z) [epsilon + c (q^r-1)/D].

Averaging U preserves this inequality because the bracket depends only on j.
Therefore any proved upper bound for the quiet moment K_j supplies an upper
bound for this residual-to-zero entry after multiplication by the bracket.
This uses a conditional affine-coset argument, not independence between output
weight and state return.

There is a second valid bound. If U meets u tiles, the nonzero old state's
evaluation polynomial has at least max(0,64-8u-3) nonzero outputs on untouched
tiles. This is a pointwise property even when a return occurs. The weighted
return is therefore also at most

    E_U[z^max(0,64-8u-3)] q^-r [epsilon+c(q^r-1)/D].

The minimum of these two expressions and the previous valid return bounds
can be used entrywise. No minimum changes the state's probability space.

## Exact small-column spectra

The dominant nonuniform refresh component is the identity branch. For its
flat-state column bound, the coupled linear map is G=F+A C. The identities
F^2=I, F A=A T, C F=T C, and C A=0 imply G^2=I. In particular, G is invertible.

The fixed G matrix has coefficients in GF256. All 64 columns have weight 64.
For each pair of columns a,b, projective representatives over GF65536 are
(1,t), for t in GF65536, and (0,1). Any cancellation in a coordinate with
b_i nonzero requires t=a_i/b_i, which belongs to GF256. Slopes outside that
subfield have the full union support of a and b. Counting the 256 subfield
slopes therefore gives the exact pair spectrum without enumerating GF65536.

The summed two-column spectra have nonzero weights 60 through 64. Their
integer counts sum to binom(64,2) q^2, including one zero word per column pair.
For a fixed target state, the identity-branch output lies in an affine coset
of this column code. The zero-coset moment bounds that column by the same
positive-Fourier argument. For a zero target, subtract only the known exact
zero-input mass q^-j.

The exact refresh decomposition is P=a I+b L+c J, where

    a=(q^3-1)/D,  b=(q-1)/(q D),  c=1-1/q.

L averages the nonzero scalar line and J averages all nonzero states.
Replacing epsilon gamma by a gamma_I+b gamma_L retains the smaller
identity-column spectrum without assuming that every scalar-line branch
has that spectrum. The implementation takes the minimum with the old
complete entry; it never subtracts a previously computed upper bound.

At the tested weak tilt near the first uncovered sparse class, this exact
column spectrum changes the failure exponent by less than 0.001 bits.
It is useful local information, but it does not resolve that class.

## Composed diagnostic

The floating provider combines the quiet moment, the
[occupied-tile flat moment](HD_RANDOMIZED_OUTER_FLAT_MOMENT.md), its D V_j
caps, and the two weighted-return bounds above. Each changed entry is the
minimum of established upper expressions. The exact small-column spectrum
is tested separately; the following diagnostic does not include it.

For each occupancy, the screen searches eleven weight tilts and seven
fractional exponents. A continuous marker search replaces the earlier
four-choice marker grid. The selected one-packet expressions remain positive:

| Active rows | Floating log2 upper expression |
| ---: | ---: |
| 16 | 2905.82 |
| 128 | 17134.91 |
| 512 | 41724.10 |
| 1024 | 53825.34 |
| 2048 | 68185.53 |
| 4096 | 58544.68 |
| 6144 | 42304.53 |

These are failed witnesses, not counterexamples and not optimized impossibility
claims. They show that the improved local moments do not by themselves close
this bounded search. Retaining more layout information before fractional
powers and restoring the outer's minimum-support restriction are separate
proof improvements worth testing next.

## Reproduction

The diagnostic helpers and their eight tests are
[column_moment](theory/hd_randomized_outer_column_moment.py),
[column_moment_tests](theory/hd_randomized_outer_column_moment_tests.py),
[refined_profile](theory/hd_randomized_outer_refined_profile.py), and
[refined_profile_tests](theory/hd_randomized_outer_refined_profile_tests.py).
Run from the theory directory:

    python -B -m unittest hd_randomized_outer_column_moment_tests hd_randomized_outer_refined_profile_tests -v
    python -B hd_randomized_outer_refined_profile.py --optimize-marker

The tests check exact integer normalization, direct GF256 projective counts,
invertibility, unchanged structural zeros, and entrywise non-regression.
The floating diagnostic is explicitly not a directed certificate.
