# Cheap state refresh and a larger certified input-weight class

**Current status:** this profile is parked after the [realized-route
window obstruction](WINDOW_DIMENSION_SEARCH.md). The class bounds below
remain valid, but the missing sparse class contains verified bad routes
independently of the four-state inner coefficients. They must not be
presented as a nearly complete whole-code certificate.

The new higher-rate candidate uses a Cauchy-form emission map and four
fresh nonzero GF256 state scalars per packet. Its local rank bound is nearly
as strong as a full random state matrix. The fixed maps also have cheaper
circuits than the prior cubic FFT implementation.

This is a changed construction, not an optimization of the same sampled
code. Its target remains rate 11/16, 10% symbol distance, and forty bits of
whole-code setup-failure margin. The whole-code claim is still open.

The new directed replay covers every inverse-inner input of weight at
least 0.11N. The prior once-per-region candidate covered weights at least
0.19N. Inputs below 0.11N are the remaining proof problem for this candidate.

## Construction and cost

Use the fixed GRS[32,22] outer over AES GF256 and one regional routing stage:
each of the 32 outer-coordinate regions independently permutes its rows.
Every routed coordinate receives an independent nonzero GF256 scalar.
These setup components are independent of the state diagonals below.

The [circuit note](CAUCHY_DIAGONAL_PACKET_COST.md) specifies fixed matrices
A of shape 16 by 4 and C of shape 4 by 16 for each of four phases. Their
properties needed here are:

- every square minor of A is nonzero;
- every set of at most four columns of C is independent;
- CA=0;
- for eta=2, every matrix I+eta C_J A_J is invertible.

For packet input X and four-symbol state S, sample a diagonal D_s with
four independent entries from GF256 without zero. The recurrence is

    R=D_s S,   Y=X+eta A R,   S'=R+C X.

The initial state is zero and the final state is discarded. No random
region-entry matrix is needed for the result here. The matrices vary
through the prescribed four phases; the residue-coordinate state does
not require a phase-basis conversion. Setup is sampled once and fixed
for every message.

The causal inner is invertible: earlier recovered inputs determine the
current state and hence the current input from Y. It is independent of
the outer route and coordinate diagonals, as required by the existing
dense-class theorem.

The direct four-point circuits give a conservative full product volume
of 169/88 GFNI512 equivalents per source record, versus 181/88 for the
previous eta-scaled candidate. This 6.63% reduction includes the outer,
input diagonals, and new state diagonals. It is not a runtime claim;
packing, XORs, routing, and memory traffic still need native measurement.

## Why diagonal randomness is sufficient for the root equations

Fix a common zero-equation history and a positive-rank state-form map.
Its image is an arbitrary r-dimensional subspace U of GF(q)^s, where
q=256 and s=4. Fix z zero-input, zero-output positions. Their emission
matrix A_Z is superregular. The fresh state diagonal D is independent
of U and of every previously imposed equation.

Let H be the dimension of the surviving image:

    H = dim(U intersect ker(A_Z D)).

This notation identifies the pre-diagonal coordinates of the image.
The state remains nonzero only when H>0. Imposing these root equations
adds r-H independent equations to the common span.

For any fixed nonzero u in U, write h for its coordinate support size.
If h<=z, the relevant columns of A_Z are independent, so A_Z D u cannot
vanish. If h>z, fix all but z diagonal entries on that support. Any z
remaining support columns form an invertible matrix. There is at most
one assignment to those entries that solves A_Z D u=0. Therefore

    Pr_D[A_Z D u=0] <= (q-1)^(-z),
    E_D[q^H-1] <= (q^r-1)/(q-1)^z.

The expectation counts nonzero vectors only in the at-most-four-dimensional
state image. It is not a first moment over all messages of the outer code.

For 1<=tau<=q-1, define a polygon by the points

    (0,0),   (q^h-1, tau^(h-r))  for 1<=h<=min(r,s-z).

Its first slope is tau^(1-r)/(q-1). Later successive slopes are

    (tau-1) tau^(h-r) / ((q-1) q^h).

They decrease by a factor tau/q. The polygon is therefore increasing
and concave, including tau=1. Jensen's inequality bounds the expected
alive rank potential by its value at

    min((q^r-1)/(q-1)^z, q^min(r,s-z)-1).

Call this value G(r,z,tau). When z=0 use G=1; when z>=s use G=0.
Set Gamma(z,tau)=max_r G(r,z,tau) for 1<=r<=4.
This is uniform over the entire preceding equation history. It does
not assert a uniform distribution on subspaces after diagonal refresh.

For example, at tau=14 the factors Gamma for z=1,2,3 are approximately
0.0714570, 0.00510510, and 0.000364431. The corresponding exact-GL bounds
are 0.0714427, 0.00510305, and 0.000364431. The implementation retains
the exact rational upper bound, not these displayed rounded values.

## Complete support transfer

Fix input and output supports I,J in one packet. Let i=|I minus J|,
a=|I intersect J|, o=|J minus I|, and z=16-i-a-o. On a nonzero-state
branch, four root positions kill the branch by the emission MDS property.
For z<4, the fresh state diagonal contributes Gamma(z,tau).

Next impose the i zero-output equations at active inputs. Their independent
coordinate diagonals each contribute tau^-1 by the alive diagonal lemma.
No state-diagonal randomness is used again in these steps.

For a zero-state checkpoint, the remaining input diagonals give the
independent pivot charge tau^-min(a,4). Without overlap, the invertible
matrix I+eta C_J A_J forbids return. A return with exactly one nonzero
output has all sixteen inputs nonzero: the input/output pairing is the
same cubic pairing as in the eta-scaled reference.

Thus the nonzero-to-nonzero transfer sums terms

    multinomial(16;z,i,a,o) x^(i+a) y^(o+a)
      * Gamma(z,tau) * tau^(-i).

The nonzero-to-zero transfer uses the same terms multiplied by
tau^-min(a,4), with a>0 and nonempty input/output, excluding output
weight one unless input weight is sixteen. Dropping output-nonzero tests
only increases the bound.

From a zero state, I=J. A nonempty zero-to-zero packet has at least five
inputs and receives tau^-4. Births receive no further rank charge.
These are the same two zero-state transfer entries as in
[the return-pivot certificate](INTERMEDIATE_RETURN_CERTIFICATE.md).

All conditional input-support branches are included. Multiplying the
positive two-state transfers therefore covers multiple excursions,
different row supports, and excursions crossing region boundaries.
Unknown formal state ranks are bounded by Gamma before composition;
they are not treated as a Markov chain with an assumed exact law.

The outer-route averaging expression is unchanged. For fixed total input
weight W, alpha=W/N, positive lambda, x>0, and 0<y<=1, define

    A(alpha)=(1-alpha)^32 + tau^-10 (1+(tau*lambda-1)alpha)^32.

The probability of a bad word in that class is at most

    (L+1)^32 A(alpha)^L e_0 T^(N/16) 1
      / (tau (x*lambda)^W y^e),

where e=ceil(N/10)-1. The return-pivot note proves this profile-independent
route bound. Its conditional diagonal argument is uniform over the new
state randomness as well. The moment lemma above supplies the additional
root charge that distinguishes this construction.

## Directed finite result and remaining work

Let L=64 ceil(2^ell/(22*64)), K=22L, and N=32L. The replay uses frozen
rational witnesses, 80-digit outward arithmetic, positive matrix powers,
and 1,300 density cells. Integer endpoints cover each listed weight exactly
once. All probabilities refer to ideal independently sampled setup.

| Requested size | Actual K | N | Covered W | Class margin at least |
|---|---:|---:|---:|---:|
| K16 | 66176 | 96256 | 10589 through 23101 | 202 bits |
| K18 | 263296 | 382976 | 42128 through 91914 | 1833 bits |
| K20 | 1048960 | 1525760 | 167834 through 366182 | 8535 bits |

The existing dense theorem applies to this independent invertible inner
and covers W>=ceil(.24N), with margins 578/3381/14745. A conservative union
bound therefore covers every W>=ceil(.11N), with margins 201/1832/8534.
These are class margins, **not the code's distance margin**.

Within the lower class, the separate [exact one-row
analysis](CAUCHY_FIRSTMOMENT_GATE.md) now has directed margins
107/105/103 bits at K16/K18/K20. It uses the exact GRS support counts and
the joint state/input-diagonal return probability, not the older
constant-state certificate. The new [two-row certificate](CAUCHY_Q2_TRANSFER.md)
adds 65/72/77 bits for exactly two active rows. The remaining gap consists
of at least three active outer rows with 0<W<.11N. Pointwise floating screens at input
densities .001, .01, .04, and .08 still give positive exponents. A separate
multirow first-moment screen also fails. Joint return bounds for two or
three inputs and a route bound preserving regional activity are next;
more margin in the already covered class is not the target.

The subsequent [joint return analysis](CAUCHY_JOINT_RETURNS.md) supplies
an audited selected-zero rank action for two and three packet inputs.
This action is distinct from its exact fixed-word probability formula.
Combining it with exact regional counts gives promising three/four-row
screens, but has not closed the larger sparse classes or changed the
whole-code status above.

The rank-moment tests exhaust every subspace of GF8^3 and every nonzero
diagonal: 50,421 subspace/diagonal pairs and 882 exact potential comparisons.
An independent reviewer checked the moment argument and Jensen step and
replayed the finite K16 result. Separate tests check the positive transfer
against an independently indexed exact sum and verify all 1,300 K16 cells.
The map reference passes 20,480 scalar/basis comparisons and 128 adjoints.

```sh
python -B research/workstreams/char2_field_spin/theory/diagonal_refresh_rank.py --check --compare
python -B research/workstreams/char2_field_spin/theory/diagonal_refresh_certificate.py
python -B research/workstreams/char2_field_spin/theory/cauchy_q1_certificate.py
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_diagonal_refresh.py
python -B research/workstreams/char2_field_spin/theory/cauchy_diagonal_packet_reference.py --check
```

All coefficients lie in GF256. A completed GF256 distance theorem would
extend to GF(2^m), m in {8,32,64,128}, by scalar extension without changing
the setup-failure margin. No such whole-code theorem or new runtime is
claimed here. Binary production and earlier proved references are unchanged.
