# Exact forced-birth coefficients: a bounded checkpoint

Extracting exactly four forced births leaves the tested global bound
strongly positive. The new central-death-interval exponent is
`+674889.164`, about 155 bits below the recorded evaluated-birth-marker
control. This is a FLOAT diagnostic, not a distance certificate or a
counterexample. The construction and full 45%/40-bit goal are unchanged.
Separately, the [selector-capacity argument](QUARTET_SELECTOR_CAPACITY.md)
closes this entire hard profile with a conservative bound below 2^-510.
The failure recorded here belongs to the moment control.

The helper is
[quartet_actual_birth_coefficients.py](quartet_actual_birth_coefficients.py).
It uses the positive actual-state/MDS count from
[the joint-count note](QUARTET_ACTUAL_MDS_JOINT.md) and the scalar field
correction in
[quartet_actual_death_profiles.py](quartet_actual_death_profiles.py).

## Objects retained by the count

The setup is sampled once and then fixed for every message. Fix a
complete candidate actual support history S before averaging its
quartet partitions and routing labels. Start with S_0 empty and leave
S_4096 free; there is no terminal reset or terminal expansion.
For update t, let Gamma_t(S_(t-1)) be the union of quartets meeting the
old support. Define

    V_t = Gamma_t(S_(t-1)) union S_t,
    B_t = S_t \ Gamma_t(S_(t-1)),
    Ls = sum_t |S_t|,       f = sum_t |B_t|,
    F = sum_t |V_t|,        G = F-f,
    D = sum_t |S_(t-1) \ S_t|.

Here f counts forced births outside the old-hit union, not all coordinate
births. These are actual-state statistics, not availability statistics.
For a candidate eligible-row set H of size h, write
C_out=R-128h. The regional count restores exactly h eligible labels in
each 8192-position region. It still relaxes other conditions needed for
an actual codeword, as described in the joint-count note.

Use positive markers s, d, z, nu, eta for Ls, D, G, f, C_out,
respectively. In the helper, `log_t` means log2(d); it is independent of
the field-bound parameter x and fractional power theta.

## Extract births without enlarging the population state

Let P_p(X,Y;nu) be the joint-count note's quartet polynomial, evaluated
at death marker d and Bernoulli eligible probability p. Its old-empty
part is (1+p*s*nu*Y)^4. Its old-live part is independent of nu and
retains the necessary MDS condition

    |H_local^c \ S_new| <= |S_old|-1.

Put P0=P_p(X,Y;0), so P0 has constant old-empty row 1. Then

    P_p = P0 + ((1+p*s*nu*Y)^4-1).

For width n divisible by four, put g=n/4 and
c_(b,j)=[w^j]((1+w)^4-1)^b. The bare one-update birth coefficient is

    K_(p,j)(u,v)
      = 1/binom(n,u) * sum_(b=0)^min(g,j) binom(g,b)c_(b,j)(p*s)^j
                         [X^u Y^(v-j)] P0^(g-b).          (1)

Coefficients outside the polynomial support are zero. Formula (1)
follows by choosing the b quartets contributing a nonconstant
old-empty term and then extracting nu^j. For n=512 and j<=4, a single
positive convolution pass need save only P0^124,...,P0^128. The physical
state remains the 513 possible old populations u.

Let Z=e_0 e_0^T. Truncated polynomial matrices are represented as
Z+sum_(j=0)^4 B_j nu^j. Under squaring, their positive coefficients obey

    B'_j = Z B_j + B_j Z + sum_(a=0)^j B_a B_(j-a).       (2)

Four squarings compose the 16 updates in one region. Degree-zero B_0
has source row zero identically zero: a path starting empty cannot
become nonempty without a forced birth. Positive-degree B_j(0,0)
entries must remain, since they count excursions that return to empty.

For each p, divide these regional positive coefficients by

    pi_h(p) = binom(8192,h)p^h(1-p)^(8192-h).

Only then take their entrywise minimum over p, separately for each
degree, in the common actual-population basis. Denote the resulting
matrices by C_j and set R(nu)=Z+sum_j C_j nu^j. To justify this minimum,
retain a formal eligible-label counter in every local polynomial,
including uncovered positions. Extracting the exact regional quota
and dividing by pi_h(p) yields the same conditional physical operator
for every p. Each adjusted positive matrix therefore bounds that
operator coefficientwise. The exact empty operator Z is not divided
by pi_h(p).

The global bare coefficient is bounded by

    M_j = [nu^j] e_0 R(nu)^256 1.                         (3)

The implementation propagates five row vectors through the 256 regions,
using v'_j=v_j Z+sum_(a=0)^j v_(j-a)C_a. This avoids global polynomial
matrix powers. The endpoint vector is all ones, and M_0=1 exactly.
Truncation at degree four is exact for these coefficients because all
birth degrees are nonnegative.

## Attach a field bound to the selected death interval

Let q=2^32, W=943718, 0<x<=1, and 0<=theta<=1. The actual-state field
bound for a fixed candidate S, after dropping its terminal factor at
most one, gives a field contribution P_S satisfying

    P_S <= min(1, Bbar_x t_x^D),
    Bbar_x = 2^-31 q^-R x^-W (1+qx)^Ls,
    t_x = q*x^2/(1+qx).

For an integer death interval I, define

    c_I = max_(D in I) {
             theta*min(0, log2(Bbar_x)+D*log2(t_x))
             -D*log2(d) }.

The corresponding fixed-profile exponent is bounded by

    log2(M_f) + log2 binom(8192,h)
      -Ls*log2(s) -(F-f)*log2(z) -(R-128h)*log2(eta) + c_I. (4)

The H-set union occurs once. Regional quota divisors are already inside
M_f, and no separate routing factor is used. Since M_f is the bare
coefficient [nu^f], (4) contains no nu^-f factor. The maximum over D
precedes summation over histories, so there is no D-union factor within
I. Bounds on disjoint intervals must be summed. A global union over
different profiles has not been applied to (4).

## Checks and the recorded control

Parent-run checks passed, including 530 exact width-4/8 physical birth
entries, all 35 width-8 quartet partitions, 255 regional checks,
eligible-quota conditioning independent of p, and positive birth-return
excursions. Floating matrices and global coefficients agree with the
independent rational references. Degree-zero propagation remains exactly
e_0; optional upward floors preserve its structural zeros. Independent
static audit passed the coefficient formula, both Z cross terms,
quota/minimum order, scale handling, and free endpoint.

The single three-tilt probe completed in 95.79 seconds under the
175-second, 2-GiB guard. It used

    (h,R,Ls,f,F) = (4096,786432,786432,4,1048576),
    x=2^-16, theta=1/16, s=2, d=1/2, z=1, eta=1/4,
    p in {1/5,1/3,1/2}, I=[34570,196607].

| Quantity | FLOAT log2 value |
|---|---:|
| Bare coefficient M_4 | 771563.0618 |
| Bare coefficient M_1 | 771518.595 |
| Exact-f=4 exponent (4), central interval | +674889.164 |
| Earlier evaluated-nu control, nu=2^-18 | +675044.0045 |

The corrected exponents differ by about 155 bits. The earlier control
uses its matched death marker and broader death scope, so this is not
an identical-marker ablation. Bare coefficient moments must not be
compared directly with a moment evaluated at nu=2^-18. The large M_1
is a diagnostic of the relaxed count; f=1 is itself incompatible with
this fixed F because F-f must be divisible by four.

Separate exact arguments cover the low and high death tails at this
profile. The present run targets the central interval and does not
close it by this moment bound. The separate selector-capacity argument
now covers the whole profile simultaneously over all death counts,
without a state-history union. Ordinary floating arithmetic is not a
directed certificate, even when the optional positive floor is used.

The result shows that removing histories with many forced activations
does not repair this control. Birth-free persistence after a small
initial activation remains, consistent with the independently proved
[four-birth obstruction](QUARTET_LOW_BIRTH_MOMENT_OBSTRUCTION.md).
That obstruction and the comparison moments are not bad codewords.
More birth-degree features are now parked. The next proof step is to
extend the selector-capacity argument to the other eligible-row counts
and profiles. The full goal remains open.
