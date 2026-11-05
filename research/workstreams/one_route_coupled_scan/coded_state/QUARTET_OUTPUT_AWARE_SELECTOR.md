# Output weight in the live-selector quotient

The actual-volume quotient can retain output weight through a positive
two-state run transfer. Its spectral envelope preserves a closed-run
penalty. Dropping that penalty gives a useful negative control: with
the current projective prefactor, this envelope can improve the remaining
R=L profile by at most 30 bits and cannot improve the R<L profile.
Thus a volume-only output-marker search cannot close the recorded gaps.
The construction and full-code goal are unchanged.

This note uses the candidate selectors, births, and field-count order in
[the live-selector quotient](QUARTET_LIVE_SELECTOR_QUOTIENT.md).
The exact algebra and small-field checks are in
[quartet_output_selector_envelope.py](quartet_output_selector_envelope.py).
No new global geometry bound is claimed here.

## Keep actual volume in the run count

Fix routing, partitions, a candidate pair (U,B), and actual state volume
L. Write f=|B|, m=L-f, and A0=sum_t|U_t intersect U_(t+1)|, with the
formal free endpoint U_(T+1)=[512]. The allowed state positions are
A=(U_t intersect U_(t+1)) union B_t, so |A|=A0+f.
Let C_A be the number of closed coordinate runs of A.

Fix every declared birth value to a nonzero field element. Mark a
nonzero state by y>0 and a nonzero output by 0<x<=1. The transfer between
zero and nonzero state types is

    M = [[1, (q-1)xy],
         [x, y(1+(q-2)x)]].                            (1)

For a run of length ell, its closed weighted count is e0 M^ell (1,x)^T
if unforced, and xy e1 M^(ell-1)(1,x)^T if its first value is fixed
nonzero. Replace the final vector by (1,1)^T for an open run. The first
forced value is fixed, so its factor is xy, not (q-1)xy.

These counts are independent of the specified nonzero birth values:
scaling every state in a run preserves both state and output supports.
Their product E(x,y) therefore has the same value on every compatible
birth fibre. The coefficient [y^L]E retains actual state volume before
the fractional power. It does not replace L by allowed volume.

Apply the raw-field rank bound while every selected quartet is still
required to be live, as in the quotient note. Only then enlarge the
amplitude family to all allowed runs. The resulting probability for
this fixed pair is at most

    min(1, gamma/(q-1) * q^(f-R) * x^(-W) * [y^L]E(x,y)). (2)

## A common envelope for all four run types

Normalize the state marker as s=(q-1)y and put

    d0=x+(1-x)/(q-1),
    M=[[1,sx],[x,s*d0]].

Let lambda be its Perron eigenvalue, let v=(1,rho)^T be a positive
right eigenvector, and let kappa be its other eigenvalue. For x>0 and
s>0, the determinant is nonnegative and the characteristic polynomial
is negative at 1, so 0<=kappa<1<lambda. Define

    d=(1+s*x^2)/lambda,     c=max(1,1/rho).

Then rho>=x, lambda/s>=x, and 0<d<=1. The closed unforced run obeys

    e0 M^ell (1,x)^T <= d*lambda^ell.                   (3)

Indeed, its normalized count is a nonnegative combination of 1 and
(kappa/lambda)^ell and decreases with ell. For a forced closed run,
the normalization by xy*rho*lambda^(ell-1) increases to

    alpha=lambda*(1-kappa)/(lambda-kappa).

The identities

    x/rho=1-kappa,
    d-alpha=kappa^2*(lambda-1)/(lambda*(lambda-kappa))>=0

give the same closed-run factor d. For open runs, (1,1)^T<=c*v gives
the corresponding envelope with factor c. There are at most512 open
runs and no closing output after the last update. Multiplying the run
bounds gives

    E(x,y) <= lambda^A0 * (xy*rho)^f * d^C_A * c^512.   (4)

Using [y^L]E<=y^(-L)E and
q^f(q-1)^m<=q^L in (2), a safe fixed-pair bound is

    min(1, K * lambda^A0 * d^C_A),
    K=2^[32(L-R)-30] * x^(-W) * s^(-m)
                      * (x*rho)^f * c^512.             (5)

The constant -30 uses the established upper bound on gamma/(q-1).
It is conservative and specific to the current q=2^32 proof.

## What must remain in the global sum

For 0<theta<=1, apply the fractional power to (5) after grouping the
amplitudes. The necessary condition A0>=m also permits a marker b>=1.
The geometry weight becomes

    beta^A0 * d^(theta*C_A),    beta=b*lambda^theta,

with external factor K^theta*b^(-m). The birth-mask and routing factors
are unchanged from the fixed-actual-L quotient; the routing factor
already includes the H-set union.

One option is to retain the closed-run statistic in the selector sum.
If C0 is the number of closed runs of the sets U_t intersect U_(t+1),
then C_A>=C0, so replacing d^C_A by d^C0 gives an upper bound independent
of B. This is additional temporal information, not a function of G and
A0 alone. Closed singleton visits of U can occur at every update and
need not produce a closed allowed run. In particular, C_A cannot be
replaced by G-A0 up to a fixed endpoint correction.

A second option retains an actual-death class D>=D0 before enlargement.
There are H_A=A0+f-L=A0-m allowed positions whose actual states are zero.
Every actual death either closes an allowed run or is followed by one
of these holes. Hence D<=C_A+H_A and

    lambda^A0*d^C_A
       <= (lambda/d)^A0 * d^(D0+m).                    (6)

Equation (6) reuses the existing 129-state overlap moment: replace its
marker by beta=b*(lambda/d)^theta and add the external factor
d^[theta(D0+m)]. This applies to the retained death class only. Bounds
for different death intervals must be combined with their proper union
budget. No death restriction may be inferred after dropping it during
amplitude enlargement.

## Why dropping the closed-run factor cannot repair the current gaps

Suppose d^C_A is discarded. Compare the resulting bound with the
theta=0 actual-volume bound using the same beta and volume marker.
Their external-factor ratio is

    [2^[32(L-R)-30] * x^(-W) * (lambda/s)^m
                       * (x*rho)^f * c^512]^theta.     (7)

Since lambda/s>=x, rho>=x, and c>=1, (7) is at least

    [2^[32(L-R)-30] * x^(L+f-W)]^theta.                 (8)

For both remaining profiles, L+f<=W and R<=L. Thus the h=2048 profile
with R=L=819200 can gain at most 30 bits from this envelope. For the
h=2560 profile with R=819200 and L=829030, every positive theta makes
the comparison worse. These statements hold for all x,s,b,theta in the
stated domains, not just a tested grid. They concern this envelope with
the closed-run factor removed, not every possible output-aware count.

There is also a direct low-death limitation. The output support lies in
supp(Q_t) union supp(Q_(t-1)), so its weight is at most L+D. At the two
profiles, D<=124518 and D<=114688 respectively makes the output-weight
condition automatic. Filtering output weight cannot remove those
subfamilies by itself.

The next useful combination therefore retains closed-run or death
information, or strengthens the consistent-selector geometry through
the regional-capacity and required-birth conditions. No volume-only
output-marker search is warranted by (7)--(8).

The parent-run exact checker passed 1008 bivariate run comparisons,
54 volume-only negative controls, 720 death/hole envelope checks, and
four algebra checks at q=2^32 in 0.021 seconds. These validate the
envelope identities and their scope; they do not certify a new global
geometry sum.
