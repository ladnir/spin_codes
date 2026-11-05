# Reducing the cost of the guarded support count

Status: analytic upper bounds, not an implemented support operator or a
whole-code certificate. These reductions change no encoder operation,
coefficient distribution, or routing rule.

The positive envelope in [Guarded regional moments](GUARDED_REGIONAL_MOMENT.md)
has G binomial terms. With G=8192, evaluating a separate regional support
operator for every term would be expensive. Two reductions avoid that
literal expansion. A third removes regional cardinality counters, at the
cost of another explicit relaxation.

## A sparse-activity class needs only two products

Retain the notation of the guarded-moment note: B=2k regions, M=pG rows,
N=BM, and fixed T subset U with U a union of whole packets. The regional
densities are alpha_j=|T_j|/M and upsilon_j=|U_j|/M. Its comparison row has
occupancies Z in T and Q in U. For t=tau^p>1, define

    H=E[t^((Z-k)_+)],
    J=E[1[Q>k] t^((Z-k)_+)].

The guarded route moment is at most G J H^(G-1). All expectations here
refer to the unconditioned comparison distribution. The guard is never
used to renormalize setup probabilities.
Fixed-support route moments use the auxiliary lane-shuffled route.
Only the complete invariant mask sum transfers to the implemented route.

Choose lambda>=t and eta>=1, and put

    c=(t-1)/(lambda-1),       gamma=c lambda^(-k),
    a_j=1+(lambda-1)alpha_j,
    b_j=1+(eta-1)upsilon_j,
    d_j=b_j+eta(lambda-1)alpha_j.

For every integer e>=1,

    t^e-1=(t-1) sum_(r=0)^(e-1) t^r
          <= c(lambda^e-1) <= c lambda^e.

For e<=0, t^(e_+)-1=0. Consequently

    H <= 1+gamma product_j a_j,
    J <= eta^(-k-1)
         [product_j b_j+gamma product_j d_j].                 (1)

The second inequality also uses 1[Q>k]<=eta^(Q-k-1).
When lambda=t, (1) is the original positive envelope. Allowing lambda>t
retains a factor tending to zero as t tends to one. At t=1, instead use
the exact identities H=1 and J=Pr[Q>k].

Now restrict the support class to |U|<=rho N. Since T subset U,
the average of the alpha_j is at most rho. Concavity and monotonicity of
log(1+(lambda-1)x) give

    product_j a_j <= [1+(lambda-1)rho]^B.

Thus, with H_rho=1+gamma[1+(lambda-1)rho]^B,

    E_route[tau^D(T) 1[D(U)>0]]
      <= G H_rho^(G-1) eta^(-k-1)
           [product_j b_j+gamma product_j d_j].               (2)

Only two separable regional products remain. Equation (2) is useful when
the scalar H_rho^(G-1) is affordable; no such numerical assertion is made
here. The bound is uniform over fixed support pairs in this restricted
class. Applying it to the complete mask sum still requires the mask-orbit
transfer and multiplicity rules of the guarded-moment note.

Retain deterministic exclusions separately. If U touches at most k
regions, the guarded moment is exactly zero. The positive envelope for J
need not preserve that zero. Dense classes still require a complementary
bound, such as the deterministic regional information-set cap.

## Compress the binomial sum by positive powers

This reduction applies to any positive product A=gamma product_j a_j
for which finite bounds A_min<=A<=A_max are known. Write g=G-1. If g=0,
the factor below is identically one. Otherwise set

    s=log A,              f(s)=g log(1+exp s),
    s_min=log A_min,      s_max=log A_max.

Partition [s_min,s_max] into m intervals. On interval i, let
ell_i(s)=c_i+kappa_i s be the chord joining its two values of f.
Since f is convex and increasing, 0<=kappa_i<=g. The chord lies above f
inside its interval. Its affine extension lies below f outside that
interval, by monotonicity of secant slopes.

Suppose every chord has excess at most epsilon on its own interval.
At each s, one chord is at least f(s), and every chord is at most
f(s)+epsilon. Exponentiating and summing proves

    (1+A)^g
      <= sum_(i=1)^m exp(c_i) A^kappa_i
      <= m exp(epsilon) (1+A)^g.                            (3)

This is a positive envelope with a quantified loss, not a truncation of
the binomial expansion. Its regional factors are a_j^kappa_i. Multiplying
by the two products in (1) leaves 2m separable terms instead of 2G terms.
The factor gamma^kappa_i belongs in each term's scalar coefficient.

The interpolation-error bound supplies a verification rule. If interval i
has length Delta_i and f''<=L_i there, its chord excess is at most
L_i Delta_i^2/8. Here

    f''(s)=g exp(s)/(1+exp(s))^2 <= g/4.                    (4)

Far from zero, the smaller local curvature permits longer intervals.
For (1), one may use A_min=gamma and A_max=gamma lambda^B.
No particular partition size or numerical speedup is asserted. A future
finite checker must verify the endpoint, curvature, and rounded affine
inequalities with directed arithmetic; ordinary floating chords are not
a certificate. Equation (3) controls loss relative to the positive
envelope, not relative to the exact route moment.

## Tangent bounds remove regional count tables

At fixed kappa>=0, the two regional factors produced above are

    F_b(alpha,upsilon)=a(alpha)^kappa b(upsilon),
    F_d(alpha,upsilon)=a(alpha)^kappa d(alpha,upsilon).

Their logarithms are concave on 0<=alpha<=upsilon<=1: each is a
nonnegative sum of logarithms of positive affine functions. Fix an anchor
(alpha_0,upsilon_0) in this domain. For either phi=log F, its tangent gives

    phi(alpha,upsilon)
      <= phi_0+A(alpha-alpha_0)+C(upsilon-upsilon_0),
    A=partial_alpha phi at the anchor,
    C=partial_upsilon phi at the anchor.

Writing a_0,b_0,d_0 for the affine factors at the anchor, the derivatives
are explicit:

    F_b: A=kappa(lambda-1)/a_0,
         C=(eta-1)/b_0;
    F_d: A=kappa(lambda-1)/a_0+eta(lambda-1)/d_0,
         C=(eta-1)/d_0.

Therefore each regional factor is bounded by

    exp(phi_0-A alpha_0-C upsilon_0)
        exp(A |T_j|/M+C |U_j|/M).                           (5)

The second factor splits across batches. For T=A, a batch contributes
its h selected state coordinates; for T=S it contributes h+w-o.
For a union U of complete active batches, it contributes b coordinates
when h+w>0 and zero otherwise. Thus (5) replaces the pair of regional
count tables by local exponential weights. Other support multiplicities,
inner-rank weights, and the output-weight restriction are unchanged.

This tangent relaxation can lose regional heterogeneity. It must not
silently replace the deterministic cap that handles the phased example.
Fix each tangent for the count being bounded, or assign it to an explicitly
covered class. Pointwise adaptation does not itself provide a separable
complete-count algorithm.

## Remaining implementation and proof obligations

Equations (2)--(5) reduce the number or shape of positive terms. They do
not supply the remaining amplitude-and-support sum. Any implementation
must retain zero initialization, impossible quiet/live transitions,
required resets, the nonzero-output guard, exact multiplicity placement,
and a free terminal. No particular finite-state transfer is asserted here.

A bounded next test is the sparse-activity class using (2), with local
rank weights and a verified positive potential. Keep the exact guarded
moment as its small-instance oracle. Only after that test should a
compressed general-activity operator replace the explicit binomial sum.
