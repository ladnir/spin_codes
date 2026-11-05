# Exact active-zero rows for the return moment

The total-zero refinement can combine moment bounds that require different
active-zero distributions. Retaining the exact active marginal imposes an
additional consistency constraint. This note derives a positive rowwise
relaxation, not the full joint-moment linear program.

Use the [V60 construction and principal event](HD_FIELD128_WIDE_STATE.md).
Fix a state dimension s in {32,36,40}, q=2^128, a good gamma, and a potential
input support S of size j in the64 packet coordinates. Sample X uniformly
from F^S. Write R=CX and Y=(I+gamma AC)X. Let A0 count zero coordinates
of Y inside S, and B0 count zeros outside S. Every measure below includes
the indicator R!=0; none is conditioned on that event.

## Exact active marginal

Principal invertibility makes Y_S uniform on F^S. The identity
C_S Y_S=M_S CX, with M_S invertible, preserves the syndrome-zero kernel.
For 0<=c<=j, the exact active factorial moment is therefore

    L_c = E[binom(A0,c) 1_(R!=0)]
        = binom(j,c) q^(-c) (1-q^(-min(s,j-c))).

Indeed, after specifying c active zeros, the remaining j-c columns of
C_S have rank min(s,j-c). Let m_a=Pr[A0=a,R!=0]. Binomial inversion gives

    m_a = sum_(c=a..j) (-1)^(c-a) binom(c,a) L_c.

These are nonnegative exact rational numbers, with m_j=0 and
sum_a m_a=L_0=1-q^(-min(s,j)). The formula also applies when j>s:
it subtracts the shortened MDS kernel rather than treating every nonzero
active vector as a nonzero syndrome. For j<=s it simplifies to
m_a=binom(j,a)(q-1)^(j-a)/q^j for a<j. At j=0 the measure is zero.
The implementation performs inversion exactly and checks the masses and
all reconstructed active factorial moments.

## Inactive moments within one exact row

Put r=64-j. V60 bounds each specified joint-zero event by p_(j,c,e),
using the original input denominator q^j. Consequently

    U_(c,e) = binom(j,c) binom(r,e) p_(j,c,e)
            >= E[binom(A0,c) binom(B0,e) 1_(R!=0)].

Fix an exact row a with m_a>0. For every c<=a, positivity gives

    E[binom(B0,e) 1_(A0=a,R!=0)] <= U_(c,e)/binom(a,c).

The same moment is at most binom(r,e)m_a because B0<=r. Thus set

    H^(a)_0 = m_a,
    H^(a)_e = min(binom(r,e)m_a,
                   min_(0<=c<=a) U_(c,e)/binom(a,c))  (1<=e<=r).

Every term in this minimum is a proved upper bound for the same finite
measure. There is no independence assumption between A0 and B0. Distinct
rows may reuse the same U budget; this weakens the relaxation, not validity.
Zero-mass rows contribute exactly zero and are omitted from evaluation.

## Positive evaluation and packet integration

Apply the generic [V57 refinement](HD_FIELD128_RETURN_TAIL_MOMENTS.md)
to H^(a), with width r, mass m_a, and the same exact 0<z<1. Its clipped
Bernstein and factorial-tail caps bound

    J_a = E[z^(r-B0) 1_(A0=a,R!=0)].

Call their minimum B_a(z). The original Bernstein cap is unnecessary in
evaluation because V57 proves that clipping cannot increase it. Then

    J_true = sum_a z^(j-a) J_a
           <= J_rows := sum_a z^(j-a) B_a(z).

Use min(J_V60,J_rows). This cannot worsen the existing return cap, but no
strict improvement is claimed. Both caps hold for every support S and
every gamma in the same principal event. All row bounds and clipping are
cached independently of z; evaluation uses precomputed exact powers.

Only the return numerator changes. Divide it by q^s-1 once in packet
entry10; quiet, birth, live, density, and the j=0 live cap are unchanged.
Keep the once-global principal failure term outside a future code union.
This changes neither the construction nor alpha, route counts, chronology,
or outer prefix. It does not permit combining bounds across state sizes.

[V64](theory/hd_field128_joint_rows_v64.py) exposes exact rational packets
and proposal-only logarithms. Its tiny tests enumerate joint zero counts
over GF16 with nonzero CA and verify all support sizes, including j>s.
No tests or numerical screens were executed when this source note was
written. Directed coefficient evaluation and full profile coverage remain
separate requirements; this is not a distance certificate.
