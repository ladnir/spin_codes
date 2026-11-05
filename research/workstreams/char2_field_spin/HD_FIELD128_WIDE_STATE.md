# Wider state without a self-orthogonal feedback premise

Larger state strengthens the return denominator and input kernel but lowers
emission distance. These are new candidates, not distance or timing claims.

## Construction and one global setup condition

Let E=GF(256) be embedded in F=GF(2^128), and put q=|F|. Fix
s in {32,36,40}. Let A interpolate from byte-labelled nodes128..127+s
to nodes0..63, and set C=A^T. Both maps have the MDS row/column ranks.
Let H_s interpolate from0..s-1 to s..2s-1. These Lagrange matrices have
nonsingular square minors.

Use the packet recurrence

    Y = X + A Q,
    Q_next = D2 H_s D1 (Q + gamma C X).

Use fresh independent nonzero full-field diagonals. Keep the outer ensemble,
paired route, initial zero and unrestricted final state. Independently sample
gamma once from Gamma16=F minus GF(2^64). Setup is fixed across messages.
Each s identifies a separate construction.

Define G to require M_D=I_s+gamma C_D A_D invertible for every subset D
of the64 output coordinates. The determinant polynomial has constant1
and degree at most s. The union bound therefore gives

    Pr[not G] <= s * 2^64 / (q - 2^64).

Charge this risk once outside a future conditional whole-code bound.
The following statements hold for each fixed gamma satisfying G.
We neither assert CA=0 nor use an involution or the translated rank lemma.
The state32 control here has feedback A128^T, not the V59 feedback A64^T.
Thus V59 numerical bounds cannot be retained across this change.

## Exact moments and the active-output bijection

Fix a potential input support S of size j, and sample X uniformly on F^S.
Put R=CX. On a return, incoming Q=gamma R is nonzero and
Y=(I+gamma AC)X. Sylvester's determinant identity makes
U_S=I_j+gamma A_S C_S invertible under G. Consequently Y_S=U_S X_S
is uniform on F^S. Moreover,

    C_S Y_S = (I_s+gamma C_S A_S) C_S X_S = M_S R.

Thus the bijection preserves the syndrome-zero kernel even when CA is
nonzero. Since 0<z<1, discarding inactive output weight proves that the
return moment is at most the birth moment on S.

Write w=1-z, t=z+(1-z)/q, u=min(s,j), and d=q^s-1. The exact quiet
and birth moments are

    K_j = sum_c C(j,c) q^(-c-min(s,j-c)) z^(j-c) w^c,
    B_j = sum_(c<j) C(j,c) q^-c z^(j-c) w^c (1-q^-min(s,j-c)).

For an independent uniform nonzero incoming state, the exact live moment is

    V_j = t^j/(1-q^-s)
          * sum_(0<=e< s, e<=64-j)
            C(64-j,e) q^-e z^(64-j-e) w^e (1-q^(-(s-e))).

MDS ranks suffice also for s>32. The elementary return caps are

    J_j <= min(B_j, (1-q^-u) z^max(0,65-s-j), d V_j).

The second cap uses emission distance65-s. The last includes the unique
return state as one summand of the sum over all nonzero incoming states.

## Dimension-only joint zeros

Fix c active zero coordinates D in S and e inactive zero coordinates Z.
Set T=S minus D, t0=j-c and u0=min(s,t0). Elimination gives

    R = M_D^-1 C_T X_T,    A_Z R = 0.

Here rank C_T=u0 and rank A_Z=min(e,s). Restricting A_Z to the u0-dimensional
image of M_D^-1 C_T gives residual rank at least

    b = max(0, min(e,s)+u0-s).

There are at most q^(t0-b) solutions and exactly q^(t0-u0) of them
have zero syndrome. The original sample space still has size q^j, so

    p_(j,c,e) = (q^(u0-b)-1) q^(t0-u0) / q^j

is a valid upper mass for the joint-zero event with R nonzero. It is
exactly zero for c=j or e>=s. No principal-nullity, reversal or
half-Cauchy assertion enters this bound.

Sum C(j,c) C(64-j,e) p_(j,c,e) by k=c+e to obtain upper factorial
moments H_k of the total number of output zeros on R nonzero. H_0 is
exactly1-q^-min(s,j). The positive expansion

    J_j <= sum_k H_k z^(64-k) (1-z)^k

can be refined by the unchanged finite-measure clipping and tail identities
in [the factorial-moment note](HD_FIELD128_RETURN_TAIL_MOMENTS.md).
Take their minimum with the elementary caps above.

## Packet interface and scope

Set rho=d/(q-1)^s. The row-major packet cap before alpha is

    (K_j, rho B_j, J_j/d, rho V_j).

At j=0 also cap the live entry by z^(65-s). Refresh density multiplies
birth and live, not return. Route multiplicities and support counts remain
outside alpha; the chronological transfer and each witness's own outer
factor remain unchanged. Each state requires fresh coefficients and full
original-weight profile sums.

[V60](theory/hd_field128_wide_state_v60.py) implements rational caps and optional
V57 refinement; FLOAT logs use only those rationals. Its unexecuted tests
include width3/state2 with CA nonzero and independent return/live enumeration.
A certificate still needs directed arithmetic, full coverage and the
once-global setup risk. No numerical gain or timing is claimed.
