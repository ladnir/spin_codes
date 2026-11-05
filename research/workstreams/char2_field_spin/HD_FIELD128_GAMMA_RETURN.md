# A fixed extension-field feedback scalar

The dense FLOAT failures identify a limitation of the current upper bounds,
not low-weight codewords. At active-group weight4096, the saved full-RS64/state32
45% search has an upper-bound logarithm of about20586.7. The corresponding
width128/state48 and state64 screens at48% give about30211.2 and29601.8.
Their selected alias entries are above the FLOAT cancellation floor.
This rules out the floor fallback as the explanation for those runs; it
does not distinguish local-cap, coefficient-alias, or common-potential loss.

One identifiable local loss occurs in the return transition. The former
rank/fiber count permits many syndromes after imposing all active output
zeros, without checking the resulting feedback consistency equation.
This note changes the feedback map so that equation is always invertible.
It then counts active and inactive output-zero constraints jointly.

## Construction and principal matrices

Let E=GF(256), F=GF(2^128), and fix the embedding E into F.
Choose gamma in F of degree16 over E. For the public profiles here, let
s be16,24,28, or30. Let A be the64-by-s Lagrange evaluation matrix from
nodes128..127+s to nodes0..63 in E, and put C=A^T.
The evaluation set is an additive F2-subspace, not a64-element subfield.
The degree bound2s-2<63 gives CA=0 by the additive-subspace power sums.

Keep the original outer ensemble and paired route. Replace the inner
recurrence by

    Y_t = X_t + A Q_t,
    Q_(t+1) = D2_t H_s D1_t (Q_t + gamma C X_t).

The initial state is zero; state persists across regions; the final state
is unrestricted. The fresh diagonals have independent nonzero entries in F.
The fixed H_s is the existing totally nonsingular Lagrange map.
The scalar gamma is fixed, not newly sampled at each step.
This is a different construction from feedback C without gamma.

For a coordinate set D of size c, write A_D for the corresponding rows.
Set M_D=A_D A_D^T and p_D(t)=det(I_c+t M_D) in E[t].
Each Lagrange row sums to1. In characteristic2, every diagonal entry of
M_D is therefore1. The matrix is symmetric.

In the determinant expansion, a cycle of length at least3 cancels with
its reversed cycle. Only fixed points and disjoint transpositions remain.
For a matching B of m transpositions, its contribution is

    t^(2m) (1+t)^(c-2m) product_{ij in B} (M_D)_(ij)^2.

Consequently p_D(t)=R_D(t)^2 for even c, and
p_D(t)=(1+t)R_D(t)^2 for odd c, where R_D(0)=1.
Sylvester's determinant identity gives

    p_D(t) = det(I_s+t A_D^T A_D),

so deg p_D<=s. Thus deg R_D<=floor(s/2)<=15.
The degree16 minimal polynomial of gamma cannot divide R_D, and gamma!=1.
It follows that p_D(gamma)!=0 for every D simultaneously.
Equivalently, every I_s+gamma C_D A_D is invertible.
This is deterministic; there is no union bound or bad-setup event.
The argument does not certify state32, whose square-factor degree can be16.

## Joint return count

Fix a potential input support S of size j and let X be uniform on F^S.
Fix D subset S of size c and E0 subset S^c of size e.
For the return contribution, put R=CX and require R!=0.
The incoming state that returns to zero is Q=gamma R, so the output is
Y=X+gamma A R. Requiring Y_D=0 and Y_E0=0 gives

    X_D = gamma A_D R,
    C_T X_T = (I_s+gamma C_D A_D) R,   T=S\D,
    A_E0 R = 0.

Put u=min(s,j-c). The MDS property gives rank(C_T)=u.
The middle matrix is invertible by the preceding lemma.
Hence R lies in the intersection of an u-dimensional subspace and
ker(A_E0). For e<s, the latter kernel has dimension s-e.
There are at most q^min(u,s-e)-1 possible nonzero syndromes.
Each syndrome has q^(j-c-u) preimages X_T, and X_D is then fixed.
For q=|F|, the joint event therefore has probability at most

    p_(j,c,e) = (q^min(u,s-e)-1) q^(j-c-u) / q^j,   e<s.

For e>=s the event is impossible. In particular p_(j,j,e)=0.
The bound also implies p_(j,c,e)<=(q^(j-c)-1)/q^j<=q^-c.
It imposes active/inactive constraints jointly; it does not assume their
independence or claim rank c+e in every case.

For 0<z<1, the positive expansion of z^wt(Y) now gives

    J_gamma <= sum_{c=0}^{j-1} sum_{e=0}^{min(s-1,n-j)}
        binom(j,c) binom(n-j,e) z^(n-c-e) (1-z)^(c+e) p_(j,c,e),

where n=64. Take the minimum with the existing basic return bounds.
The code evaluates the displayed expression directly, without subtraction
of nearly equal quantities. Terms c=j are identically zero, not numerically
clipped. The former joint count is at least this new count term by term.

## Unchanged caps and proof obligations

The zero/nonzero syndrome split is unchanged because gamma!=0. Thus the
quiet moment K and birth moment T-K are unchanged. The live moment V
averages over a uniform nonzero incoming state and is also unchanged.
The same refresh density multiplies only transitions into nonzero state.
The local matrix remains

    [[K, B(T-K)], [J_gamma/(q^s-1), B V]].

The inherited J<=T-K cap is still valid. The map I+gamma AC is invertible
and fixes ker C. Its restriction to F^S is a j-dimensional linear code,
whose weight moment is at most T; subtracting the fixed quiet subcode gives
the claim. The cap J<=(q^s-1)V also remains valid by inclusion in the
sum over all nonzero incoming states. The emission-distance cap is unchanged.

The new operator still requires a fresh coefficient computation, all-profile
support summation, and directed all-weight coverage. Old positive or negative
numerical bounds do not transfer. No diagnosis here treats an upper-bound
term as an attainable codeword event.

## Cost and current evidence

In the existing sparse tower F=E[Y]/(Y^16+Y^5+Y^3+0x32), gamma=Y has the
required degree. The independent polynomial oracle's irreducibility test
was replayed successfully. Multiplication by Y is a byte shift and sparse
reduction. The new recurrence requires s such fixed-scalar operations per
packet, before adding the feedback to the state. It does not require new
general-field diagonal products or random keys.

The scalar cannot simply be absorbed into a single D1 diagonal: the state
and feedback terms would need different multipliers. No encoder timing
claim follows from the operation count.

`theory/hd_field128_full_rs_gamma_screen.py` supplies exact local formulas,
positive log-domain local caps, and a compatible FLOAT alias coefficient
provider. Its tiny tests check all principal subsets and joint zero counts
for an independent GF64/GF8 toy, the actual CA identities, exact/log agreement,
and the alias interface. No K20 distance gate was run while implementing it.

The next discriminating check is to replace only the local return cap at
saved marks and marker banks. Root should compare state24 at r256/40% and
r4096/45% before optimizing again. State30 is a separately named candidate
under the same determinant lemma. A gain would support the return-loss
diagnosis; a remaining positive bound would require profiling the other
relaxations rather than claiming the construction fails.
