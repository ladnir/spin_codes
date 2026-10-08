# Exact joint transitions for the GF256 alternative

This note audits a different inner from the native-field three-shear
candidate. Work over F=GF(q), with q=256 in the proposed implementation.
Each step independently samples a uniform G in SL(2,F), or in GL(2,F),
and uses

    R = G Q,   Y = X + A R,   Q_next = R + C X.

Here A has row A_j=(1,j), C=A^T, and the t lane labels are distinct.
Each outer coordinate receives its own independently uniform nonzero
scalar before the inner. For a fixed source word and routing, active
X_j are therefore independent uniform elements of F minus zero.
These scalars and all G are mutually independent. If Q is nonzero,
R is uniform over F^2 minus zero, even conditional on the previous steps.
If Q=0, R=0. Thus a two-state description is exact for a fixed input mask,
provided it retains the **joint** output-weight and next-state transition.

## The nonzero-solution count is exact

For distinct labels U, write n=|U|. Let N_U(h) count tuples of nonzero
field elements x_j satisfying sum_j (1,j)^T x_j=h. Then

    q^2 N_U(0) = (q-1)^n + (-1)^n (q^2-1-nq(q-1)).

For h nonzero, let e=1 if h lies on one of the lines spanned by (1,j)^T
for j in U, and e=0 otherwise. Then

    q^2 N_U(h) = (q-1)^n + (-1)^n (-1+nq-eq^2).

These formulas also handle U empty. Additive-character orthogonality gives
a short proof: each nonzero dual vector annihilates at most one column.
On the dual line perpendicular to that column, the character sum is q-1
if h lies on the column line and -1 otherwise. All other factors are -1.

## Output weight and return are not independent

Fix an active input support S of size h and a nonzero R. Put
ell_j=A_j R and b=number of j outside S with ell_j nonzero. The return
generating polynomial, conditional on R, is

    J_(S,R)(z) = z^b/(q-1)^h
        sum_(T subset S, all ell_j nonzero on T)
        z^(h-|T|) (1-z)^|T|
        N_(S minus T)((I + sum_(j in T) C_j A_j) R).

This follows by expanding each active output factor as
z+(1-z)*1[X_j=ell_j]. Average over nonzero R to obtain J_S.
Scaling R by a nonzero scalar leaves each summand unchanged, so it is
enough to average over q+1 projective representatives.

For incoming zero state, the output weight is exactly h. The return
probability is N_S(0)/(q-1)^h, so that row of the transfer matrix is
immediate. For incoming nonzero state, J_S gives the return entry;
subtract it from the unconditional output polynomial for the other entry.

The joint polynomial really depends on S, not just h. For S={0}, a
return has output weight t-1; for S={1}, a return has output weight t.
Both returns have probability 1/(q^2-1). Indeed a return forces R=C_j X_j,
giving Y_j=j^2 X_j and Y_l=(1+lj)X_j for inactive l.

The unconditional output-weight marginal does depend only on h. With
f(z)=z+(1-z)/(q-1), it is

    ((q+1-t) z^(t-h) f(z)^h
      + h z^(t-h+1) f(z)^(h-1)
      + (t-h) z^(t-h-1) f(z)^h) / (q+1),

omitting a term when its leading integer coefficient is zero.

## Uniform regional routing permits mask averaging

For a fixed outer word, a uniform regional permutation makes its active
positions a uniform w-subset of that region. Conditional on the step
occupancies h_i, the masks in different steps are independent uniform
h_i-subsets. The independent coordinate scalars remove the remaining
nonzero-value correlations. Consequently one may average the exact
mask-specific transfer W_S(z) over masks of each size; one may not
replace each W_S by an identical weight-only law before this averaging.

For a region of length L divisible by t, its exact two-state transfer is

    [u^w] (sum_(S subset {0,...,t-1}) u^|S| W_S(z))^(L/t) / binom(L,w).

Here [u^w] extracts the coefficient of u^w and the power is a matrix
power. This conditions on exactly w active positions, not a Bernoulli
approximation. Combining the eight regions still needs the outer's
joint column-weight information or a sound bound on it.

For fast exact mask averaging, the matrix sum over a subset T is

    sum_(j in T) C_j A_j = [[|T| mod2, s], [s, s^2]],
    s = XOR_(j in T) j.

Thus its action depends only on subset-size parity and label XOR.
Counts by subset size and XOR, with specified labels excluded, can be
computed by a small exact dynamic program when labels form a binary
subspace. This can reduce the cost of the subset expansion.

## Scope

The formulas above do not apply unchanged to the small-set native-field
three shears. They require the stated uniform nonzero orbit and uniform
nonzero coordinate scalars. They also concern base-field words over
GF256. If the whole generator has GF256 entries, scalar extension
preserves its minimum symbol distance, allowing a proof over GF256
followed by use over GF(2^32), GF(2^64), or GF(2^128). Directly treating
an arbitrary extension-field state as a uniform GF256 two-vector would
be incorrect. The representation and subfield embedding must be declared.

`test_joint_transfer.py` independently enumerates all GF4 supports and
small-field nonzero-solution counts, verifies the subset expansion,
checks the geometry counterexample, and checks that uniform SL2 has the
required nonzero orbit. These are algebra tests, not a distance certificate.
