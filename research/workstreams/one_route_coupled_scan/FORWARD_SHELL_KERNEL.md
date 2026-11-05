# Forward support bounds for one-diagonal feedback

This note bounds the actual forward transition of
[one-diagonal feedback](SINGLE_DIAGONAL.md). It does not replace its partial
inverse distribution by uniform whole-batch feedback. The bounds below are
local statements; no potential contraction or whole-code certificate is claimed.

## Fixed-input experiment

Let F have characteristic two and size q. Fix a batch width b and a
superregular b by b matrix D: every square minor of D is nonsingular.
Fix the entire route and an outer word x before exposing any inner keys.
Condition on the past keys. At the current step, the vectors V,X in F^b
are therefore fixed. Sample the current key r from F^b and compute

    Y = X + r .* V,
    Z_next = V + Y,
    V_next = D Z_next.                                      (1)

D is the next batch's fixed change of state coordinates. The terminal
batch has no required outgoing state; bounds using V_next may omit that
batch or introduce a fictitious invertible D only for analysis.

Put A=supp(V), h=|A|. For h>0, Y_A is uniform on F^h and Y outside A
equals the fixed vector X outside A. All probabilities in this note are
conditional on the fixed word, route, and past. Their bounds are uniform
over that conditioning. Consequently they compose by conditional
expectation; the states need not form an independent process.

Choose an output-weight marker 0<z<=1. Omitting the fixed weight outside A
only increases z^wt(Y). This relaxation is used below.

## A support-shell kernel

Fix M subset A of size m as the exact support of Y_A, and fix an outgoing
support H' of size h'>0. After setting Y on A minus M to zero, the free
coordinates are Y_M. The equations V_next outside H'=0 have coefficient
matrix D[H'^c,M]. Its rank is min(m,b-h'). Thus at most

    q^max(0,m-b+h')                                         (2)

assignments satisfy those affine equations. Requiring the selected input
and outgoing coordinates to be nonzero can only reduce this count.
Since Y_A is uniform, summing M and H' gives

    E_r[z^wt(Y) 1[wt(V_next)=h']]
      <= K(h,h';z)
      := q^(-h) binom(b,h')
         sum_(m=0)^h binom(h,m) q^max(0,m-b+h') z^m.          (3)

The binomial factor in (3) can lose substantially. For example, m=0 leaves
one deterministic outgoing vector, not one vector for every support H'.

### Count outgoing zeros without enumerating their full support

Fix M and put d=b-h'. For j<=min(m,d), sum over all j-subsets of outgoing
zero coordinates. Each such system has at most q^(m-j) solutions.
Every assignment with exactly d outgoing zeros occurs binom(d,j) times.
The number of these assignments is therefore at most

    B_q(b;m,d)
      := min_(0<=j<=min(m,d))
           [binom(b,j)/binom(d,j)] q^(m-j).                  (4)

Here j=0 contributes q^m, including d=0 or m=0. Equation (4) counts an
affine image; no zero offset is assumed. It also permits some Y_M entries
to vanish, which is a safe overcount. A second kernel is

    K_sharp(h,h';z)
      := q^(-h) sum_(m=0)^h binom(h,m)
           B_q(b;m,b-h') z^m.                               (5)

For q>=b, the minimum in (4) occurs at j=min(m,d). Indeed the ratio of
successive expressions is (b-j)/(q(d-j))<=1. Either kernel, or their
pointwise minimum, upper-bounds the actual live transition.

### Exact reset and activation cases

Since D is invertible, V_next=0 exactly when Y=V. For h>0,

    E_r[z^wt(Y) 1[V_next=0]]
      = 1[X outside A=0] q^(-h) z^h.                        (6)

At h=0, Y=X and V_next=DX deterministically. If X=0 the state remains
zero. Otherwise superregularity gives

    wt(V_next) >= b+1-wt(X).                                (7)

For the stratified direct-sum MDS outer, let U be a whole-batch activity
union and d_out=dim(C_out intersect F^U). An outer row can be active in
this shortened space only if its occupancy exceeds the constituent
dimension. Each such row contributes at least one to d_out. A batch
contained in one region has at most one coordinate from each row, so any
X from this space satisfies wt(X)<=min(b,d_out). This gives the activation
bound b+1-min(b,d_out) for a nonzero activation. It does not assert that
every subsequent state has that weight.

## A direct potential bound without outgoing shell overcount

The affine weighted sum has an additional property which avoids summing
incompatible outgoing supports. For any F-linear map G:F^h -> F^b,
any a in F^b, and 0<=z,s<=1,

    sum_(u in F^h) z^wt(u) s^wt(Gu+a)
      <= sum_(u in F^h) z^wt(u) s^wt(Gu).                   (8)

Interpret 0^0 as 1 at the endpoints. To prove (8), use additive characters
on F. The unnormalized Fourier transform of a single-coordinate Hamming
function f_z is 1+(q-1)z at the trivial character and 1-z at every other
character. Both values are nonnegative. The transforms of the product
functions f_z(u)=z^wt(u) and g_s(v)=s^wt(v) are therefore nonnegative.
Fourier expansion gives

    sum_u f_z(u) g_s(Gu+a)
      = q^(-b) sum_chi hat(g_s)(chi)
                       hat(f_z)(G^*chi) chi(a).             (9)

Characters in characteristic two have values +1 or -1. Replacing each
chi(a) by 1 can only increase (9), and yields its value at a=0.
The adjoint G^* here acts on additive characters; no trace pairing is
assumed implicitly.

In (1), set G=D[:,A]. The next state is Gu+a with
a=D(V+X_outside_A), where X_outside_A is zero on A. Thus (8) applies to
the actual affine offset. It does not discard that offset by assumption.

### Exact homogeneous split enumerator

The graph code {(u,Gu):u in F^h} is an MDS [b+h,h] code. To verify this,
select any h rows of [I_h;G], eliminate the selected identity rows, and
use the corresponding square minor of D. Its minimum distance is b+1.

For a fixed support of size ell>0 in this code, let a_ell count words
with exactly that support. MDS shortening and inclusion-exclusion give

    a_ell = sum_(i=0)^ell (-1)^i binom(ell,i)
                              q^max(0,ell-b-i).             (10)

In particular a_ell=0 for 1<=ell<=b. For ell>=b+1, equivalently,

    a_ell = (q-1) sum_(i=0)^(ell-b-1)
                  (-1)^i binom(ell-1,i) q^(ell-b-1-i).       (11)

These are exact nonnegative integers. An implementation should evaluate
the signed sums exactly before applying any one-sided numerical rounding.

The split enumerator is consequently

    F_h(z,s) := sum_u z^wt(u) s^wt(Gu)
      = 1 + sum_(m=1)^h sum_(j=b+1-m)^b
              binom(h,m) binom(b,j) a_(m+j) z^m s^j.         (12)

It depends on h,b,q, not on the chosen support A or superregular D.

### Potential interface

Fix 0<=alpha<=1, set s=q^(-alpha), and define for 0<=j<=b

    phi(j)=q^(alpha(b-j))=s^(j-b).                          (13)

In particular phi(0)=q^(alpha b); zero is not assigned the live-state
normalization 1. Combining the uniform Y_A law, (8), and (12) gives

    E_r[z^wt(Y) phi(wt(V_next))]
      <= q^(alpha b-h) F_h(z,q^(-alpha)).                   (14)

After division by phi(h), the conditional ratio is at most

    rho_h(z,alpha)
      := q^(-(1-alpha)h) F_h(z,q^(-alpha)).                 (15)

Equations (14)--(15) include the reset contribution with its actual
phi(0) weight. Removing reset paths only decreases the left side.
Alternatively a continuation can keep (6) separately. Bounds (3), (5),
and (14) are alternative upper bounds, not independent factors to multiply.

Homogenization has a specific loss. The term 1 in (12) represents
u=0 and Gu=0. In the actual live transition, Y_A=0 and V_next=0 cannot
occur together: V+X_outside_A is nonzero on A. The resulting contribution
q^(-(1-alpha)h) can dominate (15) when h is small. Equation (8) does not
justify subtracting the homogeneous term 1 to obtain a live-only bound.
Any improvement using this forbidden origin needs a separate argument.

A uniform check rho_h<=rho would supply a one-step weighted continuation
bound for every live state. Such a check has not been performed here.
Activations still require (7), and a whole-code proof must count outer
projective lines, route events, excursion placements, and a free terminal.
No whole-code claim follows from a fixed-input transition alone.

## Bounded exact regression source

`theory/test_forward_shell.py` uses GF(16), width two, and a directly
constructed Cauchy matrix. It checks the two shell bounds, exact reset
probabilities, activation weights, the homogeneous split enumerator,
the affine Fourier inequality, and the potential interface using exact
fractions. Its two-batch fixtures include every second input for selected
first inputs. The source writes no files and is not a distance certificate.

Run from this workstream with `python -B theory/test_forward_shell.py`.
The root's serial replay passed under a 30-second timeout and a 256-MiB
memory cap. It checked 5 Cauchy minors, 12 fixed-input fixtures with 3072
key assignments, 60 joint-shell inequalities, 108 potential inequalities,
3 exact split enumerators, 162 affine-Fourier inequalities, 256 activation
cases, and 512 two-batch input pairs with 131072 key assignments. These are
exact local regression checks, not a numerical full-size contraction test
or a minimum-distance certificate. The source wrote no files.
