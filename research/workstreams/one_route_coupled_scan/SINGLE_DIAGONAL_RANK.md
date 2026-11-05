# Delayed syndrome constraints for one-diagonal feedback

Status: an analytic counting interface with passed exact regression tests.
There is no finite distance certificate or native timing claim here.
The construction remains the proposed [one-diagonal inner](SINGLE_DIAGONAL.md).

The previous rank bound charges only the immediately preceding output
batch. It misses constraints on older output variables when the syndrome
persists through zero-output batches. Tracking the reachable syndrome
space recovers those constraints without changing the encoder.

## Exact dimension of the linearized output set

Fix the entire output-support sequence W_t and inverse-randomness-support
sequence H_t. These are exact supports in the event being counted. Write
w_t=|W_t| and h_t=|H_t|. Fix the deterministic moment matrices C_t from
ADDITIVE_MDS.md; their columns use globally distinct evaluation points.
The state reconstructed from an output satisfies

    S_0=0,    V_t=C_t^(-1)S_t,    S_(t+1)=S_t+C_tY_t.

For counting, introduce free variables at the coordinates in W_t and
impose V_t[j]=0 for j outside H_t. Temporarily omit all required nonzero
conditions. Let L_t be this linear space of output prefixes before
processing batch t. Define

    m_t=dim L_t,
    U_t={S_t: output prefix in L_t},    d_t=dim U_t.

For a coordinate mask J, let C_t[J] denote the span of its moment columns.
Set

    E_t=U_t intersect C_t[H_t],        e_t=dim E_t,
    nu_t=dim(E_t intersect C_t[W_t]).

Imposing the current V-zero equations cuts L_t by exactly d_t-e_t
dimensions. Appending the w_t output variables gives the exact recurrences

    U_(t+1)=E_t+C_t[W_t],
    d_(t+1)=e_t+w_t-nu_t,
    m_(t+1)=m_t-d_t+e_t+w_t.                              (1)

The fiber argument is important. The map L_t -> U_t is surjective, with
every fiber having q^(m_t-d_t) elements. Restricting its image to E_t
leaves q^(m_t-d_t+e_t) prefixes, after which the new variables contribute
q^w_t. The per-step factor q^(w_t-d_t+e_t) may be fractional; the total
q^m_t is always an integer. This is deterministic dimension accounting,
not a probability assigned to another random coefficient.

Equation (1) counts the linearized masks exactly. It upper-bounds the
number of outputs realizing the exact masks. If H_t is empty, e_t=0.
If H_t is nonempty but e_t=0, no prefix can realize that exact state
support, and the mask sequence can be discarded. When e_t>0, the
linearized space still contains zero syndromes. Keeping those histories
is a safe overcount, not a statement that all surviving states are live.

The final state is free: no extra equation is imposed after the last
output. If an interval is instead required to end at zero syndrome,
its terminal constraint reduces m_end by exactly d_end.

## A finite dimension-only relaxation

After batch t-1, U_t lies in the span of the preceding coset's columns
indexed by H_(t-1) union W_(t-1). Write its cardinality as ell. At the
current batch write

    h=|H_t|,    w=|W_t|,    o=|H_t intersect W_t|.

The previous and current cosets are disjoint. Every set of at most b
global moment columns is independent. Consequently the two containing
column spaces have intersection dimension (ell+h-b)_+, and

    0<=e<=min(d,(ell+h-b)_+),
    max(0,e+o-h)<=nu<=min(e,o).                            (2)

The second line follows because E_t lies in C_t[H_t], and
C_t[H_t] intersect C_t[W_t] has dimension o. The outgoing relaxation
state initially has the containing size

    (d',ell')=(e+w-nu,h+w-o).                             (3)

Use e=0 when h=0. For h>0, only e>=1 can contribute to an exact support.
These are necessary conditions; they do not assert that every allowed
dimension choice is attainable. Initial state is (d,ell)=(0,0).

There is an additional known subspace, not just a containing space.
If k=|W_(t-1)|, then U_t contains the k-dimensional column space
C_(t-1)[W_(t-1)]. Quotienting U_t by that subspace gives

    e <= (d-k)+(k+h-b)_+ = d-min(k,b-h).                  (4)

Indeed the kernel of the quotient map on U_t intersect C_t[H_t] is the
intersection of two disjoint-coset column spaces, of dimension (k+h-b)_+.
Its image has dimension at most d-k. Combine (4) with (2).
This retains the original fresh-variable rank charge rather than replacing
it by the weaker containing-space bound.

If nu=e, then E_t is contained in C_t[W_t]. In that case U_(t+1) equals
C_t[W_t], and ell' can be tightened to w. Equivalently this occurs when
d'=w. Otherwise retain ell'=h+w-o. The outgoing known core has dimension
k'=w in both cases. A practical relaxed state is therefore (d,ell,k),
starting at (0,0,0).

There are exactly

    binom(b,h) binom(h,o) binom(b-h,w-o)                   (5)

mask pairs with the given h,w,o. A positive continuation bound may sum
over these classes and maximize over the allowed e,nu in (2). It must
not sum over e,nu as if they were independently chosen setup outcomes.
For example, with positive markers z,u, a continuation f is bounded by

    T(f)(d,ell,k) = sum_(h,w,o) (5) z^w u^h
      max_(e,nu allowed) q^(w-d+e) f(d',ell',w),           (6)

where (2), (4), and the refined rule for ell' define the allowed choices.

An empty maximum contributes zero. Formula (6) is uniform in the actual
incoming syndrome subspace and mask positions. It is not an exact Markov
chain. Coefficientwise dimension/weight tables can retain more information,
but a full numerical evaluation of this operator is not needed for the
amortized bound below.

## Amortized rank over an excursion

Consider one excursion of L batches. Batch zero activates a zero state,
and batches 1,...,L-1 enter with nonzero state. Let w=sum_t w_t,
H=sum_t h_t, and V=sum_t nu_t. The exact-mask conditions give
h_0=e_0=o_0=nu_0=0. For every later batch that can realize its exact mask,
e_t>=1, so (2) gives

    h_(t-1)+w_(t-1)-o_(t-1)+h_t >= b+e_t.

Sum these inequalities, retaining both endpoint terms:

    2H >= b(L-1)-w+h_last+w_last
          +sum_t e_t+sum_t o_t-o_last.                   (7)

Since nu_t<=e_t,o_t,

    sum_t e_t+sum_t o_t-o_last >= 2V-nu_last.

Also h_last>=e_last and d_end=e_last+w_last-nu_last.
Therefore

    2(H-V) >= b(L-1)-w+d_end.                            (8)

Telescoping (1) yields m_end=V+d_end. A closed excursion has output
dimension m_closed=V after imposing its terminal zero. A free-terminal
excursion has dimension m_free=V+d_end. Equation (8) proves

    m_closed-H <= (w-b(L-1))/2,
    m_free-H   <= (w-b(L-1)+b)/2.                         (9)

The first inequality could retain the additional -d_end/2 improvement.
The displayed weaker version is sufficient for composition.

For c excursions occupying v batches, put a=b(v-c), and let E be the
dimension of all linearized output assignments after the c-1 required
resets. Closed excursions have disjoint output variables. Thus (9) gives

    E-H <= (w-a+b)/2.                                    (10)

Every eligible closed excursion has at least b+1 nonzero outputs. The
global MDS reset equations consequently contribute b independent ranks
per closed excursion, giving the complementary bound

    E <= min(w-b(c-1), H+(w-a+b)/2).                     (11)

Zero-dimensional branches cannot contain a nonzero output. Classes with
an impossible closed excursion must be removed before using the first
term of (11); negative dimensions are not assigned to those classes.

The previously problematic periodic pattern illustrates the improvement.
For h=(b,b,k), w=(b,0,0), with 1<=k<=b, a full output batch restores the
reachable syndrome space to dimension b. The later h=k observation
then charges b-k ranks, even though its immediately preceding output
is zero. The local fresh-variable formula charged no rank there.
This is a comparison of bounds, not a constructed low-weight codeword.

## One observation batch sees a Cauchy matrix

The moment interface gives more structure than two arbitrary MDS spaces.
Let s_b be the vanishing polynomial of the base additive subspace, with
constant derivative kappa!=0. For an observation coordinate x=a_t+beta_j
and a source coordinate z=a_s+beta_k from an earlier batch s<t, the
impulse coefficient is

    (C_t^(-1)C_s)[j,k]
      = (s_b(a_t)+s_b(a_s))/(kappa*(x+z)).                (12)

This follows from the symmetric additive-coset evaluation formula.
The numerator is nonzero because the cosets are disjoint. For a fixed
observation batch t it depends on the source column, but not the row.
Thus the matrix from any collection of earlier output coordinates to
any coordinates of V_t is a column-scaled Cauchy matrix. All its square
minors are nonzero: observation points are distinct, source points are
distinct, and the two sets are disjoint.

Choose observation times 1<=t_1<...<t_m<n, and set t_0=0. At time t_i,
use only the output variables from batches t_(i-1),...,t_i-1. Earlier
selected observation equations contain none of these variables, because
V_t is formed before Y_t is appended. Equation (12) supplies

    r_i=min(b-h_(t_i), sum_(s=t_(i-1))^(t_i-1) w_s)

independent pivots on those variables. The selected rows and columns
form a block-triangular nonsingular minor. Hence the full zero-equation
rank is at least

    R_partition=max_(selected times) sum_i r_i.           (13)

This delays a charge until an observation batch can use it. It does not
charge the same variable twice. The maximum can be computed over temporal
partitions; even a fixed bounded selection gives a valid rank lower bound.
It depends only on h_t,w_t, so it preserves the mask symmetry used below.

For a required closed terminal, one may add a final full-syndrome
observation. Its fresh-column rank is min(b, fresh output count) by the
global Vandermonde property. Do not add that observation at a free terminal.
In particular, arbitrary minors using rows from several observation times
have not been proved nonzero. Their causal zeros and structured numerator
prevent treating the entire impulse matrix as one generic Cauchy matrix.

## Combining with the actual packet route

For a fixed exact output/state mask pair, let A be the H random inverse
coordinates. The affine inverse law bounds the expected number of valid
outputs by q^(E+D_out(A)-H). Dividing by q-1 counts projective lines:
the masks and every imposed zero equation are invariant under common
nonzero scaling. For 0<nu<=1 and tau=q^nu, a valid fractional bound is

    (q-1)^(-nu) tau^(E-H) E_route[tau^D_out(A)].           (14)

To obtain (14), cap the conditional projective count at one before taking
its nu-th power. Apply (11) and E<=w-R_partition before discarding H.
The earlier complete-packet shortening
formula must not be applied directly to this partial-packet A.

The dimension relaxation above depends only on h,w,o,d,ell,k, not coordinate
labels. Its complete mask family is invariant under jointly permuting
the W and H coordinates within each routed packet. Therefore its summed
route moment can be symmetrized by independently averaging those local
permutations. This does not add a shuffle to the encoder. It is an
identity for the relaxed mask sum, not for the realizable masks alone.
Any resulting moment bound must retain the shared packet-count law;
symmetrization does not make packet lanes independent.

This supplies a bounded next certificate interface: combine the excursion
count with these dimension bounds, a proved symmetrized partial-packet moment, and the actual
output/state mask multiplicities. Keep quiet stretches, activation costs,
and the free terminal. No margin follows from the rank inequality alone.

## Exact regression gate

`theory/test_syndrome_rank.py` compares (1) against dense GF(16) constraint
matrices, checks the dimension/core relaxation and (9), and verifies
the single-observation Cauchy minors and selected-time rank bound. It also
enumerates all two-batch outputs for width two and compares exact support
counts with the linearized bound. The finite tests do not certify K20.

Run serially with `python -B`. Root's bounded replay passes 1026 dense
fixtures, including 291 fractional-step cases and 277 eligible excursion
fixtures. It verifies 120 Cauchy entries, 2428 single-observation minors,
and all 65536 width-two/two-batch outputs across 29 exact support classes.
The tests also compare selected-time free and closed rank bounds with
dense elimination. No files are generated.
The amortized inequality with only D_out(A)<=H/2 still leaves a positive
dense exponent above roughly one-third distance. It should not be tuned
as if it already supplied the desired 35% margin. The next bounded gate
should combine the fresh/core constraint, delayed charges, and the actual
partial-packet moment. Preserve the exact subspace calculation as the
small-instance oracle, not an unrestricted large-field state-space search.
