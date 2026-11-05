# Two-tap proof: exact partial-state ranks and a closed trajectory class

Scope: the two-tap construction in [DESIGN.md](DESIGN.md), with F of size
q=2^128, s=64, b=128, N=2K, and K=2^20. The fixed outer and regional route
are unchanged. These results do not yet establish the full minimum distance.

## Partial state gives two shorter RS spaces

Fix an incoming state with support I of size r>0. Let g_I be the product
of (X-z) over the points z in Z outside I. The interpolation polynomials
of DQ and EQ both vanish on Z outside I. Write them as g_I*u and g_I*v,
where u and v have degree less than r.

On the first half of Y, feedback evaluates
g_I*(u+(1+alpha)*v). On the second half it evaluates g_I*(u+alpha*v).
The change from (u,v) to these two degree-<r polynomials is invertible.
The weights g_I are nonzero on Y. Thus the unconditioned feedback space
V_I is the direct sum of two generalized RS[64,r,65-r] codes on disjoint
coordinate halves. The actual pair sampling excludes affine hyperplanes
inside this space; it does not change its span.

If z_L and z_R output coordinates are prescribed, the corresponding
linear observation on V_I has exact rank

    rho_out = min(r,z_L) + min(r,z_R).

This applies to prescribed zeros or any other prescribed values. For a
fixed input and state, the probability of satisfying consistent output
constraints is at most q^(-rho_out)*(q/(q-1))^r. This is a conditional
local statement, not an all-message union bound.

## Add next-state constraints without assuming independence

Also prescribe t next-state coordinates whose indices lie in I. Put
a=min(r,z_L), c=min(r,z_R). The joint observation rank is at least

    rho_joint >= a+c+min(t, r-min(a,c)).

To prove the bound, first impose the output constraints homogeneously.
The remaining polynomial spaces in the two halves have dimensions r-a
and r-c. Evaluation on the selected t points of Z has ranks min(t,r-a)
and min(t,r-c), respectively. Indeed, the imposed roots lie in Y, disjoint
from Z, so each restriction is another shortened RS evaluation map.
The next-state observation is a nonzero weighted sum of these evaluations.
Its rank is at least the larger of the two ranks. Add this to a+c.

Coordinates outside I have no current diagonal randomness. Their next-state
values are determined by BX and must be accounted for separately. The rank
bound need not be exact: geometry can make the two remaining images less
aligned. The checker finds such cases, so a later refinement may retain
actual root sets rather than only their cardinalities.

## Separate state randomness from output randomness

The two-tap distribution has an equivalent, simpler sampling interface.
Write H=I+D+E. For each coordinate, the change of variables (d,e) to
(h=1+d+e,e) is a bijection from the allowed pairs to F* times F.
Thus H has independent nonzero entries, E has independent unrestricted
entries, and H and E are independent. With K0=A+A2, the recurrence becomes

    Q_t = B X_t + H_t Q_{t-1},
    Y_t = X_t + A (H_t+I) Q_{t-1} + K0 E_t Q_{t-1}.

The output still determines Q_t=Q_{t-1}+B Y_t. This is a reparameterization
of exactly the same setup distribution, not a new heuristic.

Fix an input and condition on the incoming state before exposing H_t,E_t.
Let z=B X_t. For coordinate j:

- If Q[j]=0, then Q_t[j]=z[j] deterministically.
- If Q[j]!=0 and z[j]=0, the next coordinate is always nonzero.
- If Q[j]!=0 and z[j]!=0, the next coordinate is zero with probability
  1/(q-1); otherwise it is nonzero.

These transitions are independent across coordinates under this conditioning.
They describe a fixed-input process. They must be combined with an all-message
count or rank argument, rather than interpreted as an adaptive-message bound.

Conditional on the whole state transition, E remains independent. For a
state support I of size r, the output is uniform in an affine translate
of im(K0[:,I]), an [128,r,129-r] generalized RS code. Indeed its polynomials
have the form (v_6+alpha)*g_I*u, deg(u)<r, whose multiplier is nonzero on Y.
Consequently prescribing z output coordinates has rank min(z,r) in this
remaining randomness. Counting output words of weight w in any such coset
gives the upper bound

    min(q^r, binom(128,w)*q^max(0,w-(128-r))).

For each support of size w, projection onto its complement has rank
min(128-w,r), which proves the bound. A coset can contain a low-weight
word even though the underlying linear code has large distance; do not
substitute the linear minimum distance for this coset count.

For a routed input batch of weight h between1 and64, superregularity of B
gives wt(B X_t)>=65-h. In particular, a single nonzero input activates all64
state coordinates from an empty state. This links sparse outer inputs to
the transition law above. Batches heavier than64 may lie in ker(B), and
cannot be treated as activations without further outer analysis.

## A global bound for almost fully supported state trajectories

For a fixed output y, reconstruct Q_0=0 and Q_t=Q_{t-1}+B y_t.
This state trajectory is independent of the diagonal setup. Define

    M(y) = sum_t (64 - wt(Q_{t-1})),
    R(y) = sum_t wt(Q_{t-1}) = K-M(y).

The first batch contributes64 to M. Zero coordinates are counted separately
at every batch; M is not the number of distinct coordinates that ever vanish.

Fix any route and let C be its K-dimensional outer image. The inverse input
has (q(q-1))^R equiprobable possibilities, because [A,A2] is invertible and
all active diagonal pairs are independent. Its affine support has dimension
2R. At most q^K of these possibilities belong to C. Therefore

    Pr_setup[inverse(y) in C] <= q^(K-2R)*(q/(q-1))^R
                              = q^(-K+2M)*(q/(q-1))^R.

This upper bound can exceed one for other classes; it is used only in the
specified deficit range. Since R<=K and 2K<q, Bernoulli's inequality gives
(1-1/q)^R >=1-K/q >1/2. The density factor is consequently less than2.

For an integer W, there are at most 2^N*q^W words of weight at most W.
Summing over every nonzero candidate output with M<=m gives

    Pr[exists nonzero output with wt(y)<=W and M(y)<=m]
        < 2^(1+N+128*(W-K+2m)).

The union is over fixed output words before setup is exposed. It therefore
covers messages chosen after setup; there is no fixed-message assumption.
The bound holds for every route, and hence also after averaging the route.

At K=2^20:

| Distance threshold | W | Deficit budget m | Failure bound for this class |
|---|---:|---:|---:|
| 49% | 1,027,604 | 2,293 | <2^-255 |
| 49.1% | 1,029,701 | 1,245 | <2^-127 |

These are exact integer bounds using 2^N instead of a numerical entropy
estimate. They leave all larger-deficit trajectories unresolved.

## Remove whole inactive regions without paying a deficit

A long zero prefix should not consume M when the outer already pays for
vanishing there. Let U be an interval of whole regions. Consider output
trajectories with y=0 and Q_{t-1}=0 in every batch outside U. A nonterminal
interval must therefore end at state zero. A terminal interval need not.

The inverse input is zero outside U. If U has c regions and L=8192*c
coordinates, its shortened outer dimension is exactly

    D(U) = 8192*max(0,c-128).

Each outer block contributes one coordinate to each region; the route does
not affect this count. The dimension follows from shortening RS[256,128].
For c<=128 the shortened outer contains no nonzero word, so this trajectory
class is empty. For c>128, D(U)=L-K.

Define M_U by summing the missing state coordinates only inside U. The
inverse support has dimension L-2M_U. Replace q^K by q^D(U) in the preceding
membership count to obtain q^(-K+2M_U) times the same density correction.
There are at most2^L*q^W supported output candidates. Using L<=N and taking
the union over all32896 region intervals costs fewer than16 additional bits.
Thus the two rows of the table become <2^-239 and <2^-111, respectively,
for the union of these aligned-interval classes with the same deficit budgets.

This also covers a nonzero state at the final output when U reaches the
end of the code. It does not cover arbitrary gaps, partial-region starts,
or persistent state deficits within U. Rounding an interval to full regions
can spend more deficit than allowed; do not silently treat rounding as free.

## Evidence and next obstruction

`check_partial_rank.py` verifies20480 toy cases over GF16, including every
nonempty state support, output-zero set, and selected active next-state
coordinates. It confirms the exact output rank and the joint lower bound.
The latter is strict in1548 cases, by at most2 ranks in this toy model.
A GF4 enumeration also checks the full-state conditional-uniform law exactly.

`check_dense_class.py` replays the two finite bounds and the interval union
using integer arithmetic only. Neither checker covers the remaining class.

The next global step must exploit the outer constraints together with the
new local rank formulas. In particular, partially supported feedback is a
pair of short RS spaces, not a coordinate subspace. Replacing its intersection
with the routed outer by the dimension bound K is too costly when deficits
are large. The earlier regional occupation moments may handle unaligned
zero prefixes, but do not by themselves control these non-coordinate spaces.
