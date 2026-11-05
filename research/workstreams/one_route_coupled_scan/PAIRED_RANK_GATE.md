# A local rank gate for one or two active inputs

The scalar active-step envelope forgets deterministic zero equations and
allows unrelated rank choices at its entrance and exit. This note retains
the span of the old state and the current inputs throughout one step.
It does not introduce a table indexed by the global message dimension.
The resulting Bellman bound is a local proof interface, not a numerical
closure or a whole-code certificate.

Use the actual recurrence from [CONSTRUCTION.md](CONSTRUCTION.md), with
width w=4, antipodal state wiring, coefficient set A=theta+GF(256), and
q=256. The argument also applies to another width with R+D invertible.
The input support I has size k=1 or 2. Fix the outer realization, full
route, and the history before exposing this step's coefficients.

## 1. Retain two local ranks and the remaining input directions

Let V be the fixed shortened message space. Keep the global zero-form
space W and every required nonzero form, exactly as in [PROOF.md](PROOF.md).
The potential is

    Phi = 1[all requirements retained] tau^(-dim W),  tau>=1.

Let Q_old be the state at the entrance of this step, and set

    U = span_F{Q_old[j]},
    Xspace = span_F{x_j : j in I},
    Uplus = U + Xspace.

The form spaces U, Xspace, and Uplus are fixed during the step. Only W
grows as zero equations are imposed. At every intermediate branch, retain

    a = dim((U+W)/W),
    b = dim((Uplus+W)/W),
    P = {j in I : x_j is not in U+W}.

Thus P records individual input forms still outside the old-state span.
It is not the input support: every x_j with j in I remains required
nonzero, whether or not j belongs to P.

At entrance, a is the incoming residual state rank r. Put ell=b-a.
For k<=2 the following restrictions hold:

- 0<=ell<=k and b<=a+k.
- If ell=0, P is empty.
- If ell=1, P is a nonempty subset of I.
- If ell=2, k=2 and P=I; the two input images are independent modulo U+W.

For a live incoming support H, at least one old state coordinate is a
required nonzero form. Therefore a>=1 throughout every retained branch,
even if the next state becomes quiet. The condition b>=1 alone is weaker:
new input directions must not allow a zero equation to kill the entire
required old state.

If no more input information is available, maximize over the permitted
entrance choices of b and P. This is a uniform relaxation of the actual
quotient, not a new sampling distribution. For fresh independent inputs,
as in the proved complete-column prefix, use exactly

    b=a+k,  P=I.

No global rank is reset at either boundary. Every independent equation
still increases the same dim W and costs the factor 1/tau.

## 2. Two possible updates for an independent zero

Every output or next-state zero equation belongs to Uplus. If its form f
is independent modulo W, imposing f=0 reduces b by exactly one.
There are two cases.

**An old-span equation.** If f belongs to U+W, it also reduces a by one:

    old(a,b,P) = (a-1,b-1,P).                            (1)

This update is forbidden when the incoming state is live and a=1.
Indeed, it would put every old state form into W. The membership mask P
does not change, because U+(W+F f)=U+W.

**A pivot from input j outside the old span.** At an active coordinate,
each zero form has the shape f=x_j+u with u in U. If j belongs to P,
that form is independent of U+W for every allowed coefficient. The update
is therefore

    new_j(a,b,P) = (a,b-1,P_new),                        (2)

where P_new is empty when ell=1, and P_new=P minus {j} when ell=2.
To verify (2), observe that U+W enlarges by exactly the direction of x_j.
At ell=1, this spans every input image. At ell=2, the other input image
remains independent. The intersection of the enlarged W with U is
unchanged, so a does not decrease.

At an active coordinate j outside P, all its zero forms lie in U+W.
Any independent zero there must instead use (1). Inactive-coordinate
zero equations also lie in U and use (1). Dependent equations add no rank
and change neither a, b, nor P.

Retain every input and old-state nonzero requirement during these updates.
An equation can also kill another input or an earlier nonzero output.
The actual branch must then be discarded. An abstract operator may ignore
such contradictions when they are not determined by (a,b,P); this only
enlarges the bound. It may not ignore the forced a=0 contradiction.

## 3. Joint outcomes at an active coordinate

At coordinate j in I, write the current quotient forms as

    x=x_j,  s=Q_old[j],  v=(R Q_old)[j],
    y=x+c s,  p=v+y,  c sampled independently from A.

The nonzero/zero status of s and v is fixed by H and R(H). Their nonzero
requirements survive all earlier equations in this step. The same
coefficient c occurs in y and p; never price them as independent tests.

Let C, Y, and Z denote the following nonnegative continuation functions:

| Label | y | p | Output factor |
|---|---|---|---|
| C | nonzero | nonzero | z |
| Y | zero | nonzero | 1 |
| Z | nonzero | zero | z |

Each function includes its output factor, the chosen next-support bit,
and all remaining within-step work. A superscript minus means evaluation
after the applicable rank update (1) or (2), not an extra factor 1/tau.
Unmarked functions use the current (a,b,P). Impossible branches have
continuation zero.

### Zero current state coordinate: s=0

Now y=x is required nonzero, and c is irrelevant.

- If v=0, only C is possible.
- If v!=0 and x=v modulo W, only Z is possible.
- If v!=0 and x is a different nonzero multiple of v, only C survives.
- If x and v are independent modulo W, the support-branch sum is bounded
  by C+tau^(-1) Z^-.

When j belongs to P and v!=0, only the independent case is possible and
its zero branch uses new_j. When j is outside P, maximize over the valid
equality, unequal-multiple, and independent cases; the last uses old and
requires a>=2. This recovers deterministic equations such as x=v rather
than silently assigning them zero rank cost in every case.

### Zero rotated state coordinate: s!=0 and v=0

Here y=p. Use N for both nonzero, including output factor z, and B for
both zero, including output factor one. The possible complete branch sums
are

    N + tau^(-1) B^-                    independent x,s,
    (1-1/q)N + (1/q)B                   dependent root in A,
    N                                  dependent root outside A.       (3)

If j belongs to P, only the first case is possible, with update new_j.
Otherwise its independent case uses old and requires a>=2. The dependent
cases add no rank. Take a maximum over complete cases, not a separate
maximum for each support outcome.

### Both state coordinates nonzero: s!=0 and v!=0

The outcome y=p=0 is impossible because y+p=v is required nonzero.
The remaining cases are those of [GLOBAL_RANK_GATE.md](GLOBAL_RANK_GATE.md),
with the local updates specified above.

If both x and x+v belong to F s modulo W, their roots are distinct.
For the subset T of the two labeled roots that belongs to A, the bound is

    (1-|T|/q) C + (1/q) sum_{L in T} L,                 (4)

where the root labels are Y and Z. Maximize over all four root subsets
if the quotient does not retain their actual membership. No rank changes.

If only the y test is dependent, put epsilon=1 when its root belongs to A
and zero otherwise. The complete bound is

    (epsilon/q)Y
      + (1-epsilon/q)[C + tau^(-1) Z^-].                (5)

At the exceptional coefficient, y=0 and p=v!=0, so the Z branch is
impossible. This explains the factor 1-epsilon/q on its rank term.

If only the p test is dependent, interchange Y and Z:

    (epsilon/q)Z
      + (1-epsilon/q)[C + tau^(-1) Y^-].                (6)

If both zero tests are independent, a safe bound is

    C + tau^(-1)(Y^- + Z^-).                           (7)

The zero branches in (7) are alternatives, each imposing one equation.
They do not jointly impose two ranks. At j in P, only (7) is possible
and both zero branches use new_j. At j outside P, the independent terms
in (5)--(7) use old and require a>=2. When a=1, only (4) is available.

The actual quotient selects one complete classification before c is
exposed. Maximize over these classifications, including epsilon or T
where needed. Do not add their bounds or independently select the largest
term from each formula.

## 4. Inactive coordinates and exposure order

For j outside I, y=c s and p=v+c s.

- If s=v=0, both outputs at this coordinate are zero.
- If s=0 and v!=0, y=0 and p!=0.
- If s!=0 and v=0, y and p are both nonzero; include one factor z.
- If s and v are nonzero, y is nonzero and p is optional.

In the last case, include the factor z and use the usual three alternatives

    F_nonzero + tau^(-1) F_zero^-       [only if a>=2],
    (1-1/q)F_nonzero + (1/q)F_zero,
    F_nonzero.                                         (8)

The independent branch uses old. Restricting it by b>=2 instead of a>=2
would allow an impossible independent relation between rank-one old-state
forms merely because an unused input supplied a second direction.

Process all w coordinates in a fixed order. A within-step state needs the
coordinate index, (a,b,P), and the partial next support. Output factors can
be multiplied locally rather than storing an output-weight counter.
At each coordinate, take the maximum of the complete applicable cases.
Earlier equations can alter the next coordinate's classification, so this
is a sequential Bellman bound.

Every unexposed coefficient remains uniform on A conditional on the
preceding exposed coefficients. Nonzero requirements are retained
indicators, not events on which the coefficient distribution is renormalized.
All case bounds are uniform over the unrecorded forms and requirements.
Consequently the nested maxima safely overestimate any actual sequence
of quotient classifications, including correlations between the two inputs.

## 5. End the step with the actual span identity

Let H_next be the completed next-state support. For every coefficient
realization, R+D is invertible and

    Q_next = (R+D)Q_old + X,
    span(Q_next, Xspace) = Uplus.                       (9)

The identity holds in the quotient by the final W of this step. If r_next
is the next-state residual rank, (9) gives r_next>=b-k. There is a second
lower bound: modulo U+W, the next-state coordinates have exactly the input
images x_j. Their image span has dimension b-a. Hence r_next>=b-a.

For a live next state, retain only ranks in

    max(1,b-k,b-a) <= r_next <= min(b,w,|H_next|).       (10)

For a quiet next state, rank zero is possible only if

    P is empty and a=b<=k.                              (11)

A terminal continuation is bounded by its maximum over the allowed ranks,
not their sum: each actual branch has one rank. If the rank range is
empty, the branch contributes zero. Conditions (10)--(11) are necessary
conditions; they do not assert that every allowed rank is attainable.

For a quiet incoming state, handle the whole step directly. It emits X,
the next state is X, and no zero equation or coefficient charge is needed.
If the residual input rank is e, the exact continuation is z^k f(e,I).
Without further information maximize over 1<=e<=k; fresh independent
inputs have e=k. This case does not enter a live-state recurrence requiring
a>=1.

## 6. Compact common envelopes and limits

For width four and k<=2, only a<=4, b<=6, and masks P within two inputs
are needed. These are local ranks, not a reset of the global rank budget.
At the next step recompute its old-state span from Q_next while retaining
the same W. Forgetting the global remaining dimension is a relaxation;
it does not change the source factor tau^(dim V-1).

The local Bellman operator is a maximum of nonnegative linear
continuations. It is monotone, positively homogeneous, and subadditive.
Thus a common two-state envelope can be obtained by applying it separately
to the quiet indicator and a verified positive live vector. Bound each
result uniformly over incoming (r,H) and permitted input memberships.
The resulting nonnegative matrix dominates every nonnegative combination
of these two basis functions.

Establish that common linear bound for each fixed input mask before
averaging the actual route. Do not move a nonlinear Bellman maximum
through the route expectation. Pair-slot collisions, matching parity,
source dimensions, and the free terminal remain separate obligations from
[TWO_GROUP_ROUTE.md](TWO_GROUP_ROUTE.md).

This operator retains deterministic ranks, coupled coefficient outcomes,
and a local augmented-span budget. It still discards many relations among
old state forms, input forms, and earlier nonzero requirements. No claim
is made that it closes the previously failed two-group classes. Bounded
small-field tests and an independent implementation audit are required
before its output can support a numerical certificate.
