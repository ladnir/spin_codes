# A global rank budget across singleton inputs

The empty-span bound remembers how many state dimensions survive inside
one gap. A stronger continuation can also remember the global message
dimension remaining after earlier zero equations. This note gives the
singleton-input interface needed for that refinement.

The interface is proved below, but no operator implementing it or numerical
certificate is supplied here. The simpler separation between the rank
parameter and coefficient-root probabilities remains a separate proof path.
Neither refinement changes the encoder.

## 1. Keep two ranks, not two independent potentials

Fix an outer realization, the complete route, and an exact source support.
Let V be its shortened F-linear message space, of dimension D. Keep one
common space W of imposed zero forms and all required nonzero forms.
Discard a branch if W contains any required nonzero form.

For any tau>=1, use the global potential

    Phi = 1[branch retained] tau^(-dim W).

The global remaining dimension and current residual state rank are

    d = D - dim W,
    r = dim((span(Q)+W)/W).

A retained branch has d>=1. Its support H is empty exactly when r=0.
Otherwise 1<=r<=min(d,|H|). The zero state is allowed even when d>0;
it does not mean that the message is zero.

An independent zero equation decreases d by one and multiplies the
potential by 1/tau. A dependent affine equation that survives at one
allowed coefficient has probability 1/q, where q=256. These are distinct
effects. Keeping their distinct factors allows tau>q. The old scalar
fixed-mask charge, if used elsewhere, must instead use 1/min(tau,q).

The existence indicator for a nonzero witness remains bounded by
tau^(D-1) times the summed final potential. This inequality is valid for
every tau>=1. The GF256 support bridges stated for tau<=q must not be
silently reused with tau>q.

## 2. The augmented span at a singleton input

Fix a step with exact input support {j}, and let x be its required nonzero
input form. Put

    U = span(Q),
    U_plus = U + F x,
    b = dim((U_plus+W)/W).

Before exposing the current coefficients, b is either r or r+1, and b<=d.
For a quiet incoming state, b=1. Every output or next-state zero equation
in this step belongs to U_plus. Therefore each independent equation
decreases both d and b by one. Dependent exceptional-root equations
decrease neither.

There is also an exact relation at the end of the step. Write A=R+D_t.
For every allowed coefficient choice, A is invertible and

    Q_next = A Q + e_j x.

Consequently

    span(Q_next,x) = span(Q,x) = U_plus.           (1)

The equality holds before and after quotienting by the final W of the
step. If b now denotes its residual dimension, the outgoing state rank
satisfies

    r_next is either b-1 or b.                     (2)

It must also satisfy 0<=r_next<=min(w,d,|H_next|), with rank zero exactly
when H_next is empty. In particular, a quiet outgoing state requires b=1.
This is stronger than allowing every rank at most r+1 independently at
each boundary.

Equation (1) follows in both directions: Q_next is in span(Q,x), and
Q=A^(-1)(Q_next-e_j x) is in span(Q_next,x). An independent zero equation
lies in U_plus, so its residual dimension drops by exactly one. Earlier
nonzero requirements can cause additional rejection, never an increase.

## 3. Joint outcomes at the active coordinate

The active coordinate uses one coefficient for both output and next
state. They must not be treated as two independent coefficient tests.
In the current quotient by W, write

    s = Q[j],       v = (R Q)[j],
    y = x + a s,    p = v + y,

where a is still uniform on A=theta+E and independent of the exposed
history. All indicated nonzero input and incoming-state forms remain
required nonzero.

Let F_C(d,b), F_Y(d,b), and F_P(d,b) denote nonnegative continuation bounds
for the following outcomes:

| Outcome | Output y | Next coordinate p |
|---|---|---|
| C | nonzero | nonzero |
| Y | zero | nonzero |
| P | nonzero | zero |

These functions include the corresponding output marker: C and P receive
one factor z, whereas Y receives none. They also include the chosen next
support coordinate and all remaining within-step work. Write

    F_B^- = F_B(d-1,b-1).

A term involving F_B^- is allowed only when b>=2. No branch renormalizes
the coefficient distribution after imposing a nonzero requirement.

### Zero current state coordinate: s=0

Now y=x is nonzero and a has no effect.

- If v=0, only C is possible.
- If v is nonzero and x=v modulo W, only P is possible, without a rank gain.
- If x and v are nonzero proportional forms with x!=v modulo W, only C
  survives. Imposing p=0 would kill a required nonzero form.
- If x and v are independent modulo W, the sum of the two support branches
  is bounded by F_C(d,b)+(1/tau)F_P^-.

Thus the last case spends one global rank on its zero branch. Dropping
this deterministic rank was one relaxation in the support-only transfer.

### Zero rotated state coordinate: s!=0 and v=0

Here p=y. There are only two outcomes: both nonzero, denoted N, or both
zero, denoted Z. Their output factors are z and one respectively.
The complete branch sum has one of the following bounds:

    F_N(d,b) + (1/tau) F_Z^-             if x,s are independent;
    (1-1/q)F_N(d,b)+(1/q)F_Z(d,b)        if x=kappa s and kappa is in A;
    F_N(d,b)                             if x=kappa s and kappa is outside A.

Take a maximum over these cases only when the quotient does not retain
which case holds. The independent case requires b>=2.

### Both state coordinates nonzero: s!=0 and v!=0

The simultaneous outcome y=p=0 is impossible, because p+y=v is required
nonzero. The other three outcomes are coupled as follows.

**Both affine zero tests are dependent.** If x and x+v both belong to F s
modulo W, let their potential roots be alpha and beta. They are distinct
because v is nonzero. For the subset T of these roots that belongs to A,
the branch sum is bounded by

    (1-|T|/q) F_C(d,b) + (1/q) sum_{B in T} F_B(d,b),       (3)

where the two possible labels are Y and P. There is no rank gain. A
support-only envelope may maximize over all four subsets T; this safely
allows some root-membership patterns that a more detailed quotient could
exclude.

**Only the y test is dependent.** Suppose x belongs to F s modulo W, but
x+v does not. Let epsilon be one if the root of y belongs to A, and zero
otherwise. Then the bound is

    (epsilon/q) F_Y(d,b)
      + (1-epsilon/q) [F_C(d,b)+(1/tau)F_P^-].              (4)

At the exceptional coefficient, y=0 identically in the quotient and
p=v is nonzero. At every other coefficient, imposing p=0 adds one rank.
The P branch is impossible at the exceptional coefficient: it would
impose v=0. This is why its factor is also 1-epsilon/q.

**Only the p test is dependent.** Interchanging Y and P gives

    (epsilon/q) F_P(d,b)
      + (1-epsilon/q) [F_C(d,b)+(1/tau)F_Y^-].              (5)

**Both affine zero tests are independent.** If neither x nor x+v belongs
to F s modulo W, each individual zero branch adds one rank. A safe bound is

    F_C(d,b) + (1/tau)F_Y^- + (1/tau)F_P^-.                (6)

Some coefficient choices may make a zero branch contradict another
required nonzero form. Ignoring these additional contradictions only
increases (6). The two zero branches are alternative support histories,
not simultaneous equations; each spends one rank, not two.

The actual quotient selects one of (3)--(6) before a is exposed. A uniform
Bellman envelope takes their maximum, including both epsilon choices and
the possible root subsets. It must not add these classifications or
maximize each zero/nonzero branch independently. Cases (4)--(6) require
b>=2. At b=1 only (3) is available.

## 4. A finite-state continuation

A possible complete-step state is (d,r,H). For D=w=4 there are 108 valid
states, including the quiet state at each d. Within an active step, also
keep the temporary augmented rank b, the partial next support, and the
output weight accumulated so far.

For a singleton step, choose a fixed coordinate-exposure order. All
inactive coordinates use the empty-input zero/nonzero rule. Their output
supports are fixed by the incoming H. The active coordinate uses the joint
cases above. Every independent zero changes (d,b) to (d-1,b-1), with a
factor 1/tau. Retain every earlier nonzero requirement throughout.

At the end, use the maximum terminal continuation over the valid ranks in
(2), rather than their sum. The actual realization has only one outgoing
rank. A support H_next with no valid rank contributes zero. A quiet
incoming state has the simpler exact update: one nonzero output, singleton
next support, state rank one, and unchanged global d.

In an empty step, the augmented rank is just the state rank. The existing
empty-span Bellman rule applies, but each independent zero now decreases
both d and r. Never reset d at a region boundary or when the state becomes
quiet. At every step r<=d; at most D-1 independent equations can survive
across the complete encoder.

These rules are uniform over all quotients having the stated ranks and
supports. For an independent case, each allowed coefficient contributes
at most the indicated rank-weighted continuation. For a dependent case,
the exceptional-root probabilities are computed in the original coefficient
space. Sequential maxima therefore give an upper bound by conditional
expectation, even though earlier equations change later classifications.

The resulting Bellman map is generally nonlinear. A route average cannot
be moved inside its maximum. One safe implementation would verify a common
positive linear envelope for every local classification before averaging
regional insertion positions. Another would average fixed-route bounds
without interchanging the maximum and expectation.

## 5. Intended use and remaining work

For one group occupying exactly 33 complete RS columns, D=4. The proposed
operator would preserve its budget of three independent equations across
all 132 occupied regions and their intervening gaps. It could also retain
the deterministic cancellation rank at a singleton input and exclude
incompatible state-rank choices at successive active boundaries.

This is a reserve refinement, not an explanation of every earlier loss.
The simpler parameter-separated empty-span envelope can be tried first.
Before using the new operator for a certificate, its joint table needs
independent small-field checks, its common route envelope needs verification,
and its free-terminal/source-count accounting must be established.

The next application beyond complete columns is a partially supported
one-group class. Its exact shortened dimension must be supplied separately;
the local operator does not prove an outer spectrum or justify replacing
the actual pair route by independent coordinates.
