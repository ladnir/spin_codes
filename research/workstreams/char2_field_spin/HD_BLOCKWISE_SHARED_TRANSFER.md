# A rank transfer for the whole joint state/feedback image

The shared-transform proof must control a family of source vectors without
assigning fresh independent amplitudes to them. This note gives one such
conditional bound. It analyzes the entire joint state/feedback image after
the current output-zero equations, before sampling the current mixer.

The result is a positive one-step rank bound and a capped return event. It
does not multiply a point probability by the number of source words. It is
not yet a closed whole-code distance proof.

## Chronology and joint image

Let F be a field of characteristic two and size q, and let the state have s
coordinates. The intended
case is q=65536 and s=4. Condition on the route, all earlier randomness, and
the linear output-zero constraints of the current packet. The current packet
output uses the old state, so those constraints do not depend on its mixer.

For the surviving source space, let

    U = image of (Q,CX) in F^s x F^s.

This is one joint image. Its two components may have arbitrary dependencies.
Repeated source columns and source directions invisible to both components
are quotiented out. Write d=dim U. A vector of U is denoted (h,c).

Sample u from F^s excluding zero. Sample v from F^s conditioned on v^T u !=1,
independently of the conditioned source history. The successor-state map is

    R_(u,v)(h,c) = h+c+u(v^T h).

This is the actual M=I+uv^T mixer followed by feedback. The condition on v
makes M invertible. Individual zero entries of u or v remain allowed.

The theorem concerns the image rank of R on all of U. A later full-support
restriction can only remove candidate return vectors. It cannot justify
treating the remaining amplitudes as independent random values.

## Geometry before the mixer

Define

    D(h,c)=h+c,                 W=image D,       r=dim W,
    E={h : (h,h) belongs to U},                  e=dim E,
    k=dim(E intersect W).

The kernel of D is the diagonal copy of E, so d=r+e.

For w in W, choose h_w such that (h_w,h_w+w) belongs to U. The choices differ
by E. Thus the following quotient maps are well-defined and linear:

    H: W -> F^s/E,       H(w)=h_w+E,
    J: W -> F^s/E,       J(w)=w+E.

Only their three ranks will be needed:

    a0=rank H,       a1=rank(H+J),       a01=rank(H,J).

They can be computed directly from a joint matrix, without enumerating w or
choosing a quotient basis. In particular, a0=rank(Q)-e and
a1=rank(CX)-e. To obtain a01, quotient both blocks of U by E: append the
independent column spaces (E,0) and (0,E), take the resulting rank, and
subtract 2e.

## Exact support of the successor rank

The successor rank can only be r-1, r, or r+1.

To prove this and identify the exceptional cases, first fix u. If u is
outside W, a return vector must satisfy D(h,c)=0 and v^T h=0. Its kernel
dimension is e-rank(v restricted to E). The rank increases to r+1 exactly
when v does not annihilate E.

If u belongs to W, write a return vector as

    D(h,c)=t u,      h=t h_u+a,      a in E.

The return equation becomes one scalar equation on E x F:

    (1+v^T h_u)t + v^T a =0.

Its kernel has dimension e, except when both coefficients vanish. Therefore
the rank decreases to r-1 exactly when

    u in W,       v annihilates E,       v^T h_u=1.

Otherwise it remains r. These statements use the same joint source space;
there is no independent choice of h and c.

## Rank-up probability

Put N_q=q^s-1. The rank-up probability is exactly

    p_up = ((q^s-q^r)(1-q^(-e))
             -(q^e-q^k)q^(-e)/(q-1)) / N_q.

There are q^s-q^r choices of u outside W. For u outside E, the conditional
probability that v annihilates E is q^(-e). For u in E excluding W, it is
q^(1-e)/(q-1): annihilating E already makes v^T u=0, so the forbidden value1
removes no vectors. There are q^e-q^k such u. Subtracting the annihilation
probabilities proves the formula, including e=0.

## Rank-down probability from three quotient ranks

Define the integer

    G=q^r-q^(r-a0)-q^(r-a1)+q^(r-a01).

It counts the vectors w in W for which H(w) is neither zero nor J(w).
This follows by inclusion-exclusion on the kernels of H and H+J. Their common
kernel is the kernel of (H,J). The zero vector is excluded automatically.

For any other w, the rank-down event is impossible. For a counted w, condition
on u=w. Among v satisfying v^T u !=1, the probability of

    v annihilates E and v^T h_u=1

is either q^(-e-1) or q^(-e)/(q-1). The first case occurs when the quotient
vectors H(w) and J(w) are independent. The second occurs when they are
dependent, including J(w)=0. They are not equal, and H(w) is nonzero.

Consequently

    G/(q^(e+1) N_q) <= p_down <= G/(q^e(q-1) N_q).

The upper bound can also be capped by 1-p_up. Call the resulting lower and
upper bounds l and h. They differ by at most the factor q/(q-1), before that
cap. At q=65536, no enumeration of65536 possible eigenvalues or field vectors
is needed.

## Positive domination and rank moments

For every nonnegative function f of successor rank,

    E[f(rank R)]
      <= h f(r-1) + (1-p_up-l)f(r) + p_up f(r+1).

A term with zero coefficient is omitted. In particular, no negative-rank
term is needed at r=0. Every displayed coefficient is nonnegative. Their sum
is 1+h-l; this is an upper operator, not an asserted probability law.

For a rank-moment mark tau>=1, a slightly tighter bound uses monotonicity:

    E[tau^(-rank R)]
      <= tau^(-r) * (1+p_up*(tau^(-1)-1)+h*(tau-1)).

These inequalities are conditional on any prior history that fixes U before
the mixer. They can therefore be used successively with conditional
expectations. For example, a future continuation bounded uniformly by f(h)
for every reachable source configuration with state-image rank h satisfies
the positive domination above.

This does not make rank alone a sufficient global state description. To
construct such a uniform continuation bound, one must retain or safely bound
the future source forms, output constraints, and outer-code dependencies.
Plugging in a continuation valid only for independent amplitudes would be
invalid. The present result supplies the mixer step, not that missing closure.

## A return event capped before counting source words

Let Return(U) mean that some nonzero (h,c) in U satisfies R(h,c)=0.
Because M is invertible, every such pair has h!=0. Thus it represents a
nonquiet return at the level of the joint image.

The event is equivalent to rank R<d. Hence

    e=0:  Pr[Return(U)] = p_down <= h;
    e=1:  Pr[Return(U)] = 1-p_up;
    e>=2: Pr[Return(U)] =1.

Each whole joint subspace is counted once. A large returned kernel does not
introduce a factor for all of its vectors or projective lines. This local cap
does not by itself control a union over many source-support histories or the
rare-route mechanism in the whole code.

The last case exposes a concrete limitation of rank-one updates. If U contains
a two-dimensional diagonal subspace (h,h), then a nonzero h in that subspace
always satisfies v^T h=0. It returns to zero for every sampled mixer. One
cannot charge a q^(-4) return probability to that whole subspace.

The statement is about unrestricted vectors in U. Required nonzero input
coordinates or earlier liveness constraints can remove those vectors. The
source-embedding test from
[HD_BLOCKWISE_SHARED_SUPPORT.md](HD_BLOCKWISE_SHARED_SUPPORT.md) is still
needed before calling this an actual admissible codeword.

## The saved shift64 cancellation line

The rank-three residual example from the support note has a one-dimensional
joint image generated by

    h=(34,232,1,0),       c=(72,27,211,244).

These vectors are noncollinear, also over GF65536. Here e=0 and r=1. The exact
return probability under the fresh mixer is

    (1-1/q)/(q^4-1).

The upper rank formula gives 1/(q^4-1), while its lower bound equals the exact
probability. This recovers the useful four-field-equation scale for this
entire source line without multiplying by its q-1 nonzero vectors.

By contrast, a one-dimensional diagonal joint image has return probability
q^(s-1)/(q^s-1), approximately1/q. The subspace geometry distinguishes these
cases before any full-code enumeration.

## Implementation and validation

[theory/hd_blockwise_shared_transfer.py](theory/hd_blockwise_shared_transfer.py)
computes the joint geometry, positive rank operator, rank moment, and capped
return bound. Its field arithmetic includes the actual
GF256[u]/(u^2+u+0x20) representation of GF65536. No subfield-only assumption is
made about source matrices.

Five tests pass. The principal test exhausts all67 subspaces of GF(2)^4 and
all529 subspaces of GF(4)^4. For every subspace, it enumerates every allowed
(u,v) pair for state dimension2 and checks the exact rank-up probability,
rank-down interval, positive domination, moment bound, and return event.
Additional tests cover a two-dimensional diagonal subspace, the saved shift64
line, genuine GF65536 coefficients, and duplicated source columns.

    python -B -m unittest hd_blockwise_shared_transfer_tests -v
    python -B hd_blockwise_shared_transfer.py

The next proof task is a bounded collection of reachable joint-image classes,
with a justified rule for joining source support, output residual planes, and
feedback images. The exact mixer operator should be applied before any
field-word count. No construction change or native experiment is required
to test that interface first.
