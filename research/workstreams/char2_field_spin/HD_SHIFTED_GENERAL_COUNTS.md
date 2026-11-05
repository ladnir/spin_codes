# General support counts for alternating adjacent pairs

The alternating route needs the joint multiplicity of both matching
histograms. Counting each matching separately loses their shared cyclic
geometry. The formulas below cover every indexed support, including
unequal numbers of active even and odd groups. They do not establish a
whole-code distance bound.

Fix a positive integer M and a support on the 2M indexed cyclic positions.
Write its bits as e_i,o_i for positions 2i,2i+1, with indices i modulo M.
Define

    E = sum_i e_i,                 O = sum_i o_i,
    a = sum_i e_i o_i,             b = sum_i o_i e_(i+1),
    r = E+O,                      c = r-a-b.

Let N_M(E,O,a,b) count these indexed supports. Rotations are counted
separately. Matching 0 uses (e_i,o_i); matching 1 uses (o_i,e_(i+1)).
Their histograms, ordered as (00,10,01,11), are

    H0 = (M-r+a, E-a, O-a, a),
    H1 = (M-r+b, O-b, E-b, b).

For a nonconstant support, c is its number of cyclic runs of ones.
The two constant supports have c=0.

## Exact count and feasibility

For a binary cycle of length m, let B(m,k,s) count indexed strings of
weight k with s runs. For 0<k<m,

    B(m,k,s) = m/s * binom(k-1,s-1) * binom(m-k-1,s-1)

when 1<=s<=min(k,m-k), and it is zero otherwise. Define
B(m,0,0)=B(m,m,0)=1; all other constant-weight cases are zero.

Fix the even-bit cycle and let s be its run count. Its transition counts
between consecutive even bits are

    (n00,n01,n10,n11) = (M-E-s,s,s,E-s).

An odd bit lies between each such transition. Let t count selected odd
bits on 11 transitions. The required numbers selected from the four
transition types are

    (x00,x01,x10,x11) = (O-a-b+t,b-t,a-t,t).

Consequently, with binomial coefficients outside their range defined as
zero, the exact joint count is

    N_M(E,O,a,b)
      = sum_s B(M,E,s)
          * sum_t binom(M-E-s,O-a-b+t)
                  binom(s,b-t) binom(s,a-t) binom(E-s,t).

Every fixed even cycle has the stated transition counts. Choosing its odd
bits independently within the four transition classes gives the product
of binomials. Each support has one even cycle and one such selection, so
the sum counts every support once.

The formula also gives a constant-time feasibility test. First require
0<=E,O,a,b<=M and nonnegative entries in H0 and H1. Define

    s0 = 1 if 0<E<M, and 0 otherwise,
    s1 = min(E,M-E),
    t0 = max(0,a+b-O),             t1 = min(a,b),
    u0 = max(a,b,s0+t0),           u1 = min(E,M-c,s1+t1).

The class is feasible exactly when u0<=u1. Indeed, the placement
constraints are equivalent to

    s0<=s<=s1,     t0<=t<=t1,     max(a,b)<=s+t<=min(E,M-c).

The attainable integer sums s+t fill an interval. For a feasible class,
s=max(s0,u0-t1) and t=u0-s give a witness.
This condition excludes classes that the two histograms and total run
count admit. For example, M=4,E=O=2,a=b=0 is impossible: a maximum
independent support on an eight-cycle must occupy one parity entirely.

The reference implementation evaluates the exact count with at most
O(M^2) binomial products. It shortens the sum by using parity exchange
and complementation, which preserve counts. Full pairs have a simpler
formula: if E=O=a or E=O=b, then

    N_M(E,O,a,b) = B(M,E,c).

## Integer upper bounds

The following four counts give sound bounds for a feasible class:

    binom(M,E) binom(M,O),
    M! / ((M-r+a)! (E-a)! (O-a)! a!),
    M! / ((M-r+b)! (E-b)! (O-b)! b!),
    B(2M,r,c).

The first fixes parity weights. The next two fix one matching histogram.
The last fixes the total weight and cyclic run count. The function
`joint_count_upper(M,E,O,a,b)` returns their minimum, handles the two
constant supports separately, and uses exact contraction for full pairs.
It returns zero for every infeasible class.

A further bound restricts the runs of the even-bit cycle. Define

    slo = max(s0,abs(a-b)),
    shi = min(E,M-E,min(E,M-c)-t0),
    S_E = sum_(s=slo)^shi B(M,E,s),
    D_E = min_(h in {a,b}) binom(E,h) binom(M-E,O-h).

For each even cycle, either matching alone permits its corresponding
conditional binomial number of odd selections. Thus N_M<=S_E D_E.
Exchange E and O for another bound. The function `joint_run_count_upper`
takes the minimum with the preceding bounds. It sums the run counts in
O(M) integer operations using

    B(M,E,s+1)/B(M,E,s) = (E-s)(M-E-s)/(s(s+1)).

At M=4096, reference values are:

| (E,O,a,b) | log2 exact count | log2 cheap upper |
| --- | ---: | ---: |
| (512,512,512,384) | 1200.816176290 | 1200.816176290 |
| (512,480,450,420) | 1383.291061098 | 1413.500493891 |
| (512,512,400,400) | 2194.415539248 | 2201.869382163 |
| (512,512,128,128) | 4331.596230503 | 4340.560296640 |

The integers are exact. These displayed logarithms are ordinary floating
evaluations and are not directed enclosures.

## Covering classes without listing all quadruples

Introduce positive markers x,y,p,q for E,O,a,b. Define the 2-by-2 matrix

    K[u,v] = x^v (1 + y p^u q^v),       u,v in {0,1}.

Its trace generating function is

    G_M(x,y,p,q) = trace(K^M)
                = sum_(E,O,a,b) N_M(E,O,a,b) x^E y^O p^a q^b.

Expanding a trace term chooses the even cycle. Choosing either summand
of each K entry chooses its intervening odd bit. This proves the
coefficient identity, including M=1 and both constant supports.

A profile box Q specifies an inclusive integer interval for each of
E,O,a,b. For positive rational markers, let m_Q be the minimum of
x^E y^O p^a q^b over the box. Each coordinate contributes its lower
endpoint when its marker is at least one, and its upper endpoint
otherwise. Positivity gives

    sum_(profile in Q) N_M(profile) <= floor(G_M(x,y,p,q)/m_Q).

The implementation evaluates this bound exactly with rational matrix
arithmetic and O(log M) matrix products. `box_count_upper` also takes
the minimum with the product of the two parity-binomial interval sums.
The routine neither lists profiles nor assumes that every point in Q
is feasible. Rational markers chosen by a floating search remain sound
when the final evaluation uses these exact operations.

This count can support an adaptive proof cover. Suppose each profile h
has a nonnegative scalar C(h), regional nonnegative operators R0(h),
R1(h), and a valid contribution bound

    N_M(h) C(h) e0 (R0(h) R1(h))^128 1,

where e0=(1,0) and 1=(1,1)^T. The scalar includes any outer-count and
distance-marking factors. For a box Q, bound C(h)<=C_Q and
Rj(h)<=Sj,Q entrywise for every feasible h in Q. The entire box then
contributes at most

    box_count_upper(M,Q) C_Q e0 (S0,Q S1,Q)^128 1.

All operator products retain their order. A direct entrywise envelope
of R0(h)R1(h) can replace the product of separate envelopes and may be
tighter. Numerical regional envelopes need directed arithmetic and
parameters valid uniformly throughout their box.

For all nonempty supports, start with these disjoint boxes:

    E=0,     1<=O<=M,  a=b=0;
    1<=E<=M, 0<=O<=M,  0<=a,b<=M.

Bisect boxes whose contribution bounds remain too large. Retain accepted
boxes and sum their bounds, so each feasible class appears exactly once.
`ProfileBox.split()` provides a disjoint integer bisection. The count
bound works even for the initial boxes; it has no O(M^4) enumeration
step. Achieving a useful distance margin depends on tight regional
envelopes. This approach does not promise that every instance admits a
small successful cover.

## Symmetries and the initial state

Three transformations preserve support counts:

    (E,O,a,b) -> (O,E,a,b),
    (E,O,a,b) -> (E,O,b,a),
    (E,O,a,b) -> (M-E,M-O,M-E-O+a,M-E-O+b).

The first is reflection across an even matching edge. The second is
reflection across a vertex. The third complements every support bit.
These identities justify reductions in exact counting only.

Exchanging E and O swaps the singleton orientations in both histograms.
It also preserves the regional operators if their construction is
invariant under exchanging the two lanes. For an h4 polynomial built
from a symmetric local occupancy kernel, this follows by bijectively
swapping the two lane occupancies in each half. That operator property
must be established for the actual regional bound being used.

Exchanging a and b exchanges the two regional operators once lane
symmetry holds. The initial-zero moment need not be invariant under
that exchange. For example, take

    A = [[1,2],[0,1]],       B = [[1,0],[3,1]].

Then e0 AB 1=9, but e0 BA 1=3. Repeating alternating blocks does not
create an operator-exchange identity. A support-count symmetry alone
therefore cannot halve the ordered (a,b) search. Complementation also
requires a separate operator theorem before use in a distance bound.

## Validation and remaining work

`hd_shifted_general_count_tests.py` compares the formula and DP with
exhaustive cycles through M=6, including every impossible quadruple.
It checks marginals and all three counting symmetries through M=8.
Additional tests cover feasibility witnesses, rational transfer values,
disjoint box coverage, box bounds, and large full-pair reference cases.

The next step is to construct directed regional envelopes for the
general histograms and use them in a disjoint cover of all nonempty
supports. Neither the exact counts nor isolated regional witnesses
establish the requested 35% distance or a 40-bit whole-code margin.
