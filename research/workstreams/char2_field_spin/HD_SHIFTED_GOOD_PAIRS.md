# Composing steps with few singleton pairs

The general single-step bound loses information when it applies the
fractional exponent before summing over an intermediate state. The loss
can remain significant after improving support counts and alias grids.
The following bound recovers that information for steps with at most
one singleton pair. It changes the proof calculation, not the encoder.

For the support profile (E,O,a,b)=(512,480,450,420), the initial FLOAT
grid-64 bound was +425.179068573 in log2 units. At the same witness,
the new bound is -394.206707880. These are screens, not directed
receipts or whole-code certificates.

## A positive two-step bound

Fix 0<z<1 and 0<alpha<=1. For each local support arrangement gamma,
let A_gamma be its nonnegative 2-by-2 h4 state-operator cap before
applying alpha. The cap includes the established density and quiet
factors. Write A^[alpha] for entrywise exponentiation.

A step contains 32 source pairs. Call a local arrangement good when
at most one pair is a singleton. Use variables x10,x01,y for the
numbers of 10,01,11 pairs. Define the matrix polynomials

    G = sum_(good gamma) monomial(gamma) A_gamma^[alpha],
    D = sum_(other gamma) monomial(gamma) A_gamma^[alpha].

Each arrangement is counted once. Equivalently, group arrangements with
equal degree and operator cap and multiply by their exact multiplicity.
The existing step polynomial is S=G+D.

For two consecutive good arrangements, compose their state operators
before applying alpha. Define

    C = sum_(good gamma,good eta)
          monomial(gamma) monomial(eta) (A_gamma A_eta)^[alpha],
    T = C + GD + DG + DD.

The four summands partition the cases for a consecutive pair of steps.
For nonnegative matrices A,B, the scalar inequality
(u+v)^alpha<=u^alpha+v^alpha gives

    (AB)^[alpha] <= A^[alpha] B^[alpha]

entrywise. Thus T bounds every two-step case. It also satisfies
T<=S^2 coefficientwise, and equality holds when alpha=1. All
coefficients of T are nonnegative. No commutation of state operators
is used.

Fix a regional histogram H=(n00,n10,n01,n11) summing to 4096.
A uniform pair permutation induces the uniform distribution on the
multinomial number of arrangements of these four types. Consequently,
a valid regional operator is

    R(H) = [x10^n10 x01^n01 y^n11] T^64
           / multinomial(4096;n00,n10,n01,n11).

This enumerates all 128 steps, grouped consecutively in pairs. The
coefficient retains the dependence caused by sampling without replacement
inside a region. It does not assume that adjacent step types are
independent. The outer comparison measure and independent regional
permutations remain those of the existing alternating-route proof.

The two matching histograms may use different positive coefficient
markers. Their regional bounds still compose in the order (R0 R1)^128,
starting in state zero. Equal histograms may reuse a deterministic
operator calculation without identifying their random route permutations.

## Reducing the good geometry

Let H(l0,l1,l2,l3) denote the local cap for the four lane occupancies.
The h4 cap depends symmetrically on these occupancies. Each step has
two halves, each containing 16 source pairs.

A clean step has full-pair counts u,v in its halves. Sorting u<=v
gives 153 groups with

    cap       H(u,u,v,v),
    full      u+v,
    count     binom(16,u) binom(16,v) * (1 if u=v else 2).

For one specified singleton orientation, let u count full pairs in
the half containing that singleton, and let v count full pairs in
the other half. The ranges 0<=u<=15 and 0<=v<=16 give 272 groups:

    cap       H(u,u+1,v,v), with the occupancies sorted,
    full      u+v,
    count     32 binom(15,u) binom(16,v).

The factor 32 selects the singleton's half and its position. Changing
its orientation exchanges two lane occupancies, so the cap is unchanged.
Thus the good polynomial has the form

    G = G0(y) + (x10+x01) G1(y).

This identity retains both singleton orientations. In particular,
independent markers for x10 and x01 remain valid. It does not replace
their two target degrees with their sum.

Only 425 groups are needed to construct C. Their 180625 ordered
pairs produce a polynomial of degree at most two in X=x10+x01 and
at most 64 in y. Expanding X^2 gives the exact factor two for one
singleton of each orientation. The four matrix products remain ordered.

The unmarked scalar multiplicities provide independent checks:

    sum_(clean,full=j) count = binom(32,j),
    sum_(one specified singleton,full=j) count = 32 binom(31,j).

The second formula vanishes at j=32.

## Coefficient evaluation and validation

The FLOAT implementation constructs D by setting the zero-singleton
and one-singleton coefficients of S to zero. It never obtains D by
subtracting nearly equal evaluations. It evaluates G and C from their
small X,y representations and forms T with its four positive terms.

On a grid of size 64, the y^64 coefficient of C must fold to index
zero before the FFT. Dropping that coefficient would change the
polynomial. The implementation folds all degrees modulo the grid.
After evaluating T, it raises T to the 64th power and extracts a positive
alias sum for the target coefficient. Positive markers cancel through
the same regional multinomial normalizer as in the unpaired bound.

The six tests in `theory/hd_shifted_good_pair_screen_tests.py` check the
exact geometry multiplicities, agreement with the
full oriented step polynomial, alpha=1 equality, coefficientwise
improvement, noncommuting matrix products, and marker orientation with
degree-64 folding. Disabling composition before alpha returns
+425.179068573 at the reference witness, reproducing the ordinary
single-step alias calculation.

The reference witness is

    z     = 0.49367570922234283,
    alpha = 0.09802545440976157,
    even marker log2 values
          = (-5.408942950886019,-6.4540781593556,-2.2088270710712097),
    odd marker log2 values
          = (-5.44207015077779,-4.8278390848501544,-2.296457984764048).

At grid 64 the new normalized alias entries are at least
1.796e-5 for matching zero and 1.028e-5 for matching one. These
floating diagnostics are above the screen's cancellation floor. They
do not replace directed rounding.

The next steps are a directed replay of this positive decomposition
and stress tests over other support geometries. A whole-code claim
still requires a disjoint cover of every feasible nonempty support
class and a bound on their total contribution.
