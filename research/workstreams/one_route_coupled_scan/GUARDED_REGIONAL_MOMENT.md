# Regional moments with an outer-support guard

Status: an analytic route-moment and complete-count interface, not a
whole-code certificate. The construction, coefficient distribution, and
physical packet route are unchanged. No numerical result is claimed here.

The regional information-set cap in [Outer compatibility](AFFINE_COMPATIBILITY.md)
handles concentrated dense supports. It cannot replace route averaging:
some valid routes admit short words for every inner setup. Sparse activity
needs the probability that its allowed input coordinates support a nonzero
outer word. The bound below retains that event and the actual regional
densities. It does not replace them by one global density.

## Objects and probability space

Use the outer and route from [Partial-packet moments](SHUFFLED_PACKET_ROUTE.md).
There are G groups of p independent MDS [B,k] rows, where k=B/2, and
M=pG coordinates per region. Thus N=MB. Each region independently permutes
the G groups. The implemented route has fixed packet-lane order.

First analyze the auxiliary distribution that additionally shuffles the
p lanes independently in every destination packet. All choices are sampled
once and fixed for all messages. The last section transfers the bound to
the appropriate mask sum for the implemented route.

Fix coordinate sets T and U independently of these route choices. Require
T subset U and require U to be a union of whole destination packets. T may
cut packets. For example, a complete output/state support profile fixes a
union U of whole active inner batches. Both the random inverse support A
and the union support S=A union W are subsets of U.

For a route pi, let D_pi(T) be the dimension of the routed outer shortened
to T. The exact MDS formula is the sum of (row occupancy minus k)_+.
If an eligible nonzero output has inverse input supported in U, injectivity
implies D_pi(U)>0. This is an existence guard, not a condition under which
the route distribution is renormalized.

Put alpha_j=|T intersect region j|/M and upsilon_j=|U intersect region j|/M.
For each region j, define an independent pair of indicator variables
(A_j,V_j) with probabilities

    Pr[(A_j,V_j)=(0,0)] = 1-upsilon_j,
    Pr[(A_j,V_j)=(0,1)] = upsilon_j-alpha_j,
    Pr[(A_j,V_j)=(1,1)] = alpha_j.

These pairs are a comparison distribution for one source row; they are not
claimed independent across rows of the actual packet route. Set

    Z=sum_j A_j,       Q=sum_j V_j,       t=tau^p,
    H=E[t^((Z-k)_+)],
    J=E[1[Q>k] t^((Z-k)_+)],                    tau>=1.       (1)

Both expectations are finite bivariate Poisson-binomial sums. Their
generating polynomial is

    product_j [(1-upsilon_j)+(upsilon_j-alpha_j)y+alpha_j xy].

## Guarded moment

Under the auxiliary shuffled route,

    E_route[tau^D(T) 1[D(U)>0]]
      <= min(H^G, G J H^(G-1)).                               (2)

Here G>=1. The guard can vanish exactly, for example when U touches at
most k regions. No division by the probability of the guard occurs.

To prove (2), record each destination packet by its selected counts
(a,u) for T and U. Because U is whole-packet and T subset U, the possible
marks form the ordered chain

    (0,0), (0,p), (1,p), ..., (p,p).

After averaging a uniform a-subset of source-row lanes, a nonnegative
coordinatewise increasing group factor remains increasing along this
chain. Couple the subsets by prefixes of one uniform ordering. At the
first transition the T subset remains empty and every U indicator changes
from zero to one.

The monotone without-replacement comparison in the linked note therefore
replaces regional assignments of these marks by independent draws, one
region at a time. It applies both to an ordinary group factor
tau^(sum_l(Z_l-k)_+) and to that factor multiplied by 1[Q>k]. For a whole
U packet, all p rows of a group have the same U occupancy Q.

After replacement, source groups are independent. For an unmarked group,
Holder's inequality across its p rows bounds its expectation by H.
For a marked group, keep the common indicator inside Holder:

    E[1[Q>k] product_l tau^((Z_l-k)_+)]
      <= product_l (E[1[Q>k] tau^(p(Z_l-k)_+)])^(1/p)
       = J.                                                   (3)

The p row marginals are identical. Since D(U)>0 exactly when at least one
group has Q>k, a union over the G marked groups gives G J H^(G-1).
Dropping the guard gives H^G. This proves (2) without treating packet
lanes as independent and without subtracting two upper bounds.

## A separable positive envelope

For eta>=1 define, separately in each region,

    a_j = 1+(t-1)alpha_j,
    b_j = 1+(eta-1)upsilon_j,
    c_j = b_j+eta(t-1)alpha_j.

The pointwise inequalities

    t^((Z-k)_+) <= 1+t^(Z-k),
    1[Q>k] <= eta^(Q-k-1)

give

    H <= H_hat := 1+t^(-k) product_j a_j,
    J <= J_hat := eta^(-k-1)
                       [product_j b_j+t^(-k) product_j c_j].  (4)

Consequently the marked bound in (2) is at most

    G eta^(-k-1) sum_(r=0)^(G-1) binom(G-1,r)t^(-kr)
      [product_j(a_j^r b_j)
        +t^(-k) product_j(a_j^r c_j)].                         (5)

This is a positive sum, not a cancellation-based approximation. At fixed
r and either bracketed term, its regional factor depends only on the two
regional cardinalities |T_j| and |U_j|. Thus a complete regional support
operator can compose these factors without storing a vector of all B
regional occupancies. Formula (5) is coarser than the exact H/J bound.
It still requires summing r and every support profile; separability is
not a certificate of numerical convergence or a runtime estimate.

## Use in the actual unshuffled encoder

Take a complete family of output/state mask pairs closed under joint
packet-lane permutations, after relaxing realizability to invariant rank
and amplitude bounds. A whole-packet U is fixed by every such permutation.
The orbit identity therefore applies to

    tau^D(T) 1[D(U)>0]

as well as to the unguarded moment. Consequently (2) bounds its summed
expectation under the implemented route. It does not bound a fixed T under
that route. All quantities in (1) depend only on the regional cardinalities
and are invariant within the mask class.

For H grouping, fix A first and sum W masks before fractionalization.
For S grouping, fix S=A union W first and sum its eligible mask pairs.
The same guarded moment applies with T=A or T=S. Keep exactly the earlier
multiplicity rules: only counts grouped at fixed T enter the fractional
power; the number of T masks stays outside. The projective correction,
inner rank bound, output-weight restriction, activation constraints, and
free terminal remain unchanged. U is fixed by the absolute per-batch
support profile, not resampled for each word.

A complementary moment or deterministic regional cap may separately
upper-bound the same guarded expectation by dropping the guard. Taking
the minimum is valid. Do not apply the guarded lemma directly to T^c:
usually T^c is not a subset of U, and its packet marks do not form the
ordered chain used in the proof. Likewise, scalar-cut U does not satisfy
the shared-Q step in (3).

## Bounded check and remaining gate

`theory/test_guarded_regional_moment.py` defines exact Fraction checks on
all 4096 labeled routes for p=2, B=4, G=2, for a small fixed list of nested
(T,U) supports. It checks the bivariate marginals, (2), (4), (5), and
guarded orbit identities. A scalar-cut U negative control disproves an
extension without the whole-packet hypothesis. The test is not executed
by this note; the root agent owns runs with a 30-second/256-MiB cap.

The remaining gate is a complete amplitude-and-support count with the
guard retained. Formula (5) supplies separable regional weights for that
calculation; it neither omits sparse activity classes nor proves that the
resulting upper bound is below one. Preserve the deterministic regional
cap as a complementary dense-support bound, not a uniform-route theorem.
