# A profile-measure bound for dense two-bit inputs

This bound applies to the same two-bit ensemble as [README.md](README.md).
It does not replace the regional route by a global permutation or change
the inner. The iid inputs introduced below are a comparison measure used
in the proof. Every conditioning factor is explicit.

## Bound an entire outer measure once

Let A(u,b) be the expected number of nonzero BCH row pairs with union size
u and intersection size b, averaged over the two independent coordinate
shuffles. Let a=u-b. Independent lane swaps make every packet string of
that profile equally likely. There are

    Q(u,b) = binomial(256,u) binomial(u,b) 2^a

such strings. Therefore the counting measure of each string is A(u,b)/Q(u,b).
The exact positive formula in `profiles.py`, evaluated with authenticated
componentwise spectrum caps, supplies a cap A_cap(u,b).

Choose comparison probabilities 0<p,r<1. In each coordinate, independently,
put probability 1-p on packet 00, p(1-r)/2 on each of 01 and 10, and pr on
11. Its profile probability is

    P(u,b) = binomial(256,u) p^u (1-p)^(256-u)
             * binomial(u,b) r^b (1-r)^(u-b).

For any set I of allowed pair profiles, define

    C_I(p,r) = max_{(u,b) in I} A_cap(u,b)/P(u,b).

The outer counting measure restricted to I is pointwise at most C_I times
the iid comparison measure. Thus, for **any nonnegative function** f of the
pair's packet string,

    sum_x mu_I(x) f(x) <= C_I E_reference[f(X)].

This is a bound on all allowed profiles together. Charging the unconditional
reference expectation separately for every profile would replace the
maximum by a sum. That older bound is valid but unnecessarily loose.

The same argument applies to the union-support-only operator. If S_u bounds
the mass of shell u, then its density against Bernoulli(p) supports is at
most max_u S_u/Pr[Binomial(256,p)=u]. A CDF upper at u can cap S_u, but a
difference of CDF uppers cannot. `measure.py` implements this optional
`full_cover.py --density-fold` interface without changing the old one.

Independent pairs permit products of these pointwise inequalities, even
when each pair has a different allowed set and comparison probabilities.
The implementation below currently evaluates equal allowed sets and equal
probabilities. It does not yet assemble a complete mixed-class cover.

## Regional placement still has to be paid for

Fix q active pairs. Under the comparison measure, their packet types are
independent across pairs and regions. A uniform regional permutation places
the q pair labels uniformly among G=4096 slots.

For one region, independently mark each of G slots with probability theta.
A marked slot receives a comparison packet; an unmarked slot receives 00.
Conditioned on exactly q marks, this is exactly the required regional
comparison distribution. The conditioning probability is

    beta_q(theta) = binomial(G,q) theta^q (1-theta)^(G-q).

Unconditionally, all packets are iid with active probability v=theta*p and
double probability r conditional on activity. Positivity lets us upper-bound
the q-mark conditional transition kernel by the unconditional kernel divided
by beta_q(theta). Pay this factor separately in each of 256 regions.
The mark count is auxiliary, not the actual number of nonzero packets.

If an epoch envelope L bounds this iid kernel for every entering state,
the resulting all-pairs-in-I bound is

    exp(lambda*209715) binomial(G,q) C_I(p,r)^q
       * e_zero L^(64*256) tau / beta_q(theta)^256.

For a scalar envelope replace the matrix moment by L^(64*256).
The restriction theta=1 is allowed only at q=G.

For frozen lambda,theta,p,r,I and one frozen envelope, the logarithm is
convex in q: apart from affine terms it contains -255 log binomial(G,q).
One frozen witness at both endpoints therefore bounds an integer occupancy
interval; multiplying by its length pays the interval union. Optimizing the
endpoints independently does not establish this interpolation.

## Scalar and two-state envelopes

Put z=exp(-lambda). For an expanded state word, let (h0,h1,h2) count the
two-bit slots with zero, one, and two set bits. For iid packets define

    f0 = (1-v) + v(1-r)z + vr z^2,
    f1 = (1-v)z + v(1-r)(1+z^2)/2 + vr z,
    f2 = (1-v)z^2 + v(1-r)z + vr.

The exact output moment is f0^h0 f1^h1 f2^h2. Maximizing over all actual
expanded-state histograms, including the zero state, gives a scalar envelope.
There are only 12 histograms in this construction. In the special case
v=3/4,r=1/3, the input bits are uniform and the moment is
((1+z)/2)^128 for every entering state.

The scalar bound forgets persistence of a nonzero state. `iid_kernel.py`
retains two coordinates: zero-state mass and total nonzero-state mass.
Let H0=f0^64, and let H(lambda) be the maximum moment over nonzero expanded
states. Let J0=E[z^|X| 1{C(X)=0}], where C is the fixed feedback map.
The zero-state row is exactly (J0,H0-J0).

Both transvections precede feedback but follow the output computation.
For a nonzero state their combined distribution is lazy with probability
alpha=1/4 and a uniform nonzero refresh with probability beta=3/4.
If a_+ bounds every nonzero feedback atom Pr[C(X)=s], Cauchy--Schwarz gives

    E[z^|X+E(s)| 1{C(X)=s}] <= sqrt(H(2lambda)*a_+).

It is also at most H(lambda). Consequently a valid second row is

    (alpha*min(H(lambda),sqrt(H(2lambda)*a_+))
       + beta*H(lambda)/(2^19-1),
     H(lambda)).

A sharper bound handles feedback after the output tilt directly. For each
two-bit expanded-state mask y, normalize the four input probabilities
pi(x) z^|x+y|. For a local character c, let rho_c be the maximum absolute
Fourier coefficient of this normalized distribution over all four y.
Lane symmetry gives the same bound rho_1 for either one-coordinate
character; write rho_2 for the two-coordinate character. Then every
feedback atom under this tilted input law is at most

    A_tilt = 2^-19 sum_characters rho_1^n1 rho_2^n2.

Indeed, the normalized tilted input still factors over the 64 packets.
Fourier inversion and the triangle inequality bound each syndrome by
the displayed character sum, for every fixed expanded-state word.
Thus the lazy zero-return contribution is also at most H(lambda) A_tilt.
`iid_kernel.py` takes the minimum of all three bounds. This is an upper
bound on the actual tilted feedback distribution, not an independence
assumption about feedback and output weight.

The second entry deliberately does not subtract an upper bound on return
to zero. These two entries may overcount mass, but cannot undercount it.

J0 and a_+ are obtained from the feedback characters. For a character with
two-bit sign histogram (n0,n1,n2), its tilted Fourier value is

    Phi_z = f0^n0 ((1-v)-vr z^2)^n1
              * ((1-v)-v(1-r)z+vr z^2)^n2.

J0 is the average of Phi_z over all 2^19 characters. A Fourier L1 bound
gives a feedback atom upper. For nonzero targets one may first subtract
the all-zero input mass (1-v)^64 from every character: this removes only
an atom at feedback zero. The implementation takes the smaller of the
two valid L1 bounds. All numerical certificates use outward arithmetic;
binary64 optimizers merely propose rational witnesses.

## Results and scope

The scalar profile-conditioned pilot closed several selected dense events,
including q=2048,u=192 (more than 94864 bits) and q=4096,u=192 (more than
33789 bits). The two-state version also closed u=192 at q=256,512,1024,
2048,4096; at q=512 its margin exceeded 17841 bits. These were homogeneous
union supports, not full occupancy covers.

With pointwise measure domination and the two-state envelope, the complete
box I={160,...,192} for every active pair gives:

| q | Outward log2 upper |
|---:|---:|
| 256 | +1082.397541 |
| 512 | -4690.929002 |
| 1024 | +17118.443396 |
| 2048 | -1013.654658 |
| 4096 | -53578.624635 |

These boxes include all unequal supports in the interval, all permitted
intersections and messages, and all pair locations. Positive entries do not
certify anything. No interpolation between these independently optimized
rows has been claimed. The unrestricted box {38,...,256} remains positive
at every tested occupancy; this branch is complementary to the complete
sparse-occupancy covers, not a replacement for them.

The last table records the earlier Cauchy--Schwarz version. The driver now
also uses the tilted-feedback bound, so reruns can improve these numbers.
Reproduce the corresponding diagnostics with:

```sh
python -B research/workstreams/permutation_locality/two_bit/profile_dense.py --density-envelope --stateful --groups 256 512 1024 2048 4096 --interval 160 192 --interval 38 256 --outward
python -B -m unittest discover -s research/workstreams/permutation_locality/two_bit -p 'test_*.py'
```

Forty-nine tests currently pass. New exhaustive checks compare the iid
moment to every six-bit input/state word, the two-state envelope to the
complete small-state chain, and the profile and support density identities
to direct enumeration. The posterior Fourier atom bound is also checked
against every input mask and syndrome on a six-bit example. The next task
is a complete mixed-support cover,
with separate treatment of rare light and heavy profiles, plus the final
all-occupancy aggregate. Full-code closure remains open.
