# Which outer randomizers support the distance proof?

Independent nonzero GF65536 coefficients are sufficient for the exact
randomized-outer measure used by the proposed proof. They also permit one
coefficient per outer row to be fixed to one without changing the random code
ensemble. This is the smallest proven implementation saving identified here.

GF256 coefficients are a valid alternative for quiet-word analysis, because
quiet distance still descends to GF256. Their much weaker charge does not
transfer the full-field dense bound for general words. Smaller coefficient
sets admit a precise intermediate guarantee, stated below.

These are distribution and proof-transfer results. They do not close the
remaining source classes or certify 40% distance for the whole construction.

## Construction and probability space

Let B=GF256 and F=GF65536. The fixed generator of the outer RS[256,128] code
has entries in B and defines an F-linear code C. There are R=8192 independent
source rows and N=256R=2^21 output positions. Put p=256-128=128.

For source row r and outer position j, sample a nonzero coefficient d_(r,j).
Multiply that outer symbol by d_(r,j), then apply the 256 independent regional
row permutations. The intended full distribution samples every coefficient
independently from F excluding zero. Setup is sampled once and reused for all
messages.

The inner map is the shared beta64 construction:

    Y_i = F_64 X_i + A Q_i,
    Q_(i+1) = M_i Q_i + C_feedback X_i.

The packet width is 64, the state dimension is four, and F_64 has eight MDS
blocks of width eight. Initial state is zero and final state is free. The
fixed matrices have entries in B. The rank-one mixers have entries in F.
Outer coefficients, regional permutations, and mixers are independent.

Here C denotes the outer code only; C_feedback denotes the packet feedback
matrix. This distinction avoids using one symbol for both linear maps.

For one randomized outer row, define the nonzero expected measure

    mu(x) = E_D[ number of nonzero messages whose encoded row equals x ].

The encoder is injective for each setup, so mu(x)=Pr_D[x belongs to D C] for
x!=0, and mu(0)=0. Independence between different messages is not assumed.

## Exact full-field randomization

For fixed a!=0, d a is uniform on F excluding zero when d has that law.
Independent coefficients give the product law for a fixed outer word.
For a specified support S of size h, shortening C to S leaves dimension
h-p when h>=p+1. Its exact-support count A_S depends only on h, because C
is MDS. Therefore

    mu(x) = A_S / (q-1)^h,       q=65536,       S=support(x).

This is the same fixed-input law supplied by a uniform invertible map on the
16 underlying bits. Field multiplication preserves F-linearity as well.
The transitivity argument in [TRANSITIVE_SYMBOL_MAPS.md](../k16_design/TRANSITIVE_SYMBOL_MAPS.md)
therefore applies directly.

A distribution over deterministic invertible maps which sends every fixed
nonzero symbol uniformly to all q-1 outputs needs at least q-1 maps in its
support. Each map can have probability at most 1/(q-1). Thus the nonzero
scalar family has the minimum support size and minimum entropy for this
exact fixed-input requirement. This does not assert minimum instruction cost.

## A weaker theorem needs only coefficient entropy

Fix any nonzero scalar distribution nu on F, sampled independently at every
coordinate. Let rho=max_c nu(c). For every nonzero x,

    mu(x) <= rho^p.

If x has support size below p+1, the probability is zero. Otherwise, shorten
C to its support S and choose h-p information coordinates I. Expose the
coefficients in I. The values d_j^(-1) x_j at those coordinates uniquely
determine a candidate shortened codeword c. If a required coordinate is zero,
membership fails. Otherwise every remaining coordinate requires the unique
coefficient d_j=x_j/c_j. The p remaining coefficients are independent, and
each required value has probability at most rho. Averaging the exposure proves
the bound.

In particular, independent uniform coefficients from any set T contained in
F excluding zero, of size m, give

    mu(x) <= m^-p,
    mu <= beta_m Law(U),       beta_m = q^256 / m^128,

where U is uniform on F^256. T need not be a subgroup. Nonuniform or biased
sampling remains valid if its actual maximum atom rho is used.

This theorem replaces exact transitivity with a weaker atom bound. It does
not claim the full-field exact support measure. A continuation that needs
uniform nonzero labels cannot silently replace that premise by this theorem.
An envelope proof using only positive functions can use beta_m instead.

The coefficient law must remain independent of the route and inner setup.
Marginal uniformity without joint independence is insufficient: a common
uniform coefficient on every row coordinate has uniform nonzero marginals,
but leaves the entire row code unchanged.

## Remove one coefficient per row exactly

For uniform coefficients from a multiplicative subgroup H of F excluding
zero, write

    e_(r,0)=1,       e_(r,j)=d_(r,j)/d_(r,0) for j>0.

F-linearity gives diag(d_r) C=diag(e_r) C. The map from d_r to
(d_(r,0),e_(r,1),...,e_(r,255)) is a bijection on H^256. Consequently the
remaining 255 coefficients are still independent uniform H elements.

One may therefore set a predetermined coordinate coefficient to one in every
row and sample only the other 255. The distribution of the entire randomized
code image is unchanged, including its correlations between messages. This
preserves every distance event, not merely a first moment.

The saving is 8192 coefficient multiplications out of 2^21, or 1/256 of the
randomization layer. The fixed-message distribution changes; equality is at
the code-image level. An arbitrary set T lacks this normalization unless its
law has the required multiplicative invariance.

For this MDS code, common scalars are the entire diagonal stabilizer. Indeed,
each prescribed support of size p+1 has a one-dimensional shortening whose
nonzero words have that exact support. A diagonal map preserving C must
preserve this line, so its entries agree on the prescribed support. Every
coordinate pair lies in one such support, so all diagonal entries agree.
The full coefficient ensemble therefore has exactly (q-1)^255 distinct row
code images, with the uniform distribution. Sampling 255 independent
coefficients attains the minimum entropy needed to reproduce that exact
diagonal code-image ensemble. Further savings require a different ensemble
or a cheaper way to apply its coefficients.

A coefficient shared across a whole row, or across the entire outer block,
only changes the message basis. Such coefficients supply no code-image
randomness. Sharing a coefficient across rows at a given outer position does
change the code, but removes the independence needed to multiply row measures.
Identical source rows give a direct counterexample to that product law.

## GF256 coefficients and quiet descent

Suppose all coordinate coefficients lie in B. For each fixed coefficient
setup and route, the outer code and the quiet constraints are defined over B.
A word is quiet when every emitted state is zero. Its constraints are

    C_feedback X_i=0       for every packet except the final packet,
    Y_i=F_64 X_i.

No mixer coefficient occurs in these equations. Write an F-valued quiet word
as X=X_0+eta X_1 in a fixed B basis of F. Both components satisfy the base-field
outer and quiet constraints, and each component output support is contained
in the original output support. A nonzero component has nonzero output.
Embedding B-valued words proves equality of the two minimum symbol distances.
Existence with at most a active rows descends; an exact active-row count need
not be preserved.

Independent uniform B coefficients excluding zero consequently permit quiet
counting over B, with row charge 255^-128. The directed projective quiet-ball
bound already proved in [HD_BLOCKWISE_SHARED_QUIET_COMPOSE.md](HD_BLOCKWISE_SHARED_QUIET_COMPOSE.md)
is

    log2 L_quiet <= 7693823.1287572013310355600912018628354165293620.

For every fixed route, a fixed quiet output determines one inverse input and
one active-row set. Thus quiet witnesses with at least a active rows have
failure probability at most

    L_quiet * 255^(-128a).

There is no extra union over active-row sets. At a=7519 the directed exponent
is below -198.4381857430. Therefore the B-valued witnesses have a certified
dense tail. The support-only quiet bound through 256 active rows also survives
coordinate scaling. The intervening B-valued classes 257 through 7518 remain
open under these two bounds.

This statement does not classify every extension-field quiet word by its
original active-row count. A descended component can have fewer active rows.
It certifies the dense class of base-field witnesses used in the descent
partition, which must eventually be combined with the other base-field classes.

Full F coefficients destroy this descent premise. For example, multiplying
the repetition code {(a,a):a in F} by diag(1,eta), with eta outside B, gives
{(a,eta a):a in F}. Its only word in B^2 is zero. General F-valued mixer
equations cause the same issue for nonquiet histories even when outer
coefficients lie in B. A field equation x_0+eta x_1=0 has nonzero F solutions
and no nonzero B solution. Quiet descent cannot be used for those histories.

## How much entropy is lost?

For general F-valued words, replacing full coefficients by B coefficients
loses

    128 log2(65535/255) = 1024.7199422968 bits

of charge per active outer row. One must still count output labels in F.
The ordinary ball at relative weight delta has leading binary exponent

    N [ h_2(delta) + delta log2(q-1) ],

where h_2 is binary entropy. For all rows active, the m-element coefficient
envelope therefore has leading exponent per output symbol

    h_2(delta) + delta log2(q-1) - (1/2)log2(m).

This is a diagnostic of that envelope, not an impossibility theorem for the
construction. At q=65536 and delta=0.4:

| Independent coefficient set | Size m | Exponent per symbol |
| --- | ---: | ---: |
| GF256 excluding zero | 255 | +3.3737650704 |
| Norm-one subgroup | 257 | +3.3681295143 |
| Largest proper multiplicative subgroup | 21845 | +0.1634340462 |
| Any 32768-element subset | 32768 | -0.1290582111 |
| Full field excluding zero | 65535 | -0.6290472042 |

The threshold is log2(m)>14.7418835778, or approximately m>27399.91.
Since 65535=3*5*17*257, even the largest proper subgroup misses this ordinary
dense gate. The arbitrary half-sized set clears it, but still requires a
general F multiplication. It supplies no identified arithmetic saving and
has a weaker proof margin. Full coefficients remain the primary candidate.

For the quiet B-valued problem, the extra four equations per width-64 packet
make the ordinary rate-one-half ball comparison too coarse. The exact quiet
enumerator above gives a useful bound at 40%. Conversely, using only B-valued
mixers to make every word descend does not itself make a rate-one-half GF256
code pass the ordinary 40% entropy gate. Quiet constraints apply only to quiet
words.

## Subgroup products and inexpensive mixing

Multiplying several independently sampled coefficients from the same subgroup
never leaves that subgroup. Repeated GF256 coefficients cannot create a full
GF65536 orbit. A fixed extra scalar only translates the orbit to one coset.

There is an exact factorization of full-field sampling. The group of order
65535 is the product of its subgroups of orders 255 and 257, whose intersection
is trivial. An independent GF256 nonzero coefficient and an independent
norm-one coefficient therefore multiply to a uniform full-field coefficient.
Sampling one representative from each of the 257 GF256 multiplicative cosets,
independently of a nonzero GF256 coefficient, gives the same law.

This is a valid sampling decomposition. Applying the two factors separately
does not establish a faster encoder. Their product can instead be precomputed
as the ordinary full-field coefficient.

In the quadratic tower eta^2+eta+mu=0, multiplication by c=c_0+c_1 eta uses
three fixed GF256 products:

    p_0=c_0 x_0,
    p_1=(mu c_1) x_1,
    p_2=(c_0+c_1)(x_0+x_1),
    c x=(p_0+p_1) + (p_0+p_2) eta.

GF256 scaling uses two products, one per component. This explains its possible
arithmetic attraction, and the proof loss it must justify. Operation counts
are not native timings. All constants can be formed at setup.

A pair of independent componentwise GF256 scalings is not transitive on
F excluding zero: a vector with one zero component retains that zero
component. Adding swaps preserves this obstruction. More general byte-mixing
families would need an exact orbit law or a proved maximum-atom bound before
they could replace the coefficient family.

## Is a local packet shuffle necessary?

The uniform-input envelope in [rs_uniform_envelope.py](../k16_design/rs_uniform_envelope.py)
is pointwise. Any fixed permutation preserves it; an independent random
permutation does so after averaging. Therefore the present regional route
needs no additional local shuffle for that transfer.

At the present granularity, each routed atom is one F symbol. Full-field
coefficients make mu(x) depend only on the number of nonzero symbols, because
the outer code is MDS. Thus even the exact expected row measure is invariant
under arbitrary permutations of its 256 positions. An independent row-local
shuffle does not change this first-moment measure.

The older binary proof has a different interface. A randomized large symbol
is split into several smaller packets. Conditional on total small-packet
weight, its support need not be uniform across the whole outer group.
For example, with three two-bit symbols and exactly two active symbols, a
two-bit support occupying one complete symbol is impossible; two bits in
different symbols can occur. The independent local shuffle in
[RS8_TRANSFER.md](../k20_two_ms/RS8_TRANSFER.md) supplies the missing uniform
small-packet support law. Removing that shuffle invalidates that shell
formula, though its uniform envelope remains valid.

GF256-only or arbitrary-subset coefficients do not give the full-field
exchangeability of labels used in the exact measure. Their weaker atom
envelope still survives independent routing.

## Larger alphabets and higher target distance

For full-field coefficients and outer rate k/n, the ordinary dense gate is

    h_2(delta) < (1-k/n-delta) log2(q-1).

For rate one half, the approximate bit-size thresholds on the alphabet are:

| Target relative distance | Required log2(q-1) |
| --- | ---: |
| 40% | greater than 9.70951 |
| 42.5% | greater than 13.11611 |
| 45% | greater than 19.85549 |
| 47.5% | greater than 39.92784 |
| 49% | greater than 99.97115 |

These are leading-exponent gates, not full finite-length certificates.
Keeping GF256 fixed coefficients requires an extension degree over GF256;
the natural field bit sizes are multiples of eight. Thus GF65536 clears the
40% ordinary dense gate, while a 45% target would first clear it at GF(2^24).
The outer and fixed FFT8 matrices remain valid under these field extensions.
The cost and the nonquiet proof for wider field arithmetic remain separate.

Increasing outer redundancy improves the gate by increasing 1-k/n, at the
cost of message rate. Increasing n at the same rate does not change this
leading dense inequality. It can affect sparse classes, routing costs, and
finite terms, which must be analyzed separately.

## What the new randomizers do not remove

The rare-route dimension witness in [HD_BLOCKWISE_SHARED_CUTS.md](HD_BLOCKWISE_SHARED_CUTS.md)
survives every nonzero diagonal coefficient setup. On its route event, all
coordinates of 640 chosen rows remain confined to the selected packet windows.
Diagonal scaling changes values but not confinement or source dimension.
Imposing quiet feedback still leaves at least 29496 field dimensions.

For B coefficients, the existing lower bound above 2^4110 on the raw expected
number of low-weight base-field projective quiet lines therefore survives.
For full F coefficients, the same dimension calculation over F is even larger.
This is a limitation of the uncapped first moment, not a lower bound near one
on distance failure. A proof must still cap or separately bound the rare-route
events before counting all words in their large kernels.

## Validation and next step

[theory/hd_randomized_outer_options.py](theory/hd_randomized_outer_options.py)
implements the scalar-law audit and the directed GF256 quiet tail. Its six
[tests](theory/hd_randomized_outer_options_tests.py) check arbitrary subsets,
nonuniform laws, entire-code normalization, the failure of a subfield law to
give the full-field envelope, common row scalars, and the descent limitation.
They also certify the dense quiet endpoint with Arb and reject the preceding
row threshold. All six pass.

    C:/Python314/python.exe -B -m unittest hd_randomized_outer_options_tests -v
    C:/Python314/python.exe -B hd_randomized_outer_options.py

The immediate proof choice is full-field coefficients with one normalized
coordinate per row. Continue the middle-density proof using the exact scalar
law. Keep GF256 coefficients as a distinct quiet-only comparison until a
nonquiet continuation supplies the missing entropy. No native implementation,
benchmark, remote run, or existing proof file was changed for this audit.
