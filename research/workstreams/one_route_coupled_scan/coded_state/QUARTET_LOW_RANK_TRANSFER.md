# A small transfer for the global footprint count

The forced-footprint upper kernel admits a positive130-state factorization.
It now composes every update and whole regions, rather than repeating the
largest one-step row. This is an exact reduction of the stated relaxed
count, not a claim that the relaxed count proves the distance target.

Use a=q^theta, output marker z=x^theta, m=((q+1)x)^theta, and
b=(1+(q+1)x)^theta as in
[the grouped-footprint bound](QUARTET_GROUPED_FOOTPRINT.md).
Use chi for a separate coverage marker. A physical boundary state records
f forced and o optional covered positions. Its uncovered positions are
implicit. Optional count o is a multiple of four.

## Positive factorization and terminal condition

First consider the relaxed kernel (5) in that note. It factors through
k=0,...,128, the number of full outgoing quartets:

    U_((f,o),k)=m^f [t^(4k-f)](1+bt)^o(1+t)^(512-f-o)
                       / binom(512,4k),

    V_(k,(f',o'))=binom(128,k) chi^(f'+o')
                  [y^f'](1+4y+6y^2+4y^3)^(128-k)
                  1[o'=4k].

U is zero if 4k<f. To make every activation after an empty footprint
entirely forced, zero the empty row of U and append a column that is one
only on that row. Append to V the exact empty-input row

    V_(*,(f',o'))=binom(512,f') chi^f' 1[o'=0].

Call the resulting factors A and B_chi. The physical kernel is A B_chi;
the reduced matrix C(chi)=B_chi A has dimension130 and nonnegative entries.
This kernel is stronger than (5), but weaker than full (7): a newly full
quartet inside an otherwise nonempty footprint is still labeled optional.

The initial reduced row is e_*. The physical terminal function m^f b^o
equals A(e_128+e_*). Thus the entire count, including the free terminal,
is exactly

    e_*^T C(chi)^4096 (e_128+e_*).                       (1)

No population classes, destinations, or terminal states are truncated.

## Direct construction from univariate polynomials

Let s=m chi and

    Q_s(t)=(1+(1+s)t)^4-s^4t^4.

For i,j in0,...,128,

    C_ij=binom(128,i) chi^(4i)/binom(512,4j)
           [t^(4j)](1+bt)^(4i)Q_s(t)^(128-i)-1[i=0].

Also C_(i,*)=1[i=0], C_(*,j)=(1+s)^(4j)-1, and C_(*,*)=1.
The implementation removes the i=0 empty contribution through a positive
polynomial recurrence, rather than subtracting nearly equal large values.
It uses only degree512 univariate polynomials and the130-square matrix.

At chi=0 the matrix is

    Z=(e_0+e_*)e_*^T,    Z^2=Z.

The positive regional operator is C(chi)^16-Z. Repeated squaring uses
Q_(2n)=Q_n^2+ZQ_n+Q_nZ, so no empty-path cancellation is needed.
The rank-one renewal identity counts paths occupying at least129
coordinate-regions. Its scalar sequences are

    alpha_k=e_*^T Q^k(e_0+e_*),
    beta_k=e_*^T Q^k(e_128+e_*).

This is the same renewal mechanism as before, with the correct reduced
initial, zero, and terminal vectors.

## Combine a fixed eligible-row set with the count

Fix a candidate row set H of size h before routing, and later test whether
it is exactly the rows whose covered degree is at least129. On compatible
histories in this event,

    R(V)=128h+|V outside H|.

Every guaranteed birth must belong to H. With independent Bernoulli
H-labels of probability p, a guaranteed covered birth contributes p.
A full optional quartet instead contributes

    r^4,    r=p+(1-p)/a,

from its covered outside-H rank. Therefore the same factorization works:
replace the partial/empty-row coverage marker chi by p chi, and the full
optional marker by r chi. A is unchanged. Sum over binom(8192,h) fixed
candidates and retain the global factor a^(-128h).

Condition the8192 labels in each region on exactly h H-slots. If
pi=Pr[Bin(8192,p)=h], the regional upper operator is

    Z+pi^(-1)(C_p^16-Z).                                (2)

Zero coverage contributes Z under every label realization; only the
positive part pays pi^(-1). All choices of p bound one common conditional
operator in the same reduced basis, so entrywise minima over p are valid.
They do not resample a code along different paths.

Aggregate eligibility can be retained with X>=1, Y<=1 by using

    forced marker pX,
    optional marker pX+(1-p)Y/a,
    global factor X^(-129h)Y^(-128(8192-h)).

This still drops the individual row-degree conditions. All candidate-set
conditioning is applied before imposing those conditions; treating H as
already conditioned to be the actual set would not give uniform labels.

## Numerical evidence and controls

`quartet_lowrank_footprint.py --checks` verifies261 factorization entries
against explicit small quartet-partition enumeration,8 composed paths,
25 eligible-rank substitutions, and6 rank-one renewal comparisons. It
also checks the positive sixteen-update composition against direct powers.

Full-width calculations use long double and logarithmic scaling. They
are not directed-rounding certificates. Selected whole-history bounds are:

| Candidate eligible rows | theta | output log2 z | Tested log2 upper |
|---:|---:|---:|---:|
| 1 | 1/8 | -1/2048 | +494 |
| 8 | 1/8 | -1/256 | +586 |
| 64 | 1/8 | -1/32 | -1313 |
| 512 | 1/8 | -1/4 | +50009 |
| 2048 | 1/16 | -1/2 | +149251 |
| 4096 | 1/16 | -1 | +458926 |

The last three use17 conditioning tilts; the first three use3. The
negative64-row proposal would need a verified replay. It is not full
row-count coverage. Fine quota tilts reduce the gap but do not remove it;
small aggregate-degree perturbations also leave the hard classes positive.

There is a mandatory structural control. The preserved invariant-plateau
bad setup has228 coordinates and3648 outer rows, with probability

    [binom(128,57)/binom(512,228)]^4095
       *binom(8192,3648)^(-256).

The relaxed unsplit theta=1/4 moment includes a contribution exceeding
2^33152 for every output marker. This is a proof-method obstruction, not
a frequent bad code. Theta=1/8 can suppress this family. The existing
high-coverage theorem using theta=1/4 is unaffected: its class excludes
this footprint. One common moment parameter must not be imposed on all
classes without checking this control.

All runs were serial, file-free on Peach, single-threaded, below the
180-second watchdog, and within the2GiB address-space bound. No benchmark
or construction change was made.

## Subsequent restriction and propagation audit

The broader [covered-footprint theorems](QUARTET_GLOBAL_COVERAGE_RUNS.md)
now remove full/empty footprints and an expansion class. See the
[rank-gap checkpoint](QUARTET_GLOBAL_RANK_GAP.md) for their direct insertion,
a new R-L>=65537 theorem below2^-63, and the forced-activation audit.
That checkpoint derives a separate129-state global control which rejects
the old1->512->1 overcount. Tested middle-class bounds remain positive.
The full45%/40-bit target remains open.
