# Quiet fractional moments conditioned on tile occupancies

Conditioning on the eight tile occupancies removes a substantial loss from
the current quiet bound. At two fixed middle-range witnesses, the improvement
is about 26,584 and 14,417 bits before exact regional coefficient extraction.
Both resulting bounds remain positive. These calculations diagnose a proof
loss; they do not certify either class or exhibit a low-distance word.

The construction and its randomness are unchanged. The calculation uses the
same residual-rank lower bound as the
[joint quiet moment](HD_RANDOMIZED_OUTER_JOINT_QUIET.md), but averages fewer
route variables before taking the fractional power.

## Keep the tile profile until after the fractional power

Use q=65536, 8192 outer rows, 256 independent regional permutations, and
32768 packets of 64 positions. Each packet has eight tiles of width eight.
The shared beta64 packet maps are F and C. For a fixed potential subset
S of one packet, sample X_S from the uniform distribution on F_q^S and
set the other input coordinates to zero. These artificial uniform amplitudes
belong to the outer comparison, not to a change in the encoder.

For 0<z<=1, define

    K_S(z) = E_X[1[CX=0] z^wt(FX)],
    T_S(z) = E_X[z^wt(FX)].

The final packet uses T_S because its outgoing state is unconstrained.
Let a=(a_1,...,a_8) be the tile profile, where a_l counts potential inputs
in tile l. Put j=sum_l a_l.

Fix output-zero counts v=(v_1,...,v_8), with each v_l in {0,...,8}, and put
d_l=max(a_l-v_l,0). The previous tile-rank argument gives output rank
sum_l min(a_l,v_l). Its residual feedback image has rank at least rho,
computed by starting at zero and applying

    rho <- 2                if rho=min(4,d_l)=1,
    rho <- max(rho,min(4,d_l)) otherwise.

This bound holds for every choice of potential and output-zero masks with
these counts. It uses distinct residual lines from distinct tiles. No new
claim about higher-dimensional intersections is needed here.

Define the profile bound

    J_a(z) = sum_(v in {0,...,8}^8)
               product_l [binom(8,v_l) z^(8-v_l)(1-z)^v_l
                          q^(-min(a_l,v_l))] * q^(-rho),
    Khat_a(z) = min(J_a(z), q^(-min(j,4))).

The second bound uses the full-rank minors of C. For every S with profile a,

    K_S(z) <= Khat_a(z).

Omitting the final factor q^(-rho) from J_a gives the exact free moment
T_a(z)=T_S(z). The fixed-profile formula has no binom(8,a_l) factor: its
potential occupancies have already been fixed.

For a uniform j-subset S, sort a into nondecreasing order. A sorted profile
has multiplicity

    m(a) = 8! / product_(h=0)^8 n_h! * product_l binom(8,a_l),

where n_h is the number of tiles with occupancy h. Thus
sum_(a:sum a=j) m(a)=binom(64,j). There are binom(16,8)=12870 sorted profiles.
Both J_a and the rank lower bound are invariant under permutation of the
tile counts, although the actual feedback subspaces need not be identical.

For a fractional exponent 0<alpha<=1, define

    Psi_j = sum_(a:sum a=j) m(a) Khat_a(z)^alpha / binom(64,j),
    Phi_j = sum_(a:sum a=j) m(a) T_a(z)^alpha / binom(64,j).

Conditioning the route on a before applying conditional Markov is valid.
The pointwise bound K_S<=Khat_a then gives

    E_S[K_S(z)^alpha | |S|=j] <= Psi_j.

The old calculation instead used a bound on E_S K_S and raised that bound
to alpha. Concavity gives

    Psi_j <= (E_a Khat_a(z))^alpha
           <= min(E_a J_a(z),q^(-min(j,4)))^alpha.

The first inequality measures the loss from averaging profiles too early.
The second also allows a syndrome cap to be applied separately to each
profile. A control calculation below isolates these two changes.

## Compose the already-fractional moments

Fix r active outer rows and let D=838860 be the output-weight threshold.
Put beta=q^256/(q-1)^128. Conditional on packet counts, the within-packet
subsets are independent and uniform. This remains true under the fixed
regional total of r potential positions.

The regional polynomials now have coefficients binom(64,j) Psi_j and
binom(64,j) Phi_j. These coefficients already contain the fractional power;
raising them to alpha again would be incorrect.

For a marker 0<p<1, put

    S_p = sum_j binom(64,j) p^j (1-p)^(64-j) Psi_j,
    T_p = sum_j binom(64,j) p^j (1-p)^(64-j) Phi_j,
    c_r = binom(8192,r) p^r (1-p)^(8192-r).

The complete r-row quiet failure class satisfies the upper bound

    B_r = binom(8192,r) [beta^r z^(-D)/(q-1)]^alpha
            * S_p^32767 T_p / c_r^256.

The outer binomial factor remains outside alpha. The projective divisor
q-1 remains inside alpha. The regional normalization is still c_r^256;
the marker does not replace the actual routing law.

Exact regional coefficient extraction uses

    A(x) = sum_j binom(64,j) Psi_j x^j,
    E(x) = sum_j binom(64,j) Phi_j x^j,
    G_r = [x^r] A(x)^128 / binom(8192,r),
    H_r = [x^r] A(x)^127 E(x) / binom(8192,r).

Replacing S_p^32767 T_p/c_r^256 by G_r^255 H_r gives a stronger valid
bound. The helper reuses the existing positive coefficient engine with
its exponent argument set to one.

## Exact local arithmetic and two directed checks

For rational z=n/d, multiply every local tile coefficient by d^8 and q^a.
The resulting coefficients are nonnegative integers. A five-component rank
vector is propagated along shared prefixes of the sorted profiles. At a
leaf, division by d^64 q^j restores the moment. Only 24,309 nonempty sorted
prefixes are needed. All local profile moments are exact rational numbers.

Arb is used only for fractional powers, coefficient products, and logarithms.
The following results were replayed at 256-bit precision, using alpha=1/32.
The marker is optimized separately for each displayed envelope.

| r | z | Old marker bound | Tile-conditioned marker bound | Tile-conditioned exact regional bound |
| ---: | ---: | ---: | ---: | ---: |
| 1536 | 1/4 | +35168.87729 | +8584.890630253987 | +6898.557283911752 |
| 4096 | 1/4096 | +76784.49991 | +62367.062195666380 | +60619.506541191010 |

All entries are base-two logarithms of upper expressions. Positive values
are uninformative probability bounds. None is an acceptance threshold or
a construction counterexample.

The control first averages the profilewise capped Khat_a and then applies
alpha. Its floating values match the old marker bounds to numerical precision.
Thus essentially all the displayed improvement comes from postponing the
profile average, not from moving the syndrome cap. Exact regional extraction
then adds about 1,686 and 1,748 bits at these two witnesses.

The rational markers used in the directed tile-conditioned replay are
197598598381/941270402265 for r=1536 and 541623685321/924641301185 for r=4096.
The respective marker enclosures are

    8584.890630253987170887448124591247163731102096313712526027218873198785962
        +/- 5.21e-70,
    62367.062195666380472990873305986295304608006037362259524590772987972190191
        +/- 5.47e-70.

A final bounded re-optimization checks that the first diagnosis is not
solely an artifact of using a witness selected for the old bound. It keeps
r=1536 and v=1, and evaluates exactly the 49 pairs

    z in {1,1/2,1/4,1/8,1/16,1/64,1/256},
    alpha in {1/64,1/32,1/16,1/8,1/4,1/2,1}.

Each pair uses its own optimized marker. The best floating expression is
+5438.724016621243 at z=1/8 and alpha=1/32, with base-two marker log-odds
-1.7944536117517722. The next two expressions are +6063.14914 and +6568.41799,
also at z=1/8, with alpha=1/64 and 1/16. No expression falls below +2000,
the preselected trigger for another directed exact-regional replay.
Consequently this check triggers no further replay or occupancy scan.
It improves the witness but does not close the class; it does not establish
global optimality over z and alpha.

## Combine with the actual input-support floor

Every nonzero RS[256,128] outer row has at least 129 nonzero symbols.
Consequently a word with r active rows has input weight at least 129r.
For an input marker v>=1, the artificial uniform comparison can retain
this constraint by multiplying the complete bound by v^(-129r alpha).

The zero input contributes exactly q^-j to both local moments. Every other
input has weight at most j. Therefore, before taking alpha, valid marked
profile moments are

    Khat_a(v,z) = q^-j + v^j (Khat_a(z)-q^-j),
    That_a(v,z) = q^-j + v^j (T_a(z)-q^-j).

The nonzero remainder is nonnegative. The implementation subtracts the
zero atom in exact rational arithmetic, not by subtracting approximate
upper bounds. The terminal profile must receive the same input marker.

At the existing witnesses, one bounded check with v=257/256 gives directed
optimized-marker logarithms +8616.685302774310 for r=1536 and
+62274.439421852213 for r=4096. The first worsens by 31.79 bits; the second
improves by 92.62 bits. Neither closes its class. These are local diagnostic
checks, not an optimization over v or a claim about every witness.

## Interface, tests, and next step

The new helper is
[theory/hd_randomized_outer_quiet_tile.py](theory/hd_randomized_outer_quiet_tile.py).
Its interfaces are:

- `profile_moments(z)` returns all sorted profiles, packet counts, exact
  multiplicities, exact rational Khat_a, and exact rational T_a.
- `fractional_profile_logs(z,alpha,input_mark=1)` returns two 65-element
  arrays of floating base-two logarithms of already-fractional moments.
- `fractional_profile_enclosures(z,alpha,input_mark=1)` returns the same
  quantities as raw Arb enclosures, without logarithms.
- `directed_endpoint(r,z,alpha,input_mark=1,exact_regional=False)` replays
  the optimized-marker bound and optionally the exact regional bound.

Five new tests check profile normalization, the known zero atom, a direct
two-tile output-mask sum, the exact z=1 syndrome law, Jensen ordering,
agreement with the previous mean moment, and directed arithmetic. The six
sparse-ladder tests also passed during the independent review of that proof.

    python -B -m unittest hd_randomized_outer_quiet_tile_tests -v
    python -B hd_randomized_outer_quiet_tile.py

The next useful step is to retain actual nonzero-input weight information
jointly with output weight, then test a bounded set of witnesses. The input
marker should be included only where its derivative or a short check predicts
a gain. Enlarging the feedback-rank state alone does not recover the
profile-averaging loss measured here.

No previous proof helper, encoder implementation, benchmark, remote host,
or commit was changed for this calculation.
