# XOR2 multigroup bound with nonzero amplitudes

The unchanged scalar route now has a directed certificate for one substantial
XOR2 message class. It covers every message with 1 through 8192 active groups,
where each active group has exactly one active RS constituent. The summed
setup-failure bound has logarithm below

    -97.40067454.

The certificate includes every group choice and every constituent-row choice.
It retains the exact RS support enumerator. It does not certify all XOR2
messages: groups with multiple active constituents, and mixtures of those
groups with single-constituent groups, require additional coverage.

## Construction and probability space

The candidate is the normalized two-diagonal column map from
[HD_XOR_COLUMN.md](HD_XOR_COLUMN.md). Over F=GF256, let H=I+11^T on eight
coordinates. Each group contains eight independent RS[32,16] constituents.
Each of its 32 columns uses D2 H D1, with independent nonzero input and
output diagonals. Normalization preserves the unnormalized matrix law.

There are 8192 groups and 256 regions. A group contributes one coordinate
to each region. Each region independently permutes its 8192 coordinates.
The unchanged scalar inner processes 128 packets per region. Each packet
has 64 inputs, 16 state coordinates, and four emission blocks. Each block
uses the established RS[16,4] emission map and its transpose as feedback.
The final length is 2^21 and the forbidden weight threshold is 734003.

All probabilities below concern the setup diagonals, regional permutations,
and independent inner setup coins. Fix a message before sampling this setup.
Condition on any resulting outer support pattern with positive probability.
The D2 diagonals then make all occupied outer amplitudes independent uniform
elements of F*. The regional counts may remain correlated.
This conditional amplitude statement requires D2. It does not apply directly
to H D1 alone.

## Packet operator

The earlier one-group calculation only encountered empty packets and packets
with one occupied input. Multiple groups require every feasible occupancy.
In particular, five nonzero inputs in one feedback block can cancel.

Let q=256, d=q^16-1, b=d/(q-1)^16, and kappa=q/(q-1). The existing two-state
envelope tracks zero state and a density cap relative to the uniform nonzero
state. The refresh charges b on transitions into its nonzero-state component.
Fix a weight mark 0<z<1.

For a selected support U of size l within one emission block, let X have
independent F* entries on U and zero entries elsewhere. Define

    e_l = Pr[C_U X_U=0].

Any four columns of the block feedback matrix C are independent. Therefore

    e_0=1,
    e_l=0                                      for 1<=l<=4,
    e_l=B(l,l-4;l)/(q-1)^l                     for 5<=l<=16.

Here B(n,k;w) is the number of MDS[n,k] words on one prescribed support of
size w. For w>=n-k+1 it equals

    (q-1) sum_(j=0)^(w-(n-k+1)) (-1)^j binom(w-1,j) q^(w-(n-k+1)-j).

The kernel of C_U is an MDS[l,l-4] code, which proves the last equality.
For example, e_5=1/(q-1)^4. Feedback cancellation requires every occupied
emission block to cancel; one block with occupancy 1 through 4 prevents it.

For packet occupancies l=(l_1,l_2,l_3,l_4), put j=sum_i l_i and
e(l)=product_i e_(l_i). Starting from zero, the emitted vector is X itself.
Thus the operator's first row is

    P_l[0,0] = z^j e(l),
    P_l[0,1] = b z^j (1-e(l)).

This row accounts for all occupied inputs and their possible feedback
cancellation. It never replaces an occupied packet by a single input.

The second row uses the established flat-state bounds. A return to zero
requires the incoming state to equal C X. At that state the emitted vector
is (I+A C)X. For one block, define its marked nonquiet-return moment

    J_l = E[z^wt((I+A C)X) 1_(C X != 0)].

Let barJ_l be the existing bound with uniform F amplitudes, including zero.
Conditioning all l amplitudes to be nonzero gives J_l<=kappa^l barJ_l.
The additional support bound is

    J_l <= (1-e_l) z^max(1,13-l).

To prove it, C X!=0 gives an RS[16,4] word A C X of weight at least 13.
Adding X changes only l positions. Moreover, I+A C is invertible because
C A=0, so the emitted vector remains nonzero. The implementation takes the
minimum of this bound, the conditioned free/quiet bound, and the existing
principal-nullity return bound.

Put K_l=z^l e_l. For four blocks, a valid packet return cap is

    J_packet = product_i (K_(l_i)+J_(l_i)) - product_i K_(l_i).

The code evaluates this expression through a positive recurrence. The flat
return entry is J_packet/d. The flat live entry is at most b kappa^j barV_l,
where barV_l is the established uniform-amplitude live cap. Empty packets
also use the existing inverse-refresh cap. These are bounds on the same
scalar inner used by the previous XOR column screen.

## Regional coefficient and a common potential

For every regional nonzero count m, average the packet operators over all
subsets of m positions in that region. Write the resulting matrix as R_m.
If F(t) sums the packet operators over all 64-position supports, marked by
their cardinality, then

    R_m = [t^m] F(t)^128 / binom(8192,m).

The matrix product follows packet order. Profile multiplicities include
every labelled support, and the coefficient retains all packet collisions.
The calculation does not impose a maximum occupancy below m.

Fix a maximum active-group count U. Suppose positive numbers v and u satisfy

    R_0 (1,v)^T <= (1,v)^T,
    R_m (1,v)^T <= u^m (1,v)^T              for every 1<=m<=U.

All inequalities are componentwise. They can be checked using two scalar
inequalities per occupancy. Define c_v=max(1,1/v).

Now fix any complete outer support pattern with regional counts
m_1,...,m_256<=U. The regional permutations remain independent after this
conditioning. Fresh inner coins and the state-envelope property give

    E[z^W | outer support]
       <= e_0^T R_(m_1) ... R_(m_256) 1
       <= c_v u^(sum_i m_i).

The second inequality follows by applying the displayed potential inequality
successively from the right. The initial state is zero, so its potential
is one. The endpoint comparison costs c_v once over the complete codeword.

This bound holds for each support pattern separately. Averaging the support
pattern therefore preserves it even when all regional counts are correlated.
The scalar potential loses some information about where support occurs;
it does not assume independent support masks.

## Exact union for one active constituent per group

Fix one active constituent in an active group. Every nonzero constituent
coordinate produces exactly seven nonzero outputs under H. If its RS weight
is w, its group output weight is exactly 7w.

Let B_w=B(32,16;w). Define the exact marked count

    S(u) = sum_(w=17)^32 binom(32,w) B_w u^(7w).

Each active group has eight possible active constituents. Consequently, for
exactly r active groups, with 1<=r<=U, the preceding potential gives

    Pr[there exists such a nonzero word of output weight <=734003]
       <= z^(-734003) c_v/(q-1) * binom(8192,r) * (8 S(u))^r.       (1)

The proof applies Markov's inequality to the number of projective message
classes with forbidden weight. Each encoder is F-linear, so every nonzero
scalar orbit contains q-1 messages of the same weight. Divide by q-1 once
for the complete message. Summing over all messages before averaging setup
includes messages whose choice depends on setup.

For an integer interval [L,U], sum the right side of (1) over all r in that
interval. The same u and v cover every occupancy from 0 through U. The
support floor 17 remains inside S(u); it is never replaced by a mean weight.

## Directed receipts and numerical screens

The helper saves 84 disjoint rational witnesses covering every integer
r=1,...,8192. It fixes rational z and defines v=2^a and u=2^s for saved rational
exponents a and s. Arb verifies all regional potential inequalities and
evaluates each interval union at 512-bit precision. No floating optimizer
is used during replay. Summing the interval receipts gives

    log2 failure < -97.4006745419751485245.

The first 29 witnesses cover r=1,...,128. The other 55 cover r=129,...,8192;
their largest interval bound is -163.16741650 on [176,191]. The widest
interval, [6895,8192], has bound -389.50883161. The largest individual
bound over all intervals is -97.49115357 on [43,48]. This is complete
coverage of the stated message class, not an interpolation between samples.

For the larger counts, interpolation between successful FLOAT marks only
selected candidate values of z. Each saved witness was then independently
replayed for every regional occupancy through its interval's upper endpoint.
The final replay reads only the literal rational witnesses. It performs
neither interpolation nor optimization.

For comparison, selected first-moment FLOAT screens are:

| Active groups r | One active constituent in every group | Two active constituents in every group |
| ---: | ---: | ---: |
| 2 | -889.2 | +224.6 |
| 4 | -1027.6 | +451.0 |
| 16 | -1560.0 | +1733.8 |
| 128 | -3547.4 | +9655.4 |
| 256 | -4829.0 | not screened |
| 512 | -6976.6 | not screened |
| 1024 | -11836.8 | not screened |
| 2048 | -29884.4 | not screened |
| 4096 | -119089.1 | not screened |
| 8192 | -763014.9 | not screened |

The two-constituent screen uses the exact two-row outer weight enumerator
from HD_XOR_COLUMN.md and replaces 8 S(u) in (1) by its 28-row-choice moment.
It retains the two-row equal-value cancellation branch. Its positive results
expose a limitation of this potential comparison. The same comparison gives
about +109.4 for one group, whereas the earlier chronological operator gives
about -300.6. Thus these positive results do not justify a construction change.

Two additional directed intervals provide independent numerical checks
outside the sparse endpoint:

| Group interval | z | log2(v) | log2(u) | Directed log2 interval union |
| --- | --- | --- | --- | ---: |
| [972,1024] | 833836287/1000000000 | -2171/32 | -41/20 | -1370.5943 |
| [1924,2048] | 167141457/250000000 | -1173/16 | -11181/5000 | -2033.4970 |

These supplementary checks overlap the complete coverage. The 84 disjoint
intervals determine the reported family union; these two checks are not
added to that union.

## Validation and next step

The implementation is [theory/hd_xor_multigroup.py](theory/hd_xor_multigroup.py).
Its [twelve tests](theory/hd_xor_multigroup_tests.py) pass. They exhaust small
field kernel laws, verify the first-row mass and five-input cancellation,
check noncommuting polynomial composition, compare the one-input operator
against the earlier screen, compare the fast proposal arithmetic with
logarithmic convolution, and replay the 29 sparse witnesses. The separate
full replay checks all 84 directed witnesses and sums their bounds. An
invalid potential is deliberately rejected. No native benchmark, remote
operation, or construction change was performed.

Run from the theory directory:

```text
python -B -m unittest hd_xor_multigroup_tests -v
python -B hd_xor_multigroup.py --replay-one-row-all
python -B hd_xor_multigroup.py --replay-large-intervals
python -B hd_xor_multigroup.py --active 1 2 --groups 2 4 16 128 --iterations 40
python -B hd_xor_multigroup.py --active 1 --groups 256 512 1024 2048 --iterations 35 --intervals
```

The next proof step should retain the eight-region column chronology for
groups with two active constituents. The nonzero packet operators developed
here can be reused. The homogeneous one-constituent family is now covered
for all active-group counts. A whole-code union must also cover different
active-constituent counts within the same message.
