# Independent audit of the wider-outer actual24 proposal

No mathematical gap was found in the declared construction or its proof
wiring. The completed `wider_actual24_allq_v1.json` proposal covers every
occupancy for a legitimate 256-to-512-bit outer group followed by the actual
24-bit byte-native inner. This is a different construction from the measured
small-outer encoder; it is not a metadata reinterpretation of that encoder.

The audit reviewed `wider_coverage.py`, `cached_wider24.py`,
`wider_fractional.py`, the frozen wider geometry, the exact outer shell
generator, and the pointwise uniform envelope. No reviewed source was changed.

## Concrete outer and geometry

One group consists of eight parallel GF16 RS[16,8] rows. At each evaluation
position, concatenate the eight four-bit row symbols into one 32-bit symbol.
Choose a basis of GF(2^32) over its GF16 subfield. In that basis, the same
base-field generator matrix defines an MDS[16,8] code over Q=2^32. Its minors
remain nonzero in the extension field. Therefore the aligned 32-bit symbol
support has the required MDS law.

At every symbol position, sample an independent invertible binary map whose
image of each fixed nonzero input is uniform over the Q-1 nonzero symbols.
Use independent maps between groups. These maps are sampled once at setup
and then reused for every message. Neither the shell calculation nor the
pointwise envelope requires independence between messages.

Split each randomized symbol into four eight-bit packets. Independently
permute all 64 packets within each group, assigning one packet to each region.
Then independently permute the 256 group slots within each of the 64 regions.

The geometry used by the sources is consequently:

| Object | Size |
|---|---:|
| Message / output | 65,536 / 131,072 bits |
| Outer groups | 256 |
| Group dimension / output | 256 / 512 bits |
| Byte regions | 64 |
| Slots per region | 256 |
| Physical steps per region | 32, each with eight byte slots |
| Total physical steps | 2,048 |
| Persistent inner state | 24 bits |

The state starts at zero, continues through all regions, and has no final
flush. The low-weight cutoff is floor(131072/10)=13107.

## Pointwise envelope and exact shell law

For any fixed MDS-symbol support H of size h>=9, shortening has dimension
h-8 over GF(Q). Projecting onto h-8 information coordinates bounds the number
of words with exact support H by (Q-1)^(h-8). Independent symbol maps assign
any fixed nonzero labeling probability (Q-1)^(-h). Thus the expected
nonzero-message measure satisfies

```
mu(x) <= (Q-1)^(-8),
mu <= beta * Uniform({0,1}^512),
beta = 2^512 / (2^32-1)^8.
```

For supports below nine symbols, mu is zero. The zero message is excluded,
so mu(0)=0; the uniform majorant may nevertheless include artificial zero
mass. Independent group setups permit products of the expected measures.
An independent fixed or random subsequent permutation preserves domination.

For q=1, the implementation uses the sharper exact expected shell counts.
A uniformly randomized nonzero 32-bit symbol has byte-support polynomial

```
P(t) = ((1+255t)^4-1)/(2^32-1).
```

If A_h is the MDS symbol-weight enumerator, the group polynomial is
sum_{h=9}^{16} A_h P(t)^h. The exact coefficients sum to 2^256-1, and their
zero coefficient is zero. Conditional on a byte-support pattern, active
labels are independent uniform nonzero bytes. The independent group packet
permutation makes the support uniform among its 64-region subsets.

## Fractional conditioning for q>=2

Fix a set S of q nonzero outer groups. Let J record its potential occupancy
in every physical step, regardless of the sampled byte values. Conditioning
on J concerns regional routing only; outer maps and inner updates remain
independent of that routing. Within each physical step, selected positions
are uniformly distributed conditional on the count.

Under the pointwise outer majorant, every potential byte is uniform including
zero. The sources first thin potential occupancies using p=255/256, then
rebase the weighted birth families. The common stochastic rebase preserves
every unpowered conditional moment. Four chronological physical operators
are multiplied before taking the entrywise fractional power. Fine-count
multiplicities are averaged afterward, outside that power.

For one region, the macro recurrence uses eight macros, each with 32 slots.
Thus its placement denominator is C(256,q), not C(512,q). For the complete
code, the group-subset union is a separate C(256,q) factor. The exponent is

```
log C(256,q)
  + alpha * (q log beta + theta*13107)
  + log(e0 R_q^64 1).
```

The subset union is outside alpha. The beta and cutoff factors are inside
alpha because clipping is applied to the conditional first-moment bound for
the fixed subset S. All code paths preserve this distinction.

The complete regional matrix is raised to the 64th power. It is not replaced
by a scalar zero-start regional moment. Hence state is not reset at any
regional boundary. The alpha-one regressions also compare full regional
matrices with the ordinary 32-step placement.

## q1 and the final occupancy sum

The q1 code recovers authenticated active-label operators, not powered or
potential-thinned operators. A region has either zero or one active byte.
The one-byte position is uniform over 256 slots, implemented by 32 physical
steps of eight slots.

For support size v, it computes the coefficient of u^v in

```
e0 (R0+u R1)^64 1,
```

divides by C(64,v), and applies the weight tilt and the trivial probability
cap one. A separate minimum over tilts for each v is valid. Multiplying by
the exact expected shell coefficients and by 256 accounts for all single
nonzero-group messages. There is no extra beta or message-count factor.

For q>=2, every integer occupancy through 256 is evaluated; no interpolation
is used. Selecting the strongest valid witness separately for each q is
legitimate. The final log-sum-exp adds the q1 contribution and all 255 other
occupancy contributions.

## Independent numerical replay

The completed first wider receipt has 40 authenticated source/file pins;
all matched during this audit. Five wider-geometry/coverage tests passed.

The audit rebuilt all eleven powered macro families from the authenticated
active24 caches. Their arrays matched the saved arrays exactly. It then
replayed every q=2,...,256 using scaled placement where safe and 64 sequential
regional row updates, rather than the producer's binary matrix powers.
At the largest tilt, detected scaled underflow triggered the existing full
logarithmic fallback. The largest per-trial margin difference was
3.15e-10 bits.

An independent q1 support-coefficient loop reproduced its saved value.
The resulting figures are:

| Quantity | Independent replay |
|---|---:|
| q1 margin | 110.78588490959346 bits |
| Weakest occupancy | q=25 |
| Weakest individual margin | 55.767198796693535 bits |
| Combined margin | 55.53832939285241 bits |

The combined margin differed from the saved result by 3.77e-13 bits. These
checks authenticate the mathematical wiring and numerical proposal; they do
not replace the separately planned outward evaluation.

## GL3 and scalar24 refresh

The proof is valid with either fresh independent uniform GL(3,GF256) updates
or fresh uniform nonzero GF(2^24) scalar updates in a fixed binary basis.
Both fix zero and send any determined nonzero entering state uniformly onto
the nonzero state space. Conditioning on a fixed message, outer setup, route,
and previous execution therefore gives identical next-state laws. Induction
gives the same emitted-word law for that fixed message.

Binary-adjoint scalar matrices have the same transitivity property. For fixed
a!=0, c maps to L_c^T a is injective because every nonzero L_(c-d)^T is
invertible. Thus the transpose implementation can realize the required law.

The two update ensembles need not have the same joint behavior on multiple
messages. This proof uses only their identical fixed-message probabilities
inside a conditional first moment. Clipping that expected message count
remains valid, without claiming equality of ensemble failure probabilities.

## Implementation obligations and next step

The declared ensemble requires the aligned eight-row symbols and independent
32-bit nonzero-transitive randomizers. Two separate 16-bit randomizers are not
a substitute. Group packet permutations and regional permutations must also
have the stated sizes and independence. A prospective kernel must preserve
the actual24 A/C maps, update freshness, chronology, and continuous state.

Subject to these explicit premises, the bound describes a legitimate fixed
random linear code ensemble. The next proof step is outward evaluation of
the selected witnesses and their sum. The next implementation step is literal
map and adjoint verification for this wider construction, not promotion of
the existing narrow-outer kernel.

## Supplement: complete v2 witness collection

`wider_actual24_allq_v2.json` keeps the same sources, maps, geometry, and
proof model. It adds the witnesses (theta,alpha)=(.2,.35) and (1.4,.5).
All 43 recorded pins matched during this supplementary audit.

The audit independently rebuilt all thirteen macro families from their
authenticated local caches; each matched its saved array exactly. It then
replayed every q=2,...,256 using scaled placement with the existing safe
fallback and 64 sequential regional row updates. The two new families'
largest margin differences were 1.82e-11 and 1.32e-10 bits. Across all
families, the largest difference remained 3.15e-10 bits.

The q1 support replay was unchanged at 110.78588490959346 bits. The weakest
occupancy moved to q=3, with independently replayed margin
68.89467510223749 bits. The combined replay was 68.89399468051678 bits,
within 3.56e-13 bits of the saved value. No new mathematical issue was found.
