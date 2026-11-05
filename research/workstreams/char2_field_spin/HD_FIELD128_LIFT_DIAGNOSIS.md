# Why the literal alphabet lift still needs a stronger bound

The repaired directed calculation covers every active-group weight from
1 through 64 at 35%, with 130.421601248370 bits of margin. It does not
cover the remaining weights. At 40%, the latest weight-64 FLOAT search
returns a positive log bound, +109.09573184619, after 48 evaluations.
Its witness is z=0.49785001765857895, alpha=0.006820573473687547;
the largest matching-profile term has (a,b)=(22,31).
These observations come from the root task's retained runs.
Neither numerical optimization nor this diagnosis proves an obstruction
to the code's actual distance, or a global optimum of the comparator.

## What increasing the alphabet does and does not buy

Use the [literal-lift ensemble and operator](HD_FIELD128_LIFT_PROOF.md).
Write b=log2(q), theta=alpha*b, and lambda=-log2(z)/b. A local term
q^(-d) z^w becomes, after the fractional power,

    (q^(-d) z^w)^alpha = 2^(-theta*(d+lambda*w)).

The leading outer and projective prefix becomes

    (q^(128*r) / (q*z^L))^alpha
        = 2^(theta*(128*r-1+lambda*L)).

Thus the leading rank terms depend on theta, not on b separately.
Support multiplicities and route probabilities remain outside alpha.
More formally, if I_ij(lambda) is the large-b rate
`-log2(M_ij(q^(-lambda);q))/b`, then
`M_ij^alpha` tends to `2^(-theta*I_ij(lambda))` at fixed theta,lambda.
The fixed GF256 ranks and route geometry survive this limit.
This is not exact alphabet independence: finite polynomial sums and
their coefficients still matter, and small local differences accumulate
over 32768 packets. No monotonicity claim about the optimized bound follows.

At the weight-64 witness, theta=0.873033404632006 and
lambda=0.007861069637329. The projective divisor supplies about 0.873 bits;
a rank-16 return factor supplies about 13.969 bits after alpha, not 2048.
The same scaled marks over GF256 would be z=0.957345406761104 and
alpha=0.109129175579001. This explains why a larger alphabet does not
remove a support-geometry bottleneck by itself.

The current witness needs 149.096 bits of improvement to reach a -40
log bound. Raising the target distance by one percentage point costs
about 143.927 bits at these fixed marks. The corresponding deficits at
45% and 49% would be about 868.731 and 1444.439 bits. These are witness
sensitivities, not optimized distance predictions or whole-code margins.

## First proof-only experiment: delay alpha across four packets

The existing [good-pair argument](HD_SHIFTED_GOOD_PAIRS.md) composes two
raw state operators before alpha when both packets have at most one
singleton pair. It still separates state paths at the next boundary.

Choose a smaller collection V of good arrangements, for example those
with at most two full pairs. Partition the existing positive two-packet
polynomial as `T2=C_VV+R`, where R contains every other case. Enumerate
four V arrangements in chronological order and form C_VVVV by applying
alpha only after their four raw matrices have been multiplied. Then

    T4 = C_VVVV + C_VV R + R C_VV + R^2 <= T2^2

coefficientwise. The inequality follows from fractional subadditivity;
all four displayed summands are positive and partition the cases.
The regional polynomial is T4^32, with the same exact slot normalizer.
Restricting V bounds the preparation cost without dropping any arrangement.
An all-good four-packet calculation can follow if this smaller test helps.

This changes neither the encoder nor its setup. It requires a new
coefficient provider and directed replay, including fair-swap counts.
Compose raw matrices, not matrices already raised to alpha.
Isolated empty stretches offer no improvement: their operator is diagonal,
so raw powers commute with entrywise alpha. A useful block must include
activation/return alternatives, possibly separated by an empty gap.

## Structural alternative: one RS[64,16] emission map

Keep the outer, route, 64-symbol packet width, 16-symbol state, and refresh
distribution. Replace the four RS[16,4] blocks by one RS[64,16] map A;
set C=A^T. A concrete candidate is the Lagrange evaluation map from
nodes 128,...,143 to nodes 0,...,63 in AES GF256, embedded in GF(2^128).
It is already constructed as a fixed map by
[`hd_blockwise_maps.maps`](theory/hd_blockwise_maps.py), but no retained
numerical result for another ensemble certifies this proposed lift.

The 64 evaluation nodes form an F2-linear additive subspace U of GF256.
They are not a GF64 subfield: GF64 does not embed in GF256. Its subspace
polynomial `P_U(X)=product_(u in U)(X-u)` has constant nonzero derivative.
The Lagrange leading-coefficient identity therefore gives
`sum_(u in U) f(u)=0` for every polynomial of degree below 63.
Products of two degree-below-16 basis polynomials have degree at most 30.
Consequently A^T A=0, so CA=0 and (I+AC)^2=I.

A nonzero emitted state has weight at least 49. Every j-column
restriction of C has rank min(16,j), and its kernel has minimum weight
17. In contrast, the four-block map permits rank saturation at four
inside a block and has emission minimum 13. These are genuine algebraic
changes, not a change of coordinates in the existing state space.

A tiny exact GF256 check reproduced CA=0 and all power sums of degrees
0 through 62. The concrete row-major A digest is
`4a8399e103c2fb2ca5824883d7da1f6d3e235c4e5873c8467701a6cde071586a`.
This check is not a distance gate or a native implementation result.

### A local operator without a 64-coordinate nullity scan

Fix a potential input set S of size j; draw X uniformly on F^S.
Let `t=z+(1-z)/q`, `d=q^16-1`, `B=d/(q-1)^16`, and rho=min(16,j).
The MDS rank identities give

    K_j = sum_(c=0)^j binom(j,c) z^(j-c)(1-z)^c
                         q^(-c-min(16,j-c)),
    T_j = t^j,
    f_n = sum_(e=0)^n binom(n,e) z^(n-e)(1-z)^e q^-min(16,e),
    V_j = T_j (q^16 f_(64-j)-1)/d.

Define the actual return moment
`J_S=E_X[1[CX!=0] z^wt((I+AC)X)]`. A first cap is

    Jbar_j = min(T_j-K_j,
                 (1-q^-rho) z^max(0,49-j),
                 R_j, d V_j),

where the following positive bound uses only MDS ranks:

    R_j = sum_(c=0)^j sum_(e=0)^min(15,64-j)
          binom(j,c) binom(64-j,e) z^(64-c-e)(1-z)^(c+e)
          * min(1-q^-rho,
                q^(max(0,j-c-16)-j) (q^min(rho,16-e)-1)).

To justify the second cap inside this sum, write R=CX. Requiring e
inactive outputs to vanish leaves at most `q^min(rho,16-e)-1` nonzero
syndromes. Fixing R and c active output zeros fixes those c inputs.
The remaining feedback equations have at most `q^max(0,j-c-16)` solutions.
For e>=16, no nonzero R survives. Separately, invertibility of I+AC and
an information set give the T_j-K_j cap; emission distance gives 49-j.

The same zero/nonzero density comparison therefore permits

    M_j = [[K_j,       B(T_j-K_j)],
           [Jbar_j/d,  B V_j]].

For j=0 the live entry can also be capped by z^49: every inverse refresh
is bijective and every nonzero state emits at least 49 symbols.
These caps depend only on total occupancy j, so they reuse the route
counts without enumerating all labelled 64-coordinate supports.
They may be loose. Computing the full principal-nullity table over 2^64
subsets is not a reasonable prerequisite; selected joint-rank refinements
should follow only if the coarse cap fails a bounded witness test.

State size, full-field diagonal products, outer full-field products, and
key bytes stay unchanged. A dense A or C has 1024 fixed coefficients
instead of 256, so the naive fixed-map work is four times larger.
An additive-FFT implementation may reduce that cost; no timing is inferred.
Historical h16 screens in [the grouped-packet note](HD_GROUPED_PACKET_PROOF.md)
were mixed and used different route classes and return caps. They do not
establish superiority of this candidate.

## Cheap mixing: distinguish a basis change from a new return map

For any fixed invertible H, simultaneous replacement

    A'=AH,       C'=H^-1 C

leaves A'C'=AC. It preserves K, T, V, and every exact J_S, as well as
the emission image and feedback ranks. Thus the uniform-density packet
operator is unchanged. Conjugating the refresh to H^-1 R H makes the
whole encoder identical after a state-coordinate change. Leaving its
distribution unchanged may change the ensemble and a specialized empty
cap, but supplies no new rank or return-weight gain in this comparator.

A meaningful alternative is one-sided mixing: use A'=A, C'=HC, or
A'=AH, C'=C. Then CA remains zero but the return map becomes I+AHC.
The quiet and uniform-live moments remain unchanged; only the return
moment needs a new bound. One fixed H16 transform adds no full-field
random products or key bytes. The earlier
[scattered-check proposal](HD_COMMON_REFRESH.md#scattering-check-coordinates-without-an-extra-mixer)
is an even cheaper lead: it changed feedback constants and a permutation,
with singleton return weights 56 through 63 instead of 15 or 16.
Its old screens and timing do not transfer to the present route or field.

## Recommended order

1. Test the positive four-packet refinement at weight 64 and the recorded
   profile/witness. This directly measures recoverable fractional-path
   loss, with zero encoding cost. Then cover all profiles at that weight.
2. Inspect the return contribution before investing in one-sided mixing
   or sharper labelled return enumerators. A return-only change cannot
   improve the unchanged quiet branch. Retain averaging/alpha order;
   an average return cap is not automatically a pointwise cap.
3. If the rank geometry remains decisive, test the stated full-width
   RS operator before implementing it. Increasing state to 24 or 32 is
   a later alternative: it strengthens the return-rank charge but adds
   refresh keys/products and weakens each fixed-length RS emission distance.

More mixing within the same outer group cannot improve the leading
q^128 density factor: it reflects dimension 128 inside length 256.
Replacing Toeplitz by a more uniform mixer changes only negligible
conditioning corrections at q=2^128 under this comparison. Useful outer
work must instead retain RS correlations in the proof or genuinely mix
encoded coordinates across groups. An invertible premix of all source
messages merely reparametrizes the same code and cannot raise its distance.

The next decision is a bounded proof comparison, not another alphabet
increase or a native implementation. No new 40%, 45%, or 49% certificate,
large numerical gate, benchmark, production change, or paper change is
claimed here.

## A possible cheaper full-field outer mixer

Identify four F symbols with one element of a degree-four extension E/F.
Multiplication by a uniform nonzero element of E sends any fixed nonzero
input uniformly over E minus zero. It is F-linear. Replacing each
Toeplitz mixer by this operation would therefore give atom bound
1/(q^4-1) and outer density q^256/(q^4-1)^32. The projective divisor stays
q-1. This improves the negligible conditioning factor, not the leading
dimension penalty or the known distance obstruction.

The possible benefit is implementation cost: four F coefficients replace
seven. A native binary-polynomial 512-bit multiplication has a 27-product
64-bit Karatsuba decomposition, before reduction. This is not a fair
27-versus-81 speedup claim. Native Toeplitz also needs 27 raw products;
its nine intermediate products may share reduction into four final
outputs. In the existing byte-compatible representation, a straightforward
64-byte polynomial product instead uses 1024-bit spaced coefficients and
81 raw products. Moreover, concatenating four current F representations
does not automatically give a compatible native E basis. Conversion and
reduction costs must be included, and omitting a required conversion can
break F-linearity. No implementation, field-basis construction, or timing
has been produced for this option.
