# Structural checks for the two-symbol batched inner

The candidate uses a GRS[8,4] outer over a characteristic-two field F.
One coordinate from each outer group goes to each of eight regions.
Each region independently permutes its L groups. Thus N=8L and K=4L.
Assume t divides both N and the region length. One physical step consumes
t symbols, for t in {8,16}.

Let A have row (1,j) at position j=0,...,t-1, with j interpreted in the
declared polynomial basis. Put C=A^T. The state starts at zero, and

    Y_i = X_i + A Q_i,
    Q_(i+1) = M_i Q_i + C X_i.

The final state is discarded. Each M_i is the product of three independently
sampled nonzero-coefficient shears described in [the local note](README.md).
The following obstructions hold for every realization of routing and M_i.

## A ceiling from the transparent subcode

There are R=N/t steps. Choose q outer groups and set every other message
block to zero. This leaves 4q message variables. Impose C X_i=0 at the
first R-1 steps, giving at most 2(R-1) homogeneous equations. If

    q = floor((R-1)/2) + 1 <= L,

some nonzero message satisfies all these equations. Every output-relevant
state remains zero. Its encoded output is therefore Y=X and has weight
at most 8q. The last update requires no constraint because its state is
discarded. For even R, the bound becomes

    d_min / N <= 4/t.

Thus t=16 cannot reach 50% distance: its ceiling is 25%. For t=8 the
ceiling is 50%. Neither result obstructs the current 10% target. This
argument applies even when C varies between steps, provided each C has
two rows. Independent input scalars do not remove the ceiling.

More generally, outer blocks with dimension k and length b, an s-symbol
state, and identity direct feedthrough have d_min at most
b*(floor(s*(R-1)/k)+1), when that many outer groups exist.
Avoid interpreting greater setup randomness as a way around this rank bound.

## The chosen maps make the inner an involution

The t labels form an additive binary space of dimension at least three.
Their sum is zero, their number is even, and sum_j j^2=(sum_j j)^2=0.
Consequently C A=0. Every step preserves its two checks: C Y_i=C X_i.

Applying the inner again with the same M_i produces the same state sequence,
starting at zero. Its output is Y_i+A Q_i=X_i. Thus I^2 is the identity.
This is an algebraic property, not by itself a low-distance witness.
The transparent-subcode bound above is the relevant resulting limitation.

## No one- or two-group transparent witness

Any two columns (1,j) and (1,k) of C are independent when j!=k.
Hence every nonzero vector in ker C has weight at least three.
An outer group contributes at most one coordinate to each region and
therefore to each physical step. A transparent witness with only one or
two active outer groups would give every step input weight at most two.
The equations C X_i=0 at the first R-1 steps would force every input
outside the last step to zero. An active GRS group has nonzero coordinates
in at least five regions, so it cannot be confined to that last step.

For three active groups, every nonzero constrained step must contain all
three groups. Their active-column supports coincide in the first seven
regions, since a coordinate contributed by only one or two groups cannot
cancel. The GRS minimum distance requires at least four such regions;
one further active coordinate may lie in the last region, whose final
step is unconstrained. At each of those four regions the three groups must
occupy one physical step. For any fixed three groups, a uniform regional
permutation gives this event probability

    (t-1)(t-2) / ((L-1)(L-2)).

Different regions are independent. A union over four of the first seven
regions therefore bounds this necessary event by 35 times the fourth
power of that probability. If C X_i=0 is also imposed at the final step,
at least five regions are required instead. This is a useful sparse obstruction
count, not a whole-code certificate: witnesses can also enter and leave
nonzero states, and larger active-group sets need separate treatment.

A different deterministic transparent pattern is harmless at 10%: assign
the same outer polynomial to every group. Each region then contains a
constant value. Constant physical-step inputs lie in ker C, but a
nonzero degree-three polynomial has at least five nonzero evaluations.
The resulting output weight is at least 5N/8.

## Which cheap randomizers change the construction?

Replacing j by alpha*j+beta for alpha!=0 only left-multiplies C by an
invertible two-by-two matrix. It leaves ker C unchanged. In particular,
an XOR offset does not make state activation harder to cancel.
Binary-affine permutations of lane labels preserve another conspicuous
pattern: equal values at the four corners of an affine square have both
sum and weighted sum zero. A separate row-scramble pass is therefore not
automatically a useful new randomness source.

In fact, arbitrary within-step permutations are already absorbed by the
uniform regional group permutations. Conjugating A and C by such a
permutation is equivalent to permuting the step input and output. The
input composition remains uniform within its region, and output permutations
preserve Hamming weight. This does not improve the ideal ensemble.

One scalar per complete outer group is also ineffective: it is an invertible
change of that group's message coordinates. The code image is unchanged.
Independent nonzero scalars on separate outer coordinates do change the
alignment distribution. They preserve the MDS outer and give fresh
independent coefficients for quotient-rank bounds. Those scalars are the
meaningful next comparison, if the unscaled candidate needs more proof slack.

## GF256 coefficients are a distinct efficient option

All requested fields contain GF256. A generator over GF256 has exactly the
same symbol distance over each extension, by scalar-extension invariance.
This is enough for a 10% target. The rate-half GF256 GV reference near 38%
is an existence benchmark, not an upper bound on achievable distance.

Using a GF256 basis for a larger field makes multiplication by a GF256
coefficient bytewise. It is then field-linear and compatible with GFNI.
This requires a declared compatible representation. An existing consumer's
polynomial basis may require conversion, which belongs in measured work.

Three shears with coefficients uniform in GF256 minus zero have vector
target atoms at most 1/255^2, versus 1/(256^2-1) for uniform GL(2,GF256).
The former bound is only about 0.0113 bits weaker. It does not make the
shear distribution uniform. The complete batched proof must handle its
joint quotient transitions; a bound on fixed vector atoms is insufficient.

Next retain t=8 and t=16 for the 10% performance comparison, but prefer
t=8 when a future high-distance target matters. Analyze the unscaled
regional construction first, and add independent column scalars only if
their proof gain justifies their arithmetic and setup cost.
