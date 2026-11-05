# Retaining the pre-mixer RS support floor

This is a proof refinement for the experimental G16 outer, not a change to
the encoder or a whole-code certificate. The outer uses sixteen RS[16,8]
rows per group, independent nonzero diagonal randomizers, and the affine
two-sum mixer H. The measured 3.398-ms candidate and the previously closed
paired-route checkpoint remain unchanged.

## The information discarded by the earlier bound

Fix the set A of nonzero constituent rows in one group, with |A|=a. Each
nonzero RS row has at least nine nonzero symbols. Across all sixteen
columns, the input to H therefore has at least 9a nonzero coordinates.
The earlier dense-group relaxation retained only nine nonempty columns.
These statements differ substantially when several rows are active.

For a prescribed nonzero row vector x of weight w, the expected number of
RS messages mapping to x after the independent nonzero diagonal is at most
255^(-8). Indeed, the shortened RS code on that support has dimension
w-8, and an information set gives at most 255^(w-8) words with that exact
support. Dividing by 255^w gives the bound. Rows with weight below nine
have zero mass. The comparison for the fixed set A therefore has density
at most 255^(-8a), and we may retain the total support constraint 9a.

## An exact joint counting kernel

For a fixed output support T, let

    Q[A,T;k] = #{x : supp(x) is contained in A,
                    wt(x)=k, supp(Hx)=T}.

These are nonnegative integers. They refine the existing exact count
N_A(T), since summing Q[A,T;k] over k gives N_A(T). The affine kernel makes
Q depend only on a=|A|, b=|A intersect T|, and w=|T|.

The new helper obtains Q by exact integer inclusion-exclusion over input
supports. It first counts inputs with support exactly S from the existing
N_U(T) for U contained in S. It then sums those counts over the size-k
subsets S of A. Numerical evaluation uses the resulting positive
coefficients, not an alternating floating-point sum.

For a fixed marker v>=1, define Q_A(T;v)=sum_k Q[A,T;k]v^k. The inequality

    1[total pre-H weight >=9a] <= v^(total pre-H weight-9a)

holds pointwise. It is inserted before the fractional exponent alpha.
After that exponent and the fixed-row-set Holder step across sixteen
columns, the per-column weight for an output mask of size w is

    W_a(w) = [sum_b C(w,b) C(16-w,a-b) Q_{a,b,w}(v)^(16 alpha)
              / 255^(8a alpha)]^(1/16) v^(-9a alpha/16).

The marker stays fixed for each arity across every state path and every
column. Row-set multiplicities are outside alpha and occur once inside
the Holder sum. This is distinct from the total-arity marker used after
fractionalization, which must not be alpha-powered. The final nonzero
diagonal still supplies independent uniform nonzero output amplitudes
conditional on their supports.

## Checks and first bounded control

The exact-count tests check the labelled ambient kernel, nonnegativity,
normalization by input weight, the marked total mass (1+255v)^a, the
unmarked Holder weights, and the input/output symmetry of the involution
H. A second reviewer checked the density comparison, marker chronology,
and once-only factors independently.

The first FLOAT control uses the saved packetwise fractional banks at
r=18,24,32, with alpha=1/5. It changes neither regional envelopes nor
their slopes. At each anchor it tests five fixed marker templates: zero,
uniform log2(v)=1 or 2, and two arity-specific templates that suppress the
admitted all-empty term. It also retains the low-total-arity restriction
exactly by positive polynomial convolution, rather than a Chernoff bound.

All fifteen controls fail. The unmarked template remains best, with
log2 bounds +6209.338583, +7998.822805, and +10360.850231 respectively.
The exact total-arity restriction does not improve these winners: their
dominant types already lie below the selected cutoffs. The calculation
took 6.206 seconds and built no regional bank. This is a negative result
for these fixed envelopes, not evidence of a low-weight codeword.

The frozen artifact is
`theory/hd_two_sum_input_weight_packetwise_float.json`. Its driver source
hash is `2dad89adad812ec6ed9f810e87c7cf55ac60b3b1fea629d0c37ff794442b64aa`;
the exact joint-kernel source hash is
`5b13c75a6bcc1a47308dd634ed6b3382f67d1a009f6c7ed419881f21aeb7d0cc`.

The next structural proof change separates the empty regional 0-to-0
entry from its nonempty envelope. The current single exponential is
forced through the empty value 1 and loses roughly 28 bits per region
in the r=18 control. The joint support weights can be reused with that
sharper regional calculation.
