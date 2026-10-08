# A half-rate K20 certificate with independent Toeplitz outer blocks

A longer outer block closes the K20 proof without increasing the inner state.
The new outer maps 32 GF16 symbols to 64 symbols using a systematic Toeplitz code.
The inner has three GF16 state symbols and width sixteen.
One structured permutation routes the outer symbols into 64 independently permuted regions.
There is no input diagonal.

The complete directed replay bounds setup failure by

    4.471426459321588e-25 < 2^-80.

This displayed decimal rounds upward.
The result concerns the new ideal ensemble defined below.
It does not change the existing GRS[16,8] certificates or establish performance of a native implementation.

## Construction and finite statement

Let F=GF(16), k=32, b=64, and L=32768.
There are L source rows u_r in F^k, so K=kL=1048576 and N=bL=2097152.
For each row r, sample 2k-1 independent coefficients

    theta_(r,d) <- F,       -(k-1) <= d <= k-1.

Define the k-by-k Toeplitz matrix D_r by

    (D_r)_(i,j) = theta_(r,i-j),       0 <= i,j < k,

and encode the row as O_r(u_r)=(u_r,D_r u_r).
The matrices D_r are independent between rows.
Their coefficients include zero, and setup does not condition on invertibility or another property of D_r.
The systematic coordinates make every O_r injective, including when D_r is singular.

Place coordinate j of every outer row in region j, for 0<=j<64.
Independently sample a permutation of the L individual symbols in each region.
Apply these regional permutations as one global routing stage.
Use each routed symbol directly as an input to the inner.
No routing stage follows the inner code.

Process the region-major stream in sixteen-symbol steps.
Label the sixteen lanes by j in F, and define

    A_j=(1,j,j^2),       C=A^T.

Initialize Q=0 in F^3.
For each step, independently sample M from GL(3,F).
Compute

    R=M Q,       Y=X+A R,       Q'=R+C X.

Carry Q across region boundaries and discard the final state.
All outer matrices, permutations, and refresh matrices are independent setup choices.
Setup fixes the resulting linear encoder E for every source word.

**Finite statement.** Over this complete setup distribution, the probability that
E has minimum symbol distance below 209716 is at most 4.471426459321588e-25.
Thus the ensemble has rate 1/2 and distance at least 10% except with probability below 2^-80.

The inner is injective: from Y and the previous Q, recover X=Y+A M Q, then update Q.
The outer is also injective because it is systematic.
All generator entries lie in GF16.
Consequently scalar extension preserves minimum symbol distance, as in [GF16_PROOF.md](GF16_PROOF.md).
The result applies over GF(2^m) for m=8,32,64,128 with a compatible subfield representation.
It requires no general multiplication in those larger fields.

## Exact averaging over a Toeplitz block

We first establish the source count used in the certificate.
The assertion averages over D; it is not a distance assertion about each sampled outer block.

**Toeplitz lemma.** Fix a nonzero u in F^k before sampling the 2k-1 diagonals of D.
Then Du is uniform in F^k.

To prove the lemma, let j_* be the largest index with u_(j_*) nonzero.
Condition on every diagonal coefficient except

    theta_(i-j_*),       0 <= i < k.

These k remaining coefficients are still independent and uniform in F.
The coefficient of theta_(r-j_*) in output coordinate i is u_(i-r+j_*), when that source index exists.
It is zero when r<i, because the source index then exceeds j_*.
For r=i, the coefficient is u_(j_*), which is nonzero.
The resulting k-by-k map is upper triangular with nonzero diagonal, and therefore invertible.
Every conditioned output is uniform in F^k.
Averaging the conditioned laws proves the lemma.

For subsets S,U of {0,...,k-1}, let N_D(S,U) count nonzero source words u satisfying

    supp(u)=S,       supp(Du)=U.

If S is nonempty, there are 15^|S| possible source words.
For each such word, the lemma gives probability 15^|U|/16^k of parity support U.
Linearity of expectation therefore gives the exact identity

    E_D N_D(S,U) = 16^(-k) 15^(|S|+|U|),       S nonempty.

For S empty, N_D(S,U)=0 by its definition.
The same positive expression consequently upper-bounds every prescribed support of total size w:

    E_D N_D(S,U) <= 16^(-k) 15^w.

For several active rows, independence of the D_r makes their expected support multiplicities multiply.
This step would fail if all rows shared one sampled Toeplitz matrix.
For example, two equal fixed source rows would then have identical parity outputs.
The tests explicitly distinguish this shared-matrix law from the required independent law.

No rank maximum is moved through a sum.
The argument sums source words and support patterns using linearity of expectation and a positive multiplicative envelope.

### The complete value measure supplies the input law

The Toeplitz lemma also identifies the complete average counting measure on outer values.
For every nonnegative function f:F^k times F^k -> R, it gives

    E_D sum_(u != 0) f(u,Du)
        = 16^(-k) sum_(u != 0, v in F^k) f(u,v)
        <= 16^(-k) sum_(u,v in F^k) f(u,v).

This identity holds for an arbitrary f, including a weight function depending on a causal inner computation.
It is stronger than the support-multiplicity identity.
For a fixed set of a active rows, independence of their D_r tensorizes the identity with factor 16^(-ka).

Fix a realization of the regional permutations and refresh matrices before applying this identity.
Choose f to mark the output weight of that fixed inner encoder.
The right side then sums over a product counting measure on all 64a outer symbols.
Within a prescribed support, its nonzero values are independent and uniform in F minus zero.
A fixed route only reorders these values.
Conditional on earlier values and states, the unprocessed nonzero values remain independent and uniform.
Averaging the independently sampled refresh matrices supplies the same local inner transfer as before.

This establishes the input-value law for the construction without input scalars.
This conclusion uses the averaged message-count measure before the inner analysis.
It does not assert that an individual fixed source word has independent systematic coordinates.
The earlier scalar-decorated variant has the same upper bound, but the strict receipt now identifies the no-diagonal ensemble.

## Composition with the regional inner

Fix an outer support pattern under the dominating product counting measure derived above.
Its routed nonzero symbols are independent and uniform in F minus zero.
Future input values remain independent of the previous inner state.
The local transfer therefore applies after averaging the outer message count, without sampling an input diagonal.

Use the local transfer from [GF16_SIZES.md](GF16_SIZES.md), with state dimension s=3.
Let T_w(z) be its entrywise upper-bound matrix for a uniform lane mask of weight w, with 0<z<1.
Its rows and columns record whether the state is zero or nonzero, and z marks output weight.
The powers of the field labels through degree four sum to zero, so C A=0.
The same MDS and invertible-return bounds therefore apply.
The exact singleton joint return moment for the stated monomial maps is

    J_1(z)=(z^13+12 z^14+3 z^15)/(16*(16^3-1)).

The implementation computes these sixteen singleton weights directly.
It reuses the existing directed local envelope with state=3.

Define

    F(v,z)=sum_(w=0)^16 binom(16,w) T_w(z) v^w,
    R_h(z)=[v^h] F(v,z)^(L/16) / binom(L,h).

The denominator is the regional placement count.
It follows from a uniform permutation of individual symbols within each region.
The product of region matrices carries the state across all boundaries.

Fix a set of a active source rows.
The Toeplitz support envelope contributes 16^(-32a) times 15 for each active outer coordinate.
In a region containing h active coordinates, there are binom(a,h) possible active-row subsets.
Define

    M_a(z)=sum_(h=0)^a binom(a,h) 15^h R_h(z).

Summing all support patterns gives the following upper bound on the a-class contribution:

    B_a(z) = binom(L,a)/15 * 16^(-32a) * z^(-209715)
               * [1,0] M_a(z)^64 [1,1]^T.

The division by 15 removes scalar multiples of each complete nonzero source word.
It is applied once, not once per active row.
The exponent 209715 is the largest bad output weight.
The positive source envelope includes impossible patterns with no systematic support; those extra terms only increase B_a.

For larger classes, fix x>0 and define

    H(x,z)=sum_(w=0)^16 binom(16,w) T_w(z) (15x)^w (1+x)^(16-w).

The binomial identity

    binom(a,h)/binom(L,h)=binom(L-h,a-h)/binom(L,a)

gives

    M_a(z)=[x^a] H(x,z)^(L/16) / binom(L,a).

Every coefficient is nonnegative.
Evaluating at positive x and multiplying 64 region bounds yields

    B_a <= binom(L,a)^(-63)/15 * 16^(-32a)
             * z^(-209715) * x^(-64a)
             * [1,0] H(x,z)^(4L) [1,1]^T.

For fixed rational x,z, the right side is log-convex in a.
The binomial coefficient is log-concave; all remaining a-dependent factors are geometric.
The larger endpoint bound therefore covers every integer class in an interval.

## Complete replay and scope

Fifteen rational sparse witnesses cover a=1 through 15.
Eleven rational endpoint witnesses cover a=16 through 32768, including the dense classes.
No separate dense-tail estimate is needed.

| Active source rows | Directed contribution, rounded upward here |
|---:|---:|
| 1 through 15 | 4.471377762564247e-25 |
| 16 through 32768 | 4.869676e-30 |
| Complete union | 4.471426459321588e-25 < 2^-80 |

The fixed rational witnesses are in `theory/rank_continuation_toeplitz_witnesses.py`.
Replay uses exact integer counts and 80-digit directed Decimal operations.
It temporarily widens the Decimal exponent range for positive matrix powers, then restores the previous limits.
Floating searches are separate from replay.

From the repository root:

```sh
python -B research/workstreams/char2_field_spin/theory/rank_continuation_toeplitz_certificate.py --certify
python -B research/workstreams/char2_field_spin/theory/rank_continuation_toeplitz_certificate.py --verify PATH_TO_RECEIPT
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_rank_continuation_toeplitz.py -v
```

The verifier regenerates every sparse class and interval endpoint.
It checks complete coverage, displayed endpoints, contributions, sums, and claimed margins.
Eighteen tests cover exhaustive GF4 Toeplitz output laws, exact support counts, the complete value-measure identity, row independence, receipt metadata and mutations, size replay, and Decimal context restoration.

### The same profile at K16 and K18

The separate size wrapper retains the no-diagonal Toeplitz[64,32] outer, 64 regions, and the same width-16, three-state inner.
Only the regional length L and bad-weight threshold change.
At K16, L=2048 and the bad-weight threshold is 13107.
At K18, L=8192 and the bad-weight threshold is 52428.
Substituting these values into the preceding bounds gives independently replayed rational certificates:

| Source symbols | Minimum distance | Failure bound, rounded upward | Integer margin |
|---:|---:|---:|---:|
| 65536 | 13108 | 2.821241822704024e-26 | 84 bits |
| 262144 | 52429 | 1.119712901898988e-25 | 82 bits |
| 1048576 | 209716 | 4.471426459321588e-25 | 80 bits |

Both smaller sizes use fifteen direct sparse classes.
Seven endpoint intervals cover the remaining K16 classes; nine cover the remaining K18 classes.
The wrapper's K20 replay agrees exactly with every contribution and endpoint from the fixed K20 engine.
This comparison checks the parameter substitution without changing the fixed K20 witness interface.

```sh
python -B research/workstreams/char2_field_spin/theory/rank_continuation_toeplitz_sizes.py --k-log 16 --certify
python -B research/workstreams/char2_field_spin/theory/rank_continuation_toeplitz_sizes.py --k-log 18 --certify
python -B research/workstreams/char2_field_spin/theory/rank_continuation_toeplitz_sizes.py --verify PATH_TO_SIZE_RECEIPT
```

The fixed size witnesses are in `theory/rank_continuation_toeplitz_size_witnesses.py`.
Generated no-diagonal receipts are kept under the ignored directory `tmp/char2-field-spin`, with names `rank_continuation_toeplitz_k16_nodiagonal_receipt.json`, `rank_continuation_toeplitz_k18_nodiagonal_receipt.json`, and `rank_continuation_toeplitz_k20_nodiagonal_receipt.json`.
The earlier scalar-decorated receipt is preserved separately.

The shorter [32,16] Toeplitz outer fails the present 40-bit one-row screen.
The [64,32] outer with two state symbols passes the sparse screen but fails around 5000 through 6200 active rows.
Those failures concern upper-bound screens, not distance counterexamples.
The three-state profile above has complete finite coverage.

An independent audit rederived the source law and the normalization, and replayed every endpoint with exact local fractions and 120-digit directed arithmetic.
Its final upper bound is below the bound reported here.
It also confirmed the simplification without input scalars.
The same independent replay subsequently checked every K16 and K18 class and endpoint, and confirmed the 84-bit and 82-bit claims.
The next step is full native encoding measurements.
Toeplitz multiplication can use a fixed block-local circuit; any local schedule must preserve the stated independently sampled matrices.
No encoding-time advantage is claimed here.
