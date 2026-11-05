# What the inner proof needs, and cheaper ways to provide it

The34.5ms implementation is not evidence of a design limit. It pays for
almost-uniform feedback on every live batch, which is stronger than the
distance event itself requires. This review separates that convenient
interface from the actual counting obligations. It also closes a cheaper
coefficient-field variant and specifies a more substantial inner candidate.
No new native timing or library promotion is claimed.

## The actual proof interface

Keep K=2^20, N=2^21, RS[256,128], one independent uniform permutation per
region, and512 state symbols. Setup is fixed for every message. Write C_pi
for the routed outer code. For the running-sum inner, fixing a candidate
output Y determines every state Q_t without revealing any feedback key.
The inverse is

    X_t = Y_t - R_t Q_t.

The inverse proof needs a bound on

    Pr_R[X(Y,R) belongs to C_pi],

summed over low-weight outputs. It does not intrinsically require every
feedback coordinate to be uniform. The forward proof instead counts outer
inputs and the probability of specified state resets. The route controls
the shortened dimensions of C_pi on unions of scan intervals.

The retained proof obtains both bounds from the stronger pointwise condition

    Pr[R_t z=f] <= eta q^(-b)  for every z!=0 and every f.

This discards all structure in the outer fiber, the state, and the prescribed
output. Two full-field diagonals around a superregular map meet this condition;
they are not logically required by minimum distance. For a weaker family,
we must retain the exceptional states or affine constraints it leaves behind.
A good fixed-message probability alone does not control a message chosen
after setup.

Three different resources must not be conflated:

1. Sufficient effective state dimension prevents an outer subspace from
   imposing many cheap resets. More random coefficients cannot repair too
   few state dimensions; RANDOM_FEEDBACK_BRIDGE.md gives explicit examples.
2. Reset equations reduce the number of low-weight output trajectories.
   Their rank depends on active coordinates and episode lengths.
3. Feedback supplies anti-concentration against the actual outer constraints.
   Requiring a fresh almost-uniform512-vector at every step is one sufficient
   method, not the only possible one.

## A concrete loss in the old proof: reset rank

Fix c episodes occupying v batches. The first r=c-1 must close; the final
episode has a free terminal state. A closed episode of length L_i satisfies
one sum-zero equation per state coordinate. Let k_i be the number of its
active coordinates, and w_i its output weight. On that fixed support,
these equations have rank exactly k_i, since their variables are disjoint.

Put R=sum_i k_i. Each closed episode is nonempty, so R>=r. More importantly,

    w_i <= L_i k_i <= 2 k_i + b(L_i-2).

Adding the free final episode gives, for total weight w,

    w <= 2R + b(v-2r).

The old bound retained only R>=r. On a two-batch episode with230 active
coordinates, it charged one equation instead of230. This is particularly
costly when there are many short episodes.

For a fixed support of weight w, the number of amplitudes is at most q^(w-R).
There are at most2^(bv) supports. Therefore the number of candidate outputs
of weight at most W, for fixed episode placements, is at most

    2^(bv) q^((W+b(v-2r))/2).

Independently, the unrestricted ambient reset equations have rank br, giving
q^(bv-br). Take the minimum of these two bounds and the old output-volume
bound. For a rectangle of classes, replace v by its upper endpoint and r
by its lower endpoint. Both new expressions are monotone in these directions.
The route moment and the forward bound remain unchanged.

This also holds for Q_(t+1)=delta Q_t+Y_t with fixed delta!=0. Reset equations
then have nonzero geometric weights, but their support ranks are identical.

## Closed consequence: GF32 coefficients suffice for46%

Re-run the proof over J=GF(2^32), rather than merely using32-bit keys in a
GF128 counting argument. The output-volume entropy term changes accordingly.
The old coarse reset count fails at45%; the stronger count above closes46%.

| Distance threshold | Rectangles | Feasible(v,c) pairs | Worst directed log2 class bound |
|---|---:|---:|---:|
|45%|289|4196352|less than -63.65608869973025|
|46%|459|4196352|less than -65.77481065587761|

Independent integer-slice coverage checks and corrupted-cover controls pass.
Both80- and120-digit directed replays pass. The union uses fewer than2^23
classes, so the respective uniform-feedback margins exceed40.656 and42.774
bits. These rounded inequalities are conservative displays of the receipts.

The same sandwich lemma applies over J. With nonzero J-valued diagonal
keys and fixed D over GF(2^16), its total inflation is

    (1-2^-32)^(-N) <= (1-N/2^32)^(-1) = 2048/2047 < 2^(1/1024).

The last inequality follows by applying Bernoulli to(1-1/2048)^1024.
Thus both profiles retain at least40 bits for the actual two-diagonal
distribution. The46% profile even retains more than42 bits under this bound.
This is a whole-code certificate for RS256/state512, not just a sparse class.

Now extend the sampled generator matrix G from J to F=GF(2^128). Write each
message as x=sum_i beta_i x_i in an F/J basis. Then Gx=sum_i beta_i Gx_i,
and support(Gx_i) is contained in support(Gx). A nonzero x has a nonzero
component, and G is injective. Consequently d_min over F is at least d_min
over J; embedding J-messages gives the reverse inequality. The two distances
are equal for every setup, so no union factor for the four lanes is needed.

The resulting code remains F-linear and accepts128-bit symbols. Its diagonal
scalings are GF32 scalars applied to four lanes, not full GF128 products.
The fixed MDS map and the outer retain their existing representations.
This is a deliberate coefficient-subfield variant, not a claim of uniform
GF128 mixing or an established49% profile. It has not been benchmarked here.

Partition hashes, in table order:

    f437e05f84714cc185af584ec16deae4083ffd4fbe26855e9e97e515cb8ef1cd
    faf9d54b427678b812fae5cfbdbcb20060b65385c66b4c8af6a5718bd00bef6e

Reproduce from coded_state, serially:

    python -B grouped_reset_certificate.py --outer 256 --group 1 --field-bits 32 --sum-rank --distance-permille 460 --precision 80 --seconds 20

Repeat at precision120, and distance450 for the other receipt. Defaults
preserve the old128-bit proof and its original partition hash. The frozen
generic proof helpers were not edited.

## Why independent tiny groups are not enough

There is an RS word supported on the last129 of the256 outer coordinates:
take the polynomial vanishing at the first127 evaluation points in one block.
After routing, it supplies one nonzero input in each of those129 regions.

Suppose the b=512 state and output coordinates are partitioned into fixed
independent groups of size g, with no communication between groups at any
time. An input in region j can affect at most g output symbols per later
batch in its group. Each region has M=8192 symbols, so its contribution is
at most g(N-jM)/b. Summing these possibly overlapping supports gives

    d_min/N <= g/(512*256) sum_(j=127)^255 (256-j)
             = 8385g/131072.

This is at most6.40%,12.80%, and25.59% for g=1,2,4, respectively. It holds
for every route and every within-group coefficient choice. Stronger local
randomness does not repair this obstruction. Groups must communicate.
It does not require them to communicate globally in one batch.

## Proposed next inner: one layer of coupled pairs per batch

Keep the outer, route, b=512, and128-bit payload symbols. Work first with
coefficients in J, while retaining F-linearity. Choose fixed
delta in GF256 other than0 or1, for example the AES element0x02.

At batch t, pair state coordinates whose indices differ in bit t mod9.
For each pair, independently sample a from J\{0} and c from

    S = J\{0,1,1+delta,delta+delta^(-1)},

then set b_key=c/a at setup. This uses b_key for a scalar key, not the state
dimension. Given pair(x,y), compute two shears:

    p = x + a*y,
    r = y + b_key*p.

Let A_t be this block-diagonal pair map. Emit and update

    Y_t = X_t + (A_t+delta I) Q_t,
    Q_(t+1) = A_t Q_t + X_t = delta Q_t + Y_t.

The final expression preserves the key-independent inverse-state interface.
The two shears use two variable products per pair, or one per emitted symbol.
Forming feedback adds one fixed GF256 scaling per symbol. There is no large
FFT, global MDS application, or second global routing pass. Pair wiring is
local to the512-symbol state. Division and excluded-key sampling are setup
operations, not encoding operations.

The two relevant pair matrices are

    A = [[1,a], [b_key,1+c]],
    A+delta I = [[1+delta,a], [b_key,1+c+delta]].

The restrictions on c make every entry and determinant of both matrices
nonzero. Indeed det(A)=1 and det(A+delta I)=(1+delta)^2+delta*c.
Hence both pair maps are superregular. For any nonzero pair z,

    wt(z)+wt(Az)>=3,       wt(z)+wt((A+delta I)z)>=3.

Consequences proved locally:

- Zero input cannot extinguish a nonzero state, because A_t is invertible.
- A nonzero state on zero input cannot emit a wholly zero batch, because
  A_t+delta I is invertible.
- One active coordinate spreads to all512 coordinates after nine successive
  empty-input batches. Each stage introduces a previously unused index bit;
  no paths collide during that first expansion, and every matrix entry is
  nonzero. This is not a claim about arbitrary multi-coordinate states.
- Fixed outputs reconstruct the state independently of all pair keys.

The exact local distribution exposes the remaining proof obligation.
For a fixed pair(x,y), feedback(f0,f1) satisfies

    f0=(1+delta)x+a*y,
    f1=b_key*(x+a*y)+(1+delta)y.

If y!=0 and p=x+a*y!=0, a prescribed feedback fixes both keys. Its probability
is at most1/((|J|-1)|S|). If y=0 but x!=0, one key remains unobserved and a
compatible point has probability1/(|J|-1). The same larger point mass occurs
when p=0. At zero state, feedback is deterministically zero.

Therefore this candidate does NOT satisfy the old uniform-batch lemma.
The missing randomness has identifiable causes: zero state coordinates and
one exceptional cancellation. For a fixed output, the former are linear
equations on its history; the latter is a specific value of a fresh key.
A global proof must count the rank of those history constraints jointly with
the key exceptions and the outer fiber. Repeated zero-state constraints can
be dependent and cannot simply be charged once per occurrence. General
states can cancel despite the single-coordinate diffusion result.

This is the next discriminating proof task, not a claimed distance certificate.
Retain the full sandwich as a certified control. If pair mixing loses too
much, add a second local layer per batch or periodic stronger mixing, and
measure the added proof strength against its explicit operation cost.

## Checks, diagnostic scope, and next actions

test_inner_requirements.py passes68 exact reset-count comparisons, four
scalar-extension examples,255 complete GF16 local feedback laws across180
allowed key pairs, and160 singleton-diffusion fixtures. The algebra above,
not these toy examples, establishes the local general statements.
Regression tests compare32/64/128-bit wrappers to the retained generic bound.
The old45% GF128 receipt still has its original hash and bound.

The weaker point-density interface itself also has slack: FLOAT covers at45%
complete with losses of1,8, or16 field dimensions per live batch, using the
new reset count over GF128. The32-dimension screen fails. These are diagnostic
covers, not directed receipts or a certificate for the pair family, whose
exceptional states need not obey any such uniform deficit.

No encoder benchmark ran during this review. Tests and screens were serial,
time-bounded, and source-only; local free space remained above20GiB. Existing
paper, production code, and prior constructions are unchanged.

Next: use the certified GF32 sandwich as the cheaper control. For the larger
structural improvement, derive a multi-batch rank/output count for the actual
pair-key distribution, beginning with repeated sparse-state histories and
reset schedules. Do not substitute average mixing speed for that count.
