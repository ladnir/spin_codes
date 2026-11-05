# Global classes defined by covered inputs

The construction remains fixed. Two further classes are now proved:
full/empty covered footprints with arbitrary internal states, and a
footprint-expansion class. Neither is the whole 45%/40-bit certificate.
The distinction matters: a partially occupied state can still cover the
entire input update, and should not automatically count as a sparse
footprint in the routing argument.

Use q=2^32, N=2097152, K=1048576, M=8192 outer rows, T=4096 updates of
width512, and W=943718. Fix candidate covered masks V_t before exposing
setup, with V_0 empty. Let F=sum|V_t|, R=sum_i min(128,|V intersect row i|),
and d=F-R. Let A_t contain V_t positions whose quartet at update t+1 is
wholly covered; set A_T=V_T and L=sum|A_t|. These are the same objects as
in [the grouped-footprint reduction](QUARTET_GROUPED_FOOTPRINT.md).
All statements use the actual once-sampled accepted setup and free terminal.

## Full/empty footprints, arbitrary state histories

Except with probability less than 2^-148 over setup, no nonzero J-message
of output weight at most W has V_t either empty or full at every update.
The states inside a full footprint are unrestricted. This strictly
extends the earlier [state-burst theorem](QUARTET_GLOBAL_BURSTS.md).

Consider one full covered run of length ell. Its first state is full:
the previous footprint is empty, so there is no incoming state and all
current coverage consists of births. If another empty footprint follows,
the last state of the run is zero. A length-one closed run is therefore
impossible. A terminal full run has no closing condition.

For a single coordinate, with output marker 0<x<=1, use

    H=[[1,qx],[x,1+qx]],    lambda=1+(q+1)x.

The weighted count for a closed run is qx(H^(ell-1))_10. The eigenvalues
of H are

    lambda_+/- = 1+qx/2 +/- x sqrt(q^2+4q)/2.

They satisfy 0<=lambda_- and lambda_+<=lambda. Hence

    qx(H^(ell-1))_10 <= x lambda^(ell-1).                (1)

For a terminal run the positive-vector bound gives

    qx e_1^T H^(ell-1) 1 <= (q+1)x lambda^(ell-1).       (2)

These bounds sum the internal state histories rather than forcing them
to stay full. Genuine histories share R(V); enlarging their state sum
after applying the joint-rank bound is legitimate.

Let r be the number of closed covered runs and epsilon indicate a
terminal full run. Put D=512r and B=512epsilon. Multiplying (1)--(2)
across coordinates and runs bounds the weighted state/output sum by

    x^D ((q+1)x)^B lambda^(F-D-B).                      (3)

Multiply (3) by (gamma/q)x^(-W)q^(-R), with gamma<2, and cap at one
before fractional averaging. Three marker choices recover the scalar
envelopes already checked in the state-burst proof:

- If F<=W, choose x=1. The field count is at most
  (gamma/q) exp(2N/q) q^(F-D-R).
- If W<F<124M, choose x=N/q. The count is at most
  (gamma/q) exp(1+N/q) q^(W-D-R) N^(F-W).
- If F>=124M, choose x=1/q. The count is at most
  (gamma/q) exp(N/q) q^(W-D-R) 2^F.

Every displayed exponential correction is below e^2. The projective
factor absorbs it: gamma e^2/q<18/q<1. Thus the old low-coverage field
cost D+(11/32)(F-W)_+ is unchanged. In the high-coverage range, retaining
2^F until F reaches K gives

    D+31F/32-W >= D+40346,       124M<=F<=K.

For F>=K, replace 2^F by 2^N and remove the deterministic excess in the
shortened rank, giving D+39322. Replacing 2^F by 2^N below K would lose
the required bound and is not used.

The geometric part also transfers unchanged. Binary footprint paths
with r deaths number binom(4097,2r+1), including both terminal choices.
They meet at most 2r+1 partial coordinate-regions, so the row-coverage
variance is at most min(64,(r+1)/2). A nonzero RS word requires coverage
in at least129 regions. The same4095 exact rational branches from the
old theorem therefore remain below exp(-112), and their union is below
2^-148. No additional setup independence is assumed.

Two independent mathematical audits checked this extension. The helper
`quartet_regional_rank_probe.py --footprint-burst-cert` replays all4095
branches and256 rational one-coordinate run comparisons.

## Expansion of the holes, read backwards

Put m_V=sum_t min(|V_t|,512-|V_t|), and let K_loss=F-L. Unlike state
minority, m_V measures intermediate input coverage directly. The grouped
bound charges K_loss because those covered positions cannot carry a
nonzero state compatible with the next footprint.

Write H_t for the complement of V_t, and Gamma_t(H_t) for the union of
quartets at update t meeting H_t. Then exactly

    K_loss=sum_(t=1)^T |Gamma_t(H_t) outside H_(t-1)|.    (4)

This is the omission cost for the reverse hole path. Its initial holes
H_T are unrestricted and its final holes H_0 fill the width. Its incoming
minority is exactly m_V. Reverse deaths are the forward coordinate
deaths of coverage, D_cov=sum_(t<T)|V_t outside V_(t+1)|.

Apply the already verified quartet expansion moment to this reverse path.
Its boundary displacement is |V_T|, giving

    E_parts[q^(-K_loss)]
      <=4^|V_T| 2^(-m_V) 2^(-30D_cov).                  (5)

The ordinary forward binary-path sum is still applicable:

    Z=sum_V 2^(-30D_cov) <2^6144 exp(35/24).

At x=1, the grouped field count is at most

    (gamma/q) q^(d-K_loss) exp(2N/q).

Apply the class inequality before dropping its route/partition-dependent
indicator. Equations (5) and the path sum prove

    Pr[exists nonzero realized footprint with
       m_V>=32d+7216] <2^-76.                           (6)

This statement does not require a low-output assumption. The correction
fits the earlier slack because 35/24+1/1024<3 ln2.

There is also an output-assisted version. Put z=max(0,L-W). Choose x=1
when L<=W and x=N/q otherwise. The grouped field count is then at most

    (gamma/q) q^(d-K_loss) (N/q)^z exp(1+N/q).

Since N/q=2^-11, the same argument proves

    Pr[exists nonzero output of weight<=W with
       m_V+11max(0,L-W)>=32d+7217] <2^-76.              (7)

Here L and z depend on partitions. This is harmless: apply the class
inequality pointwise, then average. There is no union over marker choices.
The exact numerical budget is

    35/24+1+1/2048 <4*(2/3)<4 ln2,
    -31+1024+6148-7217=-76.

Independent audit agrees. The flag `--footprint-expansion-cert` in
`quartet_global_coverage_probe.py` reuses the513 exact quartet comparisons
and505 free-terminal counts, and checks4096 reverse-hole geometries.

## The remaining event

Together with the earlier rank and high-coverage classes, a remaining bad
history has genuinely partial coverage, F<53N/64, R<1009256, and

    m_V+11max(0,L-W) <32d+7217.                          (8)

For a fixed candidate eligible-row set H of size h, on its actual
classification event d=|V intersect H|-128h. Thus (8) can be retained by
markers for footprint minority, covered H-slots, and L inside the new
global transfer. Also h<=7884 follows already from R>=128h. These are
necessary conditions for the remainder, not a proof that it is empty.

Next restrict the global calculation to this remainder instead of again
paying for the classes above. The full45%/40-bit certificate stays open.
