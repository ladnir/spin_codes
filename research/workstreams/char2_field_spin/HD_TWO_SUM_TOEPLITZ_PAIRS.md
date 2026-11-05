# Post-mixer Toeplitz pairs: exact law and first local gate

The proposed g16 outer replaces its output diagonal by eight independent
invertible Toeplitz2 blocks. The small exact checks confirm a useful
spreading law. They do not establish the 35% distance target or predict
native timing. No regional proof bank, optimizer, or native job ran here.

## Sampled map and exact pair law

Keep sixteen RS[16,8] rows, the affine two-sum matrix H, and the input
diagonal D_in. At each column independently, partition the sixteen H
outputs into eight fixed pairs. Sample one uniform invertible Toeplitz2
matrix for each pair. Apply these matrices after H D_in, replacing D_out.
The proposed sixteen-column schedule uses every nonzero XOR matching on
four-bit coordinate labels, then repeats direction one.

Let q be a power of two and fix x!=0 in F_q^2. For

\[
 T=\begin{pmatrix}a&b\\c&a\end{pmatrix}
\]

sampled uniformly among invertible matrices of this form, the exact law is

\[
 \mathcal L(Tx)
 =(1-q^{-2})\operatorname{Unif}(F_q^2\setminus\{0\})
 +q^{-2}\operatorname{Unif}(\operatorname{span}(x)\setminus\{0\}).
\]

There are q^2(q-1) accepted matrices. If x_0 is nonzero, divide x and y
by x_0 and write the normalized input as (1,t). Fixing y then leaves q
possible triples on the nonzero line through x, and q-1
otherwise. Indeed, after solving the two output equations, the determinant
is y_0^2+b(y_1+t*y_0). The input (0,1) follows directly. Thus the atom is
1/[q(q-1)] on that line and 1/q^2 at every other nonzero output.
Exhaustive GF2 and GF4 tests check every input/output atom, including zero.

The input-diagonal gauge remains valid. Its common nonzero scalar can be
absorbed into all eight Toeplitz blocks. Conditional on that scalar, this
is a bijection of each independently uniform accepted-matrix set.

Relative to the current mixer, the circuit adds eight variable products
and 24 binary XORs per column. Its total is 47 rather than 39 field-product
sites, including H's eight fixed products. Prepared K20 keys grow by
1 MiB. Fixed pair wiring needs no additional payload pass. These are
circuit counts, not a measured runtime estimate.

## Support law versus amplitude law

Use masks 0,1,2,3 for a pair's labelled support. Zero stays zero. For a
nonzero incoming mask t, the exact outgoing support masses are

\[
 P_t(u)=\frac{(q-1)^{|u|}}{q^2}+\frac{\mathbf1_{u=t}}{q^2},
 \qquad u\ne0.
\]

The lazy term preserves the input amplitude ratio when t=3. Consequently,
these normalized support probabilities cannot be combined directly with
an operator that requires independent nonzero input amplitudes. A GF4
test exhibits the failure: input (1,1) assigns too much mass to output
(1,1) for the support-averaged iid law to dominate it.

The helper supplies a separate pointwise iid-NZ majorant. Replace the
lazy coefficient by (q-1)/q^2 when t=3; keep it 1/q^2 for singleton t.
Divide each resulting support mass by (q-1)^|u| to obtain its per-vector
bound. This majorant is no larger than the earlier uniform atom cap
1/[q(q-1)] on every nonzero output. Its total mass is one for singleton
input and 1+(q-2)/q^2 for double input. The exact and majorant modes are
explicitly separate in the API and artifact.

## Positive contraction with the existing H kernel

Fix a constituent set A throughout the word. For one column, let
Q_{a,b,w}[k] be the existing exact count with input support contained in A,
input weight k, and H-output support T of weight w and overlap b=|A intersect T|.
Here a=|A|. For v>=1, set

\[
 N_{a,b,w}(v)=\sum_k Q_{a,b,w}[k]v^k.
\]

For each matching pair, enumerate its four H-output support masks. Weight
each by the chosen pair transition evaluated at the two output markers.
Track only total H-output weight w and overlap b. Multiplying the eight
positive pair polynomials gives coefficients c_{w,b}; the contraction is

\[
 \sum_{w,b}c_{w,b}N_{a,b,w}(v).
\]

This is exact for support-only observables in `exact-support` mode and an
upper bound in `iid-majorant` mode. It retains labelled output markers and
the fixed matching. An independent g8 labelled-mask enumeration agrees
exactly. At unit output markers, the exact result is (1+255v)^a. Full
arity at v=1 remains uniform under the invertible pair maps.

The auxiliary `pair_mask_count` keeps only which matching pairs are
nonzero after H. Its 256 possible g16 pair masks partition the input
space. This coarser count alone does not retain the lazy branch's within-pair
support orientation; the full contraction does.

The helper does not apply the whole-word factor 255^(-8a), RS floors,
fractional powers, source counts, or inner operators. A future controller
must place those operations at their established stages. In particular,
raising a scalar contraction to alpha is not a substitute for the required
support-resolved fractional calculation.

## Local gate results

First fix the pre-Toeplitz support to coordinates 0..7 in every column.
Use the complete XOR matching schedule. This conditions on a support word;
it does not bound how often the RS/moment constraints produce that word.
The table gives log2 of the output-weight moment over the pair randomness.

| Output marker z | Independent diagonal | Exact Toeplitz support law | Safe iid majorant |
|---:|---:|---:|---:|
| 1/2 | -128.000000 | -190.924979 | -190.747774 |
| 49/50 | -3.730732 | -5.574076 | -5.395522 |

Thus the safe comparison gains 62.748 bits on this fixed word at z=1/2,
but only 1.665 bits near z=0.98. These are local moment gains, not setup
failure margins.

The second control averages uniform input on a fixed coordinate space A
for one column. It uses matching direction eight and z=1/2. Values below
include normalization by 256^a, not the RS membership normalization.

| A | a | Independent diagonal | Exact pair law | Safe iid majorant |
|---:|---:|---:|---:|---:|
| 0x0001 | 1 | -7.994397 | -7.994040 | -7.993854 |
| 0x0003 | 2 | -14.955787 | -14.954583 | -14.931756 |
| 0x000f | 4 | -15.828007 | -15.904567 | -15.860782 |
| 0x00ff | 8 | -15.905346 | -15.910007 | -15.866053 |
| 0xffff | 16 | -15.910007 | -15.910007 | -15.866053 |

The arity-eight average gain is smaller than the iid comparison cost.
Full arity's exact law is unchanged. The candidate therefore needs a
calculation retaining the quiet or pre-H-weight structure, not a claim
that every dense column improves.

## Required next gate

For a word quiet in all sixteen columns, every active RS row contributes
at least nine nonzeros. If A has size a, the XOR schedule gives

\[
 \sum_j e_j(A)={a\choose2}+e_1(A),\qquad e_1(A)\le\lfloor a/2\rfloor,
\]

where e_j counts matching edges wholly in A. Its number B of nonzero
pre-Toeplitz pairs therefore satisfies

\[
 B\ge\max\left(\lceil9a/2\rceil,
 9a-{a\choose2}-\lfloor a/2\rfloor\right).
\]

This floor applies only to the fully quiet class. A general nonquiet
word needs its actual H kernel. Neither the floor nor the kernel permits
resampling A between columns. A controller that uses per-column arity
Hölder must check that it has not discarded the schedule's benefit.

Before a full proof grid, require a support-resolved fractional control
for a quiet dense class and an independent sparse r1/r2 control. Report
the exact-support diagnostic separately from the iid-NZ-safe bound. The
whole-code acceptance threshold remains log2 failure<=-40 after all class
unions; these local numbers do not meet or test that threshold. There is
no justified native performance claim until that proof gate is useful.

## Checks and pins

Nine tests pass in 0.093 seconds. They cover exhaustive pair atoms,
pointwise iid domination, the lazy-ratio counterexample, positive labelled
contraction, marked normalization, full-arity uniformity, matching incidence,
invalid scopes, and exact zeros. The local gate completed in approximately
one second. Logs are binary64 displays of exact rational calculations,
not directed distance receipts.

From `theory`:

```text
python -B -m unittest hd_two_sum_toeplitz_pairs_tests
python -B hd_two_sum_toeplitz_pairs.py --output NEW_OUTPUT.json
```

The first gate used helper SHA-256
`1e89275e076f8b9f2c1b321e28d8b0ad22c7fff85bcc00729a6202c43953bddd`.
The artifact `hd_two_sum_toeplitz_pairs_local_gate.json` has SHA-256
`3074d3a0cfe3f863e977ebcea9c3dd34930eb0be4e6025eaafd006a23ab95238`.
Its original CLI printed a hash of the LF string before Windows converted
line endings. The file hash above is authoritative. The helper now writes
explicit LF and hashes the actual file bytes; its arithmetic is unchanged.
Existing proof providers, native code, and failed controls remain unchanged.
