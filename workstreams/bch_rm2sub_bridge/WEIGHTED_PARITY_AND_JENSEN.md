# Weighted parity and a sharper first-moment denominator

The full second-moment argument must retain the dependence between messages
and regions. This note supplies a weighted parity bound and a sharper
first-moment lower bound for the same fixed construction. Neither proves
the requested setup-failure upper bound or a contradictory lower bound.

Use the fixed BCH [256,128] code and inner t128_s15. Let T_80 be its outer
weight-80 shell, and set Q:=2620 and L:=8192. A candidate message has Q
occupied rows in T_80; every other row is zero. Its input weight is 209600.
As before, zero B-syndrome in every epoch keeps the state zero and makes
the output equal the permuted input, for every choice of multipliers.

## Polynomial weights preserve the parity mixing estimate

Fix any two occupied supports of size Q. Independently sample the row words
of the two messages from T_80. Each row has a shared coordinate permutation;
these permutations and the word pairs are independent between row positions.
Let E denote the event that both messages have even parity in all 256 regions.

A row variable contains both words and their shared permutation at that
position. For row positions occupied in just one support, it contains that
one word and its permutation. Consider variables

\[
z_a:=\sum_i w_{a,i}g_{a,i}(R_i),\qquad
w_{a,i}\ge0,\quad \sum_i w_{a,i}=1,\quad 0\le g_{a,i}\le1.
\]

Here R_i is the row variable, and each function g_{a,i} depends only on
that row. Examples include normalized region counts J_r/Q, normalized
overlap counts in the intersecting support, and averages of any of the
four coordinate-pair indicators. Zero overlap counts can instead be
represented by the identically zero function.

Let f(z)=sum_alpha c_alpha z^alpha be a real polynomial of total degree
at most d<Q. Define its coefficient norm C_f:=sum_alpha |c_alpha| in this
specified representation. It need not have nonnegative coefficients.

For T_80, the previous exact character bounds give constants

\[
a=3/8,\qquad b=p_{\rm close}+(1-p_{\rm close})23/64,
\qquad a^2\le b<1.
\]

The probability p_close is bounded using the actual shell's intersection
caps, not a random-code approximation. Set D:=2^256-2 and

\[
\epsilon_d:=D a^{Q-d}+\frac{D^2}{4}b^{Q-d}.
\]

Then the actual shared-row experiment satisfies

\[
\left|\mathbb E[f\,1_E]-2^{-510}\mathbb E f\right|
 \le 2^{-510}C_f\epsilon_d. \tag{1}
\]

To prove (1), expand a degree-d monomial into its row averages. Each term
marks at most d row positions. At every unmarked position, the ordinary
character factors remain unchanged. At a marked position, its additional
[0,1]-valued factor bounds the absolute expectation by one.

If exactly one message character is nontrivial, at least Q-d of its rows
remain unmarked. Their product has magnitude at most a^(Q-d). If both
characters are nontrivial, let q_1,q_2 be the numbers of their remaining
rows and let r be their remaining support overlap. Independence between
rows and the previous shared-row character bound give

\[
a^{q_1+q_2-2r}b^r\le b^{\min(q_1,q_2)}\le b^{Q-d}.
\]

The first inequality uses a^2<=b. The expansion weights of each monomial
sum to one, so these bounds survive averaging over marked row choices.
Sum over coefficients using C_f. Finally, Fourier inversion has four
trivial character pairs, 4D pairs with one nontrivial character, and D^2
pairs with both nontrivial. The trivial terms give the reference term in (1).

For E f>0, equation (1) yields relative error at most
(C_f/E f)epsilon_d. `certify_weighted_parity.py` checks all constants with
exact rational arithmetic. In particular,

\[
d\le512,\quad C_f/\mathbb E f\le2^{2048}
\quad\Longrightarrow\quad
\left|\frac{\mathbb E[f1_E]}{2^{-510}\mathbb E f}-1\right|<2^{-500}.
\]

The receipt gives additional degree budgets. Exact toy tests enumerate
shared row permutations and check the marked-row Fourier bound. A separate
test checks the remaining-support inequality for unequal remaining sizes.
This theorem applies to the stated polynomial class; the actual joint
kernel probability has not yet been enclosed by such a polynomial.

## Exact parity conditioning for one region count

For the first moment, only one auxiliary message is sampled. After its row
permutations, the occupied rows are independent uniform weight-80 subsets
of 256 coordinates. This follows for every selected outer row word and
does not require knowing the shell's complete spectrum.

Fix a region and write J for its input count. Then J has distribution
Bin(Q,p), with p:=5/16. Conditional on J=j and on the identities of its
j one-bits, the remaining 255 coordinates of those j rows are independent
uniform weight-79 subsets. The other Q-j occupied rows have independent
uniform weight-80 subsets of the remaining coordinates.

Let K_h^(n)(w) be the binary Krawtchouk polynomial. Exact evaluation gives

\[
\max_{1\le h\le254}\frac{|K_h^{(255)}(79)|}{\binom{255}h}=97/255,
\qquad
\max_{1\le h\le254}\frac{|K_h^{(255)}(80)|}{\binom{255}h}=19/51.
\]

Write E_1 for all 256 column parities being even. If j is odd, E_1 is
impossible. If j is even, the two trivial Fourier characters on the other
255 coordinates both contribute one. Consequently

\[
\Pr[E_1\mid J=j]=2^{-254}(1+\eta_j),\qquad
|\eta_j|\le(2^{255}-2)(97/255)^Q=:\eta<2^{-3300}. \tag{2}
\]

The bound is independent of the identities of the one-bit rows, so it also
holds conditional only on J=j. `test_conditional_column_parity.py` checks
the full conditional Fourier identity against exhaustive toy slice counts.

## Jensen lower bound for the actual zero-state probability

Let beta_j be the exact probability that a uniform weight-j region vector
passes all 64 epoch kernel checks. These values come from the retained
kernel polynomial, not iid inputs. Define the even interval
I:={64,66,...,1536} and let G:=E_1 intersect {J_r in I for every region r}.
Let b_j:=Pr[Bin(Q,5/16)=j] and set

\[
t:=\sum_{j\notin I,\ j\text{ even}}b_j,\qquad
g:=1-512(1+\eta)t.
\]

Because Q is even, the original parity argument gives Pr[E_1]>=2^-255.
Equation (2) and a union bound therefore give

\[
\Pr[G]\ge2^{-255}g. \tag{3}
\]

For j in I, beta_j>0. Define f(j):=-log(beta_j) on I and f(j):=0 outside I.
This is a nonnegative function. Let S:=sum_r f(J_r), and define

\[
c:=\sum_{j\in I}b_j[-\log\beta_j].
\]

Equation (2) gives E[S1_G]<=E[S1_{E_1}]<=256*2^-254*(1+eta)c.
Combining this inequality with (3) yields

\[
\mathbb E[S\mid G]\le\frac{512(1+\eta)c}{g}.
\]

On G, the actual conditional zero-state probability over independent
region permutations is exp(-S). Jensen's inequality thus proves

\[
\Pr[\text{zero-state in every epoch}]
 \ge2^{-255}g\exp\left(-\frac{512(1+\eta)c}{g}\right). \tag{4}
\]

This calculation does not assume independent region counts. The conditional
expectation bound uses their single marginals, and Jensen retains their
joint law. It also does not use the polynomial lemma (1); equation (2)
handles the needed one-count conditioning directly.

`zero_state_jensen_lower.py` evaluates (2)--(4) outward at 256 bits. The
512-bit replay passed. The zero-state probability lower has diagnostic
log_2 value -244964.3269344031. Multiplying by binom(8192,2620)a_80^2620,
with the retained exact shell-size lower a_80, gives an actual bad-message
expectation greater than 2^19635. The diagnostic log_2 lower is
19635.093667001463, improving the earlier exact-weight-80 lower by
41.80401552073644 bits. This sharper denominator is useful for the still
unproved second-moment ratio.

## Why a worst-case count cap is not a useful replacement

There are admissible setups with much larger bad-message clusters than
the retained expectation lower. The fixed B columns have a kernel quartet
at coordinate indices {0,4,65,69}. Its 32 additive cosets partition all
128 epoch positions into disjoint kernel quartets; every coset is checked
directly, so no unverified symmetry of B is required.

Choose identity row and region permutations. Across the 64 row blocks,
select 655 of these disjoint quartets. In each quartet, repeat one arbitrary
word from T_80 in its four rows. Every column then has support equal to a
union of kernel quartets within each epoch. All B-syndromes vanish.

This injective family has exactly 2620 occupied rows and output weight
209600. Its size is at least a_80^655>2^64299. `certify_aligned_setup_cluster.py`
checks the partition and writes the occupied row positions. Its replay passed.

Thus any universal upper bound M on this candidate count must exceed
2^64299. Even the sharper retained expectation lower, divided by this
witnessed cluster lower, is below 2^-44000. The elementary bound
Pr[Z>0]>=E Z/M cannot produce a useful conclusion from that retained lower
and a universal cap. This does not limit a better expectation calculation
or a probability-sensitive treatment of exceptional setups.

The aligned setup has positive but extremely small probability. Its existence
does not contradict a target failure probability below 2^-40. It explains
why concentration or a second moment is still necessary here.

## Relation to existing methods and next step

Rathi studies second moments of weight enumerators for regular LDPC
ensembles. The resulting concentration statements are asymptotic and require
additional conditions on polynomial solutions. They are not a finite-length
certificate for this structured BCH/RM2Sub ensemble. The paper supports the
choice of a second-moment analysis, not a transfer of its numerical bounds.
[Rathi, preprint](https://arxiv.org/pdf/cs/0512066).

Frolov gives a syndrome-entropy upper bound on minimum distance for generalized
LDPC codes. Its theorem is asymptotic and requires constituent and constant-weight
code bounds. It does not directly provide the requested finite probability
statement. No result from that theorem is used in the certificates above.
[Frolov, Theorem 1](https://arxiv.org/html/1502.06874).

Next, construct and verify polynomial bounds for the normalized joint kernel
factor, or control its exponential moments with an explicit remainder. The
weighted parity budget now permits low-degree approximations without dropping
the shared parity dependence. Atypical types and higher-order dependence
remain part of the full obligation. Complete positive coverage is unchanged.

```powershell
python -B workstreams/bch_rm2sub_bridge/certify_weighted_parity.py
python -B workstreams/bch_rm2sub_bridge/test_weighted_parity.py
python -B workstreams/bch_rm2sub_bridge/zero_state_jensen_lower.py --verify
python -B workstreams/bch_rm2sub_bridge/test_conditional_column_parity.py
python -B workstreams/bch_rm2sub_bridge/certify_aligned_setup_cluster.py
```
