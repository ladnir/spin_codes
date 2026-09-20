# Asymptotic scaling targets for SPIN codes

## Scope

This document specifies parameter schedules and theorem targets.  It does not
promote finite diagnostics to asymptotic theorems.  Throughout, \(N\) is the
intermediate and output length, \(K_N\) is the message length, and
\(R_N:=K_N/N\).

For a target \(D_N=\lfloor\delta N\rfloor+1\), let

\[
  \mu_N(\delta):=\mathbb E[Z_{D_N}].
\]

The setup-failure probability is at most \(\mu_N(\delta)\).  Define the
first-moment margin

\[
  \Lambda_N(\delta):=-\log_2\mu_N(\delta)
\]

when \(\mu_N(\delta)>0\).  A linear-distance theorem needs only
\(\mu_N(\delta)=o(1)\).  A positive linear margin is the stronger statement
\(\Lambda_N(\delta)=\Omega(N)\).

## A reusable block schedule

Let the outer layer contain

\[
  L_N:=N/B_N
\]

local blocks, each mapping \(K_{B_N}\) bits to \(B_N\) bits.  Use admissible
lengths with \(B_N\mid N\).  The exponential-spectrum schedule is

\[
  B_N=\lceil c_B\log_2 N\rceil_{\mathcal B},
  \qquad K_{B_N}=R_{\mathrm{out}}B_N+O(1),
\]

where \(\lceil\cdot\rceil_{\mathcal B}\) rounds to an allowed block size.
Thus \(L_N=\Theta(N/\log N)\) and \(R_N\to R_{\mathrm{out}}\).

The logarithmic block size applies when the weighted contribution of one
outer block is exponentially small in \(B_N\).  A nonzero message supported
on one outer block creates only \(\Theta(B_N)\) shuffled nonzero positions.
An inner late-placement event can then cost only

\[
  2^{-\Theta(B_N)}=N^{-\Theta(1)}.
\]

There are \(\Theta(N/B_N)\) possible active blocks.  Constant \(B_N\) cannot
make this union bound vanish under a fixed per-block contraction.  Taking
\(B_N=\Theta(\log N)\) gives a tunable polynomial suppression.

More generally, let \(M_B\) upper-bound the complete one-active-block
contribution.  The required condition is

\[
  \frac{N}{B_N}M_{B_N}=o(1).
\]

If \(M_B\le B^{-p+o(1)}\), a valid alternative is

\[
  B_N=N^{\beta+o(1)},
  \qquad \beta>\frac1{p+1}.
\]

Thus the theory permits a sublinear block when logarithmic blocks are not
supported by the outer spectrum.

## Sparse and many-block margins

Let \(Q(x)\) be the number of active outer blocks in a nonzero message \(x\).
For a fixed \(q\), suppose the total first-moment contribution satisfies

\[
  \mu_{N,q}(\delta)
  \le \binom{L_N}{q}\,N^a\,2^{-\lambda B_Nq}
\]

for constants \(a\ge0\) and \(\lambda>0\).

For \(q=1\),

\[
  \mu_{N,1}(\delta)
  \le N^{a+1-\lambda c_B+o(1)}.
\]

Hence \(c_B>(a+1)/\lambda\) gives

\[
  \Lambda_{N,1}(\delta)=\Theta(\log N).
\]

This one-block class usually determines the global margin.  Increasing the
number of active blocks does not change that global conclusion unless the
one-block class is excluded by the outer construction.

If \(q=\alpha L_N\) for a constant \(\alpha>0\), the same local estimate gives
an exponent of order \(N\), because \(B_Nq=\Theta(N)\).  A sharper type or
large-deviation analysis should therefore target

\[
  \mu_{N,q}(\delta)\le2^{-\Omega(N)}
\]

for many-block classes.  Intermediate \(q\) gives an expected margin of order
\(qB_N\), minus the block-selection entropy.

The correct summary is:

- one-block and bounded-block classes: \(\Theta(\log N)\) margin;
- \(q\to\infty\) sparse classes: typically \(\Theta(q\log N)\), after the
  support-count cost;
- a positive fraction of active blocks: \(\Theta(N)\) margin.

A finite plot showing a larger bulk margin is consistent with this scaling,
but it does not prove it.

## Accumulator SPIN schedule

### Construction target

Use independent random binary outer blocks of length

\[
  B_N=\lceil c_A\log_2 N\rceil_{\mathcal B}
\]

and fixed rate \(R_A\in(0,1)\).  Apply one uniform permutation of all \(N\)
outer coordinates, followed by the length-\(N\) accumulator.

The accumulator has no growing state parameter.  Its exact transfer law is

\[
  p_h^{\mathrm{Acc}}(D_N)
  =\frac{1}{\binom Nh}
    \sum_{j\le D_N}
    \binom{j-1}{\lceil h/2\rceil-1}
    \binom{N-j}{h-\lceil h/2\rceil}.
\]

For \(D_N=\lfloor\delta N\rfloor\) and \(\delta<1/2\), the exact enumerator is
the upper tail of a hypergeometric random variable.  Maclaurin's inequality
and a Chernoff bound give, uniformly for every \(1\le h\le N\),

\[
  p_h^{\mathrm{Acc}}(D_N)
  \le \rho_A^h,
  \qquad \rho_A:=2\sqrt{\delta(1-\delta)}.
\]

This improves the earlier contraction \(\sqrt{4e\delta}^{\,h}\).

### Theorem and constants

For a fixed local outer ensemble, define its nonzero weight generating mass

\[
  G_B(z):=
  \mathbb E\!\left[\sum_{c\in C_B\setminus\{0\}}z^{\operatorname{wt}(c)}\right].
\]

For uniform random rate-\(R_A\) injections, the exact local law gives

\[
  G_B(\rho_A)\le2^{1-\lambda_AB},
  \qquad
  \lambda_A
  :=1-R_A-\log_2(1+\rho_A).
\]

Independence of the local outer blocks gives

\[
  \mu_N(\delta)
  \le (1+G_{B_N}(\rho_A))^{L_N}-1.
\]

It follows that \(\mu_N(\delta)=o(1)\) whenever

\[
  \lambda_A>0,
  \qquad
  c_A\lambda_A\ge1.
\]

The resulting family has rate approaching \(R_A\), positive relative
distance \(\delta\), and a global \(\Theta(\log N)\) margin unless a stronger
sparse-message property is proved.

At rate one half, every \(c_A>2\) gives a positive distance.  The exact
tradeoff is

\[
  \delta\le\delta_{c_A}
  :=\frac{1-\sqrt{1-(2^{1/2-1/c_A}-1)^2}}2.
\]

Thus \(c_A=3\) gives \(\delta\le0.0037634\), \(c_A=8\) gives
\(\delta\le0.0225363\), and \(c_A=17\) gives
\(\delta\le0.0330838\).  The limiting distance as \(c_A\to\infty\) is
\(0.0449101\).  ACCUMULATOR_SPIN_ASYMPTOTIC.md gives the complete ensemble,
proof, admissible-length rule, and arbitrary-length extension.

A pure \(L\)-by-\(B\) bit transpose cannot replace the uniform permutation at
rate one half.  Two adjacent random half-rate block subspaces intersect
nontrivially with probability bounded away from zero.  A common word in such
a pair produces accumulator output weight at most \(B\).  Across
\(\lfloor L/2\rfloor\) independent pairs, this gives

\[
  \Pr[d_{\min}>B]
  \le(2/3)^{\lfloor L/2\rfloor}.
\]

Thus every pure-transpose schedule with \(B_N=o(N)\) has vanishing relative
distance in probability.  Independent region shuffles break the persistent
adjacency used by this obstruction and require a separate column-profile
transfer theorem.  BIT_TRANSPOSE_ACCUMULATOR_OBSTRUCTION.md gives the exact
scope and proof.

For transpose plus independent region shuffles, that transfer is now explicit.
Conditioned on \(q\) active blocks, one region is a mixture over
\(a\sim\operatorname{Bin}(q,1/2)\) of two-state accumulator kernels.  Raising
the averaged kernel to the \(B_N\)-th power gives a finite first-moment bound.
The fixed-\(q\) continuum saddle at rate one half and distance \(0.02\)
requires constants increasing from \(3.2835\) at \(q=1\) to \(6.9265\) at
\(q=256\).  The values converge numerically toward the uniform-interleaver
constant \(6.9514\).  REGION_SHUFFLED_TRANSPOSE_ACCUMULATOR.md proves the
finite transfer and records the remaining uniformity gap for
\(q\to\infty\).

## Random SPIN schedule

### Construction target

Use independent random outer blocks with

\[
  B_N=\lceil c_R\log_2N\rceil_{\mathcal B}.
\]

Apply one uniform permutation on all \(N\) outer coordinates.  Use a random
time-varying recursive inner with step or memory parameters

\[
  t_N=\Theta(\log N),
  \qquad s_N=\Theta(\log N),
\]

unless a proof shows that one parameter can remain constant.  For the scalar
dense recursion in `innerDenseScalar.tex`, the relevant schedule is

\[
  m_N=\lceil\gamma\log_2N\rceil,
\]

and its encoder performs \(\Theta(N\log N)\) bit operations.

### Theorem target

Prove three complementary bounds for some target \(\delta_R>0\):

1. **Tiny weights.** For \(1\le h\le h_0\), where
   \(h_0=\Theta(\log N)\), prove explicit transfer bounds
   \[
     p_h^{\mathrm{in}}(D_N)\le T_{N,h}
   \]
   and show directly that
   \(\sum_{h\le h_0}A_h^{\mathrm{out}}T_{N,h}=o(1)\).
2. **Intermediate contraction.** There are \(\rho_R<1\) and
   \(\eta_0>0\) such that
   \[
     p_h^{\mathrm{in}}(D_N)\le N^{O(1)}\rho_R^h
     \quad(h_0<h<\eta_0N),
   \]
   together with an outer generating bound that makes the sum over this
   range vanish.
3. **Linear-weight gap.** For every \(\omega\) in the remaining compact
   weight range,
   \[
     A_{\lfloor\omega N\rfloor}^{\mathrm{out}}
     \le N^{O(1)}2^{N\Phi_R(\omega)},
     \qquad
     p_{\lfloor\omega N\rfloor}^{\mathrm{in}}(D_N)
     \le N^{O(1)}2^{-NE_R(\omega)},
   \]
   and
   \[
     \sup_\omega(\Phi_R(\omega)-E_R(\omega))<0.
   \]

These bounds give a \(\Theta(\log N)\) global margin and an
\(\Theta(N)\) margin on linear-weight classes.

### Current status

For uniform random rate-half block injections, the outer analysis is exact:

\[
  G_B(z)
  =
  (2^{B/2}-1)\frac{(1+z)^B-1}{2^B-1}.
\]

Thus \(B_N=\Theta(\log N)\) closes the outer sparse bound for every
\(z<\sqrt2-1\).  RANDOM_SPIN_PROOF_AUDIT.md gives the probability space,
proof, and linear-weight exponent.

The non-wrapping convolution has an exact bivariate transfer matrix.  Its
limiting Perron root gives a closed convolution exponent.  The parameterized
sparse proof is valid whenever

\[
  \gamma\log_2(1/b_*)>1,
  \qquad
  c_R\log_2\!\left(\frac{\sqrt2}{1+\rho_*}\right)
  >2\gamma+\frac32,
\]

with \(b_*\) and \(\rho_*\) defined in RANDOM_SPIN_PROOF_AUDIT.md.  The
certified choice

\[
  \delta_R=0.11002,
  \quad
  \eta_*=1/5000,
  \quad
  \gamma=51/50,
  \quad
  c_R=17
\]

gives strict margins greater than \(0.0186697\) and \(0.0576243\),
respectively.  The combined linear exponent is strictly negative for every
fixed normalized input weight \(\eta>0\); its tightest certified margin is
\(H_2(0.11002)-1/2<-2.3718\times10^{-5}\).  Consequently the proof-model
ensemble has rate one half, relative distance greater than \(0.11002\), and
setup failure probability \(o(1)\) for
\(B_N\ge17\log_2N\) and
\(m_N=\lceil(51/50)\log_2N\rceil\).  The earlier
\((\delta_R,c_R,\gamma)=(0.09,42,2)\) theorem remains as a more conservative
baseline.

The random dense source material still provides useful finite context, but
its compiled \(0.109\) integration is not a completed proof.  It uses a
systematic random sliding outer, reverses an inequality in the
constant-geometric corollary, and identifies input runs with distinct
live-state episodes.  The exact transfer matrix replaces those steps.

The legacy packet-transpose/random-step construction reports a complete
outward occupation certificate at message length \(2^{20}\), output length
\(2^{21}\), \(B=1024\), packet width 4, state size 15, and distance 188744.
Its reported first-moment margin is about 195.2839 bits.  The supporting
packing lemma is conditional on replacing each active outer word by
independent uniform bits and charging the stated conditioning factor.  Its
factored packet interleaver is not the uniform global interleaver in the
Random SPIN definition.  This certificate is strong finite evidence for the
random-convolution transfer mechanism, but it does not discharge the named
Random SPIN theorem target.

The BCH/bit-shuffle RandomStepConv result near relative distance 0.11002 uses
a modeled outer spectrum, nearest-binary64 optimization, and incomplete
parity tracking in the middle and dense occupations.  Its reported 40.9510-bit
margin is diagnostic.

## Structured SPIN schedule

The frozen instance

\[
  (B,t,s)=(256,128,19),
  \qquad (K,N)=(2^{20},2^{21})
\]

is a finite construction.  Repeating its fixed local parameters does not by
itself define an asymptotic family with vanishing failure probability.

### Conditional scaling target

The strongest target chooses

\[
  B_N=\lceil c_S\log_2N\rceil_{\mathcal B}.
\]

This schedule requires exponential sparse suppression.  A BA outer with
polynomial weighted spectrum instead uses

\[
  B_N=N^{\beta+o(1)}
\]

with \(\beta\) determined by the sparse exponent.  For the repeated EBCH128
constituent followed by two accumulators, the proved boundary path has
exponent 16.  A matching upper bound would permit every \(\beta>1/17\).

The random-outer baseline gives a sufficient structured-outer test.  Define

\[
  g_B:=
  \max_{h:\,\overline A_{B,h}>0}
  \log_2\frac{\overline A_{B,h}}
  {(2^{K_B}-1)\binom Bh/(2^B-1)}.
\]

After independent local coordinate permutations, every nonnegative
region-profile transfer at occupation \(q\) costs at most \(2^{g_Bq}\)
relative to the random outer paired with the same structured inner.  Hence
\(g_B=o(B)\) preserves the limiting random-outer sparse exponent.  If
\(g_B=\epsilon B+o(B)\), the exponent loses \(\epsilon\).  This comparison
does not transfer the random-convolution theorem to RM2Sub; the inner kernel
must be held fixed.

Let \(t_N\) and \(s_N\) be the smallest schedules for which the structured
inner-transfer hypotheses in `STRUCTURED_SPIN_THEOREM_TARGET.md` can be proved.
The theory workstream does not assume their values.  Candidate regimes are:

- \(t_N,s_N=O(1)\), if uniform transfer bounds persist as \(B_N\) grows;
- \(t_N=\Theta(\log N)\) and \(s_N=O(1)\);
- \(t_N,s_N=\Theta(\log N)\).

The linear-time audit determines the cost consequences.  A proof should not
select the cheapest regime before establishing its transfer bound.

Under the conditional sparse-profile and bulk-exponent hypotheses stated in
the Structured SPIN target, choose \(c_S\) above the explicit sparse threshold.
Then

\[
  \Pr[d_{\min}<\delta_SN]=N^{-\Omega(1)}+2^{-\Omega(N)}=o(1).
\]

The global first-moment margin is \(\Theta(\log N)\), while the bulk margin is
\(\Theta(N)\).

### Random-outer RM2Sub intermediate family

The first structured-inner proof model keeps the audited
\(t=128,s=19\) RM2Sub constituent fixed.  It allows the even outer length
\(B_N\) to grow and requires \(128\mid L_N\).  Each local outer is an
independent uniform rate-half injection.  This choice preserves the inner
recurrence and avoids introducing a family of RM2Sub constituents.

For one active block, an active transposed region contains one uniformly
placed impulse. The finite diagnostic uses two-state epoch envelopes, but
that envelope pays a punctured-state factor once per zero epoch and is not
asymptotically iterable with fixed \(s\). The asymptotic proof instead uses
the exact \(2^s\)-state epoch kernel and reduces it after taking the
long-region limit.

Set \(z=e^{-\theta/L_N}\), \(M=2^{19}-1\), and
\(p=2^{18}/M\). Define

\[
  \gamma=p\theta,\qquad
  a=e^{-\gamma},\qquad
  f=\frac{1-e^{-\gamma}}{\gamma}.
\]

The limiting inactive and one-impulse region matrices are

\[
  K_0=
  \begin{pmatrix}1&0\\0&a\end{pmatrix},
  \qquad
  K_1=
  \begin{pmatrix}
    0&f\\
    f/M&(1-1/M)a
  \end{pmatrix}.
\]

The Perron root of \(K_0+K_1\) gives the exact occupation-one exponent.
At \(\delta=0.11002\), binary64 optimization gives
\[
  \eta_{19}=0.2787705115\ldots,
  \qquad
  c_{\mathrm{threshold}}=3.587179987\ldots.
\]
An exact rational single-tilt inequality certifies \(c=18/5\), with
\[
  \Pr[Z_{\lfloor0.11002N\rfloor,1}>0]
  \le N^{-0.0026858\ldots+o(1)}.
\]
This closes occupation one only. It does not establish minimum distance for
the complete code.

Occupation two is also closed. If \(A(\gamma)\) and \(H(\gamma)\) denote the
two-impulse order-statistic averages, its limiting matrix is

\[
  K_2=
  \begin{pmatrix}
    A/M&(1-1/M)H\\
    (1-1/M)H/M&H/M+(1-1/M)^2a
  \end{pmatrix}.
\]

The exact same-epoch collision probability is \(127/(L_N-1)\), so it does
not change the continuum exponent. At \(\delta=0.11002\), binary64
optimization gives the occupation-two threshold
\[
  c>3.587809926\ldots.
\]
An exact rational inequality again certifies \(c=18/5\). Under that schedule,
the occupation-two bad-event probability is at most
\[
  N^{-0.0049419\ldots+o(1)}.
\]
Thus \(c=18/5\) closes occupations one and two. Occupation two has the
slightly larger optimized threshold.

For every fixed occupation \(Q\), let \(K_a(\theta)\) be the ordered
two-state reset integral for a region with \(a\) impulses. The complete
region matrix is
\[
  T_Q(\theta)=\sum_{a=0}^Q\binom Qa K_a(\theta).
\]
RM2SUB_FIXED_OCCUPATION_CONTINUUM.md proves
\[
  \log_2\mathbb E Z_{\lfloor\delta N\rfloor,Q}
  \le Q\log_2L-\eta_Q(\delta)B+o_Q(B)
\]
for every fixed \(Q\). It also gives an exact coefficient formula using
Dirichlet spacings and beta transforms.

The resulting common-tilt threshold is not constant in the current
diagnostic. At \(\delta=0.11\), it rises from \(3.58666\) at \(Q=1\) to
\(3.60010\) at \(Q=16\), \(4.10383\) at \(Q=128\), and \(4.55763\) at
\(Q=512\). Thus the exact \(c=18/5\) certificates for \(Q=1,2\) do not
extend automatically to all fixed occupations. This is a limitation of the
current sufficient condition, not a distance counterexample.

A different common tilt gives a rigorous uniform constant. The beta identity

\[
  \mathbb E[(1-X+X/y)^{-(Q+1)}]=y^k,
  \qquad X\sim\operatorname{Beta}(k,Q+1-k),
\]

converts the coefficient formula into a Perron bound. Choosing
\(y=1/2\) and \(\tau=11/10\) proves, for every fixed \(Q\ge3\),

\[
  \frac{\eta_Q(0.11)}Q
  >0.11604381273505275.
\]

Exact rational arithmetic therefore certifies

\[
  B_N=9\log_2N+O(1),
  \qquad
  \Pr[Z_{\lfloor0.11N\rfloor,Q}>0]
  \le N^{-0.04439431461547Q+o_Q(1)}.
\]

Together with the sharper receipts for \(Q=1,2\), one logarithmic schedule
closes every fixed occupation. RM2SUB_UNIFORM_FIXED_OCCUPATION.md gives the
proof and exact receipt.

The finite Bernoulli lift in RM2SUB_DENSE_OCCUPATION.md closes every growing
occupation. For \(Q/L\le10^{-4}\), an exact four-state Bernstein certificate
gives at least \(0.11715QB\) of natural-log exponent before the support
factor. Since \(B=9\log_2N+O(1)\), the support costs at most
\((\ln2/9+o(1))QB\), leaving a strict margin. For
\(10^{-4}\le Q/L\le1\), 31 outward interval boxes give a uniform negative
exponent on the \(N\) scale.

For the finite reduction, if \(W_0(z),W_1(z)\) are valid epoch envelopes,
define

\[
  R_0(z)=W_0(z)^{L_N/128},
  \qquad
  R_1(z)=\frac{128}{L_N}
  \sum_r W_0(z)^rW_1(z)W_0(z)^{L_N/128-1-r}.
\]

Coefficient extraction from \((R_0+uR_1)^{B_N}\) gives the finite
one-active random-outer reduction. The exact full-state continuum proves a
bound

\[
  M_{B_N,L_N}(\lfloor\delta N\rfloor)
  \le N^a2^{-\lambda_1B_N}
\]

with \(\lambda_1>0\). RM2SUB_ONE_ACTIVE_CONTINUUM.md gives the proof and
constant. RANDOM_OUTER_RM2SUB_RAMP.md gives the probability space, finite
theorem, admissible family, and staged extension to larger occupations.

The all-active endpoint is also exact. Relaxing every active random row to a
uniform binary row makes the complete inner input uniform, and inner
bijectivity preserves that law. Its rate-half exponent is
\(H_2(\delta)-1/2+o(1)\). The Bernoulli transfer now supplies the uniform
interpolation to this endpoint.

### Current status

The random-outer RM2Sub family now has a complete asymptotic distance
certificate at rate one half and relative distance \(0.11\). Its outer block
length is \(9\log_2N+O(1)\), and direct random-outer multiplication costs
\(O(N\log N)\).

The one-sampled Golay--BA-3/RM2Sub-S19 family now has a complete structured
certificate at rate one half and relative distance \(0.11\). Its schedule is

\[
  B_N=\frac{39}{4}\log_2N+O(1).
\]

A certified concave BA spectrum majorant and a 621-box joint RM2Sub interval
cover every positive occupation. Exact likelihood transfer covers growing
sparse occupation and the first two fixed occupations. A weight-coupled
coefficient transfer covers every fixed \(Q\ge3\). The family retains linear
ordinary and transposed encoding work.

This result does not certify the frozen ParityFanout-31x33 outer. Its modeled
\(B=256\) spectrum still has no asymptotic force. The companion
distance-\(0.101\), \(B=\Theta((\ln N)^2)\) theorem remains a conservative
reference for the same Golay--BA outer and RM2Sub inner.

## Parameter selection at a requested finite length

A selector should enumerate only admissible tuples

\[
  (B,t,s,\mathsf{route},D)
\]

and evaluate the finite first-moment sum for each tuple.  It should return the
tuple, the certified bound, and the wrapper loss.  The asymptotic schedules
guide this search; they do not replace the finite calculation.

## Proof obligations common to all three families

1. Give a construction manifest for every admissible \(N\).
2. Define every setup distribution and independence relation.
3. Prove the sparse classes uniformly, including one active block.
4. Prove every unfinished middle- and linear-weight partition with no
   uncovered transition range.  The Random SPIN proof-model partition is
   complete; this obligation remains for the other proposed families.
5. Convert numerical maximization over continuous parameters to interval or
   exact arithmetic.
6. Apply the arbitrary-length wrapper only after an admissible instance has a
   certified \((D,\varepsilon)\) pair.
