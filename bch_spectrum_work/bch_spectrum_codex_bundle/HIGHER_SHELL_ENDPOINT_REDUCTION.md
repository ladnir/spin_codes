# Direct tests for the weight-40 and weight-42 caps

Updated: 2026-09-04

The current occupation-one application bound needs two additional shell caps:

\[
A_{40}(C)\le 5\cdot10^{13},\qquad A_{42}(C)\le5\cdot10^{15}.
\]

Both caps now have conditional statistical acceptance after fixed-budget
tests, full replay, and independent checks; neither is a deterministic spectrum
theorem. See RANDOM_MODEL_NARROW_CERTIFICATE.md for the completed results.
This note derives an exact counting experiment for each cap without assuming
a random BCH spectrum.
It supports both statistical testing and a deterministic counting problem.

## Code and locator representation

Use the field F=GF(256), represented by modulus 0x14d. Coordinates of an
extended binary word are indexed by F. For a support T, write
\(p_j(T)=\sum_{x\in T}x^j\) for j>0, and let p_0(T) be its parity.

Let P be the even extension of the primitive narrow-sense BCH(255,131,37)
code. Its supports satisfy p_0=p_1=...=p_36=0. The map
\(p_{37}:P\to F\) has kernel Q, the extended BCH(255,123,39) code.
The quotient has binary dimension eight, so this map is surjective.
The existing exact generator calculation in `code/bch_quotient.py` verifies
these dimensions and the new cyclotomic coset containing exponent 37.

Fix a five-dimensional binary subspace S of F and define
\(C=\{c\in P:p_{37}(c)\in S\}\). This is a [256,128] member of the
previously studied sandwich Q<=C<=P. One may use S={0,...,31} in this field
representation. All such S give the same weight enumerator: scaling coordinates
by a nonzero a preserves weight and multiplies p_37 by a^37. Since
gcd(37,255)=1, each nonzero syndrome has the same number of words of each weight.

Translations x->x+b preserve every p_j for j<=37 on P, except that the
already-zero lower syndromes remain zero. Indeed, expand (x+b)^j and use
p_0=...=p_36=0. Thus translations preserve C and act transitively on coordinates.
Writing C' for C punctured at coordinate zero, double-counting incidences gives

\[
A_{w-1}(C')=\frac{w}{256}A_w(C) \quad\text{for even }w.
\tag{1}
\]

Here a punctured word of odd weight w-1 has a unique even extension of weight w.
Transitivity suffices for (1); a two-design hypothesis is unnecessary.

Fix w in {40,42}. Set

\[
r=(w-38)/2,\quad d=w-1=37+2r,\quad m=18+r,\quad n=18+2r=w-20.
\]

For a weight-d support T in C', its locator polynomial is
\(\Lambda(z)=\prod_{x\in T}(1+xz)\).
Newton identities and p_1=...=p_36=0 give

\[
\Lambda(z)=g(z^2)+z^{37}h(z^2),\qquad
\deg g\le m,\quad g(0)=1,\quad \deg h=r,\quad h(0)=\sigma=p_{37}(T)\in S.
\tag{2}
\]

The degree of h is exactly r because Lambda has odd degree d. For nonzero
u=z^2, we have z=u^128 and z^37=u^146. Consequently its d distinct roots
correspond to the d agreements

\[
g(u)=u^{146}h(u),\qquad u\in F^*.
\tag{3}
\]

Conversely, a pair satisfying the degree bounds and constant conditions in (2)
with exactly d agreements defines a degree-d locator with d distinct nonzero
roots. The inverse roots form a support with p_1=...=p_36=0 and p_37=sigma.
This follows by applying Newton identities in reverse; even power sums are
squares of lower power sums. Thus each such endpoint is a word in C'.

## Reconstruction cannot miss a genuine endpoint

Fix sigma and an n-element subset B of F*. Solve (3) on B for the n unknowns
g_1,...,g_m,h_1,...,h_r. The constants g_0=1 and h_0=sigma are fixed.

**Lemma.** If B is contained in the agreement set of a genuine weight-d locator
in (2), this n-by-n linear system is nonsingular.

**Proof.** Fix such a locator and a homogeneous solution (delta g, delta h).
Both perturbations have constant coefficient zero. The polynomial

\[
F(u)=h(u)\,\delta g(u)+\delta h(u)\,g(u)
\]

vanishes on B and at zero. Its degree is at most r+m=n. These n+1 distinct
zeros force F=0.

The polynomials g and h are coprime. Otherwise a common nonconstant factor
a(u) would make a(z^2) divide Lambda. Over the perfect field F, a(z^2) is a
square of a nonconstant polynomial, contradicting the distinct roots of Lambda.
It follows from F=0 that h divides delta h. Since deg h=r and deg delta h<=r,
we have delta h=lambda h for some scalar lambda, and then delta g=lambda g.
The equation delta g(0)=0 and g(0)=1 force lambda=0. The homogeneous kernel
is trivial. This argument also covers sigma=0. QED.

An algorithm may therefore treat every singular system as a non-hit. It must
not resample until a nonsingular system appears: that would change the law of B.
On a nonsingular system, reconstruct g,h and count all agreements on F*.
Return a hit exactly when that count is d.

## Exact incidence identity and sampling plans

In one trial, sample sigma from S and independently sample B from the set of
all n-element subsets of F*. Execute the reconstruction algorithm above.
Every weight-d word contributes exactly binom(d,n) distinct (sigma,B) hits.
The lemma excludes collisions between these contributions. Equation (1) gives

\[
\Pr[\mathrm{hit}]
=\frac{A_{w-1}(C')\binom{w-1}{w-20}}{32\binom{255}{w-20}}
=A_w(C)\frac{w\binom{w-1}{w-20}}{8192\binom{255}{w-20}}.
\tag{4}
\]

This is an exact identity for a fixed code, not a random-code heuristic.
It does not itself bound the number of hits.

For an integer cap U, a violating shell has A_w>=U+1. Substitute U+1 into
(4) to obtain p_bad. Fix a sample count N before sampling, and accept the cap
only if all N independent trials have zero hits. Under a violating shell,
the probability of acceptance is at most (1-p_bad)^N.

| Weight | Tested cap U | p_bad, approximate | Fixed trials for error <=2^-41 |
|---:|---:|---:|---:|
| 40 | 50,000,000,000,000 | 6.5087514681e-8 | 436,628,033 |
| 42 | 5,000,000,000,000,000 | 2.0382014963e-7 | 139,431,904 |

The trial budgets are certified with rational enclosures of logarithms.
Using U+1 is conservative; no additional divisibility restriction is assumed.
Both full experiments completed with zero hits and matching independent replay.
All 1,150 retained Python audit records passed. Fixed-seed structural checks
below are separate from these experiments.

These error bounds concern erroneous acceptance of a fixed shell-count claim.
They are distinct from the probability of a bad SPIN setup. Combined with the
existing weight-38 test at error 2^-40, two tests at error 2^-41 each give a
familywise union bound of 2^-39 under ideal IID randomness. The joint cap claim,
accepted only if all three tests pass, has false-accept error at most 2^-40 by
the intersection-union argument in RANDOM_MODEL_NARROW_CERTIFICATE.md.
The actual implementation used retained CNG pseudorandom bytes, requiring the same
randomness qualification as the weight-38 certificate. No deterministic upper
bound or posterior probability of correctness follows from zero hits.

## Small reconstruction systems and checks

For an efficient implementation, let P_B(u)=product_(b in B)(u+b), and define
R_j(u) as the remainder of u^(146+j) modulo P_B for j=0,...,r. Then

\[
g=\sigma R_0+\sum_{j=1}^r h_jR_j.
\]

Because deg g<n, the r constraints g_0=1 and g_(m+1)=...=g_(n-1)=0 suffice.
Weight 40 requires one scalar solve; weight 42 requires a 2-by-2 solve.
R_0 comes from interpolation on B; each later remainder comes from one multiply
by u followed by monic reduction. Thus a full 20-by-20 or 22-by-22 solve is
unnecessary in the hot path.

Reproduce the exact checks with:

    python -B code/verify_higher_shell_endpoint_reduction.py

The receipt is `generated/higher_shell_endpoint_reduction_checks.json`.
It includes the following checks:

- Complete enumeration over GF(16), with base locator exponent 3, at degrees
  5 and 7: all 7,280 and 48,048 (sigma,B) queries match independently enumerated
  supports. Each of the 168 and 435 endpoints is counted exactly 10 and 21 times.
- Independent extended-support counts 448 and 870 verify the puncturing factor.
- Thirty-two retained BCH witnesses cover both weights and both the zero and
  nonzero syndrome cases. All 2,048 selected endpoint bases recover the original
  locator with both reduced and full linear systems.
- Another 512 fixed-seed arbitrary queries cross-check the two linear solvers,
  including their singular-system decisions.

The general correctness claim rests on the proof above, not on the examples.
The examples test field conventions, the implementation, and the incidence law.

For the deterministic route, (4) replaces the unknown spectrum by a finite
count of successful queries. A uniform upper bound on that count would prove
the desired caps. The nonsingularity lemma alone supplies no such sharp count.
The fast sampler and independent replay passed structural tests, small
preflights, and the separate full fixed-budget experiments. See
HIGHER_SHELL_PREFLIGHT.md and generated/higher_endpoint_certificate_20260904/certificate.json.
