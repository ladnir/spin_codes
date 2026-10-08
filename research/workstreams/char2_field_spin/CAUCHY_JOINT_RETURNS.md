# Joint returns with two or three packet inputs

The Cauchy diagonal-refresh inner admits a stronger joint return bound than charging only the remaining input diagonals. The state diagonals also contribute independent pivots. Their contribution depends on the incoming state support; the injection matrix is not superregular.

This note gives two separate interfaces: an exact fixed-word return/output law and a conditional rank bound for existential witnesses. Neither interface changes the encoder or establishes a whole-code distance certificate.

## Packet and reveal order

Use the [Cauchy residue-coordinate maps](CAUCHY_DIAGONAL_PACKET_COST.md), with width sixteen, state dimension four, \(q=256\), and \(\eta=2\). Fix the incoming state \(S\), the input support \(I\), and its nonzero input values before revealing any current packet diagonal. Write \(w=|I|\), and let \(X_I\) include the fresh input diagonals. For fixed GF256 input values, its coordinates are independent samples from \(\mathbb F_{256}^*\).

The packet is

\[
 R=D_sS,\qquad Y=X+\eta AR,\qquad S'=R+CX.
\]

The four entries of \(D_s\) are independent nonzero GF256 scalars. Let \(H=\operatorname{supp}(S)\) and \(h=|H|>0\). A return satisfies \(D_sS=C_I X_I\). Therefore

\[
 \operatorname{supp}(C_I X_I)=H.
\]

For each such \(X_I\), the return equation fixes exactly one value of each of the \(h\) active state diagonals. The other state diagonals are unrestricted. On the return, define

\[
 B_I:=E_I+\eta AC_I,\qquad Y=B_I X_I,
 \qquad (B_I)_{j,i}=\mathbf1_{j=i}+\eta(u_j+u_i)^3.
\]

Here \(E_I\) inserts the \(w\) active inputs into their packet coordinates. All formulas use the actual labels in each of the four phases.

## Exact fixed-word moment

Let \(\mathcal P_{I,H}\) contain the projective vectors \([x]\in\mathbb P^{w-1}(\mathbb F_{256})\) whose input coordinates are all nonzero and satisfy \(\operatorname{supp}(C_Ix)=H\). For \(0<y\le1\), the exact conditional law is

\[
 \boxed{
 \mathbb E[y^{\operatorname{wt}(Y)}\mathbf1_{S'=0}]
 =255^{-(w+h-1)}
   \sum_{[x]\in\mathcal P_{I,H}}y^{\operatorname{wt}(B_Ix)}.
 }
\]

Each projective vector represents 255 input assignments, with the same supports and output weight. This explains the denominator. The law depends on the fixed incoming state only through \(H\).

This formula consumes the input and state diagonals jointly. It can replace another whole-packet return bound, or be minimized with another valid bound. It cannot be multiplied into a bound that has already conditioned on these same diagonals.

### Exact enumeration without scanning all vectors

For \(w=2\), every nonzero linear functional vanishes at one projective point. Group proportional input, state, and output functionals, then count those points and their complement.

For \(w=3\), the functionals define at most 23 distinct projective lines: three input forms, four state forms, and sixteen output forms. Count three disjoint strata:

1. intersections of at least two distinct lines;
2. points lying on exactly one line;
3. points lying on no line.

There are at most \(\binom{23}{2}\) intersection points. A line contains \(q+1\) points, and the plane contains \(q^2+q+1\). These exact counts produce the complete histogram, including zero input exclusions and the exact state support. The method does not sample field vectors.

## Actual injection exceptions

The enumeration covers all four phases and every input placement: 480 pairs and 2,240 triples. Every two- or three-input return has output weight at least thirteen. The following exceptional supports are real.

| Input weight | Exceptional incoming support size | Number of \((\text{phase},I,H)\) exceptions | Return probability | Output weight |
|---:|---:|---:|---:|---:|
| 2 | 2 | 16 | \(255^{-3}\) | exactly 16 |
| 3 | 1 | 16 | \(255^{-3}\) | exactly 16 |

Each exceptional class has eight cases in phase zero and eight in phase two. Each supported \(H\) has exactly two exceptional input placements in its phase. Phases one and three have none. Every other feasible pair \((I,H)\) satisfies

\[
 h+\operatorname{rank}(C_{H^c,I})\ge4.
\]

For example, phase-zero input indices \(I=(0,6)\) have labels \((0,224)\). The projective input \(x=(1,171)\) gives

\[
 C_Ix=(82,0,0,248),\qquad \operatorname{wt}(B_Ix)=16.
\]

For phase-zero indices \((0,9,15)\), the input \((1,75,24)\) gives state target \((82,0,0,0)\), again with output weight sixteen. Thus a uniform \(255^{-4}\) return-probability claim for arbitrary fixed placements is false.

If the input placement is uniform conditional on \(w\), independently of the fixed incoming state, averaging precedes the maximum over \(H\). The following simpler fixed-word bounds then hold:

\[
 \overline R_2(y)\le255^{-4}\max\{y^{13},\tfrac{17}{4}y^{16}\},
 \qquad
 \overline R_3(y)\le255^{-4}y^{13}.
\]

For two inputs, an exceptional support has two exceptional placements among 120, giving \(255/60=17/4\). For three inputs, it has two among 560, giving \(255/280<1\). The reference can evaluate the exact averaged polynomials instead of these coarse bounds.

## Conditional rank bound for return and selected output zeros

The rank argument does not assume a fixed field-valued message. Fix one common span \(\mathcal W\) of previously imposed equations. Required input forms \(\ell_i\) remain nonzero in its quotient. Split the incoming witness-state support into nonempty sets \(H\), and **sum** their contributions. A maximum over \(H\) would not cover different witnesses in the same quotient.

For one branch \(H\), first impose \(S_{H^c}=0\). Discard it if a required input form or a required \(S_j\), \(j\in H\), vanishes. Omitting the rank charge for these initial equations is conservative.

Let \(Z\) be a chosen set of output coordinates required to be zero. Define

\[
 K_{I,H,Z}:=\begin{pmatrix}C_{H^c,I}\\(B_I)_{Z,I}\end{pmatrix},
 \qquad r_{I,H,Z}:=\operatorname{rank}K_{I,H,Z}.
\]

The branch is impossible if \(\ker K\) forces any input coordinate to zero, or if some row \(C_{j,I}\), \(j\in H\), vanishes on \(\ker K\). The latter would force a required incoming state coordinate to zero on return.

For every remaining branch, choose \(r\) independent columns of \(K\). Reveal all nonpivot input diagonals. Row reduction expresses \(K D_{\rm in}\ell=0\) as \(r\) equations whose pivot diagonal occurs in its own required nonzero input form. The alive diagonal lemma charges \(\tau^{-r}\), for \(1\le\tau\le255\).

After all input diagonals are known, impose the \(h\) active state-return equations. Their state diagonals are still independent and fresh. Each equation has one required nonzero state form and contributes another factor \(\tau^{-1}\). Consequently, the conditional surviving rank potential contracts by at most

\[
 \boxed{\tau^{-(h+r_{I,H,Z})}.}
\]

This chronology uses each diagonal once. Additional input/state dependencies in the current quotient can kill branches; they do not invalidate the pivot charge.

For a word with output zero set \(Z(Y)\), the positive identity

\[
 y^{\operatorname{wt}(Y)}
 =\sum_{Z\subseteq Z(Y)}y^{16-|Z|}(1-y)^{|Z|}
\]

supplies the selected-zero marks. Thus define the local rank-return upper action

\[
 \mathcal R_I(y,\tau)=
 \sum_{\substack{H\ne\varnothing,\ Z\subseteq[16]\\
                         (I,H,Z)\ \mathrm{feasible}}}
 y^{16-|Z|}(1-y)^{|Z|}\tau^{-(|H|+r_{I,H,Z})}.
\]

The reference sums \(H\) and \(Z\), then maximizes over input placements and phases. Its optional placement average is valid when the encompassing transfer sums or averages uniform placements. It does not average phases unless that sampling law is supplied.

Feasibility is recovered from the exact projective strata by taking every submask of each output-zero mask, and deduplicating \((H,Z)\). No projective multiplicity enters this rank expression. Independently, feasibility can be checked by row ranks: none of the required coordinate forms or active state-target forms may belong to the row span of \(K\). These criteria agree. There are at most seven required nonzero forms, fewer than 256, so their proper hyperplanes cannot cover a nonzero kernel. The same criterion remains valid over extension fields.

Only \(|Z|\le3\) can be feasible, by the complete return-output enumeration. Exact output patterns also lift to extension fields: fewer than 256 required nonzero forms are involved. Hence this finite feasibility check is not limited to one-byte consumer records.

### A compact two-input result

Put \(a=1-y\). The worst-phase placement-averaged two-input action has the closed form

\[
 \overline{\mathcal R}_2(y,\tau)=
 y^{16}\tau^{-4}\max\{5,(\tau+73)/15\}
 +\tau^{-5}(16y^{15}a+13y^{14}a^2+4y^{13}a^3).
\]

At \(\tau=35\), \(y=1\), it is \(7.2\tau^{-4}\), about 170 times smaller than the old \(\tau^{-2}\) return term. This is a local comparison, not a whole-code failure margin.

## Coarse exact-support alternative

A separate upper cover imposes return before considering exact output supports. Put \(r=\operatorname{rank}C_{H^c,I}\) and \(d=w-r\). After return, the input-form image lies in a space of dimension at most \(d\).

Each realizable exact output-zero pattern is the closure of a span generated by its zero output forms. A canonical basis of rank \(k<d\) can be chosen from sixteen forms. There are at most \(\binom{16}{k}\) such spans, each adding rank \(k\). If \(m(I,H)\) is the minimum return-output weight from the exact histogram, a valid continuation bound is

\[
 \sum_{H\ \mathrm{feasible}}
 \tau^{-(h+r)}y^{m(I,H)}
 \sum_{k=0}^{d-1}\binom{16}{k}\tau^{-k}.
\]

The extra sum accounts for different surviving witnesses with different output supports. Bare \(y^{13}\) would omit that cover. This alternative can replace the return cover with a uniform continuation bound; it is not asserted to dominate each selected-zero branch individually.

The conditional proof and both weighted covers received an independent mathematical audit. A whole-code application must still retain one common equation span through all packets and include its own source/route counting.

## APIs and replay

The new file [cauchy_joint_return.py](theory/cauchy_joint_return.py) exports:

```python
joint_return_moment(phase, inputs, state_support, y)  # exact fixed-word law
uniform_return_bound(w, y, phase=None, average_inputs=False)
rank_return_moment(w, y, tau, phase=None, average_inputs=False)
coarse_rank_return_moment(w, y, tau, phase=None, average_inputs=False)
```

The last three optional arguments are keyword-only. Integer or Fraction marks give exact rational values. Floating marks are exploratory. Cached coefficient profiles permit repeated evaluations without repeating the geometric enumeration.

```sh
python -B research/workstreams/char2_field_spin/theory/cauchy_joint_return.py --weights 2 3
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p test_cauchy_joint_return.py
```

Eight tests passed. They check all 2,720 actual pair/triple placements, eight independent projective-vector enumerations, four independent selected-zero rank/feasibility enumerations, and the one-input regression. Small-field tests exhaust the joint diagonal experiment and 81 conditional alive-rank inequalities with dependent forms. None is a benchmark.

Next, use the audited rank-return profiles in the regional sparse transfer and replay any positive whole-class result with directed arithmetic. The low-input-weight class remains open until that separate composition succeeds.
