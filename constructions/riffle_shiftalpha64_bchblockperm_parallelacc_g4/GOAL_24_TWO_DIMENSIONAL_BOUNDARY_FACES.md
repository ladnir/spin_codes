# Goal 24: close the two-dimensional boundary faces

## Result

The two two-dimensional faces that occur in the weight-16 outer shell have
closed exact enumerators. The faces use packet weights $\{0,2,4\}$ and
$\{0,1,3\}$.

Together with Goal 23, the exact formulas cover six of the fifteen boundary
types. Their share of the exact shell rises from 43% at binary length 48 to
65% at binary length 136. Only faces of dimensions three and four remain.

## Composition convention

For nonnegative integers $m$ and $r$, define

\[
C(m,r):=
\begin{cases}
\binom{m+r-1}{r-1},&r>0,\\
1,&r=0\text{ and }m=0,\\
0,&r=0\text{ and }m>0.
\end{cases}
\tag{1}
\]

Thus $C(m,r)$ counts weak compositions of $m$ into $r$ parts, including the
case with no parts.

## Even face

Fix a histogram $(h_0,0,h_2,0,h_4)$. If $h_2$ is odd, its boundary
coefficient is zero. Suppose $h_2=2b$.

Every positive excursion has one of two forms:

\[
(4,4)
\qquad\text{or}\qquad
(2,4^j,2),\quad j\ge0.
\tag{2}
\]

The first excursion visits state weights $0,4,0$. The second visits state
weight two between its endpoint packets. Each weight-four transition at
state weight two replaces the support by its complement.

These are all possibilities. From state zero, a positive transition reaches
state two or four. State four can only return to zero. State two can either
return to zero with weight two or remain at weight two with weight four.

Let $a$ be the number of excursions of the first kind. The remaining
$j=h_4-2a$ weight-four packets occur inside the $b$ excursions of the second
kind. The exact boundary coefficient is

\[
\boxed{
A^{\mathrm{even}}(h_0,h_2,h_4)
=
6^b
\sum_{a=0}^{\lfloor h_4/2\rfloor}
\binom{a+b}{a}
C(h_4-2a,b)
\binom{h_0+a+b}{a+b}.
}
\tag{3}
\]

To justify (3), fix $a$. The factor $6^b$ chooses the first state support in
each weight-two excursion. The next factor orders the two excursion kinds.
The composition count distributes internal weight-four packets. The final
factor distributes zero packets among the zero-state gaps.

## Odd face

Fix a histogram $(h_0,h_1,0,h_3,0)$. Its boundary coefficient is zero unless
both $h_1$ and $h_3$ are even. Write

\[
h_1=2b,
\qquad
h_3=2c.
\]

Every positive excursion has one of the forms

\[
(3,3)
\qquad\text{or}\qquad
(1,3^{2t},1),\quad t\ge0.
\tag{4}
\]

For the second form, each pair of internal weight-three packets moves from
state weight one to two and back. Its multiplicity is $3\cdot2=6$.

These forms exhaust the odd face. From state zero, a positive transition
reaches state one or three. State three can only return with weight three.
State one can return with weight one or move to state two with weight three.
State two can only move back to state one with weight three.

Let $a$ be the number of excursions of the first kind. The remaining
$t=c-a$ internal pairs are distributed among the $b$ excursions of the
second kind. The exact coefficient is

\[
\boxed{
A^{\mathrm{odd}}(h_0,h_1,h_3)
=
\sum_{a=0}^{c}
4^{a+b}6^{c-a}
\binom{a+b}{a}
C(c-a,b)
\binom{h_0+a+b}{a+b}.
}
\tag{5}
\]

The power of four selects the outbound support of every excursion. The other
three factors count internal pairs, excursion order, and zero-packet gaps.

## Conditional probabilities

For either face, divide the exact coefficient by the histogram sequence
count from Goal 16. The denominators are

\[
T^{\mathrm{even}}
=
\frac{N!}{h_0!h_2!h_4!}6^{h_2}
\tag{6}
\]

and

\[
T^{\mathrm{odd}}
=
\frac{N!}{h_0!h_1!h_3!}4^{h_1+h_3}.
\tag{7}
\]

On the boundary $H=2D$, the ratios of (3) and (6), or of (5) and (7), equal
the exact conditional probability of output weight at most $D$.

## Proof-gym coverage

The union of Goals 23 and 24 covers these fractions of the exact
$H=16,D=8$ shell:

| Binary length | Exact types | Fraction covered exactly | Hybrid saddle error on complete shell |
|---:|---:|---:|---:|
| 48 | 6 of 15 | 43.0% | 0.1922 bits |
| 80 | 6 of 15 | 56.2% | 0.1640 bits |
| 112 | 6 of 15 | 62.1% | 0.1483 bits |
| 136 | 6 of 15 | 64.8% | 0.1399 bits |

The hybrid column inserts exact values for the closed faces and numerical
Gaussian saddle estimates for the other nine types. It is diagnostic, not a
proved upper bound.

## Verification

Two independent integer recurrences verify the formulas through packet
length 48. Each audit checks 20,825 histograms. All 41,650 comparisons pass.

Run

```powershell
python scripts/verify_riffle_boundary_even_face_formula.py --maximum-length 48
python scripts/verify_riffle_boundary_odd_face_formula.py --maximum-length 48
```

The receipts are `receipts/goal24_boundary_even_face_formula_audit.json` and
`receipts/goal24_boundary_odd_face_formula_audit.json`.

## Consequence

The low-dimensional boundary cases no longer need a uniform local-limit
theorem. The remaining shell types have active faces of dimensions three and
four. Their saddle points are interior after absent packet weights are
removed.

The next goal should derive a certified Fourier bound for one
three-dimensional face. The face supported on weights $\{0,1,2,3\}$ is the
best first target because it contains the largest remaining type at binary
length 48.
