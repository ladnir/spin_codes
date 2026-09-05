# Goal 04: lifted-orbit gate

## Result

The lifted-orbit viewpoint gives a finite alternative to the response-spectrum
problem, but the first useful theorem remains open.

Exact primary and independent computations determine the observation-code
distances through seven nodes:

\[
(d_1,\ldots,d_7)=(1,1,1,3,4,5,6).
\]

They also prove \(d_8\ge6\). Exact witnesses give

\[
d_8\le37,
\qquad
d_9\le44.
\]

These witnesses refute every universal gap-block argument through nine nodes.
Ten nodes is the first surviving window. The required statement is

\[
d_{10}\ge59.
\tag{G10}
\]

The best exact ten-node witness found has weight 76. Bounded SAT and MILP
probes did not prove (G10) or find a counterexample. Thus Goal 04 replaces the
large full-response spectrum question with one smaller coding theorem, but it
does not yet close the early-start row.

## Authenticated recurrence

Let \(A=\operatorname{Acc}\), let \(P\) be the systematic BCH parity map, and
let \(T=P\circ A\). During a zero-input node, let \(a_t\) and \(b_t\) be the
stored states of the first and second zero-state laps. The retained output is
\(q_t\). Direct substitution in the node recurrence gives

\[
q_t=A^2a_t+Ab_t,
\qquad
\begin{pmatrix}a_{t+1}\\b_{t+1}\end{pmatrix}
=
\begin{pmatrix}T&0\\TA&T\end{pmatrix}
\begin{pmatrix}a_t\\b_t\end{pmatrix}.
\]

The state-transition matrix is invertible. For \(m\ge1\), define

\[
\mathcal C_m
:=
\{(q_0,\ldots,q_{m-1}):(a_0,b_0)\in\mathbb F_2^{128}\}.
\]

The code has dimension 64 for \(m=1\) and dimension 128 for every \(m\ge2\).
In particular, the first two outputs determine the complete lifted state.

## Exact small-window census

The primary computation constructs a parity-check matrix for each
\(\mathcal C_m\). It exhausts every syndrome relation of weights one through
four. For \(m\in\{6,7,8\}\), it also exhausts every weight-five relation by a
pair-versus-triple meet-in-the-middle search.

The independent computation uses a different representation. It reduces the
first two output blocks to systematic form and constructs the coordinate
syndromes directly from the remaining blocks. It reproduces every rank,
minimum-distance witness, and weight-five exclusion.

The resulting exact statements are:

| \(m\) | Exact result |
|---:|---:|
| 1 | \(d_1=1\) |
| 2 | \(d_2=1\) |
| 3 | \(d_3=1\) |
| 4 | \(d_4=3\) |
| 5 | \(d_5=4\) |
| 6 | \(d_6=5\) |
| 7 | \(d_7=6\) |
| 8 | \(d_8\ge6\) |

The weight-six decision for \(\mathcal C_8\) was not certified.

## Exact low-weight witnesses

The seven-node minimum word has node weights

\[
(0,2,0,2,0,1,1).
\]

The trajectory becomes ordinary high weight at the next node. Thus the word
is a low-weight transient, not a long low-density orbit.

A separate exact solver witness for nine nodes has weights

\[
(0,2,0,10,0,5,0,26,1),
\]

whose total is 44. Direct recurrence replay authenticates its initial state
and every output word. Coordinate descent produces an eight-node word of
weight 37. Each witness is sufficient to refute the corresponding proposed
lower bound; neither search proves minimum distance.

## Gap accounting

Use the boundary-shift representation after nodes zero and one. The remaining
32,770 nodes contain at most 33 nodes occupied by the authenticated
support-33 word. Therefore, at least 32,737 zero-input nodes lie in at most 34
gaps.

First split off trajectories whose lifted state reaches zero. For a fixed
drive and a fixed time, at most one wrapped state causes this event. The
existing return-state union bound charges this subset below \(2^{-44}\).
On the complement, the lifted state is nonzero in every complete zero-input
block.

For a window of \(m\) nodes, the number of complete blocks is at least

\[
B_m
:=
\left\lceil
\frac{32737-34(m-1)}{m}
\right\rceil.
\]

A universal block proof needs

\[
d_m B_m>188{,}765.
\]

The relevant values are:

| \(m\) | \(B_m\) | Required \(d_m\) | Certified or witnessed |
|---:|---:|---:|---:|
| 4 | 8,159 | 24 | \(d_4=3\) |
| 5 | 6,521 | 29 | \(d_5=4\) |
| 6 | 5,428 | 35 | \(d_6=5\) |
| 7 | 4,648 | 41 | \(d_7=6\) |
| 8 | 4,063 | 47 | \(d_8\le37\) |
| 9 | 3,608 | 53 | \(d_9\le44\) |
| 10 | 3,244 | 59 | \(d_{10}\le76\) |

Every row through nine is refuted. The ten-node row is the first row not
refuted by an exact witness.

## Ten-node frontier

The remaining finite theorem is (G10). A direct SAT encoding of
\(d_{10}\le58\) did not finish within its bounded run. A separate 120-second
HiGHS formulation found an exact word of weight 172 and retained the trivial
dual bound zero. Deterministic coordinate descent found the stronger
weight-76 word.

These solver outcomes are diagnostic. They neither prove \(d_{10}\ge59\) nor
refute it. A theorem-facing continuation needs an independently checkable
lower-bound certificate for the binary \([640,128]\) code
\(\mathcal C_{10}\).

## Conclusion

The boundary-shift intuition yields a coherent proof route. Exact returns are
cheap to charge, and nonreturning trajectories reduce to lifted observation
codes on zero-input gaps.

The old short-window proof does not transfer directly. Structured transients
defeat every window through nine nodes. The route remains possible only
through the new ten-node distance target (G10). Until that target is proved,
the lifted approach is an alternative finite obstruction, not a replacement
for the response-spectrum proof.
