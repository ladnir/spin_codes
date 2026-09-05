# Goal 04: lifted-orbit gate

## Objective

Determine whether the 128-state lifted zero-input recurrence can replace the
open early-start response-spectrum argument.

Let \(A=\operatorname{Acc}\), \(P\) be the systematic BCH parity map, and
\(T=P\circ A\). On a node with zero original drive, the two zero-state laps
have stored states \(a_t,b_t\in\mathbb F_2^{64}\) and retained output
\(q_t\in\mathbb F_2^{64}\):

\[
q_t=A^2a_t+Ab_t,
\qquad
\begin{pmatrix}a_{t+1}\\b_{t+1}\end{pmatrix}
=
\begin{pmatrix}T&0\\TA&T\end{pmatrix}
\begin{pmatrix}a_t\\b_t\end{pmatrix}.
\tag{LO}
\]

For \(m\ge1\), define the lifted observation code

\[
\mathcal C_m
:=
\left\{(q_0,\ldots,q_{m-1}):
(a_0,b_0)\in\mathbb F_2^{128}\right\}.
\]

## Proof gate

The goal must complete the following finite tasks.

1. Reconstruct (LO) from the authenticated node recurrence.
2. Determine the rank and exact minimum distance of the first useful
   observation codes.
3. Determine the exact zero-output and low-output return constraints exposed
   by those codes.
4. Compare each certified block density with
   \(188{,}766/32{,}772\).
5. Decide whether a universal block argument can control the zero-input gaps
   between occupied nodes.

If a short-window distance has sufficient density, the next goal may combine
it with the global packet-permutation gap law. If low-weight lifted
trajectories defeat every useful short window, record an exact counterexample
and stop this route.

## Scope

This goal studies deterministic zero-input gaps. It does not assume that
different gaps are independent. It does not claim that a local distance bound
already proves the support-33 row.
