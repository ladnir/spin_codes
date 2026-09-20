# Goal 05: refutation of the paired nine-node return bound

## Result

The proposed return bound is false. Let

\[
s=(a,b)
\]

with

\[
a=\mathtt{89cf7abe4fbc2365},
\qquad
b=\mathtt{87452995c56be123}.
\]

This lifted state is nonzero. Its first 18 zero-input node weights are

\[
(0,7,0,20,0,10,0,10,0,
  5,0,18,0,9,0,11,0,7).
\]

Therefore

\[
W_9(s)=47
\]

and

\[
W_9(R^9s)=50.
\]

Hence

\[
W_9(s)+W_9(R^9s)=97<106.
\]

The state refutes the universal paired return-cost claim.

## Relation to the low list

Let

\[
\mathcal L_{52}
:=
\{c\in\mathcal C_9:\operatorname{wt}(c)\le52\}.
\]

The first nine-node word produced by (s) belongs to

\[
\mathcal L_{52}
\]

because its weight is 47. Its shifted word also belongs to

\[
\mathcal L_{52}
\]

because its weight is 50. Thus the shift map does not force the required
recovery cost on every low word.

## Independent verification

The independent verifier reconstructs the systematic parity map from

\[
\mathtt{f4845518b9582a1f}.
\]

It evaluates accumulation one bit at a time. It then applies the lifted
recurrence for 18 nodes. The verifier does not import the search recurrence,
its lookup tables, or the certificate generator.

The replay reproduces all 18 node weights and both nine-node sums. The
verification receipt records the node outputs, the state after nine nodes,
and hashes of the verifier and raw counterexample.

## Consequence

The earlier gap argument required every complete 18-node block to contribute
at least 106 output bits. This counterexample contributes only 97 bits.
Therefore that gap argument cannot close the support-33 row as stated.

The counterexample refutes the paired 18-node lemma. It does not refute the
full construction. A replacement proof may use a longer return window,
several consecutive windows, or the global permutation structure.
