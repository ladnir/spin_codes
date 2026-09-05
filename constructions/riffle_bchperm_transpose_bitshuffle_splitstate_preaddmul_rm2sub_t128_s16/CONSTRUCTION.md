# Construction

Let \(V:=\mathbb F_2^7\). Fix eight linearly independent homogeneous
quadratic functions

\[
q_1,\ldots,q_8:V\longrightarrow\mathbb F_2.
\]

For each \(x\in V\), define

\[
v_x:=(1,x,q_1(x),\ldots,q_8(x))\in\mathbb F_2^{16}.
\]

The columns \((v_x)_{x\in V}\) define a matrix
\(G\in\mathbb F_2^{16\times128}\). Define

\[
A(Q):=G^\mathsf TQ,
\qquad
B(X):=GX.
\]

The selection receipt specifies the eight quadratic functions. It also
records the 128 columns and 16 generator words.

An epoch receives \(X_i\in\mathbb F_2^{128}\) and
\(Q_i\in\mathbb F_2^{16}\). Setup fixes
\(\alpha_i\in\operatorname{GF}(2^{16})^*\). The epoch emits

\[
Y_i:=X_i+A(Q_i)
\]

and updates the state by

\[
Q_{i+1}:=\alpha_iQ_i+B(X_i).
\]

Every entry of \(GG^\mathsf T\) is the sum over \(V\) of a polynomial of
degree at most four. Such a sum is zero because \(4<7\). Hence

\[
BA=GG^\mathsf T=0.
\]

The inverse epoch uses the same identity as the parent PreAddMul
construction.
