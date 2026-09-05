# Construction

This construction inherits the outer code and permutation layers from
`riffle_fieldcheckpoint_accumulate_t32_s64_k32`.

The inner state is (q\in\mathbb F_2^{64}). The encoder processes 32 input
bits per implementation step. For bit position (p), it updates lane
\(\ell(p):=p\bmod64\) by

\[
 q_{\ell(p)}\gets q_{\ell(p)}+x_p,
 \qquad
 y_p\gets q_{\ell(p)}.
\]

An epoch contains eight implementation steps, or 256 input bits. At every
boundary between consecutive epochs, setup samples an independent

\[
 a_e\gets\mathrm{GF}(2^{64})^*.
\]

The forward encoder replaces the state by (L_{a_e}^{T}q), where (L_{a_e})
is multiplication by (a_e) in a fixed polynomial basis. The transposed
evaluator applies (L_{a_e}).

Each 8192-bit transposed region contains 32 epochs. The complete inner word
contains 8192 epochs. Changing the epoch length defines a different
candidate.
