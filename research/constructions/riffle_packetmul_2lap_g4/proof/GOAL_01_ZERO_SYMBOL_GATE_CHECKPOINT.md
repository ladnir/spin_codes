# Goal 01 checkpoint: the PacketMul zero-symbol gate

## Result

Goal 01 is complete. The PacketMul interface reduces the authenticated
support-33 terminal-zero row to one mixed-character statement.

> **Zero-symbol lemma.** For every nonzero terminal character
> \(\chi\in\mathbb F_2^{64}\), at most \(5M/8\) packet cells \(p\) satisfy
> \(a_\chi(p)=0\), where \(M=524{,}352\).

This lemma is not proved here. If it holds, then the aggregate terminal-zero
probability for the 26 authenticated support-33 outer words is less than
\(2^{-40.980718082061}\). The target is therefore below \(2^{-40}\) with more
than \(0.980718082061\) bits of margin.

No searched character violates the lemma. That negative search result is
diagnostic, not a proof.

## Terminal contribution map

Let \(T=P\circ A\) be the zero-input 64-bit state recurrence. Here \(A\) is
the accumulator map and \(P\) is the systematic BCH state-half map. Index a
packet cell by a slot \(s\in\{0,\ldots,15\}\) and a terminal-distance exponent
\(t\in\{1,\ldots,32{,}772\}\). For packet value \(y\in\mathbb F_{16}\), define

\[
z_y(s,t)=T^t\bigl(y\,2^{4s}\bigr).
\]

Reversing the physical node order only relabels the uniformly sampled
terminal-distance exponent. It does not change the multiset used below.

The primary implementation reconstructs \(P\) from the committed BCH
generator, reconstructs \(A\), and evaluates this recurrence. The independent
implementation imports no construction-map helper. It performs a separate
row reduction for \(P\), evaluates \(A\) bit by bit, and implements field
arithmetic by schoolbook polynomial reduction. Both implementations obtain
the basis-contribution digest

```text
3980852be768c5aecc6092fd2e9890bc56e331598059889aa0484b270824ecad
```

The digest covers the four basis values, all 16 slots, and all 32,772
exponents in the stated order.

Linearity of \(T\) gives

\[
z_y(p)=\bigoplus_{b:y_b=1}z_{2^b}(p).
\]

For a fixed character \(\chi\), the four parities
\(\langle\chi,z_{2^b}(p)\rangle\) specify a binary linear functional of
\(y\). The trace pairing on \(\mathbb F_{16}\), with modulus
\(X^4+X+1\), is nondegenerate. Consequently, there is one coefficient
\(a_\chi(p)\) such that

\[
\langle\chi,z_y(p)\rangle
=\operatorname{Tr}(a_\chi(p)y)
\quad\text{for every }y\in\mathbb F_{16}.
\]

For a nonzero source packet, multiplication by an independent uniform element
of \(\mathbb F_{16}^{\times}\) makes \(y\) uniform over the 15 nonzero field
elements. Therefore, with

\[
q_\chi=\frac{1}{M}\bigl|\{p:a_\chi(p)=0\}\bigr|,
\]

the normalized PacketMul character coefficient is exactly

\[
\beta_\chi
=q_\chi-\frac{1-q_\chi}{15}
=\frac{16q_\chi-1}{15}.
\]

## Cross-value collisions and Parseval

The relevant sample map is larger than the fixed-value map used for the
paused candidate:

\[
\Phi(p,y)=z_y(p),
\qquad (p,y)\in[M]\times\mathbb F_{16}^{\times}.
\]

It has \(15M=7{,}865{,}280\) labeled inputs. Exact enumeration gives the
following multiplicities.

| State multiplicity | Number of states |
|---:|---:|
| 1 | 7,472,040 |
| 2 | 196,620 |

No state has larger multiplicity. Thus \(\Phi\) is not injective across
packet values, even though every fixed-value map is injective. The number of
ordered colliding pairs is

\[
7{,}472{,}040+4(196{,}620)=8{,}258{,}520.
\]

Let \(m_x\) be the multiplicity of state \(x\). Character orthogonality gives
the exact normalized Parseval identity

\[
\sum_{\chi\in\mathbb F_2^{64}}\beta_\chi^2
=\frac{2^{64}\sum_xm_x^2}{(15M)^2}
=\frac{2479537839642119241728}{1006878735}.
\]

The independent implementation reconstructs the same collision histogram
and ordered collision count.

## Exact pure-component spectra

The transpose recurrence has pure irreducible components of degrees
\(1,2,4,9,10,18,20\). Every nonzero character in each component was
exhausted. For each character, Walsh inversion uses

\[
\bigl|\{p:a_\chi(p)=0\}\bigr|
=\frac{1}{16}\left(M+
\sum_{y\in\mathbb F_{16}^{\times}}
\sum_p(-1)^{\langle\chi,z_y(p)\rangle}\right).
\]

| Degree | Maximum zero count | Maximum \(q_\chi\) | Maximum \(|\beta_\chi|\) |
|---:|---:|---:|---:|
| 1 | 32,772 | \(1/16\) | \(0\) |
| 2 | 32,772 | \(1/16\) | \(0\) |
| 4 | 32,775 | \(10925/174784\) | \(1/122895\) |
| 9 | 33,882 | \(5647/87392\) | \(37/16386\) |
| 10 | 32,304 | \(673/10924\) | \(41/40965\) |
| 18 | 33,017 | \(33017/524352\) | \(79/122895\) |
| 20 | 33,242 | \(16621/262176\) | \(47/49158\) |

The character maximizing \(q_\chi\) need not maximize \(|\beta_\chi|\).
The independent implementation replays all ten distinct reported maximizers,
including their complete 16-symbol coefficient histograms.

## Mixed-character refutation search

The search evaluates 221 distinct seed characters. The seeds include:

- every pure-component maximizer;
- all 127 nonempty XOR combinations of the seven pure \(q\)-maximizers;
- mixed witnesses committed for the paused DP-2Lap candidate; and
- 64 pseudorandom restarts from seed `20260821`.

It performs 114 strict coordinate-ascent searches. When a start has a
component above the degree-1 and degree-2 subspace, the search rejects any
move that would remove every such component.

The largest exact value found is attained by

\[
\chi=\mathtt{0xb827b1995b22a37b}.
\]

For this character,

\[
|\{p:a_\chi(p)=0\}|=43{,}697,
\qquad
q_\chi=\frac{43{,}697}{524{,}352},
\qquad
|\beta_\chi|=\frac{2185}{98316}.
\]

The independent implementation replays this coefficient histogram and zero
count. The search does not certify a maximum over all \(2^{64}-1\) nonzero
characters.

## Support-33 terminal-zero gate

First sample the 33 labeled packet cells independently and uniformly. The
packet multipliers are independent, so the 33 contributions are iid under
this experiment. Let \(D\) be the event that the cells are distinct. Then

\[
\Pr[D]=\frac{(M)_{33}}{M^{33}}.
\]

Conditioned on \(D\), the cells have the ordered-injection law induced by the
packet permutation. Thus any event under the permutation law has probability
at most its iid probability divided by \(\Pr[D]\).

Suppose \(|\beta_\chi|\le B\) for every nonzero \(\chi\). Fourier inversion,
followed by the exact Parseval identity for two of the 33 powers, gives

\[
\Pr_{\mathrm{iid}}[L=0]
\le 2^{-64}\left[1+
\left(\sum_\chi\beta_\chi^2-1\right)B^{31}\right].
\]

All 26 authenticated words have packet support 33. The common sufficient
condition for their aggregate is therefore

\[
\frac{26}{\Pr[D]}2^{-64}
\left[1+
\left(\sum_\chi\beta_\chi^2-1\right)B^{31}\right]
\le2^{-40}.
\]

Solving the equality gives the rigorous decimal bracket

\[
0.613302417667203
\le B_{\mathrm{req}}
<0.613302417667204.
\]

Equivalently, it would suffice to prove
\(q_\chi\lesssim0.637471016563003\). The cleaner lemma
\(q_\chi\le5/8\) implies

\[
-\frac1{15}\le\beta_\chi\le\frac35,
\qquad |\beta_\chi|\le\frac35,
\]

and yields the aggregate interval

\[
2^{-40.980718082062}
\le U
\le2^{-40.980718082061}.
\]

## Evidence boundary and next obligation

The contribution map, trace interface, collision histogram, Parseval identity,
pure-component spectra, distinct-cell conditioning, and numerical threshold
are exact. The mixed search reports exact statistics for its characters, but
its coverage is diagnostic.

The remaining obligation for this ledger row is exactly the zero-symbol lemma
stated at the beginning. Proving that lemma closes only the authenticated
support-33 terminal-zero row. It does not prove the full construction, cover
other placement strata, or certify implementation performance.

## Receipts

- `../receipts/goal01_zero_symbol_primary.json`
- `../receipts/goal01_zero_symbol_independent.json`
- `../receipts/goal01_mixed_character_search.json`
- `../receipts/goal01_terminal_zero_gate.json`

