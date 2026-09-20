# Goal 01 report: exact accumulator stress test

## Result

The accumulator gives a clean interface between the constituent and the packet
permutation. Low output weight depends on an ordered prefix-state histogram,
not only on packet support or input bit weight.

Exact reduced experiments show that four-bit grouping can be mild for an
authenticated BCH profile and severe for repeated packet values. The
authenticated support-11 profile loses 4.317 probability bits relative to a
full bit permutation. Eight equal `0xf` packets lose 25.750 bits.

These results are `EXACT_REDUCED_INSTANCE` evidence. They make no full-size
distance claim.

## Exact constituent identity

Fix a sequence of \(n\) four-bit packets. Suppose its nonzero packets occur at

\[
1\le p_1<\cdots<p_h\le n
\]

and have ordered values \(v_1,\ldots,v_h\in\mathbb F_2^4\setminus\{0\}\).
Define

\[
q_j:=\sum_{i=1}^j v_i,
\qquad
a_j:=\operatorname{wt}(q_j),
\qquad
p_{h+1}:=n+1.
\]

The accumulator remains in state \(q_j\) from position \(p_j\) through
position \(p_{j+1}-1\). Therefore

\[
\operatorname{wt}(y)
=\sum_{j=1}^h a_j(p_{j+1}-p_j).
\tag{1}
\]

Equation (1) identifies the only low-weight mechanism. The active values must
be ordered so that their prefix XORs have low bit weight, and the long gaps
must occur during those low-weight states.

For a fixed active-value order, define

\[
r_0:=p_1-1,
\qquad
r_j:=p_{j+1}-p_j-1\quad(1\le j\le h).
\]

These integers are nonnegative and sum to \(n-h\). The exact output-weight
generating function over all support placements is

\[
z^{\sum_j a_j}
[t^{n-h}]
\frac{1}{1-t}
\prod_{j=1}^h\frac{1}{1-tz^{a_j}}.
\tag{2}
\]

Consequently, the placement spectrum for a fixed order depends only on the
five-bin histogram

\[
m_r:=\#\{j:a_j=r\},
\qquad 0\le r\le4.
\tag{3}
\]

The histogram in (3) is the compressed interface needed from the packet
permutation. A proof does not need the complete accumulator state trajectory
after this histogram is known.

## Scalar consistency gate

At packet width one, the exact input-output enumerator is

\[
A_{u,w}
=
\binom{n-w}{\lfloor u/2\rfloor}
\binom{w-1}{\lceil u/2\rceil-1}.
\tag{4}
\]

The packet dynamic program agrees with (4) on all three committed scalar test
cases. An independent top-down implementation also reproduces every four-bit
probability below.

## Reduced exact experiment

The experiment uses 64 packet positions and binary length 256. It counts
outputs of weight at most 23. Each row fixes the complete multiset of nonzero
packet values and averages over its uniform packet permutation.

| Packet multiset | Input bit weight | Packet permutation | Full bit permutation | Grouping loss |
|---|---:|---:|---:|---:|
| Authenticated support-11 BCH word | 26 | \(2^{-24.011606}\) | \(2^{-28.328729}\) | 4.317122 bits |
| Two copies of each basis packet | 8 | \(2^{-7.255703}\) | \(2^{-8.498834}\) | 1.243131 bits |
| Four selected equal-value pairs | 16 | \(2^{-10.587668}\) | \(2^{-18.804794}\) | 8.217126 bits |
| Eight copies of `0xf` | 32 | \(2^{-10.905112}\) | \(2^{-36.655192}\) | 25.750081 bits |

The last row exposes the packet-clumping obstruction. All four accumulator
lanes turn on and off together. A full bit permutation destroys that alignment,
but a four-bit packet permutation preserves it.

The authenticated BCH row is more favorable. Its varied packet values reduce
the grouping loss to about four bits in this instance. Thus packet width four
does not automatically defeat accumulator contraction. The outer spectrum must
exclude or sufficiently suppress the repeated-value profiles that resemble the
last two rows.

## Scaling boundary

The exact dynamic program indexes the remaining multiplicity of every nonzero
packet value. The support-11 authenticated profile finishes quickly. A direct
attempt on the complete support-33 triple did not finish within 90 seconds and
was stopped.

This scaling failure does not concern the accumulator recurrence. It comes from
enumerating packet-value orders. A full proof should therefore propagate a
bound on the five-bin histogram in (3), rather than enumerate every value
order or every constituent state.

## Proof and refutation obligations

For a fixed outer word, the uniform packet permutation samples an ordered
multiset of active packet values and its gaps. Equation (2) handles the gaps
exactly after the prefix-state histogram is fixed. The remaining proof
obligation is a tail bound for that histogram under a random ordering of the
outer word's packet multiset.

The matching refutation search should find outer words whose packet multiset
admits many orderings with large \(m_0\) or \(m_1\). Such orderings spend many
positions in states of weight zero or one and therefore produce unusually low
accumulator output.

## Conclusion

The stress test validates the proposed decomposition. The constituent side
collapses to an exact five-bin interface. The unresolved part lies entirely in
the ordered packet-value spectrum supplied by the outer word and packet
permutation.

The next bounded goal is to derive a transfer operator on
\((c_1,\ldots,c_{15},q,m_0,\ldots,m_4)\) that bounds the histogram tail without
enumerating all full trajectories. Apply it first to the 26 authenticated
support-33 words.
