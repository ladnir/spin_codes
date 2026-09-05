# Research TODOs

## Two-sided regular expanders

The prime-field EC certificate already uses two-sided regular incidence.
Investigate whether fixing both degrees improves the other ensembles:

- binary Expand--Accumulate;
- binary wrapped Expand--Convolute;
- prime-field Expand--Accumulate.

Compare ensembles at the same rate and total number of expander edges.  The
analysis must account for both effects of right-degree conditioning: reduced
low-support collisions and the entropy or conditioning cost in the union
bound.  Do not assume that two-sided regularity improves every support range.

The binary exact regional enumerator and certificate are recorded in
`notes/BINARY_TWO_SIDED_REGULAR.md`.  The EA candidate `62/31` fails at
message weight two.  Wrapped EC with degrees `10/5` and memory `15` is
certified with `32.4990` failure bits.  Its decimal markers are frozen in
`results/binary_biregular_ec_rate_half_d10_m15_gv.json`; a standard-library
checker validates the schema and support partition independently of Arb.

Wrapped EC with degrees `6/3` and memory `79` is also certified at the same
cutoff.  Its failure probability is below `3.558444e-7`, or `21.4222` bits.
The verifier uses layered binary64 matrices with independent exponent bands
for exact supports `1..383`; Arb evaluates the remaining scalar bounds.

The intermediate two-sided profiles are now certified as well.  Degrees
`14/7` with memory `6` give `26.0875` failure bits, and degrees `18/9` with
memory `4` give `21.6746` bits.  The degree-`14` profile minimizes the simple
work proxy `d_L+m`, with value `20`, among the four certified profiles.

At rate one half, `8/4` is not a fallback because even binary right degree
puts the all-one message in the kernel.  See
`notes/BINARY_DEGREE_SIX_OUTER.md`.

Memory `78` does not close with the present exact-block partition and marker
method.  The next proof question is whether a sharper partition or output-tail
bound can certify it.  Separately, replace the simple `d_L+m` proxy with a
realistic implementation cost model for all four certified profiles.
