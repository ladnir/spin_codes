# Riffle FrozenRandom512-P2-RandomStepConv g=4 sigma=20

Status: ACTIVE_MODEL_EXPLORATION

This proof-gym model replaces RM(4,9) with random linear constituents. It does
not propose a dense random encoder as the eventual implementation code.

The 16384 data symbols are partitioned into 4096 groups of four. During
setup, each group receives an independent linear injection

\[
E_i:\mathbb F_2^{256}\longrightarrow\mathbb F_2^{512}.
\]

Each injection is sampled uniformly. Setup freezes every injection, so
encoding is deterministic after setup. No local distance floor is imposed.
The random constituent typically has minimum distance near 58.

The independence between groups is a proof-gym choice. It makes the ensemble
first moment factor exactly. A later constituent can replace this ensemble if
it supplies a comparable fixed-code spectrum bound.

The baseline field outer uses two parity symbols. Each parity symbol is
encoded by the exact extended BCH \([128,64,22]\) constituent. The data and
parity outputs contain 524352 four-bit packets in total. A global packet
permutation precedes the one-lap random step convolution with 20 state bits.

The exact ensemble-average weight spectrum is

\[
\mathbb E[A_w]
=\frac{2^{256}-1}{2^{512}-1}\binom{512}{w}.
\]

The corresponding expected packet-support spectrum is

\[
\mathbb E[A_s^{\mathrm{pkt}}]
=\frac{2^{256}-1}{2^{512}-1}\binom{128}{s}15^s.
\]

The initial analysis deliberately studies the ideal constituent structure.
It does not claim an efficient implementation of a dense random encoder.

Retaining the exact two-parity BCH spectrum closes the complete one-active-
group contribution at `-91.83` bits. The earlier positive support-95 row was
an artifact of dropping both determined parities. The next sparse case has
at least two active data groups, where global parity cancellation is possible.

The first global pointwise scan is also negative. Its largest sampled exponent
is `-84.40` at total support 360, with occupation mode three. A complete
support cover would convert this diagnostic into an end-to-end first-moment
bound.

See `REPORT.md` and `GOAL_01_RANDOM_SPECTRUM_BASELINE.md`.
