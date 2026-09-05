# Initial full-stripe exploration

## Question

The global packet permutation admits a bulk late-start family at packet width
eight. The S-stripe family forces every active outer block to contribute the
same number of candidate packets to every inner region. This report asks
whether that constraint restores a 40-bit first-moment margin at 9% relative
distance.

The message contains \(2^{20}\) bits. Each outer constituent is an independent
uniform random rate-half injection. The inner encoder makes one pass and
discards its terminal state.

## Evaluator

Condition on \(a\), the number of active outer blocks. In the full-stripe
instance, each stripe contains one candidate packet from every active block.
The stripe permutation places these \(a\) candidates uniformly among the
\(L\) outer-block positions.

For each output-weight tilt, the evaluator computes the four entries of the
stripe transition matrix exactly. A log-domain dynamic program extracts every
candidate-count coefficient from zero through \(L\). The evaluator then
raises each two-state transition matrix to the \(B/g\) stripe count and sums
all occupations \(1\le a\le L\).

A fixed nonzero input to a random outer injection has a uniform nonzero
output. The calculation relaxes this output to independent uniform bits and
divides by \(1-2^{-B}\) for every active block. This change is an upper bound.
The only remaining analytic relaxation is the Chernoff bound for total output
weight.

The global-permutation control gives

\[
\lambda(4096,8,40,S=1)=-43959.07.
\]

The earlier global-support evaluator gave approximately \(-43958.12\). This
agreement checks the occupation model before striping.

## What full striping changes

At \(B=4096\) and \(\sigma=40\), exact full striping improves the margin to

\[
\lambda(4096,8,40,S=512)=-36132.37.
\]

Thus striping recovers about 7,827 bits but does not close the deficit.

The dominant tilted path has 99 active outer blocks. Its expected
stripe-level transitions are approximately:

| Transition | Expected count |
|---|---:|
| off to live | 256.00 |
| live to off | 255.00 |
| off to off | 0.78 |
| live to live | 0.23 |

The path therefore alternates between off-start and live-start stripes. An
off-start stripe activates near its end. The next stripe carries the episode
near its beginning and then terminates. One short live episode straddles each
pair of stripe boundaries. Full striping forces about one reset per pair, not
one reset per stripe.

This alternating-boundary family is the new canonical obstruction. The old
single late suffix is no longer dominant.

## Outer-size curve at sigma 40

| Outer length \(B\) | Full stripes \(S=B/8\) | Dominant occupation | Margin \(\lambda\) |
|---:|---:|---:|---:|
| 4096 | 512 | 99 | -36132.37 |
| 8192 | 1024 | 50 | -27570.36 |
| 16384 | 2048 | 25 | -10260.14 |
| 32768 | 4096 | 1 | +10278.20 |

Increasing \(B\) creates more stripe boundaries while reducing the number of
outer blocks. The alternating family must pay for proportionally more state
resets. At \(B=32768\), the bulk family falls below a one-active-block case.

The structured permutation also uses less setup entropy than a uniform global
packet permutation. The global permutation of 262144 packets uses about
4.34 million bits. Full striping uses about 1.98, 1.72, 1.47, and 1.21 million
bits for the four displayed outer lengths, respectively.

## First observed contours

The two viable boundary pairs use an output-tilt grid with spacing 0.01. The
\(B=8192\) high-memory probe uses spacing 0.02.

| Outer length \(B\) | Cell below contour | Cell above contour |
|---:|---:|---:|
| 8192 | \(\sigma=95:\lambda=-903.15\) | none observed; one-block floor |
| 16384 | \(\sigma=50:\lambda=-18.71\) | \(\sigma=51:\lambda=1005.29\) |
| 32768 | \(\sigma=27:\lambda=-1325.98\) | \(\sigma=28:\lambda=725.78\) |

Negative values mean that this upper-bound method does not prove the target.
They do not constitute refutations. The two positive cells have much more
than 40 bits of numerical slack.

## Current status

The full-stripe construction has moved the problem from a bulk late-start
floor to a tunable reset tradeoff. The cells \((B,\sigma)=(16384,51)\) and
\((32768,28)\) are candidates for a closed ensemble proof.

The next proof step should formalize the regional transition polynomial and
replace floating-point log arithmetic with outward-rounded interval
arithmetic. The finite output-tilt grid already gives valid Chernoff choices;
it does not require proving that the grid contains the analytic optimum.
