# Goal 01: accumulator stress test

## Question

Can the current proof architecture cleanly separate an exactly analyzable
recursive constituent from the packet permutation and outer-code constraints?

## Target

For Riffle ParallelAcc g=4:

1. derive the exact output-weight formula for a fixed ordered packet support;
2. recover the standard scalar accumulator enumerator at packet width one;
3. exhaust the packet-width-four distribution for reduced fixed multisets;
4. compare packet permutation with full bit permutation at equal input weight;
5. state the remaining deterministic-interleaver obligation.

The reduced experiment uses 64 four-bit packet positions. Its failure threshold
is \(\lfloor 0.09\cdot 4\cdot64\rfloor=23\) bits. One tested multiset is the
support-11 local BCH word underlying an authenticated support-33 outer word.

## Evidence labels

Exact integer dynamic programming produces `EXACT_REDUCED_INSTANCE` evidence.
The experiment does not prove or refute the full-size construction.

## Pass condition

Goal 01 passes if an independent scalar formula agrees with exhaustive dynamic
programming at width one, and the packet-width-four results expose one explicit
support statistic required from the riffle layer.
