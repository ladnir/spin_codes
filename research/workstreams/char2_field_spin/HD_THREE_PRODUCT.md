# A three-product two-state inner

One coefficient can be removed from the four-product inner while preserving
two genuine state variables and a useful support bound. This changes the
setup distribution. It is a proof-side candidate, not a measured replacement
for the 8.3-ms encoder or a new distance certificate.

At each position independently sample nonzero d1,d2,f from F=GF65536.
Starting with s1=s2=0, compute

    a = s1+d1*x;
    y = a+f*s2;
    b = s2+d2*x;
    (s1,s2) = (a,b).

This sets c=d1 in the preceding inner. The local input/state map has
determinant f*d2, so it is always invertible. The complete input/output
matrix is triangular with diagonal d1. No coefficient rejection or
conditioning factor is required.

Both state streams remain weighted XOR prefixes. Batching needs three
field products per input: d1*x, d2*x, and f times the old second state.
The existing inner needs four. An implementation would form an exclusive
second prefix for the observation; it must not substitute the updated
second state. This product count does not predict the full encoder time.

## Exact support law

Write p=1/(q-1), where q=|F|. A support state records which state coordinates
are nonzero. Fix the incoming values and the current input before sampling
the fresh coefficients.

For zero input the state is unchanged. The output is zero for state00,
nonzero for states10 and01, and zero with probability p for state11.

For nonzero input, a and b are independent. If an old coordinate is zero,
its new coordinate is nonzero. If an old coordinate is nonzero, its new
coordinate is zero with probability p. Conditional on those new support
flags, the output law is:

- If the old s2 is zero, y=a exactly.
- If the old s2 is nonzero and a=0, y is nonzero.
- If the old s2 and a are both nonzero, y is zero with probability p.

These statements follow from uniform nonzero d1,d2,f and remain valid for
arbitrary incoming nonzero values. Thus they compose even when routed input
values repeat. `theory/hd_three_product.py` implements the resulting exact
four-state transition probabilities.

The new law differs materially from the earlier conditioned ensemble.
For example, a live first state, zero second state, and live input produce
a zero output together with a zero new first state with probability p.
The earlier ensemble assigns a smaller probability to that joint event.
Its rank-pivot proof must not be reused: the output equation and new-first-
state equation can now coincide.

## A shared two-state upper bound

For a positive output mark v<=1, set b=v+(1-v)p. Aggregate the three
nonzero state supports into one live-vector state. The following positive
matrices upper-bound the weighted transitions, with row/column order
zero, live:

    T0 = [[1, 0],     T1 = [[0, v],
          [0, b]]          [v*p^2, b]].

T0 handles zero input and T1 handles nonzero input. A return of both states
requires two independent cancellations, has probability at most p^2, and
has nonzero output. The total live-to-live output moment is at most b:
the output-zero probability is at most p. Birth from the zero state has
nonzero output and a live successor.

The same matrices dominate the original conditioned four-product inner.
Each actual live support class satisfies both entrywise bounds, so every
nonnegative continuation preserves the inequality. Consequently, a complete
certificate based solely on this upper operator could also cover the
three-product variant. A certificate using the exact four-state law would
need to be evaluated again.

The envelope deliberately overcounts some joint mass. It does not infer
independence of a return and an output event. It also does not itself bound
the source components or supply an all-support certificate.

## Validation and next step

Five exact tests pass: every GF4 local value and coefficient tuple,
local invertibility, all two-step setups and messages, transition mass,
the shared envelope for both distributions, and an inverse-tail rank-two
fixture. Thus the extra state has not collapsed to a scalar inner.

```text
python -B research/workstreams/char2_field_spin/theory/hd_three_product_tests.py
```

First evaluate the common upper operator in the component proof. If it
gives useful all-support evidence, implement this distribution as a separate
cost control and measure the saved product. No native code or benchmark
was added for this candidate in this round.
