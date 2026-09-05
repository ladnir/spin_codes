# Conditioning-bound attempt

Fix a group-occupancy profile ((k_u)).  In one outer-coordinate region, the
exact lane law is the iid law conditioned on distinct lanes inside every
group.  For every nonnegative transfer entry,

\[
 R_{(k_u)}(z)
 \le \left(\prod_u p_{k_u}^{-1}\right)R_a^{\mathrm{iid}}(z),
 \qquad
 a:=\sum_u k_u.
\]

The 256 lane permutations for distinct outer coordinates are independent.
The complete inner moment therefore incurs the profile factor

\[
 \prod_u p_{k_u}^{-256}.
\]

Let (S:=8192/q).  Summing this factor over every active-block subset of size
(a) gives

\[
 \Gamma_a=
 \frac{[x^a]
 \left(\sum_{k=0}^{q}{q\choose k}p_k^{-256}x^k\right)^S}
 {{8192\choose a}}.
\]

The reduction is rigorous.  Applying it to the floating-point iid certificate
does not close 9% distance.  For (q=32), the adjusted bound first falls below
40 bits of margin at occupation 12.  At occupation 64, the surcharge is
21,155.75 bits, while the iid margin is 4,294.29 bits.

Smaller groups do not repair the bound.  The aggregate margins through
occupation 64 are negative for every tested value:

| (q) | Aggregate margin |
|---:|---:|
| 4 | -9,288.83 bits |
| 8 | -13,053.19 bits |
| 16 | -15,424.81 bits |
| 32 | -16,861.46 bits |

Rare dense groups dominate because the conditioning penalty is raised to the
256th power.  This is a failure of the general conditioning bound.  It is not
evidence of low-distance words.

The receipts are in `receipts/q*_conditioning_low64_delta09.json`.

