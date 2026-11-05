# Birth constraints and a structural rank diagnostic

The [active-row count](RANDOM_PAIRING_ACTIVE_ROWS.md) still fails at
intermediate turnover. This pass keeps the construction fixed and asks
whether omitted birth constraints explain that gap. They improve some
bounds, but the tested longest-suffix profiles remain open. No new
distance certificate or native timing is claimed.

## Retain nonzero inputs at empty-pair births

Use M=8192, R suffix regions, L=16R updates, state weight256, and B swaps.
Let U be the enlarged union of consecutive supports, with A=L*256+B+256
positions. There are B marked births after initial activation.

A birth in an empty incoming pair requires X=Q_next nonzero. Thus, if its
outer row is inactive, that birth must instead be a crossing endpoint.
The fair-coin comparison from the capped proof assigns probability one
half to this necessary crossing. At other outside-U coordinates, an
inactive row still requires the crossing coin to be zero.

Fix an exact active row set S. Introduce positive markers t and v for
total U membership and total birth membership. Represent the regional
route by independent three-category draws, then condition on both counts.
The categories are ordinary U, marked birth, and outside U. Sequential
binomial conditioning costs at most (M+1)^(2R). Before conditioning, rows
are independent. Arithmetic-geometric mean bounds across regions use

    p=A/(RM), beta=B/(RM),
    a=(p-beta)*r/t + beta*r/(t*v) + 1-p,
    b=(p-beta)/t + beta/(2*t*v) + (1-p)/2.

For 1<r<=q=2^32, the resulting all-row-set upper bound is

    H rho^T (M+1)^(2R)
    r^[(1+log2 C)/32-B-z] t^A v^B (r^-128*a^R+b^R)^M.

Here T=L-1, z=max(0,L*256+B-W), W=964689, and C=binom(T*256-B,z).
H is the previous exact-B history-count upper bound. The same rho
dominates the full inactive crossing set, including marked births;
there is no second matching-density charge for retaining these coins.

## What the selected screens show

All numbers below are FLOAT log2 upper bounds, not failure lower bounds.
The optimizer proposes positive weights; it does not certify optimality.

| R | B | Previous active-row bound | With birth compatibility |
| --- | --- | --- | --- |
| 256 | 32768 | +54115 | +56544 |
| 256 | 65536 | +200645 | +201830 |
| 256 | 131072 | +297971 | +296463 |
| 192 | 65536 | -63543 | -75233 |
| 192 | 131072 | +208414 | +189273 |

The extra conditioning charge can exceed the gain. The old bound remains
available; the new expression does not replace it when it is weaker.

For R256,B131072, the old optimized row-set sum is concentrated near6218
active rows. Its tilted U-occupancy fractions are approximately0.703 on
active rows and0.120 on inactive rows. This is a diagnostic of the moment,
not a sampled codeword or a typical setup. Most state-union positions
already avoid inactive rows, so the birth condition removes little mass.

An explicitly optimistic diagnostic requires every birth AND death to
land on an active row, whether or not the real pair arithmetic requires
it. Its formula uses beta=2B/(RM), deletes the marked term from b, and
uses v^(2B). This is NOT an upper bound for the actual code. Its selected
R256 values remain positive: +53876, +195564, +282518 at the three B values
above. This suggests that more tuning of this particular marked-row
moment is unlikely to resolve the main gap. It does not prove that no
better weights, different bounds, or joint pair constraints can do so.

## A separate forward-count check

To avoid committing exclusively to the capped count, a second helper
sums input-support patterns and state transitions directly. Fix s active
outer rows. For a specified input support of size a, the number of words
with exactly those active rows is at most q^(a-128s). If any active row
has fewer than129 positions, its actual contribution is zero; allowing
the positive formal weight only enlarges the count.

Set p=s/M. Under a Bernoulli surrogate for active-row labels, use marker t
to retain the actual Rs active labels in the suffix. Define

    a0=1-p+p/t, a1=p*q/t, f=a0+a1.

These weights sum zero and nonzero input-support choices, not field
amplitudes sampled independently of the outer. For an incoming support
of weight w and k crossing pairs, a uniform-image comparison gives the
one-step generating polynomial

    (f/q)^(w+k)
    [x*y+(1+(q-2)*y)*z]^w
    [1+(q-1)*y*z]^k
    [a0+a1*y*z]^(512-w-k).

Here x marks deaths, y output weight, and z next-state weight. Average over
the exact matching law P_w(k). The accepted-map density factor less than2
is paid globally, not per step. The source-label conditioning cost is
(M+1)^R, and its conserved marker contributes t^(Rs).

For the half-full family, extract or upper-bound the coefficient of z^256
at each subsequent step. Initial activation contributes
binom(512,256)*a1^256*a0^256*y^256. Finally include
binom(M,s)*q^(-128s), the output threshold marker, and optionally x^-B.
The helper relaxes q-1 and q-2 to q. All chosen weights are positive.

The selected forward bounds are weaker on the difficult classes. With a
turnover marker, R256,s6218,B131072 still gives about+1.81 million bits,
versus about+0.30 million for the capped sum over all row sets. The first
unconserved source-label proposal was weaker still; it was corrected to
retain Rs. These checks do not justify replacing the capped proof.

## A coefficient-independent test for a real obstruction

The next diagnostic should distinguish a loose probability bound from a
support pattern that necessarily leaves a nonzero message. The following
rank count provides such a sufficient test, without changing the code.

Fix the actual route, pairings, and any proposed state support masks H_t,
starting at H_0=empty. The masks need not be exact supports of a witness.
Let F_t be all coordinates of pairs meeting H_t. Set

    V_t=F_t union H_(t+1),
    J_t=F_t minus (H_t union H_(t+1)).

Use the true initial union here, not the padded initial batch from the
probability bound. Let D count all deaths, A_state=sum_t |H_(t+1)|, and
J=sum_t |J_t|. Let c_i count row i's routed coordinates in V. The RS
outer shortened to V has dimension d=sum_i(c_i-128)_+.

Introduce d input variables and A_state state-amplitude variables.
The recurrence has potentially nonzero equations only on V. Counting
coordinates gives

    sum_t |V_t| = A_state+J+D.

Selecting z output-zero coordinates adds at most z equations. Hence the
solution space has dimension at least

    gap = d-J-D-z.

If gap>0, a nonzero message exists for EVERY realization of the local
arithmetic with this pairing locality. Indeed, a zero input forces every
state to zero by the recurrence, so a nonzero solution cannot consist
only of auxiliary state variables. Its output is supported on the union
of consecutive masks, minus the selected zeros. A positive gap with this
remaining support at most W would therefore be an actual distance
obstruction, not merely a failed first-moment bound.

No full-size positive gap has been found or searched for in this pass.
Next run a bounded search for this gap using both route AND pairings.
The earlier route-only search did not include the pairings. A negative
search result will not prove distance, but the attained gaps and support
geometry can identify which constraints a sharper proof must retain.

## Checks and scope

Exact tests pass12 marked-row moment comparisons and10 tiny uniform-image
transfer comparisons. The latter use q=2 and4 and explicitly do not test
the accepted SL2 distribution. A separate structural-rank test passes3888
equation-count identities and669 positive-gap witness checks using small
repetition outers and invertible local GF2 gates. That argument does not
require the accepted-map distribution, but the toy examples are not
full-size counterexamples.

All numerical jobs ran serially under the20GiB guard, in under2s each,
without raw files. No construction, native kernel, paper, prior certificate,
or frozen helper was changed. Whole-code closure and timing remain open.
