# A crossing burst defeats the exact paid-set moment

The regional transfer refinements do not address the final obstruction.
The exact moment itself is too large for the current proof interface.
This note proves that statement under the actual regional quotas, without
histogram envelopes, artificial terminal resets, or independent-bit models.
It does not prove that the encoder has low minimum distance.

Keep the unchanged quartet construction and the middle class r=2100,
c=111572. There are256 regions of8192 inputs, each containing16 updates
of width512. The original inactive set has R=537600 positions.

For a paid set P, let H(P) maximize lambda*A-tau*D over its allowed state
supports, with the actual empty initial state and free final state. The
existing proof interface attempts to bound

    U = binom(8192,2100) * b^(-B(c)) * Z_c(b),
    Z_c(b) = E sum_(P subset I, |P|=c) b^H(P),  b>=1.

The expectation uses the actual routing and quartet partitions. The field
constraints give, for mu,nu,eta>=0,

    lambda=39/512+mu+nu, tau=381/512+eta,
    B(c)=(1+mu+eta)c+713983mu+743225nu-255463eta-111127.

## A feasible support across a regional boundary

Pair consecutive regions. Require the last seven updates of the first
region and the first seven updates of the second to have all inputs
allowed. This gives a fourteen-update full-input burst across their common
boundary. Select all512 state coordinates in its first thirteen updates.
Omit its final update and charge the resulting512 deaths.

Every selected state coordinate has its own input and its next quartet
inside the allowed burst. Thus the support is feasible for every quartet
partition. Repeat this construction in all128 region pairs. The total is

    A*=13*512*128=851968, D*=512*128=65536,
    H(P)>=S:=lambda*A*-tau*D*.

There is no state reset at a regional boundary. All128 deaths are charged,
including the last burst's death before the actual end of the stream.

For c_j payments in one region, its unpaid set has2100-c_j positions.
The expected number of paid sets satisfying the prescribed seven-update
condition is

    g(c_j)=binom(2100,c_j)
             * binom(4608,2100-c_j)/binom(8192,2100-c_j).

Assign435 payments to44 regions and436 to212 regions. These assignments
are disjoint as paid-set classes. Independence of the regional permutations
therefore gives event weight

    G=binom(256,44) * g(435)^44 * g(436)^212.

This counts paid sets under the actual distribution. It is not a probability
of finding a codeword, and G need not be a probability in general.

## No scalar moment or current field multipliers repair the expression

Put A0=binom(R,c). Since H>=0, Z_c(b)>=A0. The burst event also gives
Z_c(b)>=G*b^S. Define

    L0=713983+c=825555,
    K0=(39/512)L0-(c-111127), M=255463-c,
    rho=L0/A*=825555/851968.

We have

    B=lambda*L0-K0-nu*(L0-743225)-eta*M,
    B-rho*S = -K0+rho*(381/512)*D*
                -nu*(L0-743225)-eta*(M-rho*D*) <= 0.

The last inequality follows from exact positive rational comparisons.
It does not impose an onset-three condition such as2lambda<=tau.

If B<=0, then b^(-B)Z_c(b)>=A0. Otherwise S>0 and
theta=B/S lies in(0,rho]. Interpolating the two lower bounds on Z gives

    b^(-B)Z_c(b) >= A0^(1-theta) G^theta
                   >= A0^(1-rho) G^rho,

where the second inequality uses G<A0. Consequently, for EVERY b>=1 and
every mu,nu,eta>=0 in this field-constraint family,

    log2 U >= 37010890379304638123/72057594037927936
             >513.62928.

This is a directed lower bound on the proposed union-bound expression.
It rules out closing that expression by exact histogram enumeration,
finer buckets, better boundary weights, different moment markers, or the
current field multipliers. It does not rule out a different proof event,
additional field constraints, better row-set counting, or the construction.

## A stronger obstruction at the currently selected weights

Use eight updates on each side of the paired boundary instead. The support
has A*=983040,D*=65536. At lambda=707/2560,tau=381/512, its score is222720,
whereas B=165556. Directed calculations give:

| Scope | Lower bound on log2 expression, before the row-set union |
| --- | ---: |
| Current weights, b=16 | greater than118288.50735 |
| Current weights, every b>=1 | greater than19615.43430 |
| All current multipliers with2lambda<=tau, every b>=1 | greater than5087.31772 |

For the third line, B/S is at most250598845/324599808. This maximum occurs
at lambda=(381/512)/2,eta=nu=0. Exact signs of the two affine coefficients
prove the maximum; it is not inferred from a parameter sweep. The seven-
update argument above is weaker numerically but removes the onset restriction.

`quartet_crossing_burst.py` verifies28 exact paid-event counts, three support
sizes, the rational multiplier inequalities, and directed logarithms using
the existing integer/rational log helper. Runtime0.46s, serial remote,
90s watchdog, less than64MiB, no generated files. The directed fractions
are retained in the reproducible source output, not rounded into premises.

A separate preliminary helper,`quartet_geometry_lower.py`, checks3130
exact one-coordinate mean identities. Its Jensen calculation sums241/256
selected count paths per region, dropping cross-boundary rewards. The
FLOAT lower expression is-14680.28881 and therefore shows no obstruction.
The direct crossing-burst argument retains precisely the boundary rewards
that this simpler test discards; the negative Jensen result is preserved,
not presented as evidence that the moment might close.

## Consequence for the ongoing proof

Do not continue refining the same all-paid-set moment as a route to closure.
The [arbitrary-pair checkpoint](QUARTET_ARBITRARY_PAIR.md) remains valid,
but its transfer improvements cannot overcome this exact obstruction.

The paid-set moment counts supersets containing payments that its maximizing
support never uses. A [closed-paid-set reduction](QUARTET_CLOSED_PAID_SETS.md)
removes that redundancy without changing the encoder. It gives a different
regional counting factor and is the next proof interface to test. The full
45%/40-bit certificate remains open; prior certified exclusions are preserved.
