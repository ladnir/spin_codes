# Scaled feedback for larger GF16 states

This checkpoint studies the existing GRS[16,8] outer and single regional
permutation, with a larger inner state. It is a proof screen, not a new
implementation or a whole-code certificate. The current certified profiles
are unchanged.

## Construction and exact algebra

Let GF16 use the polynomial basis modulo X^4+X+1. For state dimension s,
let A have rows (1,j,...,j^(s-1)), one for each j in GF16, and take
C = gamma A^T with gamma=2. An inner step remains

    R = M Q,   Y = X + A R,   Q_next = R + C X,

where M is independent and uniform in GL(s,16). Independent nonzero input
diagonals and the sixteen independent regional symbol permutations are
unchanged. The state does not reset between regions.

For 1 <= s <= 15, the (i,j) entry of A^T A is one exactly when i+j=15,
and zero otherwise. The nonzero part is a collection of disjoint swap
matrices. Therefore I+gamma A^T A is invertible whenever gamma != 1:
each two-dimensional block has determinant 1+gamma^2 != 0. The identity
relating I+AC and I+CA then gives invertibility of B=I+AC. State16 is
excluded because degree30 contributes another nonzero power sum.

Scaling C does not change its kernel or the MDS birth counts. It does not
change the emission distribution after a uniform nonzero refresh. It does
change the output on return to zero. The model recomputes the singleton
return polynomial directly from the sixteen columns of B and retains the
invertible-B return bound for larger supports.

## Numerical screen

At K=2^20, rate 1/2, and distance 10%, the sharper regional coefficient
screen gives the following per-class margins. All use the exact outer
support multiplicities through the positive eta/rho envelope, with
rho=871/1000. These are floating upper-bound evaluations, not directed
whole-code certificates and not estimates of actual bad-code probability.

| State | 128 active outer rows | 256 active outer rows |
| ---: | ---: | ---: |
| 8 | -98.58 bits | -116.56 bits |
| 9 | -66.63 bits | -84.66 bits |
| 10 | -34.81 bits | -53.00 bits |
| 12 | 29.29 bits | 11.17 bits |
| 13, additional screen | 61.47 bits | 43.51 bits |
| 14, additional screen | 93.41 bits | 75.38 bits |

Thus the requested states through twelve do not establish 40 bits with
this proof envelope. States thirteen and fourteen show that enlarging the
state can overcome the sampled middle-class obstruction, but no coverage
of every class or summed 40-bit result is claimed. In particular, the
state13 value at a single class leaves little room for the union bound.
Performance has not been measured for any scaled-feedback profile.

The cheaper coefficient-Chernoff relaxation is substantially weaker here:
at 256 active rows it gives approximately -371.65 bits for states nine
through twelve. It cannot replace regional coefficient extraction in this
part of the proof.

## Checks and reproduction

`theory/large_state_screen.py` keeps this model separate from the preserved
certificates. It includes directed local envelopes and directed replay for
individual classes, but its search commands do not certify the whole code.

Eight tests in `theory/test_large_state_screen.py` pass. They check every
state from one through fifteen and every nonzero GF16 feedback scale,
direct singleton return outputs, agreement with the unchanged gamma=1
model, exact rational enclosure of the local envelopes, and exhaustive
input/nonzero-state enumeration for the analogous GF4 construction with
nonzero CA. The tests include invalid-parameter guards.

From the theory directory:

    python -m unittest test_large_state_screen -v
    python large_state_screen.py --state 8 9 10 12 --active 128 256 --mode sparse

The next useful comparison is a cheaper construction with a complete
certificate, rather than implementing a large state on the strength of
these isolated positive screens.
