# Memory and termination argument for the BMS connection

## Outcome

The memory/distance claim holds with an explicit terminated streaming
realization. There is no distance loss from buffering: the original
codeword is a prefix of the new codeword. Output length increases by one
outer block, O(log N).

Accumulator SPIN gives logarithmic outer memory and a one-bit,
time-invariant inner. Structured SPIN has constant recursive data state;
its position-dependent mixer schedule requires either time-varying
transitions or an O(log N) internal position counter. The latter still
gives logarithmic outer memory and sublinear inner memory.

This answers the terminated-encoder version of the BMS question. It is
not an assertion that the unmodified map exactly matches their displayed
unterminated, fixed-output-count interface, nor a historical priority claim.

## Explicit streaming outer

Fix a realized block-diagonal outer O with M blocks, input block size
k=b/2, output block size b, message length K=Mk, and N=Mb=2K.
Write its block maps as G_j and its input blocks as x_j.

Keep a previous-output buffer P of b bits and storage for the current
block calculation. Initially P=0. During input block j:

1. Read one message bit at each transition.
2. Emit the next two bits of P at each transition.
3. Accumulate or buffer the current block input.
4. After k transitions, finish G_j(x_j) and use it as the next P.

The outer emits the all-zero buffer while reading x_1, G_1(x_1) while
reading x_2, and so on. Append k fixed zeros as termination input and
emit G_M(x_M) while reading them. The resulting output is exactly

    O_stream(x || 0^k) = 0^b || O(x).

The terminated domain has K free message bits, not K+k. The additional
inputs are fixed, not fresh message coordinates. This distinction is
essential to the rate and injectivity claims.

For arbitrary fixed linear block maps, one implementation maintains a
b-bit partial sum, XORing the appropriate generator vector for each
incoming 1 bit. With P, this uses at most 2b data bits, plus counters.
It takes O(b) bit work per input bit, agreeing with the dense outer cost.
For the structured outer, buffering k input bits and evaluating the
O(b)-work constituent at the block boundary uses O(b) scratch and retains
O(b) total work per block. The streaming interface imposes no worst-case
constant-computation latency on a transition; aggregate encoding work is
what the paper measures.

A phase counter and block-position counter need O(log N) bits. For a
fixed realized code, local matrices/permutations are part of the encoder
description. This is the usual distinction between state size and
transition-description size; it does not assert that total code storage
is O(log N). A counter makes a position-dependent outer schedule into
one finite-state transition rule with O(log N) extra state bits.

## One interleaver and preservation of distance

For the original permutation Pi on N coordinates, define a permutation
Pi_ext on N+b coordinates by

    Pi_ext(0^b || z) = Pi(z) || 0^b.

This specifies the coordinate permutation independently of z. Continue
the original causal inner on the b extra zero bits. For a mixer
schedule, preserve every original mixer and choose any fixed
continuation (for example identity maps). A partial final inner group can be
serialized and stopped after the required output bits.

Then for every message x,

    E_ext(x) = E(x) || tail(x).

All maps remain linear. Consequently E_ext is injective whenever E is,
and d_min(E_ext) >= d_min(E). This is a prefix identity, not an argument
that arbitrary shortening or puncturing preserves distance.

Writing N'=N+b gives

    rate(E_ext) = N / (2(N+b)) = 1/2 - O(log N/N),
    d_min(E_ext)/N' >= delta N/(N+b) = delta - O(log N/N)

whenever d_min(E)>delta N and b=Theta(log N). The original good-setup
event implies the extended one; no new failure event is introduced.
Interleaving is still one permutation, though Pi_ext is deliberately
not uniform on all (N+b)! permutations. The BMS existence question does
not require a uniform interleaver.

For an accumulator, tail(x) is simply b copies of the accumulator's
state after its original N inputs. Its recursive state is still one bit.

## Serializing the structured inner

The native recurrence is

    Y_i = X_i + A Q_i,
    Q_(i+1) = M_i Q_i + C X_i,

with t=128 and s=19. Hold Q_i fixed while reading the t input bits.
On input bit x_j, emit x_j+(A Q_i)_j and accumulate x_j C(e_j) into
an s-bit register P. After t inputs, set Q=M_i Q+P and clear P.

This is exactly the native map, with one input and one output per
transition, two s-bit data registers, and a constant-size within-group
phase counter. Temporary fixed-width arithmetic registers do not change
the O(1) state or work-per-bit bound. No growing block buffer is needed.

The mixer M_i depends on the group index. With an external
time-dependent schedule, the data state remains constant. With the
position counter included in a time-invariant encoder's state, total
memory is O(log N). Thus “19-bit recursive state” is accurate;
“19 bits of all-inclusive autonomous machine state” is not.

## Exact scope relative to BMS

BMS Section I-A defines constituent encoders with one input bit and a
constant number of output bits per transition. It then notes that
termination inputs are often appended and sets that detail aside.
The construction above states those inputs and their effect explicitly.
It is therefore appropriate to say “a positive answer with standard
termination” or “in the terminated-encoder formulation.” Do not claim
that their literal all-inputs-free, unterminated displayed map is
identical to this terminated realization.

BMS's explicit definition of time-varying automata occurs in the
repeat-convolute discussion. The accumulator-based positive result
does not need a time-varying inner. For the structured result, the
counter construction avoids relying on an unstated extension of their
serial-concatenation definition; the distinction between recursive
data state and total control state remains visible.

No result here provides efficient decoding, a bound on total stored
setup data, or an exact-GV limit beyond the existing SPIN theorems.

## Priority and next step

A bounded search found the BMS question and related serial-concatenation
work, but does not establish that this is the first positive answer.
The truncated Block-Accumulate architecture also requires its existing
attribution. Audit the BA paper and subsequent growing-memory literature
before elevating this to a historical open-problem-resolution claim.

The mathematical model correspondence is now explicit. The remaining
editorial decision is whether the short terminated-encoder qualification
belongs in the main sentence or a technical note. Contributions have not
been changed.
