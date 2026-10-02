"""An exact uniform-input envelope for independently randomized MDS symbols.

Let C be a linear [n,k] MDS code over GF(Q), with Q a power of two. At every
symbol position sample an independent invertible binary map that sends each
fixed nonzero symbol uniformly onto the Q-1 nonzero symbols. Uniform GL maps
satisfy this premise. Define mu(x) as the expected number of nonzero messages
whose randomized codeword equals x. Setup is sampled once; no independence
between different messages is needed.

For a fixed support S of size h, shortening C to S gives dimension h-(n-k)
when h>=d=n-k+1. Choose that many information coordinates of the shortened
code. Projection onto them is injective, and an exact-support word is nonzero
at every projected coordinate. Thus at most (Q-1)^(h-(n-k)) words have exact
support S. Independent symbol randomization sends each such word to a fixed
x with probability (Q-1)^-h. Consequently

    mu(x) <= (Q-1)^-(n-k)   for every nonzero x.

For 0<h<d the shortened code is zero, so mu(x)=0. The zero message is excluded
and invertible maps preserve zero, so mu(0)=0 as well. If U is uniform on all
n*log2(Q) binary output coordinates, this proves measure domination

    mu <= beta * Law(U),   beta = Q^n / (Q-1)^(n-k).

The right side deliberately includes artificial zero-vector mass. An
independent subsequent permutation preserves the domination. Products across
groups require independent group setups. This is a first-moment count bound,
not a distance certificate or a bound on every sampled code's enumerator.

Splitting U into R packets of b bits gives independent packet activity
p=(2^b-1)/2^b and uniform nonzero labels conditional on activity. Therefore
the expected nonzero-message support shells satisfy

    A_v <= C(R,v)*(2^b-1)^v / (Q-1)^(n-k),  v>0.

Use beta and activity for the full IID comparison, or shell_caps() for the
sharper zero-excluded shell envelope. Their total masses intentionally differ.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import comb

from rs_outer import expected_group_support_counts, mds_symbol_weight_counts


@dataclass(frozen=True)
class UniformInputEnvelope:
    """Exact scalar constants and packet-shell caps for one MDS outer group.

    Defaults describe four parallel GF256 RS[8,4] rows with independent GL32
    symbol maps. The other K16 candidate uses n=16, k=8, packet_bits=4 and
    packets_per_symbol=4 (four GF16 rows with independent GL16 maps). Generic
    values describe an MDS ensemble over Q, not a claim that a particular
    parallel base-field implementation exists for those parameters.
    """

    n: int = 8
    k: int = 4
    packet_bits: int = 4
    packets_per_symbol: int = 8

    def __post_init__(self):
        for value in (self.n, self.k, self.packet_bits, self.packets_per_symbol):
            if type(value) is not int or value < 1:
                raise ValueError('positive integer MDS and packet parameters required')
        if not 1 <= self.k <= self.n <= self.q:
            raise ValueError('require 1 <= k <= n <= Q for an evaluation RS code')

    @property
    def symbol_bits(self) -> int:
        return self.packet_bits * self.packets_per_symbol

    @property
    def q(self) -> int:
        return 1 << self.symbol_bits

    @property
    def regions(self) -> int:
        return self.n * self.packets_per_symbol

    @property
    def message_bits(self) -> int:
        return self.k * self.symbol_bits

    @property
    def output_bits(self) -> int:
        return self.n * self.symbol_bits

    @property
    def minimum_symbol_weight(self) -> int:
        return self.n - self.k + 1

    @property
    def pointwise_nonzero_cap(self) -> Fraction:
        return Fraction(1, (self.q - 1) ** (self.n - self.k))

    @property
    def beta(self) -> Fraction:
        """Coefficient multiplying the probability law of the FULL uniform input."""
        return (1 << self.output_bits) * self.pointwise_nonzero_cap

    @property
    def activity(self) -> Fraction:
        """Packet activity in that IID comparison, including its possible zero vector."""
        return Fraction((1 << self.packet_bits) - 1, 1 << self.packet_bits)

    def shell_caps(self, *, include_zero: bool = False) -> tuple[Fraction, ...]:
        """Caps indexed 0..R; the default sets shell zero to its exact value zero.

        With include_zero=True these are precisely beta times Binomial(R,p),
        hence sum to beta. Otherwise their sum is beta-pointwise_nonzero_cap.
        They are upper bounds, not the exact expected support counts.
        """
        if type(include_zero) is not bool:
            raise ValueError('include_zero must be Boolean')
        density = self.pointwise_nonzero_cap
        labels = (1 << self.packet_bits) - 1
        return tuple(density * comb(self.regions, v) * labels ** v
                     if v or include_zero else Fraction(0)
                     for v in range(self.regions + 1))

    def exact_pointwise_symbol_counts(self) -> tuple[Fraction, ...]:
        """Expected multiplicity of each particular vector, indexed by symbol weight.

        MDS symmetry gives A_h/C(n,h) exact-support words for every fixed
        symbol support. GL randomization makes each nonzero labeling equally
        likely. Entry zero excludes the zero message.
        """
        counts = mds_symbol_weight_counts(self.q, self.n, self.k)
        return (Fraction(0),) + tuple(
            Fraction(counts[h], comb(self.n, h) * (self.q - 1) ** h)
            for h in range(1, self.n + 1))

    def verify_shell_domination(self) -> None:
        """Check the lemma against the independent exact shell-count generator."""
        exact = expected_group_support_counts(self.n, self.k,
            self.packet_bits, self.packets_per_symbol)
        caps = self.shell_caps()
        if any(value > cap for value, cap in zip(exact, caps)):
            raise ArithmeticError('uniform envelope does not dominate exact expected shells')
        symbol_counts = self.exact_pointwise_symbol_counts()
        if any(value > self.pointwise_nonzero_cap for value in symbol_counts):
            raise ArithmeticError('uniform envelope does not dominate pointwise symbol counts')

    def metadata(self) -> dict:
        """JSON-ready premises and exact rational constants; no numerical endpoint."""
        return dict(schema='mds-uniform-input-envelope-1', n=self.n, k=self.k,
            alphabet_size=self.q, symbol_bits=self.symbol_bits,
            packet_bits=self.packet_bits, packets_per_symbol=self.packets_per_symbol,
            regions=self.regions, message_bits=self.message_bits, output_bits=self.output_bits,
            minimum_symbol_weight=self.minimum_symbol_weight,
            pointwise_nonzero_cap=str(self.pointwise_nonzero_cap),
            beta=str(self.beta), iid_packet_activity=str(self.activity),
            actual_zero_message_excluded=True, iid_envelope_includes_artificial_zero=True,
            expected_count_bound=True, deterministic_enumerator_bound=False,
            whole_code_certificate=False,
            premises=['linear MDS[n,k] over the declared alphabet',
                'independent invertible binary symbol maps with uniform nonzero fixed-input images',
                'independent group setups when multiplying group measures',
                'any later permutation is independent of the randomized symbols'])
