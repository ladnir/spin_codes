"""Exact algebra prototype and operation counts, NOT an encoding benchmark.

Uniform feedback is unchanged. Novel-basis syndrome coordinates replace dense
global Vandermonde updates by one fixed transpose FFT and a Gray translation.
Run with -B. No output files, field tables, or unbounded searches.
"""
from collections import Counter
import itertools
import random
import time


class BinaryField:
    def __init__(self, degree, modulus):
        self.degree, self.modulus = degree, modulus

    def mul(self, a, b):
        out = 0
        while b:
            if b & 1:
                out ^= a
            a <<= 1
            if a >> self.degree:
                a ^= self.modulus
            b >>= 1
        return out

    def delta(self, x):
        return self.mul(x, x) ^ x

    def inverse(self, x):
        assert x
        out, power = 1, (1 << self.degree)-2
        while power:
            if power & 1:
                out = self.mul(out, x)
            x = self.mul(x, x)
            power >>= 1
        return out

    def chain(self):
        # Invert x -> x^2+x on its image, using binary column elimination.
        pivots = {}
        for i in range(self.degree):
            value, preimage = self.delta(1 << i), 1 << i
            while value:
                bit = value.bit_length() - 1
                if bit not in pivots:
                    pivots[bit] = value, preimage
                    break
                v, p = pivots[bit]
                value ^= v
                preimage ^= p
        result = [1]
        for _ in range(1, self.degree):
            target, x = result[-1], 0
            while target:
                v, p = pivots[target.bit_length() - 1]
                target ^= v
                x ^= p
            assert self.delta(x) == result[-1]
            result.append(x)
        # The chain starting at 1 is linearly independent (apply Delta powers).
        return result


def xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def point(bits, chain):
    return xor(x for i, x in enumerate(chain) if bits >> i & 1)


def novel_values(field, h, x):
    """p_j(x), with p_j=product Delta^i(x) over the set bits of j."""
    result = [1]
    for _ in range(h):
        result += [field.mul(v, x) for v in result]
        x = field.delta(x)
    return result


def translate(field, values, delta):
    """Lower-triangular moment translation L(delta), not polynomial shift L^T."""
    result = list(values)
    stride = 1
    while stride < len(result):
        if delta:
            for base in range(0, len(result), 2 * stride):
                for j in range(stride):
                    result[base + stride + j] ^= field.mul(delta, result[base + j])
        delta = field.delta(delta)
        stride *= 2
    return result


def transpose_evaluate(field, chain, values, offset=0):
    """Transpose of novel-basis evaluation on offset+V_h; exact scalar oracle."""
    if len(values) == 1:
        return list(values)
    h = len(values).bit_length() - 1
    half = len(values) // 2
    left = transpose_evaluate(field, chain, values[:half], offset)
    right = transpose_evaluate(field, chain, values[half:], offset ^ chain[h-1])
    c = offset
    for _ in range(h-1):
        c = field.delta(c)
    low = [a ^ b for a, b in zip(left, right)]
    high = [b ^ field.mul(c, a) for a, b in zip(low, right)]
    return low + high


def dense_syndrome(field, chain, values, offset):
    h = len(values).bit_length() - 1
    result = [0] * len(values)
    for j, value in enumerate(values):
        column = novel_values(field, h, offset ^ point(j, chain))
        for i, coefficient in enumerate(column):
            result[i] ^= field.mul(coefficient, value)
    return result


def coefficient_bits(cantor_index):
    """Minimum subfield degree for a nonzero combination with this top index."""
    return 1 << cantor_index.bit_length()


def operation_counts(h, log_n=21):
    """Exact scalar coefficient-product sites; omit 0 and 1 products."""
    b, n = 1 << h, 1 << (log_n-h)
    fixed, translation, direct = Counter(), Counter(), Counter()
    for i in range(h):
        # Fixed FFT: one node for each upper-bit prefix; 2^i butterflies/node.
        for prefix in range(1, 1 << (h-1-i)):
            top = prefix.bit_length()  # gamma_1,...,gamma_(h-1-i)
            fixed[coefficient_bits(top)] += (1 << i) * n
        for k in range(log_n-h+1):
            count = n // (1 << k) - n // (1 << (k+1))
            translation[coefficient_bits(h+k-i)] += count * (b//2)
        # Direct time-varying transpose FFT: nonzero quotient coset dominates
        # every within-batch prefix, hence all b/2 products at this level.
        for t in range(1, n+1):
            gray = t ^ (t >> 1)
            top = h + gray.bit_length()-1-i
            direct[coefficient_bits(top)] += b//2
    return dict(batch=b, batches=n, dense_products=(1 << log_n)*b,
                direct_fft=dict(sorted(direct.items())),
                fixed_fft=dict(sorted(fixed.items())),
                gray_translation=dict(sorted(translation.items())),
                combined=dict(sorted((fixed+translation).items())),
                total_products=sum(fixed.values())+sum(translation.values()),
                batches_with_GF32_transition=n//(1 << (16-h)))


def normalized_counts(field, chain, h, levels, log_n=21):
    """Diagonal state scaling: charge its reset cost, not just saved multiplies."""
    b, n = 1 << h, 1 << (log_n-h)
    multipliers = [field.inverse(chain[h-i]) if i in levels else 1 for i in range(h)]
    scales = [1]
    for multiplier in multipliers:
        scales += [field.mul(x, multiplier) for x in scales]
    classes = {0: 0, 1: 0}

    def category(value):
        if value not in classes:
            square = value
            for degree in range(1, 33):
                square = field.mul(square, square)
                if degree in (1, 2, 4, 8, 16, 32) and square == value:
                    classes[value] = degree
                    break
        return classes[value]

    count = Counter(operation_counts(h, log_n)['fixed_fft'])
    for scale in scales:
        if category(scale):
            count[category(scale)] += n
    for k in range(log_n-h+1):
        times = n//(1 << k)-n//(1 << (k+1))
        for i in range(h):
            coefficient = field.mul(chain[h+k-i], multipliers[i])
            if category(coefficient):
                count[category(coefficient)] += times*(b//2)
    return dict(levels=list(levels), combined=dict(sorted(count.items())),
                total_products=sum(count.values()))


def polynomial_mod(a, b):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a


def irreducible_power_two(field):
    # Degree is a power of two: only prime divisor is 2 (Rabin criterion).
    x = 2
    for i in range(1, field.degree+1):
        x = field.mul(x, x)
        if i == field.degree//2:
            a, b = field.modulus, x ^ 2
            while b:
                a, b = b, polynomial_mod(a, b)
            assert a == 1
    assert x == 2


def main():
    started = time.monotonic()
    deadline = started + 45
    rng = random.Random(20261022)
    fixtures = Counter()
    for degree, modulus, hs in ((4, 0x13, (1, 2)), (32, 0x10000008D, (3, 4, 8, 9))):
        field = BinaryField(degree, modulus)
        irreducible_power_two(field)
        chain = field.chain()
        for h in hs:
            b = 1 << h
            # Independent dense polynomial evaluations check transform and L orientation.
            for _ in range(3 if h < 5 else 1):
                values = [rng.randrange(1 << degree) for _ in range(b)]
                offset = rng.randrange(1 << degree)
                fixed = transpose_evaluate(field, chain, values)
                shifted = transpose_evaluate(field, chain, values, offset)
                assert translate(field, fixed, offset) == shifted
                assert translate(field, translate(field, values, offset), offset) == values
                if h < 5:
                    assert shifted == dense_syndrome(field, chain, values, offset)
                else:
                    # Check 4 independently chosen matrix columns at actual dimensions.
                    for j in (0, 1, b//2, b-1):
                        unit = [0]*b
                        unit[j] = 1
                        assert transpose_evaluate(field, chain, unit, offset) == novel_values(
                            field, h, offset ^ point(j, chain))
                fixtures['transform_and_translation'] += 1
                assert time.monotonic() < deadline
            # Full state recursion equals the global syndrome in moving coordinates.
            state = [0]*b
            syndrome = [0]*b
            previous = 0
            steps = min(12, (1 << (degree-h))-1)
            for t in range(steps):
                gray = (t+1) ^ ((t+1) >> 1)
                offset = point(gray << h, chain)
                output = [rng.randrange(1 << degree) for _ in range(b)]
                state = [a ^ c for a, c in zip(translate(field, state, offset ^ previous),
                         transpose_evaluate(field, chain, output))]
                contribution = transpose_evaluate(field, chain, output, offset)
                syndrome = [a ^ c for a, c in zip(syndrome, contribution)]
                assert state == translate(field, syndrome, offset)
                previous = offset
                fixtures['stream_state'] += 1
                assert time.monotonic() < deadline

    # The tempting one-translation recurrence has a genuine weight-two reset:
    # L(delta)e_(b-1)=e_(b-1). Removing C0 is NOT a coordinate change.
    field = BinaryField(4, 0x13)
    chain = field.chain()
    for h in (1, 2):
        b = 1 << h
        unit = [0]*(b-1)+[1]
        for delta in range(16):
            assert translate(field, unit, delta) == unit
            fixtures['removed_transform_weight_two'] += 1
    # Independent small-field rank routine; novel moments are MDS across time.
    from check_maps import Field
    rank_field = Field(4, 0x13)
    for h in (1, 2):
        b = 1 << h
        columns = [novel_values(field, h, x) for x in range(12)]
        for size in range(1, b+2):
            for chosen in itertools.combinations(range(12), size):
                matrix = [[columns[j][i] for j in chosen] for i in range(b)]
                assert rank_field.rank(matrix) == min(size, b)
                fixtures['global_reset_rank'] += 1
        # Couple arbitrary feedback keys under the invertible state change.
        for _ in range(16):
            state, syndrome, previous = [0]*b, [0]*b, 0
            for t in range((16//b)-1):
                matrix = [[rng.randrange(16) for _ in range(b)] for _ in range(b)]
                message = [rng.randrange(16) for _ in range(b)]
                output = [a ^ c for a, c in zip(message, rank_field.apply(matrix, syndrome))]
                # T'_t = T_t L(a_(t-1)); each full uniform matrix stays full uniform.
                recovered = [a ^ c for a, c in zip(output, rank_field.apply(
                    matrix, translate(field, state, previous)))]
                assert recovered == message
                offset = point(((t+1) ^ ((t+1) >> 1)) << h, chain)
                state = [a ^ c for a, c in zip(translate(field, state, offset ^ previous),
                         transpose_evaluate(field, chain, output))]
                syndrome = [a ^ c for a, c in zip(syndrome,
                            dense_syndrome(field, chain, output, offset))]
                assert state == translate(field, syndrome, offset)
                previous = offset
            fixtures['coupled_stream_inverse'] += 1
    assert time.monotonic() < deadline
    print(dict(kind='exact_algebra_checks_NOT_finite_certificate', checks=dict(fixtures),
               seconds=round(time.monotonic()-started, 3)), flush=True)
    for outer, h in ((256, 9), (512, 8)):
        print(dict(kind='operation_counts_NOT_timings', outer=outer, **operation_counts(h)), flush=True)
    field = BinaryField(32, 0x10000008D)
    chain = field.chain()
    for outer, h in ((256, 9), (512, 8)):
        for levels in (range(h), range(max(0, h-7), h)):
            # Exact coupling under the diagonal coordinate scaling.
            inv_norm = [field.inverse(chain[h-i]) if i in levels else 1 for i in range(h)]
            scales = [1]
            for value in inv_norm:
                scales += [field.mul(a, value) for a in scales]
            values = [rng.randrange(1 << 32) for _ in scales]
            delta = chain[h] ^ chain[h+1]
            expected = [field.mul(a, c) for a, c in zip(translate(field, values, delta), scales)]
            actual = [field.mul(a, c) for a, c in zip(values, scales)]
            coeff = delta
            for i in range(h):
                stride = 1 << i
                factor = field.mul(coeff, inv_norm[i])
                for base in range(0, len(actual), 2*stride):
                    for j in range(stride):
                        actual[base+stride+j] ^= field.mul(factor, actual[base+j])
                coeff = field.delta(coeff)
            assert actual == expected
            print(dict(kind='normalized_operation_counts_NOT_timings', outer=outer,
                       **normalized_counts(field, chain, h, levels)), flush=True)
    assert time.monotonic() < deadline
    print(dict(kind='all_exact_checks_PASS', normalized_coordinate_checks=4,
               total_seconds=round(time.monotonic()-started, 3)), flush=True)


if __name__ == '__main__':
    main()
