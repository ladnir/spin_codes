"""Exact local checks, not a minimum-distance certificate. No generated files."""
import itertools
import random
import time


class Field:
    def __init__(self, degree, modulus):
        self.size = 1 << degree
        def multiply(a, b):
            c = 0
            while b:
                if b & 1:
                    c ^= a
                a <<= 1
                if a & self.size:
                    a ^= modulus
                b >>= 1
            return c
        self.table = [[multiply(a, b) for b in range(self.size)]
                      for a in range(self.size)]
        self.inverse = [0] + [next(b for b in range(1, self.size)
                                  if self.table[a][b] == 1)
                              for a in range(1, self.size)]
        basis = [1]
        for _ in range(1, degree):
            basis.append(next(a for a in range(self.size)
                              if self.table[a][a] ^ a == basis[-1]))
        self.points = [self.xor(basis[j] for j in range(degree) if i >> j & 1)
                       for i in range(self.size)]
        assert len(set(self.points)) == self.size

    @staticmethod
    def xor(values):
        result = 0
        for value in values:
            result ^= value
        return result

    def delta(self, x, rounds):
        for _ in range(rounds):
            x = self.table[x][x] ^ x
        return x

    def fft(self, values, offset=0, inverse=False):
        if len(values) == 1:
            return values[:]
        half = len(values) // 2
        c = self.delta(self.points[offset], half.bit_length() - 1)
        if inverse:
            lo = self.fft(values[:half], offset, True)
            hi = self.fft(values[half:], offset ^ half, True)
            high = [a ^ b for a, b in zip(lo, hi)]
            low = [a ^ self.table[c][b] for a, b in zip(lo, high)]
            return low + high
        lo = [a ^ self.table[c][b] for a, b in zip(values[:half], values[half:])]
        hi = [a ^ b for a, b in zip(lo, values[half:])]
        return self.fft(lo, offset) + self.fft(hi, offset ^ half)

    def interpolation_matrix(self, source, destination):
        rows = []
        for x in destination:
            row = []
            for j, xj in enumerate(source):
                numerator = denominator = 1
                for k, xk in enumerate(source):
                    if j != k:
                        numerator = self.table[numerator][x ^ xk]
                        denominator = self.table[denominator][xj ^ xk]
                row.append(self.table[numerator][self.inverse[denominator]])
            rows.append(row)
        return rows

    def apply(self, matrix, vector):
        return [self.xor(self.table[a][b] for a, b in zip(row, vector))
                for row in matrix]

    def rank(self, matrix):
        a = [row[:] for row in matrix]
        pivot = 0
        for col in range(len(a[0])):
            selected = next((i for i in range(pivot, len(a)) if a[i][col]), None)
            if selected is None:
                continue
            a[pivot], a[selected] = a[selected], a[pivot]
            scale = self.inverse[a[pivot][col]]
            a[pivot] = [self.table[scale][x] for x in a[pivot]]
            for i in range(pivot + 1, len(a)):
                scale = a[i][col]
                a[i] = [x ^ self.table[scale][y] for x, y in zip(a[i], a[pivot])]
            pivot += 1
            if pivot == len(a):
                break
        return pivot


class Maps:
    def __init__(self, field, state):
        self.f, self.s, self.b = field, state, 2 * state
        self.alpha = field.delta(field.points[self.b], state.bit_length() - 1)
        assert self.alpha not in (0, 1)
        y, z = field.points[:self.b], field.points[self.b:self.b + state]
        self.a_matrix = field.interpolation_matrix(z, y)
        self.b_matrix = field.interpolation_matrix(y, z)

    def a(self, q):
        coeff = self.f.fft(q, self.b, True)
        return self.f.fft(coeff + [0] * self.s)

    def bmap(self, y):
        coeff = self.f.fft(y, 0, True)
        folded = [a ^ self.f.table[self.alpha][b]
                  for a, b in zip(coeff[:self.s], coeff[self.s:])]
        return self.f.fft(folded, self.b)

    def kernel(self, v):
        # f=(s_log(s)+alpha)*v; its novel coefficients are [alpha*v, v].
        return self.f.fft([self.f.table[self.alpha][x] for x in v] + v)

    def a2(self, q):
        coeff = self.f.fft(q, self.b, True)
        return self.f.fft([self.f.table[1 ^ self.alpha][x] for x in coeff] + coeff)


def weight(v):
    return sum(x != 0 for x in v)


def check(degree, modulus, state, exhaustive):
    f = Field(degree, modulus)
    maps = Maps(f, state)
    rng = random.Random(20261021)
    for i in range(maps.b):
        y = [int(j == i) for j in range(maps.b)]
        assert maps.bmap(y) == f.apply(maps.b_matrix, y)
    for i in range(state):
        q = [int(j == i) for j in range(state)]
        assert maps.a(q) == f.apply(maps.a_matrix, q)
        assert maps.bmap(maps.a(q)) == q
        assert maps.bmap(maps.a2(q)) == q
    for _ in range(32):
        q = [rng.randrange(f.size) for _ in range(state)]
        y = [rng.randrange(f.size) for _ in range(maps.b)]
        assert maps.a(q) == f.apply(maps.a_matrix, q)
        assert maps.bmap(y) == f.apply(maps.b_matrix, y)
    print(f"GF(2^{degree}), b={maps.b}, s={state}: FFT/Lagrange and BA=I PASS", flush=True)
    a2_columns = [maps.a2([int(j == i) for j in range(state)]) for i in range(state)]
    combined = [row + [col[i] for col in a2_columns]
                for i, row in enumerate(maps.a_matrix)]
    assert f.rank(combined) == maps.b
    print("Complementary [A,A2] full rank and B*A2=I PASS", flush=True)
    if exhaustive:
        minima = [maps.b, maps.b, maps.b + state]
        for q in itertools.product(range(f.size), repeat=state):
            if not any(q):
                continue
            aq = maps.a(list(q))
            kv = maps.kernel(list(q))
            assert not any(maps.bmap(kv))
            minima = [min(minima[0], weight(aq)), min(minima[1], weight(kv)),
                      min(minima[2], weight(aq) + weight(q))]
        assert minima == [state + 1, state + 1, 2 * state + 1]
        minors = 0
        for r in range(1, state + 1):
            for rows in itertools.combinations(range(state), r):
                for cols in itertools.combinations(range(maps.b), r):
                    assert f.rank([[maps.b_matrix[i][j] for j in cols] for i in rows]) == r
                    minors += 1
        print(f"Exhaustive {f.size ** state - 1} states/kernel words: minima={minima}; {minors} B minors PASS", flush=True)
    for persistent in (False, True):
        q = [0] * state
        inputs, outputs, keys = [], [], []
        for _ in range(40):
            d = [rng.randrange(2, f.size) for _ in range(state)]
            x = [rng.randrange(f.size) for _ in range(maps.b)]
            feedback = maps.a([f.table[a][b] for a, b in zip(d, q)])
            y = [a ^ b for a, b in zip(x, feedback)]
            b_y = maps.bmap(y)
            q = [a ^ b for a, b in zip(q, b_y)] if persistent else b_y
            inputs.append(x); outputs.append(y); keys.append(d)
        q = [0] * state
        for x, y, d in zip(inputs, outputs, keys):
            feedback = maps.a([f.table[a][b] for a, b in zip(d, q)])
            assert x == [a ^ b for a, b in zip(y, feedback)]
            b_y = maps.bmap(y)
            q = [a ^ b for a, b in zip(q, b_y)] if persistent else b_y
    # Explicit live-state, zero-output transition; local, not an outer codeword.
    q = [1] + [0] * (state - 1)
    x = maps.a([2] + [0] * (state - 1))
    y = [a ^ b for a, b in zip(x, maps.a([2] + [0] * (state - 1)))]
    assert not any(y) and not any(maps.bmap(y))
    assert [a ^ b for a, b in zip(q, maps.bmap(y))] == q
    print("Both recurrence inverses and zero-output persistence fixture PASS", flush=True)
    # Two taps: exact support preservation on empty input, for sampled fixtures.
    for _ in range(32):
        q = [rng.randrange(f.size) for _ in range(state)]
        d = [rng.randrange(f.size) for _ in range(state)]
        e = [rng.choice([x for x in range(f.size) if x != (1 ^ di)]) for di in d]
        aq = maps.a([f.table[a][b] for a, b in zip(d, q)])
        a2q = maps.a2([f.table[a][b] for a, b in zip(e, q)])
        feedback = [a ^ b for a, b in zip(aq, a2q)]
        next_q = [a ^ b for a, b in zip(q, maps.bmap(feedback))]
        expected = [f.table[1 ^ a ^ b][c] for a, b, c in zip(d, e, q)]
        assert next_q == expected
        assert [bool(x) for x in next_q] == [bool(x) for x in q]
    print("Two-tap empty-input state support preservation PASS", flush=True)


if __name__ == "__main__":
    start = time.monotonic()
    check(4, 0x13, 4, True)
    check(8, 0x11b, 64, False)
    print(f"Local checks completed in {time.monotonic() - start:.3f}s; NO whole-code certificate.")
