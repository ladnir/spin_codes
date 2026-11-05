"""Quiet-constraint rank bounds and bounded diagnostics; not a full certificate."""
from collections import deque
from fractions import Fraction
import math


def structural_rank(edges, source_dimensions, batch_capacity):
    """Integral max flow; edges[i][j] counts available coordinates, not bits."""
    left, right = len(source_dimensions), len(batch_capacity)
    total = left + right + 2
    source, sink = total - 2, total - 1
    residual = [[0] * total for _ in range(total)]
    for i, dimension in enumerate(source_dimensions):
        residual[source][i] = dimension
        for j, count in enumerate(edges[i]):
            residual[i][left + j] = count
    for j, capacity in enumerate(batch_capacity):
        residual[left + j][sink] = capacity
    value = 0
    while True:
        parent = [-1] * total
        parent[source] = source
        queue = deque([source])
        while queue and parent[sink] < 0:
            u = queue.popleft()
            for v, capacity in enumerate(residual[u]):
                if capacity and parent[v] < 0:
                    parent[v] = u
                    queue.append(v)
        if parent[sink] < 0:
            return value
        v, amount = sink, sum(source_dimensions)
        while v != source:
            amount = min(amount, residual[parent[v]][v])
            v = parent[v]
        v = sink
        while v != source:
            u = parent[v]
            residual[u][v] -= amount
            residual[v][u] += amount
            v = u
        value += amount


def nullity_moment(dimension, rank, tau, field_size):
    assert 0 <= rank <= dimension and tau >= 1 and field_size > 1
    tau = Fraction(tau)
    return tau ** (dimension-rank) * (1+(tau-1)/(field_size-1)) ** rank


def binomial_tail(n, threshold, denominator):
    return sum((Fraction(math.comb(n,j)*(denominator-1)**(n-j), denominator**n)
                for j in range(threshold,n+1)), Fraction(0))


def optimistic_mixed_screen():
    # Even perfect quiet rank cannot repair this particular forward envelope.
    # No route average, placement counts, or rank losses are charged here.
    N, b, s, n = 1 << 21, 128, 64, 16384
    W = 491*N//1000
    v, j, c, resets = n//2, n//2, 1, 0
    quiet_weight = (s+1)*j
    deficit = max(0,b*v-191*c-128*resets+quiet_weight-W)
    witness = (deficit+1)//2
    price = 64*resets + Fraction(107,128)*witness
    residual = N//2 - s*j
    exponent = residual - 1 - price + Fraction(1,64)
    assert residual==524288 and witness==275582 and exponent > 0
    print('Optimistic forward-envelope obstruction:',
          f'v={v}, j={j}, remaining dimension>={residual},',
          f'T={price}, q-exponent={exponent} > 0.', flush=True)
    print('This is a bound failure, not a low-distance codeword. '
          'Test output-side counting next for this dense class.', flush=True)


if __name__ == '__main__':
    optimistic_mixed_screen()
