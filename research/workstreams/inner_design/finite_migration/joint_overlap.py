"""Exact joint-map distance audit and a narrower Bernoulli overlap cap."""
from functools import lru_cache
from itertools import combinations

from flint import arb
import activation_density
import short_dense

model = activation_density.model


def basis(rows):
    pivots = {}
    for row in rows:
        assert type(row) is int and row >= 0
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return pivots


def remainder(word, pivots):
    for pivot in sorted(pivots, reverse=True):
        if word >> pivot & 1:
            word ^= pivots[pivot]
    return word


def distance_audit(a_rows, b_rows, length, radius=3):
    """Unique syndromes through radius r exclude nonzero weights <= 2r."""
    assert type(length) is int and length > 0
    assert type(radius) is int and 1 <= radius <= length
    assert all(0 <= row < 1 << length for row in a_rows + b_rows)
    a_basis, b_basis = basis(a_rows), basis(b_rows)
    joint = basis(a_rows + b_rows)
    assert len(a_basis) == len(a_rows) and len(b_basis) == len(b_rows)
    assert len(joint) == len(a_rows) + len(b_rows), 'Map images intersect'
    assert remainder((1 << length) - 1, b_basis) == 0, 'All-one word absent from B image'
    syndromes = [remainder(1 << j, joint) for j in range(length)]
    seen = {0: 0}
    for weight in range(1, radius + 1):
        for support in combinations(range(length), weight):
            syndrome, mask = 0, 0
            for j in support:
                syndrome ^= syndromes[j]
                mask |= 1 << j
            assert syndrome not in seen, ('Low-weight joint word', mask ^ seen.get(syndrome, 0))
            seen[syndrome] = mask
    return dict(length=length, rank=len(joint), minimum_distance_lower=2 * radius + 1,
                checked_syndromes=len(seen), disjoint_images=True, all_one_in_b=True)


@lru_cache(maxsize=1)
def audit():
    a, b, *_ = activation_density.ladder.candidate.maps('weight5_seed0')
    generators = model.independent.g.tv.fixed.maps.generators
    return distance_audit(generators(a, 19), generators(b, 19), 128)


def overlap_ends(v, w, length, distance):
    # q != 0 excludes both Aq+B^Ta=0 and Aq+B^Ta=1 because the
    # images are disjoint and 1 belongs to the B image.
    lower = max(0, v + w - length, (v + w - length + distance + 1) // 2)
    upper = min(v, w, (v + w - distance) // 2)
    assert lower <= upper
    return lower, upper


class Engine(activation_density.Engine):
    def __init__(self, exponent):
        super().__init__(exponent)
        self.joint_audit = audit()

    def bernoulli(self, theta, lam):
        assert theta > 0 and theta < 1 and lam > 0
        # Preserve the zero-state row, whose exact kernel contribution does
        # not depend on any A/B overlap estimate.
        result = list(super().bernoulli(theta, lam))
        z = (-lam).exp()
        g0, g1 = 1 - theta + theta * z, theta + (1 - theta) * z
        rho0 = min(arb(1), model.up(abs(1 - 2 * theta * z / g0)))
        rho1 = min(arb(1), model.up(abs(1 - 2 * (1 - theta) * z / g1)))
        distance = self.joint_audit['minimum_distance_lower']
        entries = []
        for v in self.levels:
            cap = arb(1)
            for w, count in self.b_spectrum.items():
                ends = overlap_ends(v, w, 128, distance)
                cap += count * max(model.up(rho0**(w-h) * rho1**h) for h in ends)
            moment = g0**(128-v) * g1**v
            entries.append(model.up(moment * (cap / (2*(self.m+1)) + arb(1)/(2*self.m))))
        for i, entry in [(1, max(entries))] + [(j+2, e) for j, e in enumerate(entries)]:
            result[i*self.n] = entry
            for j, w in enumerate(self.levels):
                result[i*self.n+j+2] = model.up(entry*self.spectrum[w])
        return tuple(map(model.up, result))


class Checker(short_dense.Checker):
    def __init__(self, exponent):
        super().__init__(exponent)
        self.engine = Engine(exponent)
        assert not self.fixed_moment.cache_info().currsize
