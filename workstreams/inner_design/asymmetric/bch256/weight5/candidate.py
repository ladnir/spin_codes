"""Candidate-specific exact maps and transfers; no frozen producer is changed."""
from collections import Counter
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import bch_model as model
import dense_search


@lru_cache(maxsize=3)
def maps(name):
    record = model.base.read(model.ASYMMETRIC.parent/'NO_CONSTANT_MAP.json')
    a = list(record['columns'])
    if name == 'balanced':
        b = list(a)
    else:
        assert name in ('weight5_seed0', 'weight5_seed1')
        pool = [sum(1 << j for j in support) for support in itertools.combinations(range(19), 5)]
        b = random.Random(int(name[-1])).sample(pool, 128)
    search = model.independent.search
    assert len(b) == len(set(b)) == 128 and search.rank(b) == 19
    assert all(0 < c < 1 << 19 for c in b)
    exact = model.independent.g.tv.fixed.maps
    sa = exact.spectrum(exact.generators(a, 19))
    sb = exact.spectrum(exact.generators(b, 19))
    assert sa == {int(w): int(n) for w, n in record['spectrum'].items()}
    kernel = exact.dual_spectrum(sb, 128, 19)
    assert sum(sb.values()) == 1 << 19 and sum(kernel.values()) == 1 << 109
    low = search.low_kernel(b)
    assert all(kernel.get(j, 0) == low[f'weight{j}'] for j in range(1, 5))
    spectrum = {w: n for w, n in sb.items() if w}
    k = [kernel.get(j, 0) for j in range(129)]
    caps = model.independent.g.fiber_caps(128, 19, spectrum, k)
    cancellation = search.low_cancellation(a, b, 19)
    assert all(cancellation[j]['nonzero_input_count'] == math.comb(128, j) for j in (1, 2))
    return a, b, {w: n for w, n in sa.items() if w}, spectrum, k, caps, cancellation


class Engine(model.Engine):
    def __init__(self, exponent, name='weight5_seed0'):
        assert exponent in (16, 18, 20)
        (self.a_columns, self.columns, self.spectrum, self.b_spectrum,
         self.kernel, self.caps, self.low) = maps(name)
        self.name = name
        self.levels = sorted(self.spectrum)
        self.n = len(self.levels)+2
        self.m = (1 << 19)-1
        self.exponent = exponent
        self.length = 1 << (exponent-7)
        self.outer_length = 256
        self.output_bits = 256*self.length
        self.cutoff = self.output_bits//10

    def identity(self):
        result = super().identity()
        result['inner'] = dict(t=128, s=19, transvection_rounds=1, feedback_name=self.name,
            expansion_columns=self.a_columns, feedback_columns=self.columns,
            feedback_sha256=hashlib.sha256(json.dumps(self.columns, separators=(',', ':')).encode()).hexdigest())
        return result

    def audit(self):
        return dict(name=self.name, dual_minimum_weight=min(self.b_spectrum),
                    kernel_low=self.kernel[:5], column_weights={str(w): n for w, n in Counter(c.bit_count() for c in self.columns).items()})


class Checker(dense_search.Checker):
    def __init__(self, exponent, name='weight5_seed0'):
        # The parent sets only outer scalars and empty caches during initialization.
        super().__init__(exponent)
        self.engine = Engine(exponent, name)
        assert self.rows == self.engine.length and self.cutoff == self.engine.cutoff
        assert not self.epoch.cache_info().currsize and not self.fixed_moment.cache_info().currsize


def sources():
    return model.sources()
