"""Length-parameterized independent-map transfer; frozen producers stay intact."""
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ASYMMETRIC = HERE.parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ASYMMETRIC))
import certify as independent
sys.path.insert(0, str(ROOT/'workstreams/bch_rm2sub_bridge'))
import certificate_search_core as outer
from occupation_three import BANDS
from flint import arb, arb_poly

base = outer.base
up = independent.up
number = independent.old.number
pack = independent.old.pack
unpack = independent.old.unpack


def exact(value):
    m, e = map(int, value.upper().man_exp())
    return F(m)*F(2)**e


def sources():
    result = outer.outer_dependencies()
    for module in list(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path:
            path = Path(path).resolve()
            if path.is_relative_to(ROOT) and path.suffix == '.py':
                result[path.relative_to(ROOT).as_posix()] = base.sha(path)
    paths = [ASYMMETRIC.parent/'NO_CONSTANT_MAP.json', ASYMMETRIC/'FEEDBACK_SCREEN.json',
             base.HERE/'MIGRATION_MANIFEST.json', base.HERE/'generated/christoffel_oa29_caps.json']
    paths += sorted((base.HERE/'generated').glob('joint_shell_*/cap.json'))
    for path in paths:
        result[path.relative_to(ROOT).as_posix()] = base.sha(path)
    return result


def authenticate(record):
    outer.authenticate(record, ROOT)


class Engine(independent.Engine):
    def __init__(self, exponent):
        assert exponent in (16, 18, 20)
        super().__init__(json.loads((ASYMMETRIC.parent/'NO_CONSTANT_MAP.json').read_text()))
        self.exponent = exponent
        self.length = 1 << (exponent-7)
        self.outer_length = 256
        self.output_bits = 256*self.length
        self.cutoff = self.output_bits//10

    def identity(self):
        return dict(message_bits=1 << self.exponent, output_bits=self.output_bits,
                    outer_rows=self.length, outer_length=256, outer_dimension=128,
                    outer='fixed-p37-syndrome-subspace-0..31-GF256-0x14D',
                    outer_manifest_sha256=base.sha(base.HERE/'MIGRATION_MANIFEST.json'),
                    cutoff=self.cutoff, distance='1/10', inner=independent.INNER,
                    setup='independent row/region permutations; independent one-round transvections; '
                          'zero start; output before update; persistent state; no flush')

    @lru_cache(maxsize=24)
    def region(self, log_lam, maximum):
        assert 0 <= maximum <= self.length
        n = self.n
        rows = self.epoch(log_lam)
        a = tuple(arb_poly([row[k]*math.comb(128, j)
                           for j, row in enumerate(rows[:min(128, maximum)+1])]) for k in range(n*n))
        current = tuple(arb_poly([int(i == j)]) for i in range(n) for j in range(n))
        remaining = self.length//128
        while remaining:
            if remaining & 1:
                current = self.poly_mul(current, a, maximum)
            remaining >>= 1
            if remaining:
                a = self.poly_mul(a, a, maximum)
        return [tuple(max(arb(0), up(p[j]/math.comb(self.length, j))) for p in current)
                for j in range(maximum+1)]

    def q1(self, log_lam, linear=False):
        zero, one = self.epoch(log_lam)[:2]
        single = independent.single
        if linear:
            rz = single.identity(self.n)
            ra = (arb(0),)*(self.n*self.n)
            for _ in range(self.length//128):
                ra = single.add(single.mul(ra, zero, self.n), single.mul(rz, one, self.n))
                rz = single.mul(rz, zero, self.n)
            region = rz, tuple(v/(self.length//128) for v in ra)
        else:
            region = single.regions(zero, one, self.n, self.length//128)
        return single.moments(*region, self.n, self.outer_length)

    def costs(self, probabilities):
        caps = outer.inputs.caps_module.caps()
        assert len(probabilities) == len(BANDS)
        result = []
        for band, p in zip(BANDS, probabilities):
            if band == (256,):
                assert p == 1 and caps[256] == 1
                result.append(arb(1))
            else:
                assert 0 < p < 1
                p = number(p)
                result.append(max(up(arb(caps[w])/math.comb(256, w)/p**w/(1-p)**(256-w)) for w in band))
        return result

    def adaptive(self, rows, probabilities, q):
        roots = [up((c.log()/256).exp()) for c in self.costs(probabilities)]
        ps = list(map(number, probabilities))
        current = rows[:q+1]
        for _ in range(q):
            current = [tuple(max(up(root*((1-p)*a[k]+p*b[k])) for root, p in zip(roots, ps))
                             for k in range(self.n*self.n)) for a, b in zip(current[:-1], current[1:])]
        return independent.terminal(current[0], self.n, 256)*math.comb(self.length, q)*len(ps)**q
