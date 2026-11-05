"""Exact projected support counts for native gamma profile cells.

Interior counts upper-bound actual supports by the minimum of three
projection sums. Singleton-free faces use the exact contracted-cycle
count. Each weight uses one-dimensional integer prefix tables; no profile
pair is enumerated. The bounded cache stores only completed tables.
"""

import ast
from collections import OrderedDict
from dataclasses import dataclass
from hashlib import sha256
from math import comb
from pathlib import Path

from hd_field128_gamma_cell_geometry import Cell


def source_pins():
    directory = Path(__file__).resolve().parent
    pending, pins = [Path(__file__).name], {}
    while pending:
        name = pending.pop()
        if name in pins:
            continue
        path = directory/name
        raw = path.read_bytes()
        pins[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw, filename=str(path))):
            modules = ([item.name for item in node.names] if isinstance(node, ast.Import)
                       else [node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in pins:
                    pending.append(target)
    return dict(sorted(pins.items()))


@dataclass(frozen=True)
class Prefix:
    first: int
    cumulative: tuple

    def interval(self, low, high):
        """Sum an inclusive integer interval, with zero outside the table."""
        low = max(low, self.first)
        high = min(high, self.first+len(self.cumulative)-2)
        if low > high:
            return 0
        return self.cumulative[high-self.first+1]-self.cumulative[low-self.first]


def _marginal_prefix(pairs, weight, budget):
    first, last = max(0, weight-pairs), weight//2
    single = weight-2*first
    value = comb(pairs, first)*comb(pairs-first, single)*(1 << single)
    cumulative, total = [0], 0
    for full in range(first, last+1):
        if not full % 128:
            budget.check()
        total += value
        cumulative.append(total)
        if full < last:
            single, empty = weight-2*full, pairs-weight+full
            value, remainder = divmod(value*single*(single-1), 4*(full+1)*(empty+1))
            if remainder:
                raise ArithmeticError('Nonintegral exact matching-marginal recurrence')
    return Prefix(first, tuple(cumulative))


def _run_prefix(length, weight, budget):
    if weight in (0, length):
        return Prefix(0, (0, 1))
    last, value, total = min(weight, length-weight), length, 0
    cumulative = [0]
    for runs in range(1, last+1):
        if not runs % 128:
            budget.check()
        total += value
        cumulative.append(total)
        if runs < last:
            value, remainder = divmod(value*(weight-runs)*(length-weight-runs), runs*(runs+1))
            if remainder:
                raise ArithmeticError('Nonintegral exact cyclic-run recurrence')
    return Prefix(1, tuple(cumulative))


class PrefixCache:
    """Bounded LRU: one interior bundle or one face table per entry.

    Each entry has O(pairs) integers. No all-weight or profile-pair table
    is retained. A caller may share this cache across a batch of cells.
    """

    def __init__(self, max_entries=4):
        if type(max_entries) is not int or not 1 <= max_entries <= 16:
            raise ValueError('Prefix cache must retain between1 and16 weight entries')
        self.max_entries = max_entries
        self._entries = OrderedDict()

    def __len__(self):
        return len(self._entries)

    def clear(self):
        self._entries.clear()

    def get(self, kind, pairs, weight, budget):
        if (kind not in ('interior', 'face') or type(pairs) is not int or pairs < 1
                or type(weight) is not int or not 0 <= weight <= (2*pairs if kind == 'interior' else pairs)):
            raise ValueError('Invalid exact prefix-table scope')
        budget.check()
        key = kind, pairs, weight
        saved = self._entries.get(key)
        if saved is not None:
            self._entries.move_to_end(key)
            return saved
        if kind == 'interior':
            saved = (_marginal_prefix(pairs, weight, budget), _run_prefix(2*pairs, weight, budget))
        else:
            saved = (_run_prefix(pairs, weight, budget),)
        budget.check()
        self._entries[key] = saved
        if len(self._entries) > self.max_entries:
            self._entries.popitem(last=False)
        return saved


_CACHE = PrefixCache()


def _count(cell, budget, cache):
    if not isinstance(cell, Cell):
        raise ValueError('Expected a native exact gamma Cell')
    if cache is None:
        cache = _CACHE
    if not isinstance(cache, PrefixCache):
        raise ValueError('Expected a bounded PrefixCache')
    budget.check()
    pairs, slices = cell.pairs, 0
    if cell.kind == 'interior':
        (rl, rh), (al, ah), (bl, bh) = cell.ranges
        even, odd, runs = 0, 0, 0
        for weight in range(max(1, rl), min(2*pairs-1, rh)+1):
            budget.check()
            # Native interior slices are rectangles: both singleton counts
            # are at least1, so a+b<=r-1 is already implied. These integer
            # endpoints are exact, not sampled real projection estimates.
            alo, ahi = max(al, 0, weight-pairs), min(ah, pairs, (weight-1)//2)
            blo, bhi = max(bl, 0, weight-pairs), min(bh, pairs, (weight-1)//2)
            if alo > ahi or blo > bhi:
                continue
            marginal, cycle = cache.get('interior', pairs, weight, budget)
            even += marginal.interval(alo, ahi)
            odd += marginal.interval(blo, bhi)
            runs += cycle.interval(weight-ahi-bhi, weight-alo-blo)
            slices += 1
        budget.check()
        return min(even, odd, runs), (even, odd, runs), slices, False
    if cell.kind in ('even_face', 'odd_face'):
        (ul, uh), (kl, kh) = cell.ranges
        total = 0
        for weight in range(max(1, ul), min(pairs-1, uh)+1):
            budget.check()
            low, high = max(1, kl), min(kh, weight, pairs-weight)
            if low > high:
                continue
            cycle, = cache.get('face', pairs, weight, budget)
            total += cycle.interval(low, high)
            slices += 1
        budget.check()
        return total, (total,), slices, True
    raise ValueError('Unsupported gamma cell kind')


def support_count(cell, budget, *, cache=None):
    """Integer upper bound on actual supports in the complete cell domain."""
    return _count(cell, budget, cache)[0]


def count_diagnostics(cell, budget, *, cache=None):
    count, projections, slices, exact = _count(cell, budget, cache)
    names = ('contracted_cycle',) if exact else ('even_marginal', 'odd_marginal', 'cycle_runs')
    return dict(support_count_upper=count, exact_support_count=exact,
                method='contracted-cycle-exact' if exact else 'minimum-of-exact-projection-sums',
                projection_sums=dict(zip(names, projections)), nonempty_integer_slices=slices,
                profile_pairs_enumerated=0, omitted_coordinate_factor=False)
