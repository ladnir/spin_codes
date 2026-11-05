"""Exact zero-column ranks for a staggered local XOR shear, before D2.

Fix original source rows A across sixteen columns. Let G_j be the16x|A|
matrix Q_j Hlambda P_j restricted to those rows, with lambda=0xbc.
Pair columns j and j'=j xor1. For a nonzero binary normal n, define
T_j={r: parity(n&r)=j mod2}. The new column is

 Y_j=G_j x_j+diag(1_T_j)G_j' x_j'.

The two target sets partition the sixteen rows, so exactly one endpoint
per row is XORed with its mate. This involutive shear costs128 XORs per
256-symbol superblock and no products or keys. Independent D2 stays AFTER
the shear. That placement preserves the usual final NZ-amplitude law.

Under the UNIFORM comparison x_j,x_j' in F256^A, each zero-column mass is
256^-rank([G_j,diag(1_T_j)G_j']). Both columns zero has mass256^(-2|A|),
because the full shear and both G_j embeddings are injective. These are
raw local comparison probabilities, NOT the actual correlated RS-word
probabilities. The independent original D1 entries only rescale columns
of these matrices and therefore do not change these ranks.

For a singleton A, all G_j coefficients are nonzero. The untouched half
first forces x_j=0, then the target half forces x_j'=0. Hence every
normal/schedule/column has zero-column rank2. A nonzero RS row occupies
at least9 original columns, hence at least5 pairs and at least10 output
columns. One active original symbol gives24 output symbols per pair;
two give at least16. These facts do not by themselves improve distance:
the conservative whole-word symbol floor80 is weaker than the old144.

Ranks over GF4 equal ranks over AES-GF256 by field extension. Symbols
2/3 below embed as0xbc/0xbd, NOT AES byte2/3. Packed Gaussian elimination
uses at most64 bits for up to32 coefficients. No amplitude enumeration,
regional provider, Renyi bound, native edit, or benchmark is performed.
"""
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from time import perf_counter

from hd_pair_scheduled_outer import column_permutations, inverse_table, SCHEDULES

DIRECTORY = Path(__file__).resolve().parent
SCHEDULE = 'input-inverse-output-xor'
PINS = {'hd_pair_scheduled_outer.py':
        '64dc50bf18171d371ae16f4484087e9d42d89cc7c9b55a7576be55b1251e4741'}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen schedule source changed: '+name)
    return {**PINS, Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def pack(coefficients):
    coefficients = tuple(coefficients)
    if len(coefficients) > 32 or any(type(x) is not int or not 0 <= x < 4 for x in coefficients):
        raise ValueError('At most32 GF4 coefficients required')
    return sum(value << (2*i) for i, value in enumerate(coefficients))


def _times_two(row, low_mask):
    low, high = row & low_mask, (row >> 1) & low_mask
    return high | ((low ^ high) << 1)


def rank_gf4(rows, width):
    """Exact packed elimination; no floating point or per-coefficient products."""
    rows = tuple(rows)
    if (type(width) is not int or not 0 <= width <= 32 or len(rows) > 32
            or any(type(row) is not int or not 0 <= row < 1 << (2*width) for row in rows)):
        raise ValueError('At most32 rows/columns of packed GF4 entries required')
    low_mask = ((1 << (2*width))-1)//3
    pivots = [0]*width
    rank = 0
    for row in rows:
        while row:
            pivot = (row.bit_length()-1)//2
            coefficient = (row >> (2*pivot)) & 3
            if coefficient != 1:
                twice = _times_two(row, low_mask)
                row = row ^ twice if coefficient == 2 else twice
            if pivots[pivot]:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                rank += 1
                break
    return rank


def _source(source):
    source = tuple(source)
    if (not 1 <= len(source) <= 16 or len(set(source)) != len(source)
            or any(type(i) is not int or not 0 <= i < 16 for i in source)):
        raise ValueError('One to sixteen distinct original row labels required')
    return source


@lru_cache(maxsize=4096)
def _column_rows(source, column, schedule):
    before, after = column_permutations(schedule, column)
    result = [0]*16
    for hrow in range(16):
        result[after[hrow]] = pack(2 ^ int(hrow == before[i]) ^ ((hrow & before[i]).bit_count() & 1)
                                  for i in source)
    return tuple(result)


def column_rows(source, column, *, schedule=SCHEDULE):
    source = _source(source)
    if type(column) is not int or not 0 <= column < 16 or schedule not in SCHEDULES:
        raise ValueError('Prescribed schedule and column0..15 required')
    return _column_rows(source, column, schedule)


def paired_rows(source, column, normal=1, *, schedule=SCHEDULE):
    """Return both output maps with variables ordered(first input A, mate A)."""
    source = _source(source)
    if type(normal) is not int or not 1 <= normal <= 15:
        raise ValueError('A nonzero binary normal1..15 is required')
    first = column_rows(source, column, schedule=schedule)
    second = column_rows(source, column ^ 1, schedule=schedule)
    shift = 2*len(source)
    left, right = [], []
    for row in range(16):
        target = ((normal & row).bit_count() & 1) == (column & 1)
        left.append(first[row] | ((second[row] << shift) if target else 0))
        right.append((second[row] << shift) | (first[row] if not target else 0))
    return tuple(left), tuple(right)


def pair_ranks(source, column, normal=1, *, schedule=SCHEDULE, verify_joint=False):
    """Exact ordered marginal ranks; joint rank follows the invertible shear."""
    source = _source(source)
    left, right = paired_rows(source, column, normal, schedule=schedule)
    a = len(source)
    ranks = rank_gf4(left, 2*a), rank_gf4(right, 2*a)
    if any(not a <= rank <= min(16, 2*a) for rank in ranks):
        raise ArithmeticError('The unchanged full-column embedding was lost')
    joint = 2*a
    if verify_joint and rank_gf4(left+right, 2*a) != joint:
        raise ArithmeticError('The full pair map is not injective')
    if a == 1 and ranks != (2, 2):
        raise ArithmeticError('The universal singleton rank proof failed')
    return dict(original_source=list(source), arity=a, column=column, mate=column ^ 1,
        normal=normal, first_zero_rank=ranks[0], second_zero_rank=ranks[1], joint_zero_rank=joint,
        first_rank_gain=ranks[0]-a, second_rank_gain=ranks[1]-a,
        probability_interpretation='P(Yj=0)=256^-first_zero_rank; P(both=0)=256^-joint_zero_rank')


def source_plan():
    """ALL labelled a1/a2; at most four fixed examples for each higher arity."""
    result = []
    inverses = inverse_table()
    for a in range(1, 17):
        if a <= 2:
            sets = tuple(combinations(range(16), a))
        else:
            even_first = tuple(range(0, 16, 2))+tuple(range(1, 16, 2))
            sets = tuple(sorted({tuple(range(a)), tuple(range(16-a, 16)),
                                 tuple(sorted(even_first[:a])), tuple(sorted(inverses[:a]))}))
        for source in sets:
            result.append(dict(arity=a, source=source, all_labelled_arity_covered=a <= 2))
    return result


def screen(*, normals=tuple(range(1, 16)), columns=tuple(range(16)),
           max_pairs=100000, max_seconds=30, progress=None):
    """Explicit bounded rank inventory. Root, not helper tests, owns real runs."""
    normals, columns = tuple(normals), tuple(columns)
    if (not normals or len(set(normals)) != len(normals)
            or any(type(n) is not int or not 1 <= n <= 15 for n in normals)
            or not columns or len(set(columns)) != len(columns)
            or any(type(j) is not int or not 0 <= j < 16 for j in columns)
            or type(max_pairs) is not int or not 1 <= max_pairs <= 100000
            or not 0 < max_seconds <= 120):
        raise ValueError('Distinct valid labels and bounded pair/time budgets required')
    plan, pins, started = source_plan(), dependencies(), perf_counter()
    requested = len(plan)*len(normals)*len(columns)
    if requested > max_pairs:
        raise ValueError('Requested labelled rank pairs exceed the explicit budget')
    rows, complete = [], True
    for item in plan:
        for normal in normals:
            # Odd/even requests differ only by swapping the two input blocks.
            cache = {}
            for column in columns:
                if perf_counter()-started > max_seconds:
                    complete = False
                    break
                even = column & ~1
                if even not in cache:
                    cache[even] = pair_ranks(item['source'], even, normal)
                row = dict(cache[even])
                if column & 1:
                    row.update(column=column, mate=even,
                        first_zero_rank=row['second_zero_rank'], second_zero_rank=row['first_zero_rank'],
                        first_rank_gain=row['second_rank_gain'], second_rank_gain=row['first_rank_gain'])
                row['all_labelled_arity_covered'] = item['all_labelled_arity_covered']
                rows.append(row)
            if not complete:
                break
        if progress:
            progress(dict(completed_source=list(item['source']), completed_pairs=len(rows), complete=complete))
        if not complete:
            break
    summaries = []
    for a in range(1, 17):
        selected = [row for row in rows if row['arity'] == a]
        if selected:
            minimum = min(min(row['first_zero_rank'], row['second_zero_rank']) for row in selected)
            summaries.append(dict(arity=a, minimum_zero_rank=minimum, minimum_rank_gain=minimum-a,
                maximum_zero_rank=max(max(row['first_zero_rank'], row['second_zero_rank']) for row in selected),
                worst_examples=[row for row in selected if min(row['first_zero_rank'],row['second_zero_rank']) == minimum][:8],
                all_labelled_arity_covered=complete and a <= 2))
    if dependencies() != pins:
        raise ValueError('Sources changed during exact rank screen')
    return dict(rows=rows, summaries=summaries, requested_pairs=requested, completed_pairs=len(rows),
        complete=complete, field_order=256, coefficient_field_order=4, lambda_aes_byte=188,
        schedule=SCHEDULE, normals=list(normals), columns=list(columns),
        all_labelled_arities=[1, 2], higher_arities_representatives_only=True,
        zero_probabilities_for_uniform_column_comparison=True, RS_word_probabilities_claimed=False,
        original_A_fixed=True, final_D2_after_shear=True, extra_XORs_per_superblock=128,
        extra_field_products=0, extra_keys=0, provider_runs=0, whole_code_certificate=False,
        reviewed_source_sha256=pins, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = screen(max_seconds=args.max_seconds,
                    progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('rows', 'reviewed_source_sha256')}))
