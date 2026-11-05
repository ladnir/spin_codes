"""Exact fourth-moment upper bounds from small parity-column witnesses.

Fix A=(0,...,11), kept throughout the eight scheduled column pairs.
The pair code before final independent D2 has dimension24 in32 positions.
Let B be its actual8-row parity matrix, including the inverse shear.
For an exact support S, put d_S=|S|-rank(B_S). An information-coordinate
projection injects full-support words into(F256*)**d_S. Consequently
n_C(S)<=255**d_S, and I4=sum_S n_C(S)**4/255**(3|S|).

Exhibit four disjoint independent column sets I_i, of sizes s_i, leaving
u columns uncovered. For every S, rank(B_S)>=max_i|S intersect I_i|.
Summing this positive upper bound over all supports gives the exact
rational ceiling
 U4=256**u * sum_(0<=k_i<=s_i) prod_i C(s_i,k_i)
                         *255**(sum_i k_i-4*max_i k_i).
There are at most9**4 size tuples, not2**32 support masks. A cruder bound
is I4<=2**32*255**u, since each S has |S|-4rank(B_S)<=u.

The bounded greedy/exchange search only constructs a witness. Failure
to cover30 columns is not an obstruction or a maximum-packing claim.
Every returned witness is checked by exact GF4 elimination. GF4 ranks
equal GF256 ranks. The parity is diag(B0,B1)*shear_inverse; applying the
forward shear to a dual generator instead would be wrong.

The zero moment is exactly1 at these parameters, so U4-1 is a valid
nonzero upper atom. The optional eight-pair activity tail retains>=5
nonzero original pairs without renormalizing the density. These are
local moment upper bounds for one fixed A, not a distance certificate.
No regional provider, native operation, or automatic production run.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
from time import perf_counter

from hd_lambda_column_shear_ranks import column_rows, rank_gf4, SCHEDULE
from hd_lambda_outer import matrix as lambda_matrix

DIRECTORY = Path(__file__).resolve().parent
SOURCE = tuple(range(12))
PINS = {
    'hd_lambda_column_shear_ranks.py': '4067013201cadef412ae647cc81e394e90b49d9403dfdd4713a739a191a96fa8',
    'hd_pair_scheduled_outer.py': '64dc50bf18171d371ae16f4484087e9d42d89cc7c9b55a7576be55b1251e4741',
    'hd_binary_syndrome_outer.py': 'b1dc2dcb6e212a28526f7c7dc7f7a55999be7e6f5f53370f78768d705f2bda96',
    'hd_lambda_outer.py': 'f0c063843e295a160a1186c6c084a07dd1c2a3745005ed68057fde59913c99e7',
    'hd_outer_pair_layers.py': '017c5ba20747ae8e25acf02393e9702a03b086038e1814ac2424cdfe49d66f60',
}


def dependencies():
    for name,digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen scheduled matrix source changed: '+name)
    # The packed rank helper constructs Hlambda directly, without
    # importing its AES matrix generator. Bind both exact definitions.
    aes = lambda_matrix(lambda_value=0xbc,q=256)
    packed = column_rows(tuple(range(16)),0,schedule='static')
    embedding = (0,1,0xbc,0xbd)
    if any(embedding[(packed[i]>>(2*j))&3] != aes[i][j]
           for i in range(16) for j in range(16)):
        raise ValueError('Packed GF4 matrix differs from the pinned AES Hlambda')
    return dict(PINS, **{Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()})


def parity_columns(source, column):
    """32 packed GF4 column vectors, with left4/right4 parity coordinates."""
    source = tuple(source)
    if (len(source) != 12 or tuple(sorted(source)) != source or len(set(source)) != 12
            or any(type(i) is not int or not 0 <= i < 16 for i in source)
            or type(column) is not int or column not in range(0,16,2)):
        raise ValueError('Twelve increasing original labels and an even column0..14 required')
    outside = tuple(i for i in range(16) if i not in source)
    # H is symmetric and involutive. Transposing Q H P[:,O] gives
    # H[P(O),Q^-1(:)], the full parity of the scheduled column code.
    left = column_rows(outside,column,schedule=SCHEDULE)
    right = column_rows(outside,column+1,schedule=SCHEDULE)
    first,second = [],[]
    for row in range(16):
        a,b = left[row],right[row] << 8
        if not row & 1:  # left is target: inverse(a,b)=(a+b,b)
            first.append(a)
            second.append(a|b)
        else:           # right is target: inverse(a,b)=(a,a+b)
            first.append(a|b)
            second.append(b)
    columns = tuple(first+second)
    if rank_gf4(columns,8) != 8:
        raise ArithmeticError('Scheduled full parity lost its eight checks')
    return columns


def verify_partition(columns, checks, colors):
    """A certificate check, independent of how the four sets were found."""
    columns,colors = tuple(columns),tuple(tuple(c) for c in colors)
    if (type(checks) is not int or not 0 <= checks <= 16
            or not 1 <= len(columns) <= 32 or len(colors) != 4
            or any(type(c) is not int or not 0 <= c < 1 << (2*checks) for c in columns)):
        raise ValueError('At most32 columns,16 checks, and exactly four color sets required')
    used = [i for color in colors for i in color]
    if (len(used) != len(set(used))
            or any(type(i) is not int or not 0 <= i < len(columns) for i in used)):
        raise ValueError('Column-index sets must be disjoint and in range')
    for color in colors:
        if rank_gf4(tuple(columns[i] for i in color),checks) != len(color):
            raise ValueError('A claimed independent color set is dependent')
    return dict(colors=[list(c) for c in colors],sizes=[len(c) for c in colors],
        covered=len(used),uncovered=[i for i in range(len(columns)) if i not in used],
        exact_independence_verified=True,maximum_coverage_claimed=False)


def independent_partition(columns, checks, *, max_attempts=32, target_uncovered=2, check_time=None):
    """Bounded deterministic greedy packing plus one-eviction exchanges.

    Only the explicitly verified sets certify a bound. The search does
    not assert an optimal matroid partition, even if the target is missed.
    """
    columns = tuple(columns)
    verify_partition(columns,checks,((),(),(),()))
    if (type(max_attempts) is not int or not 1 <= max_attempts <= 64
            or type(target_uncovered) is not int or not 0 <= target_uncovered <= len(columns)):
        raise ValueError('Use1..64 attempts and an in-range uncovered target')
    best,attempts = ([],[],[],[]),0
    def independent(indices):
        return rank_gf4(tuple(columns[i] for i in indices),checks) == len(indices)
    for attempt in range(max_attempts):
        if check_time is not None:
            check_time()
        attempts += 1
        offset,direction = attempt//2,1 if attempt%2 == 0 else -1
        order = tuple((offset+direction*i)%len(columns) for i in range(len(columns)))
        colors = [[],[],[],[]]
        pending = []
        for index in order:
            for color in sorted(range(4),key=lambda c:(len(colors[c]),c)):
                if len(colors[color]) < checks and independent(colors[color]+[index]):
                    colors[color].append(index)
                    break
            else:
                pending.append(index)
        # Each successful exchange adds a previously uncolored index.
        # Recheck both changed sets; no general augmenting-path theorem
        # is needed for the validity of the returned witness.
        for index in pending:
            if check_time is not None:
                check_time()
            inserted = False
            for color in range(4):
                if len(colors[color]) < checks and independent(colors[color]+[index]):
                    colors[color].append(index)
                    inserted = True
                    break
                for position,evicted in enumerate(tuple(colors[color])):
                    replacement = colors[color][:position]+colors[color][position+1:]+[index]
                    if not independent(replacement):
                        continue
                    for target in range(4):
                        if target != color and len(colors[target]) < checks and independent(colors[target]+[evicted]):
                            colors[color] = replacement
                            colors[target].append(evicted)
                            inserted = True
                            break
                    if inserted:
                        break
                if inserted:
                    break
        if sum(map(len,colors)) > sum(map(len,best)):
            best = tuple(list(c) for c in colors)
        if len(columns)-sum(map(len,best)) <= target_uncovered:
            break
    result = verify_partition(columns,checks,best)
    result.update(attempts=attempts,target_uncovered=target_uncovered,
        target_reached=len(result['uncovered']) <= target_uncovered,
        method='bounded-greedy-one-eviction-not-maximum-matroid-union')
    return result


@lru_cache(maxsize=512)
def _size_sum(sizes,q):
    # Group equal exponents before rational arithmetic. There are at
    # most33 exponent bins and at most9**4 tuples in the a12 case.
    histogram = Counter()
    for counts in product(*(range(size+1) for size in sizes)):
        multiplicity = 1
        for size,k in zip(sizes,counts):
            multiplicity *= comb(size,k)
        histogram[sum(counts)-4*max(counts)] += multiplicity
    return sum((F(count)*F(q-1)**power for power,count in histogram.items()),F(0))


def moment_upper_from_sizes(width,checks,sizes,*,q=256):
    """Generic p4 bound, conditional on a verified four-set witness."""
    sizes = tuple(sizes)
    if (type(width) is not int or not 1 <= width <= 32
            or type(checks) is not int or not 0 <= checks <= min(width,16)
            or len(sizes) != 4 or any(type(s) is not int or not 0 <= s <= min(checks,8) for s in sizes)
            or sum(sizes)>width or type(q) is not int or q not in (4,256)):
        raise ValueError('Bounded width/checks, four independent-set sizes<=8, and q4/256 required')
    uncovered = width-sum(sizes)
    return F(q)**(4*checks-width+uncovered)*_size_sum(tuple(sorted(sizes)),q)


def pair_moment_upper(source,column,*,max_attempts=32,check_time=None):
    columns = parity_columns(source,column)
    witness = independent_partition(columns,8,max_attempts=max_attempts,
                                    target_uncovered=2,check_time=check_time)
    # Verify once more independently of the search's result fields.
    verify_partition(columns,8,witness['colors'])
    upper = moment_upper_from_sizes(32,8,witness['sizes'])
    zero = F(1)
    if upper < zero:
        raise ArithmeticError('Moment ceiling cannot exclude its exact zero atom')
    return dict(source=list(source),column_pair=[column,column+1],parity_columns=list(columns),
        witness=witness,fourth_moment_upper=str(upper),zero_moment=str(zero),
        nonzero_fourth_moment_upper=str(upper-zero),
        cruder_fourth_moment_upper=str(F(2**32)*F(255)**len(witness['uncovered'])))


def activity_tail_upper(nonzero_atoms,minimum=5):
    """Exact positive tail with zero atom1; fixed A across eight pairs."""
    values = tuple(F(value) for value in nonzero_atoms)
    if len(values) != 8 or any(value<0 for value in values) or minimum != 5:
        raise ValueError('Eight nonnegative atoms and the prescribed five-pair floor required')
    coefficients = [F(1)]
    for value in values:
        following = [F(0)]*(len(coefficients)+1)
        for degree,coefficient in enumerate(coefficients):
            following[degree] += coefficient
            following[degree+1] += value*coefficient
        coefficients = following
    return sum(coefficients[minimum:],F(0))


def receipt(*,max_attempts=32,max_seconds=30,progress=None):
    """Root-owned exact witness test for ONE prescribed A, no provider."""
    if type(max_seconds) not in (int,float) or not 0 < max_seconds <= 60:
        raise ValueError('A bounded budget in(0,60] seconds is required')
    pins,started = dependencies(),perf_counter()
    def check_time():
        if perf_counter()-started > max_seconds:
            raise TimeoutError('Fixed-A pair-moment witness budget exceeded')
    rows,complete,reason = [],True,None
    try:
        for column in range(0,16,2):
            check_time()
            row = pair_moment_upper(SOURCE,column,max_attempts=max_attempts,check_time=check_time)
            rows.append(row)
            if progress:
                progress(dict(column_pair=row['column_pair'],witness=row['witness'],
                    fourth_moment_upper=row['fourth_moment_upper'],elapsed_seconds=perf_counter()-started))
    except TimeoutError as error:
        complete,reason = False,str(error)
    tail = activity_tail_upper(row['nonzero_fourth_moment_upper'] for row in rows) if complete else None
    if dependencies() != pins:
        raise ValueError('Frozen matrix sources changed during witness calculation')
    return dict(source=list(SOURCE),arity=12,all_labelled_rowsets_covered=False,
        original_source_fixed_across_pairs=True,schedule=SCHEDULE,normal=1,
        lambda_aes_byte=188,field_order=256,coefficient_field_order=4,
        pair_width=32,pair_dimension=24,pair_codimension=8,moment_order=4,
        complete=complete,stop_reason=reason,pairs=rows,
        fixed_A_restricted_word_moment_upper=str(tail) if tail is not None else None,
        minimum_nonzero_pairs=5,restricted_density_not_renormalized=True,
        exact_integer_and_rational_arithmetic=True,local_moment_upper_bound=True,
        maximum_packing_claimed=False,regional_provider_runs=0,
        rs_factor_included=False,source_count_included=False,rowset_count_included=False,
        failure_probability_bound=False,whole_code_certificate=False,native_operations=0,
        reviewed_source_sha256=pins,elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-attempts',type=int,default=32)
    parser.add_argument('--max-seconds',type=float,default=30)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = receipt(max_attempts=args.max_attempts,max_seconds=args.max_seconds,
                     progress=lambda row: print(json.dumps(row),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('pairs','reviewed_source_sha256')}))
