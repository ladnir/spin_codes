"""All-labelled-A12 fourth-moment witness inventory and directed gate.

For each original twelve-row set A, keep A fixed through eight scheduled
column pairs. The frozen cover helper constructs the actual eight-row
parity of each 32-coordinate pair code. Four disjoint independent column
sets certify an exact rational upper U_A,j for its fourth density moment.
The zero moment is1. Thus the positive degree>=5 tail of
prod_j(1+v*(U_A,j-1)) upper-bounds the restricted whole-group moment J_A.
Every nonzero RS row occupies>=9 columns, hence>=5 original pairs.

The raw group weight is beta*sum_A J_A^(1/4), where
beta=256^(16*12)/255^(8*12). The sum contains all1820 original A once;
there is no additional binomial factor. Whole-group Minkowski keeps A
fixed rather than resampling it between pairs. For r18, reuse the frozen
full99-profile UNIFORM reference and its reviewed method_expression at
p4, alpha1/5, theta4/15, and mark691/700. That expression charges source
counts outside theta and one projective factor inside theta.

Inventory construction is exact and independent of Arb precision.
Evaluation at256/384 bits rechecks every saved parity and independence
witness, but does not repeat the search. Neither mode invokes a provider.
A completed evaluation is an upper bound for this homogeneous r18/A12
class under the new shear construction, not a whole-code certificate.
Final D2 follows the shear; the comparison requires full independent D1
and permits only one common gauge per original column pair.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

import hd_lambda_shear_moment_cover as cover
import hd_lambda_shear_renyi_ideal as ideal

DIRECTORY = Path(__file__).resolve().parent
SOURCES = tuple(combinations(range(16),12))
PINS = {
    'hd_lambda_shear_moment_cover.py': '19c0ffb208349383d9f45bc91127c5fbb940ecc3d1ba4ddb3daa07dc1e24fb61',
    'hd_lambda_shear_renyi_ideal.py': '7c2d5e9b720cd6cae09212de0fbf11ea0c19352f023706d08bd2919b5477c7a5',
}
SCOPE = dict(schema='shear-a12-fourth-moment-inventory-v1',arity=12,
    labelled_rowsets=1820,pair_count=8,field_order=256,coefficient_field_order=4,
    pair_width=32,pair_dimension=24,pair_codimension=8,moment_order=4,
    schedule='input-inverse-output-xor',normal=1,lambda_aes_byte=188,
    original_A_fixed_across_pairs=True,minimum_active_pairs=5,
    exact_integer_and_rational_arithmetic=True,maximum_packing_claimed=False)


def dependencies():
    pins = {}
    for sources in (cover.dependencies(),ideal.dependencies(),PINS):
        for name,digest in sources.items():
            if name in pins and pins[name] != digest:
                raise ValueError('Conflicting source provenance: '+name)
            if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Frozen dependency changed: '+name)
            pins[name] = digest
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return pins


def _timer(max_seconds):
    if type(max_seconds) not in (int,float) or not 0 < max_seconds <= 600:
        raise ValueError('A finite budget in(0,600] seconds is required')
    started = perf_counter()
    def check():
        if perf_counter()-started > max_seconds:
            raise TimeoutError('All-A12 moment job exceeded its time budget')
    return started,check


def source_record(source,*,max_attempts=32,check_time=None):
    """One exact record; callers determine whether every A completed."""
    pairs = []
    for column in range(0,16,2):
        if check_time is not None:
            check_time()
        row = cover.pair_moment_upper(source,column,max_attempts=max_attempts,check_time=check_time)
        pairs.append(dict(column_pair=row['column_pair'],parity_columns=row['parity_columns'],
            colors=row['witness']['colors'],sizes=row['witness']['sizes'],
            uncovered=row['witness']['uncovered'],search_attempts=row['witness']['attempts'],
            fourth_moment_upper=row['fourth_moment_upper'],
            nonzero_fourth_moment_upper=row['nonzero_fourth_moment_upper']))
    tail = cover.activity_tail_upper(F(row['nonzero_fourth_moment_upper']) for row in pairs)
    return dict(source=list(source),pairs=pairs,restricted_word_moment_upper=str(tail))


def inventory(*,max_attempts=32,max_seconds=120,progress=None):
    """Root-owned bounded search. A timeout yields a non-certifying partial."""
    if type(max_attempts) is not int or not 1 <= max_attempts <= 64:
        raise ValueError('Use1..64 bounded packing attempts')
    pins = dependencies()
    started,check = _timer(max_seconds)
    rows,complete,reason = [],True,None
    try:
        for source in SOURCES:
            check()
            rows.append(source_record(source,max_attempts=max_attempts,check_time=check))
            if progress and (len(rows)%128 == 0 or len(rows)==len(SOURCES)):
                progress(dict(completed_rowsets=len(rows),total_rowsets=len(SOURCES),
                    elapsed_seconds=perf_counter()-started))
    except TimeoutError as error:
        complete,reason = False,str(error)
    if dependencies() != pins:
        raise ValueError('Sources changed during inventory construction')
    return dict(SCOPE,complete=complete,stop_reason=reason,completed_rowsets=len(rows),
        all_labelled_rowsets_covered=complete and len(rows)==1820,rows=rows,
        max_attempts=max_attempts,reviewed_source_sha256=pins,
        failure_probability_bound=False,whole_code_certificate=False,
        fresh_regional_banks=0,native_operations=0,elapsed_seconds=perf_counter()-started)


def validate_record(row,source,*,check_time=None):
    """Reconstruct parity, verify each set, and recompute exact rational data."""
    if (row.get('source') != list(source) or any(type(i) is not int for i in row['source'])
            or len(row.get('pairs',())) != 8):
        raise ValueError('Wrong original source or missing column pairs')
    atoms = []
    for index,pair in enumerate(row['pairs']):
        if check_time is not None:
            check_time()
        column = 2*index
        if pair.get('column_pair') != [column,column+1]:
            raise ValueError('Column-pair chronology changed')
        actual = cover.parity_columns(source,column)
        if pair.get('parity_columns') != list(actual):
            raise ValueError('Saved parity does not match the actual scheduled shear')
        witness = cover.verify_partition(actual,8,pair.get('colors',()))
        if pair.get('sizes') != witness['sizes'] or pair.get('uncovered') != witness['uncovered']:
            raise ValueError('Saved independent-set sizes/uncovered labels changed')
        upper = cover.moment_upper_from_sizes(32,8,witness['sizes'])
        if upper < 1 or F(pair['fourth_moment_upper']) != upper or F(pair['nonzero_fourth_moment_upper']) != upper-1:
            raise ValueError('Saved exact fourth-moment arithmetic changed')
        atoms.append(upper-1)
    tail = cover.activity_tail_upper(atoms)
    if F(row['restricted_word_moment_upper']) != tail or tail <= 0:
        raise ValueError('Saved positive activity tail changed')
    return tail


def validate_inventory(data,*,check_time=None,progress=None):
    for key,value in dict(SCOPE,complete=True,all_labelled_rowsets_covered=True,
                          completed_rowsets=1820,reviewed_source_sha256=dependencies()).items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError('Wrong exact inventory scope/provenance: '+key)
    rows = data.get('rows')
    if type(rows) is not list or len(rows) != 1820:
        raise ValueError('All1820 labelled original sources are required')
    tails = []
    for row,source in zip(rows,SOURCES):
        tails.append(validate_record(row,source,check_time=check_time))
        if progress and len(tails)%128 == 0:
            progress(dict(verified_rowsets=len(tails),total_rowsets=1820))
    return tuple(tails)


def checked_inventory(path,digest,*,check_time=None,progress=None):
    path = Path(path).resolve()
    if type(digest) is not str or len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
        raise ValueError('An explicit lowercase SHA256 inventory pin is required')
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != digest:
        raise ValueError('Exact inventory hash does not match its supplied pin')
    tails = validate_inventory(json.loads(raw),check_time=check_time,progress=progress)
    return tails,dict(file=str(path),sha256=digest,verified_rowsets=1820,verified_pair_witnesses=14560)


def raw_group_upper(tails):
    """Arb context required; rowset multiplicity occurs only in this sum."""
    tails = tuple(F(value) for value in tails)
    if not tails or len(tails)>1820 or any(value<=0 for value in tails):
        raise ValueError('One to1820 positive exact word-moment upper bounds required')
    multiplicities = Counter(tails)
    roots = arb(0)
    exponent = ideal.A.uniform.rational(F(1,4))
    for tail,count in sorted(multiplicities.items()):
        roots += count*ideal.A.uniform.rational(tail)**exponent
    beta = F(256)**192/F(255)**96
    raw = ideal.A.uniform.rational(beta)*roots
    if not raw > 0 or not raw.is_finite():
        raise ArithmeticError('Invalid directed raw group upper bound')
    return raw,dict(beta=str(beta),rowset_terms=len(tails),distinct_tail_values=len(multiplicities),
        fourth_root_sum_enclosure=str(roots),raw_group_log2_enclosure=str(ideal.A.log2(raw)))


def certificate(path,digest,*,precision=256,max_seconds=120,progress=None):
    """Root-owned exact-inventory replay at one matching reference precision."""
    if type(precision) is not int or precision not in (256,384):
        raise ValueError('Use256/384-bit Arb with the matching cached reference')
    pins = dependencies()
    started,check = _timer(max_seconds)
    tails,inventory_pin = checked_inventory(path,digest,check_time=check,progress=progress)
    with ctx.workprec(precision):
        check()
        reference,reference_pin = ideal.checked_reference(precision)
        raw,details = raw_group_upper(tails)
        value = ideal.method_expression(raw,4,reference)
        lower = -ideal.A.uniform.rational_ceiling(-value)
        upper = ideal.A.uniform.rational_ceiling(value)
        reference_text,value_text = str(reference),str(value)
    check()
    if dependencies() != pins:
        raise ValueError('Sources changed during directed evaluation')
    if sha256(Path(inventory_pin['file']).read_bytes()).hexdigest() != digest:
        raise ValueError('Exact inventory changed during directed evaluation')
    if sha256((DIRECTORY/reference_pin['file']).read_bytes()).hexdigest() != reference_pin['sha256']:
        raise ValueError('Cached uniform reference changed during directed evaluation')
    return dict(details,complete=True,construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        active_groups=18,homogeneous_arity=12,all_labelled_original_rowsets=True,
        original_A_fixed_across_pairs=True,source_profiles=99,moment_order=4,
        uniform_alpha='1/5',outer_theta='4/15',mark='691/700',forbidden_weight=734003,
        pair_count=8,minimum_active_pairs=5,exact_original_row_floors_retained=False,
        restricted_density_not_renormalized=True,source_count_outside_theta=True,
        rowset_sum_inside_theta=True,additional_rowset_binomial_factor=False,
        projective_normalizations=1,additional_beta_factor=False,
        full_independent_D1_comparison=True,permitted_gauges_per_pair=1,
        final_independent_D2_after_shear=True,arithmetic_enclosures_directed=True,
        precision_bits=precision,cached_uniform_log2_enclosure=reference_text,
        failure_upper_log2_enclosure=value_text,failure_upper_log2_lower=str(lower),
        failure_upper_log2_upper=str(upper),selected_class_40_bits=bool(upper<=-40),
        failure_probability_upper_bound=True,mixed_arities_covered=False,
        other_group_counts_covered=False,whole_code_certificate=False,
        exact_inventory=inventory_pin,cached_reference=reference_pin,
        reviewed_source_sha256=pins,fresh_regional_banks=0,packing_searches=0,
        native_operations=0,elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('inventory','certificate'))
    parser.add_argument('--max-attempts',type=int,default=32)
    parser.add_argument('--max-seconds',type=float,default=120)
    parser.add_argument('--precision',type=int,choices=(256,384),default=256)
    parser.add_argument('--inventory',type=Path)
    parser.add_argument('--inventory-sha256')
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    progress = lambda row: print(json.dumps(row),flush=True)
    if args.mode == 'inventory':
        result = inventory(max_attempts=args.max_attempts,max_seconds=args.max_seconds,progress=progress)
    else:
        if args.inventory is None or args.inventory_sha256 is None:
            parser.error('Certificate mode requires an inventory and its explicit SHA256')
        result = certificate(args.inventory,args.inventory_sha256,precision=args.precision,
                             max_seconds=args.max_seconds,progress=progress)
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps({key:value for key,value in result.items()
        if key not in ('rows','reviewed_source_sha256')}))
