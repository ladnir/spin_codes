"""Exact parity-cover moments for every labelled arity11..16 row set.

Fix an original A of size a. Each adjacent-column shear code has length32,
dimension2a, and b=32-2a parity checks. For an exact output support S,
shortening and an information-coordinate projection give
 n_C(S)<=(q-1)^(|S|-rank B_S).
Exhibit p disjoint independent parity-column sets of sizes s_i, p=2,3,4.
Their rank lower bound yields the exact rational upper
 U_p=q^(p*b-32+u) sum_k prod_i C(s_i,k_i)
                    (q-1)^(sum k_i-p max k_i), u=32-sum s_i.
The zero atom is exactly q^(p*b-32), subtracted once to obtain a
nonzero upper atom. At a16 the code is uniform and U_p=1 exactly.

The positive degree>=5 tail of eight pair atoms keeps the same original
A throughout. Any nonzero RS row occupies>=9 original columns, hence
>=5 pairs. The raw arity weight is
 B_a=256^(16a)/255^(8a) * sum_{|A|=a} tail_A^(1/p).
There is no additional row-set multiplicity. For each FIXED p, form the
18th power of sum_{a=11}^16 B_a x^(a-11). Its degree d covers total
arity198+d. Apply the reviewed uniform-reference expression, with one
projective factor and source counts outside theta. Only afterwards may
we select the smallest completed p-bound for each total arity.

The inventory retains explicit labelled parity/witness data. A failed
packing search still gives a valid verified witness, not an obstruction.
The directed evaluator rechecks witnesses without repeating searches.
This covers only r18 with every group arity11..16 under the scheduled
shear construction, never the whole code. No provider or native job.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
import json
from math import comb
from pathlib import Path
from time import perf_counter

from flint import arb,arb_poly,ctx

import hd_lambda_shear_moment_all_a12 as old
from hd_lambda_column_shear_ranks import column_rows,rank_gf4,SCHEDULE

DIRECTORY = Path(__file__).resolve().parent
ARITIES,ORDERS = tuple(range(11,17)),(2,3,4)
SOURCES = tuple(source for a in ARITIES for source in combinations(range(16),a))
PINS = {'hd_lambda_shear_moment_all_a12.py':
    '3544a89d48894cf6c82ec478bfac4d96f3840df9f54d49b948e1c7d2f999651c'}
SCOPE = dict(schema='shear-neighbor-moments-v1',arities=list(ARITIES),moment_orders=list(ORDERS),
    labelled_rowsets=6885,pair_count=8,field_order=256,coefficient_field_order=4,
    pair_width=32,schedule=SCHEDULE,normal=1,lambda_aes_byte=188,
    original_A_fixed_across_pairs=True,minimum_active_pairs=5,
    exact_integer_and_rational_arithmetic=True,maximum_packing_claimed=False)


def dependencies():
    pins = old.dependencies()
    for name,digest in PINS.items():
        if name in pins and pins[name] != digest:
            raise ValueError('Conflicting source pin: '+name)
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen source changed: '+name)
        pins[name] = digest
    pins[Path(__file__).name] = sha256(Path(__file__).read_bytes()).hexdigest()
    return pins


def parity_columns(source,column):
    source = tuple(source)
    if (len(source) not in ARITIES or tuple(sorted(source)) != source or len(set(source)) != len(source)
            or any(type(i) is not int or not 0<=i<16 for i in source)
            or type(column) is not int or column not in range(0,16,2)):
        raise ValueError('Increasing arity11..16 source and an even column0..14 required')
    outside = tuple(i for i in range(16) if i not in source)
    b = 2*len(outside)
    if not outside:
        return (0,)*32
    left,right = column_rows(outside,column),column_rows(outside,column+1)
    first,second = [],[]
    shift = 2*len(outside)
    for i in range(16):
        x,y = left[i],right[i]<<shift
        first.append(x if not i&1 else x|y)
        second.append(x|y if not i&1 else y)
    result = tuple(first+second)
    if rank_gf4(result,b) != b:
        raise ArithmeticError('Actual scheduled pair parity lost rank')
    return result


def verify_partition(columns,checks,colors,order):
    columns,colors = tuple(columns),tuple(tuple(c) for c in colors)
    if (type(order) is not int or order not in ORDERS or len(colors)!=order
            or type(checks) is not int or not 0<=checks<=10 or not 1<=len(columns)<=32
            or any(type(x) is not int or not 0<=x<1<<(2*checks) for x in columns)):
        raise ValueError('Bounded parity and exactly p=2/3/4 color sets required')
    used = [i for color in colors for i in color]
    if (len(set(used))!=len(used)
            or any(type(i) is not int or not 0<=i<len(columns) for i in used)):
        raise ValueError('Witness sets must be disjoint labelled column subsets')
    if any(rank_gf4(tuple(columns[i] for i in color),checks)!=len(color) for color in colors):
        raise ValueError('Dependent witness color set')
    return dict(colors=[list(c) for c in colors],sizes=[len(c) for c in colors],
        uncovered=[i for i in range(len(columns)) if i not in used],
        exact_independence_verified=True,maximum_coverage_claimed=False)


def _extend(basis,row,width):
    """Incremental packed GF4 basis for search; final rank is reverified."""
    low_mask=((1<<(2*width))-1)//3
    while row:
        pivot=(row.bit_length()-1)//2
        coefficient=(row>>(2*pivot))&3
        if coefficient!=1:
            low,high=row&low_mask,(row>>1)&low_mask
            twice=high|((low^high)<<1)
            row=row^twice if coefficient==2 else twice
        if basis[pivot]:
            row^=basis[pivot]
        else:
            result=list(basis)
            result[pivot]=row
            return tuple(result)
    return None


def independent_partition(columns,checks,order,*,max_attempts=32,check_time=None):
    columns=tuple(columns)
    verify_partition(columns,checks,[()]*order,order)
    if type(max_attempts) is not int or not 1<=max_attempts<=64:
        raise ValueError('Use1..64 bounded search attempts')
    target=min(len(columns),order*checks)
    best=[[] for _ in range(order)]
    attempts=0
    def basis_for(indices):
        basis=(0,)*checks
        for index in indices:
            basis=_extend(basis,columns[index],checks)
            if basis is None:
                return None
        return basis
    for attempt in range(max_attempts):
        if check_time is not None: check_time()
        attempts+=1
        offset,direction=attempt//2,1 if attempt%2==0 else -1
        colors=[[] for _ in range(order)]
        bases=[(0,)*checks for _ in range(order)]
        pending=[]
        for step in range(len(columns)):
            index=(offset+direction*step)%len(columns)
            for color in sorted(range(order),key=lambda c:(len(colors[c]),c)):
                updated=_extend(bases[color],columns[index],checks)
                if updated is not None:
                    colors[color].append(index)
                    bases[color]=updated
                    break
            else:
                pending.append(index)
        if sum(map(len,colors))<target:
            for index in pending:
                if check_time is not None: check_time()
                inserted=False
                for color in range(order):
                    updated=_extend(bases[color],columns[index],checks)
                    if updated is not None:
                        colors[color].append(index)
                        bases[color]=updated
                        inserted=True
                        break
                    for position,evicted in enumerate(tuple(colors[color])):
                        replacement=colors[color][:position]+colors[color][position+1:]+[index]
                        updated=basis_for(replacement)
                        if updated is None: continue
                        for other in range(order):
                            moved=None if other==color else _extend(bases[other],columns[evicted],checks)
                            if moved is not None:
                                colors[color],bases[color]=replacement,updated
                                colors[other].append(evicted)
                                bases[other]=moved
                                inserted=True
                                break
                        if inserted: break
                    if inserted: break
        if sum(map(len,colors))>sum(map(len,best)):
            best=[list(c) for c in colors]
        if sum(map(len,best))==target: break
    result=verify_partition(columns,checks,best,order)
    result.update(search_attempts=attempts,capacity_target=target,
        capacity_target_reached=sum(result['sizes'])==target)
    return result


@lru_cache(maxsize=2048)
def _size_sum(sizes,q):
    p=len(sizes)
    histogram=Counter()
    for counts in product(*(range(s+1) for s in sizes)):
        multiplicity=1
        for s,k in zip(sizes,counts): multiplicity*=comb(s,k)
        histogram[sum(counts)-p*max(counts)]+=multiplicity
    return sum((F(n)*F(q-1)**e for e,n in histogram.items()),F(0))


def moment_upper_from_sizes(width,checks,sizes,*,q=256):
    sizes=tuple(sizes)
    p=len(sizes)
    if (type(width) is not int or not 1<=width<=32 or type(checks) is not int
            or not 0<=checks<=min(width,10) or p not in ORDERS or sum(sizes)>width
            or any(type(s) is not int or not 0<=s<=checks for s in sizes)
            or type(q) is not int or q not in (4,256)):
        raise ValueError('Width<=32,checks<=10,p2/3/4, independent sizes and q4/256 required')
    zero=F(q)**(p*checks-width)
    upper=zero*F(q)**(width-sum(sizes))*_size_sum(tuple(sorted(sizes)),q)
    if upper<zero: raise ArithmeticError('Moment bound excluded the zero atom')
    return zero,upper-zero,upper


def activity_tail(zero,nonzero_atoms):
    zero,atoms=F(zero),tuple(map(F,nonzero_atoms))
    if zero<=0 or len(atoms)!=8 or any(n<0 for n in atoms):
        raise ValueError('A positive zero atom and eight nonnegative atoms required')
    coefficients=[F(1)]
    for atom in atoms:
        following=[F(0)]*(len(coefficients)+1)
        for k,value in enumerate(coefficients):
            following[k]+=zero*value
            following[k+1]+=atom*value
        coefficients=following
    return sum(coefficients[5:],F(0))


def source_record(source,*,max_attempts=32,check_time=None):
    a=len(source)
    b=32-2*a
    pairs=[]
    atoms={p:[] for p in ORDERS}
    for column in range(0,16,2):
        if check_time is not None: check_time()
        columns=parity_columns(source,column)
        orders={}
        for p in ORDERS:
            witness=independent_partition(columns,b,p,max_attempts=max_attempts,check_time=check_time)
            zero,nonzero,upper=moment_upper_from_sizes(32,b,witness['sizes'])
            orders[str(p)]=dict(witness,zero_moment=str(zero),nonzero_moment_upper=str(nonzero),
                               moment_upper=str(upper))
            atoms[p].append(nonzero)
        pairs.append(dict(column_pair=[column,column+1],parity_columns=list(columns),orders=orders))
    tails={str(p):str(activity_tail(F(256)**(p*b-32),atoms[p])) for p in ORDERS}
    return dict(source=list(source),arity=a,pair_codimension=b,pairs=pairs,
                restricted_word_moment_uppers=tails,uniform_pair_code=(a==16))


def inventory(*,max_attempts=32,max_seconds=120,progress=None):
    if type(max_attempts) is not int or not 1<=max_attempts<=64:
        raise ValueError('Use1..64 bounded search attempts')
    pins=dependencies()
    started,check=old._timer(max_seconds)
    rows,complete,reason=[],True,None
    try:
        for source in SOURCES:
            check()
            rows.append(source_record(source,max_attempts=max_attempts,check_time=check))
            if progress and (len(rows)%128==0 or len(rows)==len(SOURCES)):
                progress(dict(completed_rowsets=len(rows),total_rowsets=6885,
                    last_arity=len(source),elapsed_seconds=perf_counter()-started))
    except TimeoutError as error:
        complete,reason=False,str(error)
    if dependencies()!=pins: raise ValueError('Sources changed during inventory')
    return dict(SCOPE,rows=rows,complete=complete,stop_reason=reason,completed_rowsets=len(rows),
        all_labelled_rowsets_covered=complete and len(rows)==6885,max_attempts=max_attempts,
        reviewed_source_sha256=pins,failure_probability_bound=False,whole_code_certificate=False,
        fresh_regional_banks=0,native_operations=0,elapsed_seconds=perf_counter()-started)


def validate_record(row,source,*,check_time=None):
    a=len(source)
    b=32-2*a
    if (row.get('source')!=list(source) or any(type(i) is not int for i in row['source'])
            or row.get('arity')!=a or row.get('pair_codimension')!=b
            or row.get('uniform_pair_code') is not (a==16) or len(row.get('pairs',()))!=8
            or set(row.get('restricted_word_moment_uppers',{}))!={'2','3','4'}):
        raise ValueError('Wrong fixed original source or pair-code dimensions')
    atoms={p:[] for p in ORDERS}
    for index,pair in enumerate(row['pairs']):
        if check_time is not None: check_time()
        columns=parity_columns(source,2*index)
        if (pair.get('column_pair')!=[2*index,2*index+1]
                or pair.get('parity_columns')!=list(columns) or set(pair.get('orders',{}))!={'2','3','4'}):
            raise ValueError('Actual parity, pair chronology, or moment orders changed')
        for p in ORDERS:
            saved=pair['orders'][str(p)]
            witness=verify_partition(columns,b,saved.get('colors',()),p)
            if saved.get('sizes')!=witness['sizes'] or saved.get('uncovered')!=witness['uncovered']:
                raise ValueError('Derived witness sizes or uncovered labels differ')
            zero,nonzero,upper=moment_upper_from_sizes(32,b,witness['sizes'])
            if (F(saved['zero_moment'])!=zero or F(saved['nonzero_moment_upper'])!=nonzero
                    or F(saved['moment_upper'])!=upper):
                raise ValueError('Exact local moment arithmetic differs')
            atoms[p].append(nonzero)
    tails={p:activity_tail(F(256)**(p*b-32),atoms[p]) for p in ORDERS}
    if any(tails[p]<=0 or F(row['restricted_word_moment_uppers'][str(p)])!=tails[p] for p in ORDERS):
        raise ValueError('Positive word activity-tail arithmetic differs')
    return tails


def validate_inventory(data,*,check_time=None,progress=None):
    for key,value in dict(SCOPE,complete=True,all_labelled_rowsets_covered=True,
                          completed_rowsets=6885,reviewed_source_sha256=dependencies()).items():
        if type(data.get(key)) is not type(value) or data[key]!=value:
            raise ValueError('Wrong complete inventory scope/provenance: '+key)
    if type(data.get('rows')) is not list or len(data['rows'])!=6885:
        raise ValueError('All6885 original labelled sources are required')
    grouped={a:{p:Counter() for p in ORDERS} for a in ARITIES}
    for index,(row,source) in enumerate(zip(data['rows'],SOURCES),1):
        tails=validate_record(row,source,check_time=check_time)
        for p in ORDERS: grouped[len(source)][p][tails[p]]+=1
        if progress and index%128==0:
            progress(dict(verified_rowsets=index,total_rowsets=6885))
    return grouped


def checked_inventory(path,digest,*,check_time=None,progress=None):
    path=Path(path).resolve()
    if type(digest) is not str or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
        raise ValueError('An explicit lowercase inventory SHA256 is required')
    raw=path.read_bytes()
    if sha256(raw).hexdigest()!=digest: raise ValueError('Inventory hash differs')
    grouped=validate_inventory(json.loads(raw),check_time=check_time,progress=progress)
    return grouped,dict(file=str(path),sha256=digest,verified_rowsets=6885,
        verified_pair_parities=55080,verified_color_witnesses=165240)


def raw_weights(grouped,order):
    if type(order) is not int or order not in ORDERS: raise ValueError('Use p2/3/4')
    weights,rows=[],[]
    exponent=old.ideal.A.uniform.rational(F(1,order))
    for a in ARITIES:
        counts=grouped[a][order]
        if (sum(counts.values())!=comb(16,a) or any(type(n) is not int or n<=0 for n in counts.values())
                or any(not isinstance(tail,F) or tail<=0 for tail in counts)):
            raise ValueError('Exact positive tails must contain every labelled source once')
        norm=sum((n*old.ideal.A.uniform.rational(tail)**exponent for tail,n in sorted(counts.items())),arb(0))
        beta=F(256)**(16*a)/F(255)**(8*a)
        value=old.ideal.A.uniform.rational(beta)*norm
        if not value>0 or not value.is_finite(): raise ArithmeticError('Invalid raw group norm')
        weights.append(value)
        rows.append(dict(arity=a,labelled_rowsets=comb(16,a),distinct_tail_values=len(counts),
            beta=str(beta),raw_group_log2_enclosure=str(old.ideal.A.log2(value))))
    return weights,rows


def _total_range(total_min,total_max):
    if (type(total_min) is not int or type(total_max) is not int
            or not 198<=total_min<=total_max<=288):
        raise ValueError('A bounded total-arity range inside198..288 is required')
    return range(total_min,total_max+1)


def controls_from_reference(grouped,reference,*,total_min=198,total_max=288,check_time=None):
    totals=_total_range(total_min,total_max)
    A=old.ideal.A
    threshold=-A.uniform.rational(F(1,5))*734003*A.log2(A.uniform.rational(F(691,700)))
    controls=[]
    for p in ORDERS:
        if check_time is not None: check_time()
        weights,diagnostics=raw_weights(grouped,p)
        coefficients=arb_poly(weights)**18
        rows=[]
        for total in totals:
            # Fixed p across all18 groups. Source counts occur in reference.
            value=old.ideal.limits.total_log_bound(coefficients[total-198],p,reference,threshold)
            rows.append(dict(total_arity=total,failure_log2_upper=str(A.uniform.rational_ceiling(value))))
        controls.append(dict(moment_order=p,outer_theta=str(F(1,5)*F(p,p-1)),
            weight_diagnostics=diagnostics,totals=rows))
    selected,union=[],arb(0)
    for index,total in enumerate(totals):
        choice=min(controls,key=lambda c:F(c['totals'][index]['failure_log2_upper']))
        upper=F(choice['totals'][index]['failure_log2_upper'])
        union+=arb(2)**A.uniform.rational(upper)
        selected.append(dict(total_arity=total,moment_order=choice['moment_order'],failure_log2_upper=str(upper)))
    upper=A.uniform.rational_ceiling(A.log2(union))
    return dict(controls=controls,selected_totals=selected,all_totals_failure_log2_upper=str(upper),
                selected_class_40_bits=bool(upper<=-40))


def certificate(path,digest,*,precision=256,max_seconds=120,total_min=198,total_max=288,progress=None):
    if type(precision) is not int or precision not in (256,384): raise ValueError('Use256/384 bits')
    _total_range(total_min,total_max)
    pins=dependencies()
    started,check=old._timer(max_seconds)
    grouped,inventory_pin=checked_inventory(path,digest,check_time=check,progress=progress)
    with ctx.workprec(precision):
        reference,reference_pin=old.ideal.checked_reference(precision)
        result=controls_from_reference(grouped,reference,total_min=total_min,total_max=total_max,check_time=check)
        reference_text=str(reference)
    check()
    if dependencies()!=pins: raise ValueError('Sources changed during directed evaluation')
    if sha256(Path(inventory_pin['file']).read_bytes()).hexdigest()!=digest:
        raise ValueError('Inventory changed during directed evaluation')
    if sha256((DIRECTORY/reference_pin['file']).read_bytes()).hexdigest()!=reference_pin['sha256']:
        raise ValueError('Uniform reference changed during directed evaluation')
    return dict(result,complete=True,construction='scheduled-lambda-bc-staggered-column-shear-final-D2',
        active_groups=18,every_group_arity_range=[11,16],total_arity_range=[total_min,total_max],orders=list(ORDERS),
        source_profiles=99,field_order=256,source_dimension=1048576,output_length=2097152,
        forbidden_weight=734003,mark='691/700',uniform_alpha='1/5',
        all_labelled_original_rowsets=True,original_A_fixed_across_pairs=True,
        pair_count=8,minimum_active_pairs=5,exact_original_row_floors_retained=False,
        restricted_density_not_renormalized=True,source_count_outside_theta=True,
        rowset_sum_inside_theta=True,additional_rowset_binomial_factor=False,
        projective_normalizations=1,additional_beta_factor=False,fixed_p_within_each_total_bound=True,
        full_independent_D1_comparison=True,permitted_gauges_per_pair=1,final_independent_D2_after_shear=True,
        arithmetic_enclosures_directed=True,precision_bits=precision,failure_probability_upper_bound=True,
        mixed_arities_covered=list(ARITIES),other_group_counts_covered=False,whole_code_certificate=False,
        exact_inventory=inventory_pin,cached_reference=reference_pin,cached_uniform_log2_enclosure=reference_text,
        reviewed_source_sha256=pins,fresh_regional_banks=0,packing_searches=0,native_operations=0,
        elapsed_seconds=perf_counter()-started)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('inventory','certificate'))
    parser.add_argument('--max-attempts',type=int,default=32)
    parser.add_argument('--max-seconds',type=float,default=120)
    parser.add_argument('--precision',type=int,choices=(256,384),default=256)
    parser.add_argument('--total-min',type=int,default=198)
    parser.add_argument('--total-max',type=int,default=288)
    parser.add_argument('--inventory',type=Path)
    parser.add_argument('--inventory-sha256')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): parser.error('Refuse to overwrite evidence')
    progress=lambda row:print(json.dumps(row),flush=True)
    if args.mode=='inventory':
        result=inventory(max_attempts=args.max_attempts,max_seconds=args.max_seconds,progress=progress)
    else:
        if args.inventory is None or args.inventory_sha256 is None: parser.error('Inventory path/SHA required')
        result=certificate(args.inventory,args.inventory_sha256,precision=args.precision,
                           max_seconds=args.max_seconds,total_min=args.total_min,total_max=args.total_max,
                           progress=progress)
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps({k:v for k,v in result.items()
        if k not in ('rows','controls','selected_totals','reviewed_source_sha256')}))
