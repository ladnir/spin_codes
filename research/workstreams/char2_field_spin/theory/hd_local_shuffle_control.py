"""FLOAT screen of independent uniform S16 output shuffles per outer column.

This is a NEW ensemble: each group's16 columns independently sample an
output permutation after D2 Hstar D1. Setup is then fixed for all messages.
It does not use the inverse input schedule. No native code is changed.

For fixed original row set A, let C_A count the outside image and Z_A the
actual zero-outside inside fiber, of dimension d. Averaging a uniform
output permutation and using atoms in[0,1] gives the positive majorant
sum_w (Z_A[w]+256^d C_A[w]1_{w>0}) e_w(t)/choose(16,w).
It has exact constant1 and mass256^a. Average independent shuffles before
products over groups; only then take the r-th moment in the common measure.
The original A remains fixed across all16 columns.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_bit_syndrome_counted_general import ENDPOINTS
from hd_bit_syndrome_scalar_allarity import prepare_geometry
from hd_inverse_column_control import INPUT, INPUT_SHA
from hd_post_refresh_quiet_rows import log_matmul, log_matrix_power, moment_atoms
from hd_two_sum_gf16_state_envelope import elementary_logs

DIRECTORY = Path(__file__).resolve().parent


def shuffled_coefficients(outside, zero_fiber, dimension, *, size=16, q=256):
    if (not outside or not zero_fiber or outside[0]!=1 or zero_fiber[0]!=1
            or any(type(v) is not int or v<0 for v in (*outside,*zero_fiber))
            or not 0<=dimension<=size or len(outside)>size+1 or len(zero_fiber)>size+1):
        raise ValueError('Exact image and zero-fiber enumerators required')
    if sum(zero_fiber)!=q**dimension:
        raise ValueError('Zero-fiber dimension disagrees')
    result = []
    for w in range(size+1):
        numerator = zero_fiber[w] if w<len(zero_fiber) else 0
        if w and w<len(outside):
            numerator += q**dimension*outside[w]
        result.append(Fraction(numerator,comb(size,w)))
    if result[0]!=1 or sum(result[w]*comb(size,w) for w in range(size+1)) != q**dimension*sum(outside):
        raise AssertionError('Mass or exact zero changed')
    return tuple(result)


@lru_cache(maxsize=17)
def count_geometry(length):
    if type(length) is not int or not 0<=length<=16:
        raise ValueError('Count length0..16 required')
    if length == 0:
        return ((0,0,0,0,0),),()
    previous,layers = count_geometry(length-1)
    counts = tuple((a,b,c,d,length-a-b-c-d) for a in range(length+1)
        for b in range(length-a+1) for c in range(length-a-b+1) for d in range(length-a-b-c+1))
    lookup = {row:i for i,row in enumerate(counts)}
    targets = []
    for label in range(5):
        indices = []
        for row in previous:
            added = list(row)
            added[label] += 1
            indices.append(lookup[tuple(added)])
        targets.append(np.asarray(indices,dtype=np.int32))
    return counts,layers+((len(counts),tuple(targets)),)


def path_measure(measure, *, length=16):
    masses,empty = moment_atoms(measure)
    atoms = (-np.inf,measure['log2_atom00'],measure['log2_atom01'],measure['log2_atom10'],0.)
    if max(atoms)>0:
        raise ValueError('Output-weight monotonicity requires every atom<=1')
    counts,layers = count_geometry(length)
    current = np.full((1,2,2),-np.inf)
    current[0,0,0] = current[0,1,1] = 0.
    for size,targets in layers:
        following = np.full((size,2,2),-np.inf)
        for label,((incoming,outgoing),mass) in enumerate(zip(ENDPOINTS,masses)):
            following[targets[label],:,outgoing] = np.logaddexp2(
                following[targets[label],:,outgoing],current[:,:,incoming]+mass)
        current = following
    elementary = np.full((len(counts),length+1),-np.inf)
    for index,row in enumerate(counts):
        marks = [atom for label,atom in enumerate(atoms) if np.isfinite(atom) for _ in range(row[label])]
        elementary[index,:len(marks)+1] = elementary_logs(np.asarray(marks))
    if not np.allclose(np.logaddexp2.reduce(current,axis=0),log_matrix_power(empty,length),atol=1e-10):
        raise ArithmeticError('Count DP does not sum to the common constant measure')
    return current,elementary,log_matrix_power(empty,length)


def column_moment(coefficients, paths, elementary, groups, *, nonzero=False):
    logs = np.asarray([log2(c.numerator)-log2(c.denominator) if c else -np.inf for c in coefficients])
    if nonzero:
        logs[0] = -np.inf
    kernel = np.logaddexp2.reduce(elementary+logs[None,:],axis=1)
    return np.logaddexp2.reduce(paths+groups*kernel[:,None,None],axis=0)


def first_active_words(full, active, empty, groups, *, columns=16, minimum=9):
    current = np.full((minimum+1,*full.shape),-np.inf)
    current[0,...,0,0] = current[0,...,1,1] = 0.
    for _ in range(columns):
        following = np.full_like(current,-np.inf)
        following[minimum] = log_matmul(current[minimum],full)
        for d in range(minimum):
            following[d] = np.logaddexp2(following[d],log_matmul(current[d],empty))
            following[d+1] = np.logaddexp2(following[d+1],log_matmul(current[d],active))
        current = following
    return (np.logaddexp2.reduce(current[minimum,...,0,:],axis=-1)
            +(groups-1)*log2(comb(columns,minimum)))


def run(arities, cache_path, *, max_seconds=240, progress=lambda x:None):
    started = perf_counter()
    arities = tuple(sorted(set(arities)))
    if not arities or any(a not in range(1,17) for a in arities):
        raise ValueError('A nonempty subset of arities1..16 required')
    raw = (DIRECTORY/INPUT).read_bytes()
    if sha256(raw).hexdigest()!=INPUT_SHA:
        raise ValueError('Common input receipt changed')
    saved = json.loads(raw)
    for name,digest in saved['reviewed_sources'].items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Frozen common dependency changed: '+name)
    def check_time():
        if perf_counter()-started>max_seconds:
            raise TimeoutError('Local-shuffle screen exceeded its budget')
    masks = tuple(sum(1<<i for i in row) for a in arities for row in combinations(range(16),a))
    geometry = prepare_geometry(source_masks=masks,max_masks=len(masks),check_time=check_time)
    classes = Counter()
    for item in geometry['classes']:
        coeff = shuffled_coefficients(item['counts_by_outputweight'],item['zero_fiber_counts'],item['d'])
        classes[item['arity'],coeff] += item['multiplicity']
    progress(dict(stage='exact_geometry',classes=len(classes),elapsed_seconds=perf_counter()-started))
    paths,elementary,empty = path_measure(saved['state_measure'])
    full,active = [],[]
    for index,(a,coeff) in enumerate(classes):
        check_time()
        full.append(column_moment(coeff,paths,elementary,saved['groups']))
        active.append(column_moment(coeff,paths,elementary,saved['groups'],nonzero=True))
    full,active = np.stack(full),np.stack(active)
    words = first_active_words(full,active,empty,saved['groups'])
    prefix,r = saved['common_prefix_log2'],saved['groups']
    per = {}
    for a in arities:
        terms = [log2(multiplicity)+word/r for ((b,_),multiplicity),word in zip(classes.items(),words) if b==a]
        per[a] = prefix+r*np.logaddexp2.reduce(terms)-8*a*r*log2(255)
    mixed = prefix+r*np.logaddexp2.reduce([(v-prefix)/r for v in per.values()])
    with Path(cache_path).open('xb') as handle:
        np.savez_compressed(handle,full=full,active=active,words=words,
            arities=np.array([a for a,_ in classes]),multiplicities=np.array(list(classes.values())))
    return dict(construction='Hstar with independent uniform S16 output shuffle per group and column',
        groups=r,source_class=saved['source_class'],source_count=saved['source_count'],mark=saved['mark'],
        arities=list(arities),per_arity_log2=per,mixed_log2=mixed,common_prefix_log2=prefix,
        exact_classes=len(classes),labelled_rowsets=len(masks),all_arities=(arities==tuple(range(1,17))),
        active_column_floor=9,uniform_first9_cost_bits=(r-1)*log2(comb(16,9)),
        coefficient_classes=[dict(arity=a,multiplicity=m,coefficients=[str(c) for c in coeff])
                             for (a,coeff),m in classes.items()],
        input_sha256=INPUT_SHA,cache=str(Path(cache_path).resolve()),
        cache_sha256=sha256(Path(cache_path).read_bytes()).hexdigest(),
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        helper_sha256={name:sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in
            ('hd_bit_syndrome_scalar_allarity.py','hd_bit_syndrome_flat_geometry.py',
             'hd_post_refresh_quiet_rows.py','hd_two_sum_gf16_state_envelope.py')},
        elapsed_seconds=perf_counter()-started,directed=False,whole_code_certificate=False,
        new_ensemble=True,native_implemented=False,native_benchmarks=0,bad_setup_charge_required=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arities',type=int,nargs='+',default=list(range(1,17)))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cache',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists() or args.cache.exists():
        parser.error('Refuse to overwrite saved evidence')
    answer = run(args.arities,args.cache,progress=lambda x:print(json.dumps(x),flush=True))
    with args.output.open('x',encoding='utf-8') as handle:
        json.dump(answer,handle,indent=2,allow_nan=False)
        handle.write('\n')
    print(json.dumps({k:answer[k] for k in ('per_arity_log2','mixed_log2','elapsed_seconds')}),flush=True)
