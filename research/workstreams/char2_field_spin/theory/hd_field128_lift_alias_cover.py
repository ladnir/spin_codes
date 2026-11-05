"""Bounded fresh full-field coverage from reusable directed alias banks.

Band mode evaluates every requested integer weight and every feasible
matching histogram. Box mode verifies a supplied disjoint clipped box
partition and replays each complete box. Old manifests supply geometry
and exact witness proposals only; no old probability ceiling is reused.
"""

import ast
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
from time import monotonic

from flint import arb,ctx

import hd_field128_lift_certificate as base
import hd_field128_lift_tail_certificate as tail
import hd_field128_lift_alias_certificate as interval
import hd_field128_lift_acb_alias_certificate as native


def source_pins():
    directory = Path(__file__).resolve().parent
    pending,seen = [Path(__file__).name],{}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        path = directory/name
        raw = path.read_bytes()
        seen[name] = sha256(raw).hexdigest()
        for node in ast.walk(ast.parse(raw,filename=str(path))):
            modules = ([a.name for a in node.names] if isinstance(node,ast.Import)
                       else [node.module] if isinstance(node,ast.ImportFrom) and node.module else [])
            for module in modules:
                target = module.split('.')[0]+'.py'
                if (directory/target).is_file() and target not in seen:
                    pending.append(target)
    return dict(sorted(seen.items()))


def read_json(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    return json.loads(raw),dict(path=str(path),sha256=sha256(raw).hexdigest())


def witness_from_artifact(artifact,*,old=False):
    """Read exact literals only; numerical results are deliberately ignored."""
    if old:
        fixed = artifact['fixed_witness']
        if fixed.get('pair_swaps') is not True or fixed.get('orientation_aggregation')!='independent-fair-swaps':
            raise ValueError('Old proposal must explicitly use independent fair pair swaps')
        row = fixed['source_exact_witness']
        markers = tuple(tuple(base.exact(x) for x in pair.split()) for pair in row['markers'])
        return dict(z=base.exact(row['z'])**16,alpha=base.exact(row['alpha'])/16,
                    markers=markers,ratio=base.exact(row['potential']),grid=fixed.get('grid',128))
    if (artifact.get('construction')!=base.CONSTRUCTION or artifact.get('field_order')!=base.FIELD_ORDER
            or artifact.get('pair_swaps') is not True or artifact.get('full_field_random_coefficients') is not True):
        raise ValueError('A literal full-field construction witness is required')
    row = artifact['witness']
    return dict(z=base.exact(row['z']),alpha=base.exact(row['alpha']),
                markers=tuple(tuple(map(base.exact,pair)) for pair in row['markers']),
                ratio=base.exact(artifact['state_potential_ratio']),grid=row.get('grid',128))


def manifest_seeds(path,first=65,last=5853):
    if type(first) is not int or type(last) is not int or not 1<=first<=last<=base.GROUPS:
        raise ValueError('Invalid seed clipping range')
    manifest,pin = read_json(path)
    entries = manifest['entries']
    if len(entries)!=37:
        raise ValueError('Expected the old37-entry manifest geometry')
    cursor,result = 1,[]
    for entry in entries:
        lo,hi = entry['weight_range']
        if lo!=cursor or not lo<=hi<=base.GROUPS:
            raise ValueError('Old manifest weight ranges are not a disjoint full partition')
        cursor = hi+1
        if max(first,lo)<=min(last,hi):
            result.append(dict(id=entry['id'],kind=entry['kind'],weight_range=[max(first,lo),min(last,hi)],
                               proposal=entry.get('proposal')))
    if cursor!=base.GROUPS+1:
        raise ValueError('Old geometry leaves a weight gap')
    return dict(source=pin,requested_weight_range=[first,last],seeds=result,
                old_numerical_ceilings_reused=False,directed=False,whole_code_certificate=False)


class RegionalCache:
    """Reuse restored residues and update histogram normalization by recurrence."""

    def __init__(self,bank):
        self.bank = bank
        self.numerator = lru_cache(maxsize=max(1024,2*bank.grid*len(bank.markers)))(self._numerator)

    def _numerator(self,index,single_residue,full_residue):
        atlas = self.bank.atlases[index]
        if getattr(atlas.upper,'ndim',None)==3:
            values = (arb(float(atlas.upper[i,single_residue,full_residue]))*atlas.factors[i] for i in range(4))
        else:
            position = single_residue*self.bank.grid+full_residue
            values = (atlas.upper[i][position]*atlas.factors[i] for i in range(4))
        return tuple(min(base.upper(value),base.upper(fallback))
                     for value,fallback in zip(values,self.bank.numerators[index]))

    def regional(self,weight,budget):
        bank = self.bank
        counts = tail.full_counts(base.PAIRS,weight)
        result = {}
        first = counts.start
        initial_single = weight-2*first
        initial_count = base.comb(base.PAIRS,first)*base.comb(base.PAIRS-first,initial_single)
        for index,(x,y) in enumerate(bank.markers):
            budget.check()
            scale = base.upper(1/(initial_count*base.rational(x)**initial_single*base.rational(y)**first))
            factor = base.rational(x*x/y)
            for full in counts:
                if not full%128:
                    budget.check()
                single = weight-2*full
                numerator = self.numerator(index,single%bank.grid,full%bank.grid)
                row = tuple(base.upper(value*scale) for value in numerator)
                previous = result.get(full)
                result[full] = row if previous is None else tuple(min(a,b) for a,b in zip(row,previous))
                if full+1<counts.stop:
                    empty = base.PAIRS-weight+full
                    scale = base.upper(scale*factor*((full+1)*(empty+1))/(single*(single-1)))
        return result


def union_ceiling(rows,precision):
    if not rows:
        return None
    with ctx.workprec(precision):
        value = sum((arb(2)**base.rational(F(row['log2_upper_rational'])) for row in rows),arb(0))
        return str(base.rational_ceiling(tail.log2arb(value)))


def checked_range(first,last,max_weights):
    if (any(type(x) is not int for x in (first,last,max_weights))
            or not 1<=first<=last<=base.GROUPS or not 1<=max_weights<=base.GROUPS
            or last-first+1>max_weights):
        raise ValueError('Weight interval exceeds its explicit resource scope')


def band_receipt(bank,first,last,ratios,*,method='runs',max_weights=256,max_seconds=600.,
                 failure_bits=40,progress=None):
    checked_range(first,last,max_weights)
    ratios = tuple(dict.fromkeys(map(base.exact,ratios)))
    if not ratios or len(ratios)>32 or min(ratios)<=0 or type(failure_bits) is not int or failure_bits<0:
        raise ValueError('Use1..32 positive exact potentials and a nonnegative integer target')
    if method not in ('runs','marginal'):
        raise ValueError('Unknown support summation method')
    pins,budget = source_pins(),base.Budget(max_seconds)
    cache,rows,failure = RegionalCache(bank),[],None
    with ctx.workprec(bank.precision):
        for weight in range(first,last+1):
            try:
                budget.check()
                regions = cache.regional(weight,budget)
                candidates = []
                for ratio in ratios:
                    support = tail.support_sum(base.PAIRS,weight,regions,ratio,method=method,budget=budget)
                    value = bank.constant+bank.slope*weight+tail.log2arb(support)
                    candidates.append((base.rational_ceiling(value),ratio,str(value)))
                ceiling,ratio,enclosure = min(candidates,key=lambda item:item[0])
                rows.append(dict(active_groups=weight,regional_histograms=len(regions),
                                 log2_upper_rational=str(ceiling),log2_upper_enclosure=enclosure,
                                 state_potential_ratio=str(ratio),covers_all_supports_at_weight=True,
                                 covers_all_even_odd_splits=True))
                if progress:
                    progress(dict(stage='complete-weight',active_groups=weight,log2_upper_rational=str(ceiling),
                                  completed_weights=len(rows),requested_weights=last-first+1))
            except (TimeoutError,ArithmeticError) as error:
                failure = dict(weight=weight,error_type=type(error).__name__,message=str(error))
                break
        ceiling = union_ceiling(rows,bank.precision)
    bank.check_pins()
    if source_pins()!=pins:
        raise RuntimeError('Coverage helper changed during replay')
    next_weight = first+len(rows)
    complete = next_weight>last
    return dict(bank.metadata(),coverage_driver_sha256=pins,requested_weight_range=[first,last],
                weight_receipts=rows,weights_evaluated=len(rows),every_integer_weight_evaluated=complete,
                completed_weight_range=[first,next_weight-1] if rows else None,
                pending_weight_ranges=[] if complete else [[next_weight,last]],
                complete_requested_range=complete,completed_union_log2_upper_rational=ceiling,
                log2_upper_rational=ceiling if complete else None,
                requested_range_meets_target=complete and F(ceiling)<=-failure_bits,
                weights_above_target=[row['active_groups'] for row in rows if F(row['log2_upper_rational'])>-failure_bits],
                failure_bits=failure_bits,support_method=method,potential_candidates=list(map(str,ratios)),
                elapsed_seconds=monotonic()-budget.started,max_seconds=max_seconds,max_weights=max_weights,
                failure=failure,whole_code_certificate=False,
                scope='completed integer weights only; no endpoint interpolation or pending-weight claim')


def checked_partition(first,last,boxes):
    """Count clipped integer profiles after proving all box interiors disjoint."""
    checked_range(first,last,base.GROUPS)
    root = ((first,last),(0,base.PAIRS),(0,base.PAIRS))
    boxes = tuple(tail.validate_box(row) for row in boxes)
    if not 1<=len(boxes)<=512:
        raise ValueError('A partition of1..512 boxes is required')
    counts = []
    for index,box in enumerate(boxes):
        if box[0][0]<first or box[0][1]>last:
            raise ValueError('A box exceeds the requested weight band')
        count = tail.profile_count(box)
        if not count:
            raise ValueError('Empty clipped boxes must be omitted')
        counts.append(count)
        for earlier in boxes[:index]:
            intersection = tuple((max(a,c),min(b,d)) for (a,b),(c,d) in zip(box,earlier))
            if all(a<=b for a,b in intersection) and tail.profile_count(intersection):
                raise ValueError('Clipped boxes overlap')
    total = tail.profile_count(root)
    if sum(counts)!=total:
        raise ValueError('Clipped boxes do not cover every profile in the requested band')
    return boxes,total


def box_partition_receipt(bank,first,last,rows,*,max_seconds=600.,failure_bits=40,progress=None):
    if type(failure_bits) is not int or failure_bits<0:
        raise ValueError('Integer nonnegative failure target required')
    boxes,total = checked_partition(first,last,[row['ranges'] for row in rows])
    indices = [tuple(row['marker_indices']) for row in rows]
    if any(len(pair)!=2 or any(type(i) is not int or not 0<=i<len(bank.markers) for i in pair) for pair in indices):
        raise ValueError('Each box needs two fixed valid marker indices')
    pins,budget = source_pins(),base.Budget(max_seconds)
    receipts,failure = [],None
    for index,(box,markers) in enumerate(zip(boxes,indices)):
        try:
            budget.check()
            remaining = max_seconds-(monotonic()-budget.started)
            if remaining<=0:
                raise TimeoutError('Box partition exhausted its time budget')
            receipt = bank.box_receipt(box,marker_indices=markers,max_seconds=remaining)
            receipts.append({key:receipt[key] for key in ('ranges','marker_indices','integer_profile_triples',
                                                          'log2_upper_rational','log2_upper_enclosure')})
            if progress:
                progress(dict(stage='complete-box',index=index,log2_upper_rational=receipt['log2_upper_rational']))
        except (TimeoutError,ArithmeticError) as error:
            failure = dict(box_index=index,error_type=type(error).__name__,message=str(error))
            break
    bank.check_pins()
    if source_pins()!=pins:
        raise RuntimeError('Coverage helper changed during box replay')
    ceiling = union_ceiling(receipts,bank.precision)
    complete = len(receipts)==len(boxes)
    return dict(bank.metadata(),coverage_driver_sha256=pins,requested_weight_range=[first,last],
                partition_profiles=total,partition_boxes=len(boxes),box_receipts=receipts,
                pending_boxes=[dict(ranges=box,marker_indices=marker) for box,marker in zip(boxes[len(receipts):],indices[len(receipts):])],
                complete_requested_range=complete,completed_union_log2_upper_rational=ceiling,
                log2_upper_rational=ceiling if complete else None,
                requested_range_meets_target=complete and F(ceiling)<=-failure_bits,
                failure_bits=failure_bits,failure=failure,whole_code_certificate=False,
                scope='completed clipped boxes only; the supplied partition is checked independently')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    seeds = sub.add_parser('seeds',help='Inspect old37manifest geometry; no numerical work')
    seeds.add_argument('--manifest',type=Path,default=Path(__file__).with_name('hd_shifted_closed_cover_manifest.json'))
    seeds.add_argument('--first',type=int,default=65)
    seeds.add_argument('--last',type=int,default=5853)
    for name in ('band','boxes'):
        task = sub.add_parser(name)
        source = task.add_mutually_exclusive_group(required=True)
        source.add_argument('--witness',type=Path,help='Full-field receipt; exact witness literals only')
        source.add_argument('--old-proposal',type=Path,help='Old exact saved witness, lifted with z^16 and alpha/16')
        task.add_argument('--first',type=int,required=True)
        task.add_argument('--last',type=int,required=True)
        task.add_argument('--ratio',action='append',help='Positive exact candidate; may repeat up to32times')
        task.add_argument('--backend',choices=('interval','acb'),default='interval')
        task.add_argument('--precision',type=int,default=256)
        task.add_argument('--bits',type=int,default=80)
        task.add_argument('--grid',type=int)
        task.add_argument('--distance',default='7/20')
        task.add_argument('--failure-bits',type=int,default=40)
        task.add_argument('--build-seconds',type=float,default=600.)
        task.add_argument('--max-seconds',type=float,default=600.)
        task.add_argument('--output',type=Path,required=True)
        if name=='band':
            task.add_argument('--max-weights',type=int,default=256)
            task.add_argument('--method',choices=('runs','marginal'),default='runs')
        else:
            task.add_argument('--partition',type=Path,required=True,help='JSON list of ranges+marker_indices rows')
    args = parser.parse_args()
    if args.command=='seeds':
        print(json.dumps(manifest_seeds(args.manifest,args.first,args.last),indent=2))
        return
    if args.output.exists():
        parser.error('Output exists; choose a new path')
    artifact,input_pin = read_json(args.old_proposal or args.witness)
    witness = witness_from_artifact(artifact,old=args.old_proposal is not None)
    ratios = tuple(map(F,args.ratio)) if args.ratio else (witness['ratio'],)
    progress = lambda row:print(json.dumps(row),flush=True)
    cls = interval.AliasBank if args.backend=='interval' else native.AliasBank
    bank = cls(witness['z'],witness['alpha'],witness['markers'],distance=F(args.distance),
               grid=args.grid or witness['grid'],precision=args.precision,bits=args.bits,
               max_seconds=args.build_seconds,progress=progress)
    if args.command=='band':
        result = band_receipt(bank,args.first,args.last,ratios,method=args.method,max_weights=args.max_weights,
                              max_seconds=args.max_seconds,failure_bits=args.failure_bits,progress=progress)
    else:
        partition,partition_pin = read_json(args.partition)
        result = box_partition_receipt(bank,args.first,args.last,partition,max_seconds=args.max_seconds,
                                       failure_bits=args.failure_bits,progress=progress)
        result['partition_artifact'] = partition_pin
    result['witness_artifact'] = input_pin
    result['witness_artifact_use'] = 'exact literals only; no numerical ceiling imported'
    with args.output.open('x',encoding='utf-8') as target:
        json.dump(result,target,indent=2,allow_nan=False)
        target.write('\n')
    print(json.dumps({key:result[key] for key in ('requested_weight_range','complete_requested_range',
                                                   'completed_union_log2_upper_rational','requested_range_meets_target','failure')}))


if __name__=='__main__':
    main()
