"""Serial bounded fresh replays of old potential-band geometry.

Each old exact witness supplies its complete marker bank. The new z and
alpha are short decimal rationals obtained by rounding z_old**16 and
alpha_old/16; these are new admissible proposals, not equalities claimed
for the exact lift. Every band is regenerated for the full-field setup.
Finished bands and witnesses are saved immediately with exclusive creates.
No old probability ceiling is read into a new probability calculation.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from time import monotonic

import hd_field128_lift_certificate as base
import hd_field128_lift_alias_cover as cover


DEFAULT_FACTORS = (F(1,2),F(3,4),F(1),F(5,4),F(3,2))


def source_pins():
    # These are the only local imports above. cover's recursive closure
    # already includes base and both directed alias implementations.
    pins = cover.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def rounded_witness(artifact,source):
    exact = cover.witness_from_artifact(artifact,old=True)
    z = F(str(float(exact['z'])))
    alpha = F(str(float(exact['alpha'])))
    if not 0<z<1 or not 0<alpha<=1:
        raise ValueError('Rounded proposal left the admissible parameter range')
    return dict(schema='field128-lift-alias-witness-only-v1',
                status='unexecuted rounded proposal; no numerical result or coverage claim',
                construction=base.CONSTRUCTION,field_order=base.FIELD_ORDER,field_bits=128,
                pair_swaps=True,full_field_random_coefficients=True,directed=False,whole_code_certificate=False,
                source=source,source_witness_role='exact proposal literals only; old numerical bounds ignored',
                transform=dict(z='Fraction(str(float(Fraction(old_z)**16)))',
                               alpha='Fraction(str(float(Fraction(old_alpha)/16)))',
                               equality_claim='none; rounded values are new arbitrary exact witnesses',
                               markers='every original positive exact marker pair, unchanged'),
                state_potential_ratio=str(exact['ratio']),
                witness=dict(z=str(z),alpha=str(alpha),markers=[list(map(str,row)) for row in exact['markers']],
                             grid=exact['grid']))


def load_jobs(manifest,first,last):
    if type(first) is not int or type(last) is not int or not 257<=first<=last<=4067:
        raise ValueError('This batch covers clipped old potential bands within257..4067 only')
    plan = cover.manifest_seeds(manifest,first,last)
    directory = Path(manifest).resolve().parent
    jobs,cursor = [],first
    for seed in plan['seeds']:
        if seed['kind']!='potential-band' or seed['weight_range'][0]!=cursor:
            raise ValueError('Requested interval must have contiguous potential-band seeds')
        proposal = seed['proposal']
        artifact,pin = cover.read_json(directory/proposal['path'])
        if pin['sha256']!=proposal['sha256']:
            raise ValueError('Pinned old proposal changed: '+proposal['path'])
        witness = rounded_witness(artifact,pin)
        jobs.append(dict(id=seed['id'],weight_range=seed['weight_range'],witness=witness))
        cursor = seed['weight_range'][1]+1
    if cursor!=last+1:
        raise ValueError('The planned batch leaves a weight gap')
    return jobs,plan['source']


def write_new(path,value):
    encoded = json.dumps(value,indent=2,allow_nan=False)+'\n'
    with Path(path).open('x',encoding='utf-8',newline='\n') as stream:
        stream.write(encoded)
    return dict(path=str(Path(path).resolve()),sha256=sha256(encoded.encode('utf-8')).hexdigest())


def missing_intervals(first,last,covered):
    result,start = [],None
    for weight in range(first,last+1):
        if weight not in covered:
            if start is None:
                start = weight
        elif start is not None:
            result.append([start,weight-1])
            start = None
    if start is not None:
        result.append([start,last])
    return result


def run_batch(manifest,output_dir,*,first=257,last=4067,backend='interval',precision=256,bits=80,
              factors=DEFAULT_FACTORS,max_bands=4,max_seconds=3600.,build_seconds=180.,band_seconds=300.,
              distance=F(7,20),failure_bits=40,method='runs',progress=None):
    factors = tuple(dict.fromkeys(map(base.exact,factors)))
    distance = base.exact(distance)
    if (not factors or len(factors)>32 or min(factors)<=0 or backend not in ('interval','acb')
            or type(max_bands) is not int or not 1<=max_bands<=64
            or type(precision) is not int or not 128<=precision<=4096
            or type(bits) is not int or not 32<=bits<=512 or not 0<distance<F(1,2)
            or type(failure_bits) is not int or failure_bits<0 or method not in ('runs','marginal')
            or any(not isfinite(x) or not 0<x<=3600 for x in (build_seconds,band_seconds))):
        raise ValueError('Invalid explicit batch parameters or resource budgets')
    budget = base.Budget(max_seconds)
    jobs,manifest_pin = load_jobs(manifest,first,last)
    pins = source_pins()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=False,exist_ok=False)
    completed,attempts,rows,failure = [],[],[],None
    sequence = 0

    def snapshot(final=False):
        nonlocal sequence
        ceiling = cover.union_ceiling(rows,precision)
        covered = {row['active_groups'] for row in rows}
        pending = missing_intervals(first,last,covered)
        complete = not pending
        summary = dict(schema='field128-lift-alias-batch-v1',construction=base.CONSTRUCTION,
                       field_order=base.FIELD_ORDER,field_bits=128,target_distance=str(distance),
                       forbidden_weight=int(distance*base.OUTPUT_LENGTH),outer_measure='invertible-toeplitz',
                       pair_swaps=True,full_field_random_coefficients=True,directed=True,whole_code_certificate=False,
                       source_sha256=pins,manifest_source=manifest_pin,old_numerical_ceilings_reused=False,
                       requested_weight_range=[first,last],planned_bands=[dict(id=job['id'],weight_range=job['weight_range']) for job in jobs],
                       completed_bands=completed.copy(),attempts=attempts.copy(),weight_receipts=rows.copy(),
                       covered_weight_count=len(covered),duplicate_weight_count=len(rows)-len(covered),
                       pending_weight_ranges=pending,complete_requested_range=complete,
                       completed_union_log2_upper_rational=ceiling,
                       log2_upper_rational=ceiling if complete else None,
                       requested_range_meets_target=complete and F(ceiling)<=-failure_bits,
                       weights_above_target=[row['active_groups'] for row in rows if F(row['log2_upper_rational'])>-failure_bits],
                       failure_bits=failure_bits,potential_factors=list(map(str,factors)),
                       precision=precision,packet_mantissa_bits=bits,backend=backend,support_method=method,
                       max_bands=max_bands,max_seconds=max_seconds,build_seconds=build_seconds,band_seconds=band_seconds,
                       elapsed_seconds=monotonic()-budget.started,failure=failure,
                       scope='completed per-weight events only; pending weights have no probability claim')
        name = 'final_summary.json' if final else f'snapshot_{sequence:03d}.json'
        saved = write_new(output_dir/name,summary)
        sequence += 1
        if progress:
            progress(dict(stage='batch-snapshot',path=saved['path'],completed_bands=len(completed),
                          covered_weight_count=len(covered),pending_weight_ranges=pending,
                          completed_union_log2_upper_rational=ceiling,failure=failure))
        return summary

    snapshot()
    cls = cover.interval.AliasBank if backend=='interval' else cover.native.AliasBank
    for index,job in enumerate(jobs[:max_bands]):
        try:
            budget.check()
            stem = f'band_{index:03d}_{job["weight_range"][0]}_{job["weight_range"][1]}'
            witness_pin = write_new(output_dir/(stem+'_witness.json'),job['witness'])
            witness = cover.witness_from_artifact(job['witness'])
            ratios = tuple(witness['ratio']*factor for factor in factors)
            attempt = dict(id=job['id'],weight_range=job['weight_range'],witness=witness_pin,status='building')
            attempts.append(attempt)
            report = (lambda row:progress(dict(row,band_id=job['id'])) if progress else None)
            remaining = min(build_seconds,max_seconds-(monotonic()-budget.started))
            if remaining<=0:
                raise TimeoutError('Batch deadline exhausted before bank construction')
            bank = cls(witness['z'],witness['alpha'],witness['markers'],distance=distance,grid=witness['grid'],
                       precision=precision,bits=bits,max_seconds=remaining,progress=report)
            remaining = min(band_seconds,max_seconds-(monotonic()-budget.started))
            if remaining<=0:
                raise TimeoutError('Batch deadline exhausted after bank construction')
            attempt['status'] = 'evaluating'
            lo,hi = job['weight_range']
            receipt = cover.band_receipt(bank,lo,hi,ratios,method=method,max_weights=hi-lo+1,
                                          max_seconds=remaining,failure_bits=failure_bits,progress=report)
            receipt['batch_witness_artifact'] = witness_pin
            receipt['old_manifest_band_id'] = job['id']
            receipt_pin = write_new(output_dir/(stem+'_receipt.json'),receipt)
            rows.extend(dict(row,band_id=job['id'],receipt_sha256=receipt_pin['sha256']) for row in receipt['weight_receipts'])
            attempt['status'] = 'complete' if receipt['complete_requested_range'] else 'partial'
            attempt['receipt'] = receipt_pin
            completed.append(dict(id=job['id'],weight_range=job['weight_range'],receipt=receipt_pin,
                                  complete_requested_range=receipt['complete_requested_range'],
                                  completed_union_log2_upper_rational=receipt['completed_union_log2_upper_rational']))
            del bank
            snapshot()
            if not receipt['complete_requested_range']:
                failure = dict(id=job['id'],reason='partial-band',detail=receipt['failure'])
                break
        except (TimeoutError,ArithmeticError,RuntimeError) as error:
            failure = dict(id=job['id'],weight_range=job['weight_range'],error_type=type(error).__name__,message=str(error))
            write_new(output_dir/f'failure_{index:03d}.json',failure)
            break
    if source_pins()!=pins:
        raise RuntimeError('A pinned batch helper changed during execution')
    result = snapshot(final=True)
    return result


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,default=Path(__file__).with_name('hd_shifted_closed_cover_manifest.json'))
    parser.add_argument('--output-dir',type=Path,required=True,help='New directory; every artifact is exclusive-create')
    parser.add_argument('--first',type=int,default=257)
    parser.add_argument('--last',type=int,default=4067)
    parser.add_argument('--backend',choices=('interval','acb'),default='interval')
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--bits',type=int,default=80)
    parser.add_argument('--factor',action='append',help='Exact potential multiplier; default1/2,3/4,1,5/4,3/2')
    parser.add_argument('--max-bands',type=int,default=4)
    parser.add_argument('--max-seconds',type=float,default=3600.)
    parser.add_argument('--build-seconds',type=float,default=180.)
    parser.add_argument('--band-seconds',type=float,default=300.)
    parser.add_argument('--distance',default='7/20')
    parser.add_argument('--failure-bits',type=int,default=40)
    parser.add_argument('--method',choices=('runs','marginal'),default='runs')
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('Output directory exists; choose a new directory')
    result = run_batch(args.manifest,args.output_dir,first=args.first,last=args.last,backend=args.backend,
                       precision=args.precision,bits=args.bits,factors=tuple(map(F,args.factor)) if args.factor else DEFAULT_FACTORS,
                       max_bands=args.max_bands,max_seconds=args.max_seconds,build_seconds=args.build_seconds,
                       band_seconds=args.band_seconds,distance=F(args.distance),failure_bits=args.failure_bits,method=args.method,
                       progress=lambda row:print(json.dumps(row),flush=True))
    print(json.dumps({key:result[key] for key in ('complete_requested_range','covered_weight_count',
                                                  'pending_weight_ranges','requested_range_meets_target','failure')}))


if __name__=='__main__':
    main()
