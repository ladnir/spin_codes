"""Directed coverage search for common-refresh candidates at 35% distance.

Search emits exact rational witnesses. Replay never searches. The reported
margin includes the union of all covered classes and the diagonal density
loss; floating-point proposal results never count as certified coverage.
"""
from fractions import Fraction as F
from flint import arb, ctx

from hd35_design_screen import FORBIDDEN_WEIGHT, LENGTH, proposal
from hd_randomized_outer_proof import Parameters, dense_tail_receipt
from hd_sumrank_paired_certificate import receipt as paired
from hd_sumrank_global_certificate import receipt as global_bound

OUTER = 128
THRESHOLD = -88


def evaluate(witness, precision, outer_length=OUTER, local_rank=4):
    kind, values = witness
    options = dict(precision=precision, threshold=THRESHOLD,
                   feedthrough='common'+str(local_rank), outer_length=outer_length,
                   forbidden_weight=FORBIDDEN_WEIGHT)
    if kind == 'paired':
        return paired(*values, **options)
    if kind == 'regional-single':
        return paired(*values, paired=False, **options)
    if kind == 'global':
        return global_bound(*values, paired=True, **options)
    if kind == 'global-single':
        return global_bound(*values, paired=False, **options)
    raise ValueError('Unknown witness kind')


def initial(outer_length=OUTER):
    rows=LENGTH//outer_length
    p = Parameters(rows=rows, outer_n=outer_length, outer_k=outer_length//2,
                   width=64, state=16, limit=FORBIDDEN_WEIGHT)
    tail = dense_tail_receipt(p, failure_bits=96)
    covered = set(range(tail['active_rows_first'], rows+1))
    return tail, covered, arb(2)**arb(tail['union_log2_upper'])


def finish(covered, union, precision, witnesses, outer_length=OUTER, local_rank=4):
    log2 = lambda x: x.log()/arb(2).log()
    loss = (LENGTH//64)*log2((arb(65536)**16-1)/arb(65535)**16)
    final = log2(union)+loss
    missing = [r for r in range(1, LENGTH//outer_length+1) if r not in covered]
    return dict(construction='common-diagonal-Cauchy-refresh',
                outer_n=outer_length, outer_k=outer_length//2, local_rank=local_rank, message_length=LENGTH//2,
                field_order=65536, packet_width=64, state=16,
                forbidden_weight=FORBIDDEN_WEIGHT, precision=precision,
                ideal_union_log2_upper=str(log2(union)), density_loss_bits=str(loss),
                union_log2_upper=str(final), uncovered_rows=missing,
                directed=True, whole_code_certificate=not missing and bool(final < -40),
                witnesses=witnesses)


def accept(item, covered, union):
    intervals = item.get('certified_intervals', [])
    if not intervals:
        raise ArithmeticError('Directed witness accepted no classes')
    for first, last in intervals:
        covered.update(range(first, last+1))
    return union+arb(2)**arb(item['union_log2_upper'])


def _evaluate_job(job):
    index,witness,precision,outer_length,local_rank=job
    return index,witness,evaluate(witness,precision,outer_length,local_rank)


def replay(witnesses, precision=256, progress=None, outer_length=OUTER, local_rank=4, workers=1, allow_rejected=False):
    if not 1<=workers<=16:raise ValueError('Use one to sixteen proof workers')
    from concurrent.futures import ProcessPoolExecutor
    with ctx.workprec(precision):
        tail, covered, union = initial(outer_length)
        rejected=[]
        jobs=[(i,w,precision,outer_length,local_rank) for i,w in enumerate(witnesses)]
        pool=ProcessPoolExecutor(max_workers=workers) if workers>1 else None
        try:
            results=pool.map(_evaluate_job,jobs) if pool else map(_evaluate_job,jobs)
            for index,witness,item in results:
                if allow_rejected and not item.get('certified_intervals'):
                    rejected.append(index)
                else:union = accept(item, covered, union)
                if progress:
                    progress(dict(index=index, witness=witness, receipt=item, covered=len(covered)))
        finally:
            if pool:pool.shutdown()
        result=finish(covered, union, precision, witnesses,outer_length,local_rank)
        result['rejected_witness_indices']=rejected
        return result


def search(precision=256, progress=None, outer_length=OUTER, local_rank=4, start=1, regional_single=False):
    with ctx.workprec(precision):
        tail, covered, union = initial(outer_length)
        witnesses = []
        rows=LENGTH//outer_length
        # A partial search must retain its missing prefix in the final result.
        if not 1<=start<=rows:raise ValueError('Invalid search start')
        while len(covered) < rows:
            first = next((r for r in range(start, rows+1) if r not in covered),None)
            if first is None:break
            regional = regional_single or first < 512
            result = proposal(first, family='common' if regional and not regional_single else 'common-single', outer_length=outer_length,
                              regional=regional, iterations=90,local_rank=local_rank)
            z = F(result['mark']).limit_denominator(1000000)
            alpha = F(result['alpha']).limit_denominator(1000000)
            marker = F(1/(1+2**(-result['marker_log_odds']))).limit_denominator(10000000)
            values = [str(z), str(alpha), str(marker)]
            if regional:
                last = min(rows, first+max(4, first//3))
                witness = ['regional-single' if regional_single else 'paired', [first, last]+values]
            else:
                witness = ['global-single', [first]+values]
            item = evaluate(witness, precision,outer_length,local_rank)
            if not any(lo <= first <= hi for lo, hi in item.get('certified_intervals', [])):
                raise ArithmeticError(f'No coverage progress at {first}: {item}')
            witnesses.append(witness)
            union = accept(item, covered, union)
            if progress:
                progress(dict(index=len(witnesses), witness=witness, receipt=item,
                              covered=len(covered), screen=result['log2_upper_screen']))
        return finish(covered, union, precision, witnesses,outer_length,local_rank)


if __name__ == '__main__':
    import argparse
    import json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--search', action='store_true')
    p.add_argument('--precision', type=int, default=256)
    p.add_argument('--outer-length', type=int, choices=(128,256), default=128)
    p.add_argument('--local-rank', type=int, choices=(1,2,4), default=4)
    p.add_argument('--start',type=int,default=1)
    p.add_argument('--regional-single',action='store_true',help='Retain exact regional coefficients of single-packet profile moments')
    args = p.parse_args()
    if not args.search:
        p.error('Use --search, or run hd35_common128_witnesses.py for pinned replay')
    def progress(item):
        receipt=item.pop('receipt')
        item.update(intervals=receipt.get('certified_intervals'),bound=receipt.get('union_log2_upper'))
        print(json.dumps(item),flush=True)
    print(json.dumps(search(args.precision, progress,args.outer_length,args.local_rank,args.start,args.regional_single)), flush=True)
