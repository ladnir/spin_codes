"""Matched-map spectrum backtests and exact conditional endpoint scores.

Numerical backtests are Q1 diagnostics. Conditional endpoint bounds use
exact rational scores and an explicit weighted-spectrum assumption.
"""
import argparse
import csv
from fractions import Fraction as F
import math
from pathlib import Path
import numpy as np
import bridge as base
import dual_track_q1 as model
import frontier_ledger as ledger


def reference_counts(block,dimension,distance,even,all_one):
    step=2 if even else 1
    weights=list(range(distance,block-distance+1 if all_one else block+1,step))
    mass=sum(math.comb(block,w) for w in weights)
    result={w:F(((1<<dimension)-1-int(all_one))*math.comb(block,w),mass) for w in weights}
    result[0]=F(1)
    if all_one:result[block]=F(1)
    assert sum(result.values(),F(0))==1<<dimension
    return result


def outer_inputs():
    theory=base.ROOT/'workstreams/finite_asymptotic_theory';small=theory/'small_k_replay/spectra'
    paths=[(8,4,True,True,small/'ebch8_4_spectrum.csv'),
           (32,8,True,True,theory/'ebch32_16_delta8_spectrum.csv'),
           (64,12,False,False,small/'xbch64_32_philips_spectrum.csv'),
           (128,22,True,True,small/'ebch128_weight_counts.json')]
    result=[]
    for block,distance,even,one,path in paths:
        if path.suffix=='.json':counts={int(w):int(n) for w,n in base.read(path)['weight_counts'].items()}
        else:
            with path.open(newline='',encoding='utf-8') as handle:counts={int(r['weight']):int(r['count']) for r in csv.DictReader(handle)}
        assert sum(counts.values())==1<<(block//2) and counts[0]==1
        assert min(w for w,n in counts.items() if w and n)==distance
        if even:assert all(w%2==0 for w,n in counts.items() if n)
        if one:assert all(counts.get(w,0)==counts.get(block-w,0) for w in range(block+1))
        result.append((block,distance,even,one,path,counts))
    return result


def coefficients(block,m,spectrum,count_arrays):
    rows=(1<<m)//(block//2);assert rows%64==0
    def values(grid):
        lam=np.exp(grid)/rows
        return math.log(rows)+np.minimum(0.,model.refresh.coefficient_logs(
            *model.refresh.epoch_logs(64,20,spectrum,lam),rows//64,block)+(block*rows//10)*lam[:,None])
    grid=np.arange(-40,81)/10;coarse=values(grid);indices=np.argmin(coarse,axis=0)
    best=coarse[indices,np.arange(block+1)]
    relevant=set()
    for counts in count_arrays:
        terms={w:math.log(n)+best[w] for w,n in counts.items() if w and n};maximum=max(terms.values())
        relevant.update(w for w,v in terms.items() if v>=maximum-30*math.log(2))
    assert all(0<indices[w]<len(grid)-1 for w in relevant),'Dominant witness reaches the coarse search edge'
    fine=np.unique(np.concatenate([grid[indices[w]]+np.arange(-30,31)/500 for w in sorted(relevant)]))
    return np.minimum(best,values(fine).min(axis=0))


def margin(value):return math.log2(value.denominator)-math.log2(value.numerator)


def conditional_endpoint(full_path,q1_path):
    full=base.read(full_path);m=full['message_exponent'];assert m==28
    assert ledger.build(m,[base.ROOT/r['file'] for r in full['certificates']])==full
    q1=base.read(q1_path)
    assert q1['message_exponent']==m and full['certificates'][0]['file']==q1_path.relative_to(base.ROOT).as_posix()
    counts=model.even_binomial()
    coefficients={int(w):base.decode(v) for w,v in q1['coefficient_upper'].items()}
    reference=sum((counts[w]*v for w,v in coefficients.items()),F(0))
    higher=base.decode(full['higher_occupancy_upper'])
    rows=[]
    for bits in (0,1,2,4,8,12,16):
        upper=(1<<bits)*reference+higher
        rows.append(dict(weighted_score_inflation_bits=bits,conditional_upper=base.encode(upper),
            conditional_margin_bits=margin(upper),higher_occupancy_penalty_bits=margin((1<<bits)*reference)-margin(upper)))
    return dict(status='EXACT_RATIONAL_CONDITIONAL_SCORE_NOT_AN_UNCONDITIONAL_CERTIFICATE',
        message_exponent=m,reference_q1_upper=base.encode(reference),reference_q1_margin_bits=margin(reference),
        higher_upper=base.encode(higher),higher_margin_bits=margin(higher),
        assumption='For the saved outward Q1 coefficients c_w, sum A_w*c_w <= 2^b * sum reference_A_w*c_w.',
        cases=rows)


def run(output):
    output=output.resolve();assert not output.exists()
    _,_,spectrum,_=model.maps.load('t64_s20');inputs=outer_inputs();results=[]
    for block,distance,even,one,path,counts in inputs+[(256,38,True,True,None,None)]:
        reference=reference_counts(block,block//2,distance,even,one)
        for m in (16,20,24,28):
            co=coefficients(block,m,spectrum,[reference]+([counts] if counts else []))
            predicted=model.aggregate(reference,co);actual=model.aggregate(counts,co) if counts else None
            ratios=[(w,float(counts[w]/reference[w])) for w in counts if w and counts[w] and w in reference] if counts else []
            terms={w:math.log(n)+co[w] for w,n in reference.items() if w and n}
            norm=float(np.logaddexp.reduce(list(terms.values())))
            results.append(dict(block_bits=block,dimension=block//2,minimum_distance=distance,message_exponent=m,
                model_q1_bits=predicted['margin_bits'],exact_spectrum_q1_bits=actual['margin_bits'] if actual else None,
                model_intercept_bits=predicted['margin_bits']+m-math.log2(block//2),
                exact_intercept_bits=actual['margin_bits']+m-math.log2(block//2) if actual else None,
                weighted_inflation_bits=predicted['margin_bits']-actual['margin_bits'] if actual else None,
                model_dominant_weight=predicted['dominant_weight'],exact_dominant_weight=actual['dominant_weight'] if actual else None,
                minimum_shell_count_ratio=dict(ratios).get(distance),
                model_first_two_shell_share=sum(math.exp(v-norm) for w,v in terms.items() if w in (distance,distance+(2 if even else 1))),
                evidence='BINARY64_Q1_ONLY_MATCHED_SELECTED_MAP'))
            print('B',block,'m',m,'model',round(predicted['margin_bits'],6),'exact',round(actual['margin_bits'],6) if actual else None,flush=True)
    full_path=base.HERE/'generated/curve_k28_full_retained_v1.json';q1_path=base.HERE/'generated/frontier_k28_q1_v1.json'
    endpoint=conditional_endpoint(full_path,q1_path)
    known=[r for r in results if r['exact_spectrum_q1_bits'] is not None]
    b128=next(r for r in results if r['block_bits']==128 and r['message_exponent']==28)
    b64=next(r for r in results if r['block_bits']==64 and r['message_exponent']==28)
    b256=next(r for r in results if r['block_bits']==256 and r['message_exponent']==28)
    linear=b128['exact_intercept_bits']+2*(b128['exact_intercept_bits']-b64['exact_intercept_bits'])
    sources=[Path(__file__),Path(model.__file__),Path(model.refresh.__file__),Path(model.maps.__file__),Path(ledger.__file__),
        full_path,q1_path,base.HERE/'generated/larger_state_inputs_v1/t64_s20_selection.json',
        *[r[4] for r in inputs],base.BCH/'SPIN_HEURISTIC_EVIDENCE.md',base.BCH/'LOCATOR_INCIDENCE_DUAL_TRACK.md']
    base.write_new(output,dict(status='MATCHED_SPECTRUM_BACKTEST_AND_CONDITIONAL_ENDPOINT',configuration='t64_s20',
        rows=results,conditional_endpoint=endpoint,
        largest_observed_weighted_inflation_bits=max(r['weighted_inflation_bits'] for r in known),
        naive_64_128_linear_intercept_forecast_at_256=linear,
        b256_reference_intercept_at_28=b256['model_intercept_bits'],
        limitations=['Backtests use known spectra but binary64 transfer bounds, not observed failure probabilities.',
            'The reference conditions only on dimension, minimum distance, parity and the all-one/complement property.',
            'B64 keeps odd weights and is a shortened constituent, not a like-for-like extended BCH rung.',
            'No statistical confidence interval or unbounded constituent-size extrapolation is claimed.',
            'K28 conditional cases assume a weighted-spectrum score bound; they do not assume componentwise bounds on every shell.',
            'The exact rational conditional score uses frozen outward tilt witnesses; backtests optimize a finer diagnostic witness grid.'],
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in sources}))
    print('Conditional endpoint',[(r['weighted_score_inflation_bits'],r['conditional_margin_bits']) for r in endpoint['cases']],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.output)
