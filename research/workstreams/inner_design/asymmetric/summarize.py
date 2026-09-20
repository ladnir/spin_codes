"""Check artifact bindings and summarize measured and binary64 results.

This is an integrity check, not an outward certificate verifier. Numerical
screens and machine timings retain separate labels in the output.
"""
import hashlib
import json
import math
from pathlib import Path
from statistics import median

import search_feedback as search
g=search.g
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def bind(record):
    for name,digest in record['source_sha256'].items():assert sha(ROOT/name)==digest,name


def main():
    sources=[Path(__file__),HERE/'README.md',HERE/'TRANSFER_ARGUMENT.md']
    for name in ('FEEDBACK_SCREEN.json','OCCUPANCIES_Q8.json','OCCUPANCIES_Q128.json','DENSE_SCREEN.json','IMPLEMENTATION.json'):
        path=HERE/name;bind(read(path));sources.append(path)
    implementation=read(HERE/'IMPLEMENTATION.json')
    for candidate in implementation['candidates']:
        for name,digest in candidate['source_sha256'].items():
            path=HERE.parent/'generated'/candidate['name']/name
            assert sha(path)==digest; sources.append(path)
    exact=read(HERE/'FEEDBACK_SCREEN.json')
    a_columns=read(HERE.parent/'NO_CONSTANT_MAP.json')['columns']
    # Reconstruct integer map audits and retain low-weight cancellation summaries.
    cancellations=[]
    for candidate in exact['candidates']:
        rebuilt=search.audit(a_columns,candidate['b_columns'],19,candidate['name'])
        for key in ('rank','ba_zero','ba_rank','low_kernel','dual_minimum_weight','kernel_minimum_weight'):
            assert rebuilt[key]==candidate[key],(candidate['name'],key)
        for key in ('dual_spectrum','kernel_spectrum'):
            assert {str(k):v for k,v in rebuilt[key].items()}==candidate[key]
        low=search.low_cancellation(a_columns,candidate['b_columns'],19)
        rows=[]
        for j,data in low.items():
            total=math.comb(128,j)-int(candidate['kernel_spectrum'].get(str(j),0))
            assert data['nonzero_input_count']==total
            assert sum(sum(v.values()) for v in data['by_weight'].values())==total
            emitted=[w for hist in data['by_weight'].values() for w,count in hist.items() if count]
            rows.append(dict(input_weight=j,nonzero_syndrome_inputs=total,
                minimum_emitted_weight=min(emitted),maximum_emitted_weight=max(emitted),
                distinct_syndrome_histograms=len(data['patterns']),
                counts_by_A_syndrome_weight=data['by_weight']))
        cancellations.append(dict(candidate=candidate['name'],results=rows))
    sparse=read(HERE/'OCCUPANCIES_Q128.json');bounds=[]
    outer_weights=[w for w,n in enumerate(g.tv.smaller_outer.spectrum()) if w and n]
    for candidate in sparse['results']:
        for row in candidate['results']:
            delta=row['distance_target'];maximum=64 if delta=='33/200' else 128
            assert row['covered_occupations']==[1,maximum]
            assert [v['occupation'] for v in row['higher']]==list(range(3,maximum+1))
            path=HERE/f'DENSE_{candidate["name"]}_{delta.replace("/","_")}.json'
            cover=read(path);bind(cover);sources.append(path)
            assert cover['candidate']==candidate['name'] and cover['distance_target']==delta
            assert cover['occupation_min']==maximum+1 and cover['occupation_max']==32768
            assert cover['unresolved_leaves']==0 and cover['leaf_count']==len(cover['selected_boxes'])
            assert cover['bands']==[[w for w in outer_weights if w!=128],[128]]
            g.check_coverage(cover['selected_boxes'],32768,maximum+1)
            sm=row['partial_union_margin_bits'];dm=cover['dense_margin_bits']
            union=min(sm,dm)-math.log2(1+2**(-abs(sm-dm)))
            bounds.append(dict(candidate=candidate['name'],distance=delta,
                all_occupancy_margin_bits_diagnostic=union,dense_margin_bits=dm,
                dense_leaves=cover['leaf_count'],exact_partition_checked=True))
    measurements=HERE/'measurements';performance=[];kernel=[]
    variants=['balanced','mixed2_3_sparse','mixed2_3_masked','greedy3_2_sparse','greedy3_2_masked']
    for variant in variants:
        paths=[measurements/f'asymmetric-{variant}-repeat{i}.jsonl' for i in (1,2,3)]
        rows=[read(p) for p in paths];sources+=paths
        for row in rows:
            assert (row['m'],row['outer_length'],row['outer_dimension'],row['trials'])==(20,128,32,101)
            assert row['inplace'] and row['layout']=='packed24' and row['tile_rows']==4096
        assert len({row['output_hash'] for row in rows})==1
        performance.append(dict(variant=variant,repeat_medians_ms=[r['median_ms'] for r in rows],
            median_of_medians_ms=median(r['median_ms'] for r in rows),output_hash=rows[0]['output_hash'],
            retained_setup_bytes=rows[0]['retained_setup_bytes'],workspace_bytes=rows[0]['workspace_bytes']))
    baseline=performance[0]['median_of_medians_ms']
    for row in performance:row['time_reduction_percent']=100*(1-row['median_of_medians_ms']/baseline)
    for candidate in ('mixed2_3','greedy3_2'):
        hashes={r['output_hash'] for r in performance if r['variant'].startswith(candidate)}
        assert len(hashes)==1
    for variant in ('balanced','mixed2_3','greedy3_2'):
        paths=[measurements/f'asymmetric-kernel-{variant}-repeat{i}.jsonl' for i in (1,2,3)]
        trials=[[json.loads(line) for line in p.read_text().splitlines()] for p in paths];sources+=paths
        for i,n in enumerate((32768,4194304)):
            rows=[r[i] for r in trials]
            assert all(r['variant']==variant and r['blocks']==n and r['trials']==101 for r in rows)
            assert len({r['output_hash'] for r in rows})==1
            kernel.append(dict(variant=variant,blocks=n,repeat_medians_ms=[r['median_ms'] for r in rows],
                median_of_medians_ms=median(r['median_ms'] for r in rows),output_hash=rows[0]['output_hash']))
    log=measurements/'asymmetric-correctness.log'
    assert '100% tests passed, 0 tests failed out of 4' in log.read_text();sources.append(log)
    sources += [measurements/'asymmetric-binaries.sha256',measurements/'asymmetric-kernel-binaries.sha256']
    sources += [HERE/name for name in ('AsymmetricInner.h','kernel_benchmark.cpp','CMakeLists.txt','run_bench.sh','run_kernel_bench.sh','test_asymmetric.py')]
    result=dict(status='BOUNDED_EXPLORATION_COMPLETE_NOT_AN_OUTWARD_CERTIFICATE',
        selected_candidate='greedy3_2_sparse',measurement_host='Peach / Ryzen 9 7950X / CPU 15',
        compiler='GCC 15.2.0, Release -O3 -DNDEBUG, znver4',trials_per_run=101,
        full_encoder=performance,standalone_inner=kernel,all_occupancy_binary64=bounds,
        exact_low_weight_cancellation=cancellations,
        source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sources})
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','exact_low_weight_cancellation')},indent=2))


if __name__=='__main__':main()
