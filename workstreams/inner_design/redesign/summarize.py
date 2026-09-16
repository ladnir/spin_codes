"""Bind the bounded redesign exploration; not a distance certificate verifier."""
import argparse
import hashlib
import json
from pathlib import Path
from statistics import median
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--partial',action='store_true');args=parser.parse_args()
    sources=[p for p in HERE.iterdir() if p.suffix in ('.py','.cpp','.h','.sh','.md') or p.name=='CMakeLists.txt']
    for name in ('MACRO_Q1.json','MACRO_IMPLEMENTATION.json','ROUTING_IMPLEMENTATION.json','HEADROOM_IMPLEMENTATION.json','LOCALITY_COUNTEREXAMPLE.json'):
        path=HERE/name;record=read(path);sources.append(path)
        for filename,digest in record['source_sha256'].items():assert sha(ROOT/filename)==digest,filename
        if name in ('MACRO_IMPLEMENTATION.json','ROUTING_IMPLEMENTATION.json'):
            for candidate in record['candidates']:
                for filename,digest in candidate['generated_sha256'].items():
                    path=HERE.parent/'generated'/candidate['name']/filename
                    assert sha(path)==digest;sources.append(path)
        elif name=='HEADROOM_IMPLEMENTATION.json':
            for candidate,files in record['generated_sha256'].items():
                for filename,digest in files.items():
                    path=HERE.parent/'generated'/candidate/filename
                    assert sha(path)==digest;sources.append(path)
    sets={
        'headroom':['asymmetric_greedy3_2_sparse','redesign_no_inner','redesign_no_inner_unrolled','redesign_outer_only'],
        'macro':['asymmetric_greedy3_2_sparse','redesign_macro256_r2','redesign_macro512_r2'],
        'routing':['asymmetric_greedy3_2_sparse','redesign_gather_pf0','redesign_gather_pf64'],
        'accumulator':['control','2','3'],
        'accumulator-bucket':['control','2']}
    results=[];missing=[]
    for suite,variants in sets.items():
        suite_rows=[]
        for variant in variants:
            paths=[HERE/'measurements'/f'redesign-{suite}-{variant}-repeat{i}.jsonl' for i in (1,2,3)]
            for p in paths:
                if not p.exists():missing.append(p.name)
            paths=[p for p in paths if p.exists()];rows=[read(p) for p in paths];sources+=paths
            if not rows:continue
            assert all(r['m']==20 and r['trials']==101 and r['inplace'] for r in rows)
            assert len({r['output_hash'] for r in rows})==1
            if suite.startswith('accumulator') and variant in ('2','3'):
                assert all(r['correctness']=='PASS' and r['rounds']==int(variant) for r in rows)
            suite_rows.append(dict(variant=variant,repeat_count=len(rows),repeat_medians_ms=[r['median_ms'] for r in rows],
                median_of_medians_ms=median(r['median_ms'] for r in rows),output_hash=rows[0]['output_hash'],
                setup_bytes=rows[0].get('retained_setup_bytes',rows[0].get('setup_bytes')),
                workspace_bytes=rows[0]['workspace_bytes']))
        control=suite_rows[0]['median_of_medians_ms']
        for row in suite_rows:row['time_reduction_percent']=100*(1-row['median_of_medians_ms']/control)
        if suite=='routing':assert len({r['output_hash'] for r in suite_rows})==1
        results.append(dict(suite=suite,rows=suite_rows))
    for suite,count in (('headroom',3),('macro',2),('routing',2)):
        p=HERE/'measurements'/f'redesign-{suite}-correctness.log'
        if p.exists():
            assert f'100% tests passed, 0 tests failed out of {count}' in p.read_text();sources.append(p)
        else:missing.append(p.name)
    sources+=list((HERE/'measurements').glob('*.sha256'))
    measured=HERE/'measurements/redesign-measured-sources.sha256'
    if measured.exists():
        for line in measured.read_text().splitlines():
            digest,name=line.split(maxsplit=1)
            assert sha(ROOT/name)==digest,('compiled source mismatch',name)
    else:missing.append(measured.name)
    # The optimized and direct accumulator schedules implement the same map.
    direct=next(r for r in results if r['suite']=='accumulator')['rows']
    bucket=next(r for r in results if r['suite']=='accumulator-bucket')['rows']
    if len(bucket)>1:assert direct[1]['output_hash']==bucket[1]['output_hash']
    if not args.partial:assert not missing,missing
    witness=read(HERE/'LOCALITY_COUNTEREXAMPLE.json')
    assert witness['outgoing_state']==0 and witness['message_hex']!='0x0'
    assert witness['full_output_weight']*200<33*witness['output_bits']
    payload=dict(status='PARTIAL_EXPLORATION' if missing else 'BOUNDED_EXPLORATION_COMPLETE_NO_NEW_WINNER',
        precision_status='NO_NEW_DISTANCE_CERTIFICATE',measured_host='Peach, Ryzen 9 7950X, CPU 15, GCC 15.2.0, Release znver4',
        message_bits=1<<20,outer=[128,32,32],results=results,missing_repeat_or_test_receipts=missing,
        locality_counterexample_weight=witness['full_output_weight'],locality_counterexample_relative_weight=witness['relative_weight'],
        source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sources})
    (HERE/'RESULT.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in payload.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
