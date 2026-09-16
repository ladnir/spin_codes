"""Validate exact-map benchmark receipts and source binding, then summarize."""
import hashlib
import json
import statistics
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
DATA=HERE/'measurements'


def main():
    source_count=0
    for line in ((DATA/'compiled-sources.sha256').read_text()+(DATA/'combined-sources.sha256').read_text()).splitlines():
        digest,name=line.split(maxsplit=1)
        path=ROOT/name.strip()
        assert hashlib.sha256(path.read_bytes()).hexdigest()==digest,name
        source_count+=1
    for name in ('correctness.log','correctness-refine.log','correctness-hier.log','correctness-pages.log','correctness-combined.log'):
        text=(DATA/name).read_text()
        assert '100% tests passed, 0 tests failed' in text,name
        assert 'm=16 PASS' in text and 'm=18 PASS' in text and 'm=20 PASS' in text,name
    final={}
    for name in ('control','huge2048','huge4096','write32'):
        rows=[json.loads((DATA/f'final-{name}-{r}.jsonl').read_text()) for r in (1,2,3)]
        assert all(r['trials']==101 and r['m']==20 and r['inplace'] and r['output_hash']=='6f9101ba915d544a' for r in rows)
        final[name]=dict(medians_ms=[r['median_ms'] for r in rows],median_ms=statistics.median(r['median_ms'] for r in rows),workspace_bytes=rows[0]['workspace_bytes'],retained_setup_bytes=rows[0]['retained_setup_bytes'])
        if name.startswith('huge'):
            for r in (1,2,3):
                log=(DATA/f'final-{name}-{r}.log').read_text()
                assert log.count('hint=0 errno=0 collapse=0 errno=0')==2
    for value in final.values():
        value['time_reduction_percent']=100*(1-value['median_ms']/final['control']['median_ms'])
    combined={}
    for name in ('asymmetric_greedy3_2_sparse','routeopt_pages_both','routeopt_pages_bothwrite'):
        rows=[json.loads((DATA/f'combined-{name}-{r}.jsonl').read_text()) for r in (1,2,3)]
        assert all(r['trials']==101 and r['m']==20 and r['inplace'] and r['output_hash']=='6f9101ba915d544a' for r in rows)
        combined[name]=dict(medians_ms=[r['median_ms'] for r in rows],median_ms=statistics.median(r['median_ms'] for r in rows))
        if name.startswith('routeopt_pages'):
            for r in (1,2,3):
                assert (DATA/f'combined-{name}-{r}.log').read_text().count('hint=0 errno=0 collapse=0 errno=0')==2
    for value in combined.values():
        value['time_reduction_percent']=100*(1-value['median_ms']/combined['asymmetric_greedy3_2_sparse']['median_ms'])
    sizes={}
    for m in (16,18):
        sizes[m]={}
        expected='9e7e3b0a16dd8537' if m==16 else '4d07a0efc0e493a8'
        for name in ('asymmetric_greedy3_2_sparse','routeopt_pages_bothwrite'):
            rows=[json.loads((DATA/f'sizes-{m}-{name}-{r}.jsonl').read_text()) for r in (1,2,3)]
            assert all(r['trials']==101 and r['m']==m and r['inplace'] and r['output_hash']==expected for r in rows)
            sizes[m][name]=dict(medians_ms=[r['median_ms'] for r in rows],median_ms=statistics.median(r['median_ms'] for r in rows))
        sizes[m]['time_reduction_percent']=100*(1-sizes[m]['routeopt_pages_bothwrite']['median_ms']/sizes[m]['asymmetric_greedy3_2_sparse']['median_ms'])
    screens=[]
    for path in sorted(DATA.glob('*.jsonl')):
        if path.name.startswith(('final-','perf-','profile-','combined-','sizes-')):continue
        row=json.loads(path.read_text())
        assert row['trials']==31 and row['output_hash']=='b180fa625fba455f',path
        screens.append(dict(receipt=path.name,**row))
    result=dict(status='EXACT_MAP_OPTIONAL_LINUX_WORKSPACE_OPTIMIZATION',
        new_distance_certificate=False,source_files_verified=source_count,
        final=final,combined=combined,smaller_sizes=sizes,screens=screens,
        limitations=['Linux huge-page availability and this host affect the gain.',
                     'Workspace preparation, including huge-page collapse, is excluded from encoding time.',
                     'Profile clocks perturb timing; hardware counters cover the entire process, including setup.',
                     'Supported implementation and prior certificate artifacts were not changed.'])
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(final=final,combined=combined),indent=2))


if __name__=='__main__':main()
