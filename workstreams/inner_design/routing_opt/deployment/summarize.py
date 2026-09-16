"""Validate deployment receipts; never substitutes for a distance certificate."""
from pathlib import Path
import hashlib
import json
import statistics

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
DATA=HERE.parent/'measurements/deployment'


def main():
    count=0
    for line in (DATA/'sources.sha256').read_text().splitlines():
        digest,name=line.split(maxsplit=1)
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
        count+=1
    log=(DATA/'correctness.log').read_text()
    assert '100% tests passed, 0 tests failed out of 12' in log
    assert '100% tests passed, 0 tests failed out of 2' in (DATA/'public-correctness.log').read_text()
    for name in ('supported_off','supported_on','asymmetric_off','asymmetric_on','asymmetric_fallback','asymmetric_portable'):
        assert f'{name}_test' in log
    rows=[]
    for inner in ('supported','asymmetric'):
        for m in (16,18,20):
            cells={}
            hashes=set()
            for state in ('off','on'):
                runs=[json.loads((DATA/f'{inner}_{state}-m{m}-r{r}.jsonl').read_text()) for r in (1,2,3)]
                assert all(r['m']==m and r['inplace'] and r['trials']==101 and r['layout']=='packed24' for r in runs)
                hashes.update(r['output_hash'] for r in runs)
                cells[state]=dict(medians_ms=[r['median_ms'] for r in runs],median_ms=statistics.median(r['median_ms'] for r in runs),workspace_bytes=runs[0]['workspace_bytes'],retained_setup_bytes=runs[0]['retained_setup_bytes'])
            assert len(hashes)==1
            assert cells['off']['workspace_bytes']==cells['on']['workspace_bytes']
            assert cells['off']['retained_setup_bytes']==cells['on']['retained_setup_bytes']
            rows.append(dict(inner=inner,m=m,output_hash=next(iter(hashes)),**cells,time_reduction_percent=100*(1-cells['on']['median_ms']/cells['off']['median_ms'])))
    result=dict(status='OPTIONAL_SIZE_DEPENDENT_ROUTING_VALIDATED',default_enabled=False,
                policy='Linux, quarter rate, K>=2^18; all other cells retain original map instantiation',
                source_files_verified=count,test_targets_passed=12,public_build_test_targets_passed=2,results=rows,
                certificate_status='No new distance claim; unchanged code maps; prior balanced certificate integrity passes using archived historical build inputs',
                limitations=['Three processes of 101 in-place calls per cell; setup and workspace preparation excluded.',
                             'Unavailable-page and non-Linux control paths simulated in dedicated tests; no second hardware platform measured.',
                             'The optional kernel is not activated at K=2^16; small timing differences between builds remain.'])
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in rows:print(row['inner'],row['m'],row['off']['median_ms'],row['on']['median_ms'],row['time_reduction_percent'])


if __name__=='__main__':main()
