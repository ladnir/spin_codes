"""Validate and summarize final-map serial performance receipts.

This report binds the generated map and implementation, not a distance claim.
The separate all-occupancy certificate supplies the mathematical statement.
"""
import hashlib
import json
from pathlib import Path
import statistics

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    directory=HERE/'measurements/balanced'
    manifest=json.loads((HERE/'BALANCED_IMPLEMENTATION.json').read_text())
    for name,digest in manifest['source_sha256'].items():assert sha(ROOT/name)==digest,name
    for candidate in manifest['candidates']:
        for name,digest in candidate['source_sha256'].items():
            assert sha(HERE/'generated'/candidate['name']/name)==digest,name
    environment=json.loads((directory/'ENVIRONMENT.json').read_text())
    for name,digest in environment['remote_source_sha256'].items():assert sha(HERE/name)==digest,name
    correctness=(directory/'balanced-correctness.log').read_text()
    assert '100% tests passed, 0 tests failed out of 2' in correctness
    for m in (16,18,20):assert correctness.count(f'm={m} PASS')==2
    records={};sources=[Path(__file__),HERE/'BALANCED_IMPLEMENTATION.json',HERE/'NO_CONSTANT_MAP.json',
                        HERE/'run_balanced_bench.sh',directory/'ENVIRONMENT.json',
                        directory/'balanced-correctness.log',directory/'balanced-binaries.sha256']
    for variant in ('baseline','sparse','masked'):
        rows=[]
        for run in (1,2,3):
            path=directory/f'balanced-{variant}-repeat{run}.jsonl';sources.append(path)
            row,=[json.loads(line) for line in path.read_text().splitlines() if line.strip()]
            rows.append(row)
        records[variant]=rows
    fields=('m','outer_length','outer_dimension','inplace','tile_rows','layout','trials','workspace_bytes')
    reference=records['baseline'][0]
    assert reference['m']==20 and reference['trials']==101 and reference['inplace'] is True
    for rows in records.values():
        for row in rows:assert all(row[k]==reference[k] for k in fields)
    assert len({r['output_hash'] for name,rows in records.items() if name!='baseline' for r in rows})==1
    assert {r['output_hash'] for r in records['baseline']}=={'26d0b43de7571c4c'}
    baseline=statistics.median(r['median_ms'] for r in records['baseline'])
    results=[]
    for name,rows in records.items():
        median=statistics.median(r['median_ms'] for r in rows)
        results.append(dict(name=name,run_medians_ms=[r['median_ms'] for r in rows],
             median_of_medians_ms=median,relative_time=median/baseline,
             paired_time_reduction_percent=[100*(1-r['median_ms']/b['median_ms']) for r,b in zip(rows,records['baseline'])],
             retained_setup_bytes=rows[0]['retained_setup_bytes'],output_hash=rows[0]['output_hash']))
    result=dict(status='SERIAL_BALANCED_MAP_PERFORMANCE',workload={k:reference[k] for k in fields},
                source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sources},results=results)
    (HERE/'BALANCED_PERFORMANCE.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
