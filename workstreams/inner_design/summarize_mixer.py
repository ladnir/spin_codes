"""Check and summarize serial timings of the same one-round mixer map."""
import hashlib
import json
from pathlib import Path
import statistics

HERE=Path(__file__).resolve().parent


def main():
    records={};sources={}
    for path in sorted((HERE/'measurements/mixer').glob('mixer-*-repeat*.jsonl')):
        row,=[json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        name=path.stem.split('-')[1];records.setdefault(name,[]).append(row)
        sources[path.relative_to(HERE).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
    assert set(records)=={'baseline','rows','sparse','masked'}
    assert all(len(rows)==3 for rows in records.values())
    fields=('m','outer_length','outer_dimension','inplace','tile_rows','layout','trials','workspace_bytes')
    reference=records['baseline'][0]
    for rows in records.values():
        for row in rows:assert all(row[k]==reference[k] for k in fields)
    assert len({row['output_hash'] for name,rows in records.items() if name!='baseline' for row in rows})==1
    assert len({row['output_hash'] for row in records['baseline']})==1
    baseline=statistics.median(row['median_ms'] for row in records['baseline']);results=[]
    for name,rows in records.items():
        median=statistics.median(row['median_ms'] for row in rows)
        results.append(dict(name=name,run_medians_ms=[r['median_ms'] for r in rows],
                            median_of_medians_ms=median,relative_time=median/baseline,
                            retained_setup_bytes=rows[0]['retained_setup_bytes']))
    result=dict(status='SERIAL_NEW_MIXER_PERFORMANCE_NOT_DISTANCE_CERTIFICATE',
                workload={k:reference[k] for k in fields},source_sha256=sources,results=results)
    (HERE/'MIXER_PERFORMANCE.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
