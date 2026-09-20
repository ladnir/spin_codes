"""Summarize serial full-encoder comparisons; no benchmark is launched."""
import hashlib
import json
from pathlib import Path
import statistics

HERE = Path(__file__).resolve().parent


def main():
    records = {}
    sources = {}
    for path in sorted((HERE/'measurements').glob('*-repeat*.jsonl')):
        row, = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        name = path.stem.rsplit('-repeat', 1)[0]
        records.setdefault(name, []).append(row)
        sources[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    assert len(records) == 5 and all(len(rows) == 3 for rows in records.values())
    fields = ('configuration', 'm', 'outer_length', 'outer_dimension', 'inplace',
              'tile_rows', 'layout', 'trials', 'output_hash')
    reference = records['baseline'][0]
    assert all(all(row[key] == reference[key] for key in fields)
               for rows in records.values() for row in rows)
    baseline = statistics.median(row['median_ms'] for row in records['baseline'])
    result = []
    for name, rows in records.items():
        median = statistics.median(row['median_ms'] for row in rows)
        result.append(dict(name=name, run_medians_ms=[row['median_ms'] for row in rows],
                           median_of_medians_ms=median, relative_time=median/baseline))
    payload = dict(status='SERIAL_BENCHMARK_SUMMARY_NOT_A_SIGNIFICANCE_TEST',
                   workload={key: reference[key] for key in fields},
                   source_sha256=sources, results=result,
                   decision='Keep unchanged baseline: no material repeatable improvement.')
    (HERE/'BASIS_PERFORMANCE.json').write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
