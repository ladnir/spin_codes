"""Matched serial timings; generated samples stay in the ignored measurements tree."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('root',type=Path)
a=p.parse_args()
out=a.root/'workstreams/spin_optimized/measurements/k16-r2'
out.mkdir(parents=True,exist_ok=False)
cases=[('old12819','build-k16-r2/base/spin_benchmark','12819'),
       ('sub12r1','build-k16-integrated/spin_benchmark','6412'),
       ('sub12r2','build-k16-r2/base/spin_benchmark','6412')]
rows=[]
for seed in (1,17):
    for repeat in range(5):
        ordered=cases if repeat%2==0 else list(reversed(cases))
        for label,binary,config in ordered:
            proc=subprocess.run([str(a.root/binary),'16','auto','101',str(seed),'0','0',config],
                                check=True,capture_output=True,text=True)
            row=json.loads(proc.stdout)
            (out/f'{label}-s{seed}-r{repeat}.json').write_text(proc.stdout)
            rows.append(dict(label=label,seed=seed,repeat=repeat,**row))
            print(label,seed,repeat,row['median_ms'],flush=True)
summary=[]
for label,binary,_ in cases:
    for seed in (1,17):
        group=[r for r in rows if r['label']==label and r['seed']==seed]
        assert len({r['output_hash'] for r in group})==1
        summary.append(dict(label=label,seed=seed,median_ms=statistics.median(r['median_ms'] for r in group),
            process_medians_ms=[r['median_ms'] for r in group],
            retained_setup_bytes=group[0]['retained_setup_bytes'],workspace_bytes=group[0]['workspace_bytes']))
result=dict(summary=summary,binary_sha256={b:hashlib.sha256((a.root/b).read_bytes()).hexdigest() for _,b,_ in cases},
    policy='serial; CPU15; 5 processes x 101 calls; 3 warmups; in-place; setup excluded; no input copy/reset')
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
