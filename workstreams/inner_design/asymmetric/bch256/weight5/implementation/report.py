"""Bind the certified maps to tested sources and summarize serial measurements."""
import hashlib
import json
from pathlib import Path
import re
import statistics
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import verify_full
model=verify_full.model
ROOT=model.ROOT
sys.path.insert(0,str(ROOT/'workstreams/bare_bch_rm2sub'))
import generate as shared


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())


def main():
    proof=verify_full.verify()
    manifest=read(HERE/'IMPLEMENTATION.json')
    model.authenticate(manifest)
    assert manifest['instance']==proof['instance']
    assert manifest['certificate_sha256']==sha(HERE.parent/'FULL_M20_VERIFIED.json')
    for name,digest in manifest['generated_sha256'].items():
        assert sha(HERE/name)==digest,name
    maps=re.findall(r'columns\{([^}]+)\};',(HERE/'AsymmetricMap.h').read_text())
    assert len(maps)==2
    parsed=[[int(c,0) for c in row.split(',')] for row in maps]
    assert parsed==[proof['instance']['inner']['expansion_columns'],proof['instance']['inner']['feedback_columns']]
    header=(ROOT/'workstreams/bare_bch_rm2sub/generated/BchCircuit.h').read_text()
    rows=[sum(int(w,16)<<(64*j) for j,w in enumerate(re.findall(r'0x([0-9a-f]+)ULL',line)))
          for line in header.splitlines() if line.startswith('{0x')]
    p,q=shared.bch.generator_polynomial(37),shared.bch.generator_polynomial(39)
    def extend(w):return w|((w.bit_count()&1)<<255)
    expected=shared.reduce_rows([extend(q<<j) for j in range(123)]+[extend(p<<j) for j in range(5)])
    assert rows==expected and len(rows)==128
    measurements=HERE/'measurements'
    for filename in ('correctness-final.log','correctness-verbose.log'):
        assert '100% tests passed, 0 tests failed out of 12' in (measurements/filename).read_text()
    assert '100% tests passed, 0 tests failed out of 4' in (measurements/'sanitizer.log').read_text()
    assert '100% tests passed, 0 tests failed out of 1' in (measurements/'sanitizer-selected.log').read_text()
    remote_sources={}
    for line in (measurements/'sources.sha256').read_text().splitlines():
        digest,name=line.split(None,1)
        name=name.strip()
        assert sha(ROOT/name)==digest,name
        remote_sources[name]=digest
    for name in ('Weight5Spin.cpp','Weight5Inner.h','AsymmetricMap.h','CMakeLists.txt','inner_identity.cpp'):
        assert (HERE/name).relative_to(ROOT).as_posix() in remote_sources
    primary=('baseline','baseline_tuned','sparse','sparse_pages','sparse_tuned','grouped_tuned')
    variants=primary+('masked_tuned','shared_tuned','sparse_pf128')
    summaries={}
    for name in variants:
        paths=[measurements/f'confirm-{name}-{repeat}.jsonl' for repeat in (1,2,3)]
        if name in primary:
            paths=[measurements/f'final-{name}-{repeat}.jsonl' for repeat in (1,2,3)]+paths
        data=[read(p) for p in paths]
        for row in data:
            assert (row['m'],row['outer_length'],row['outer_dimension'],row['trials'],row['tile_rows'],row['layout'],row['inplace'])==(20,256,128,101,2048,'packed24',True)
            assert 0<row['p10_ms']<=row['median_ms']<=row['p90_ms']
        assert len({r['output_hash'] for r in data})==1
        assert len({r['retained_setup_bytes'] for r in data})==1
        assert all(r['workspace_bytes']==41943040 for r in data)
        flags=(measurements/f'flags-{name}.txt').read_text()
        assert '-O3' in flags and '-march=znver4' in flags
        summaries[name]=dict(median_ms=statistics.median(r['median_ms'] for r in data),processes=len(data),
            run_medians_ms=[r['median_ms'] for r in data],retained_setup_bytes=data[0]['retained_setup_bytes'],
            workspace_bytes=data[0]['workspace_bytes'],output_hash=data[0]['output_hash'])
    assert summaries['baseline']['output_hash']==summaries['baseline_tuned']['output_hash']
    assert len({summaries[n]['output_hash'] for n in variants[2:]})==1
    tune=[read(path) for path in measurements.glob('tune-*.jsonl')]
    assert len(tune)==40
    for configuration in {row['configuration'] for row in tune}:
        assert len({r['output_hash'] for r in tune if r['configuration']==configuration})==1
    selected=summaries['sparse_pages']['median_ms']
    source_paths=[HERE/'IMPLEMENTATION.json',HERE/'generate.py',Path(__file__),HERE/'run.sh',HERE/'tune.sh',HERE/'final.sh',HERE/'confirm.sh',HERE/'screen_final.sh']
    source_paths+=list(measurements.glob('*'))
    result=dict(status='CERTIFIED_WEIGHT5_IMPLEMENTATION_AND_SERIAL_PERFORMANCE_BOUND',
        selected='sparse_pages',summaries=summaries,proof_margin_bits=proof['margin_bits'],
        reduction_vs_supported_percent=100*(1-selected/summaries['baseline']['median_ms']),
        reduction_vs_routing_tuned_supported_percent=100*(1-selected/summaries['baseline_tuned']['median_ms']),
        exact_A_B_columns_bound=True,exact_outer_rows_bound=True,correctness_tests=12,sanitizer_tests=5,
        setup_and_workspace_construction_excluded=True,default_changed=False,
        source_sha256={**remote_sources,**{p.relative_to(ROOT).as_posix():sha(p) for p in source_paths if p.is_file()}})
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
