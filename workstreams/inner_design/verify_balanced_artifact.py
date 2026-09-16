"""Bind the full certificate, exact map, tested source, and timing receipts.

This fast integrity check does not replace the full higher-precision replay.
It requires that replay's receipt and checks exact dyadic unions independently.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re

import generate_balanced as generated
import general_occupancies as g

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def check_dyadic_union(terms,upper,target):
    # Common power-of-two denominator: avoid repeated million-bit rational
    # gcds when very small dense bounds are added to the sparse contribution.
    exponent=min([int(v['exponent']) for v in terms]+[int(upper['exponent']),-target])
    assert all(int(v['mantissa'])>=0 for v in terms)
    total=sum(int(v['mantissa'])<<(int(v['exponent'])-exponent) for v in terms)
    ceiling=int(upper['mantissa'])<<(int(upper['exponent'])-exponent)
    assert total<=ceiling<(1<<(-target-exponent))


def sources(record,root,historical_build=False):
    for name,digest in record['source_sha256'].items():
        normalized=name.replace('\\','/')
        path=root/normalized
        if historical_build and normalized in ('workstreams/bare_bch_rm2sub/Spin.cpp',
                                               'workstreams/bare_bch_rm2sub/CMakeLists.txt') and sha(path)!=digest:
            # Historical build inputs only. Keep each original manifest digest;
            # live candidates and numerical inputs never use a substitution.
            path=HERE/'source_snapshots'/digest/path.name
        assert sha(path)==digest,name


def partition(bands,weights):
    flattened=[w for band in bands for w in band]
    assert sorted(flattened)==sorted(weights) and len(flattened)==len(set(flattened))


def main():
    certificate_path=HERE/'NO_CONSTANT_MARGIN_CERTIFICATE.json'
    replay_path=HERE/'NO_CONSTANT_MARGIN_CERTIFICATE_REPLAY.json'
    certificate=read(certificate_path);replay=read(replay_path)
    assert certificate['status']=='OUTWARD_ALL_OCCUPANCY_CERTIFICATE'
    assert certificate['precision_bits']==256
    assert (certificate['message_bits'],certificate['output_bits'])==(1<<20,1<<22)
    assert certificate['inner']==dict(t=128,s=19,transvection_rounds=1,map='NO_CONSTANT_MAP.json')
    assert replay['status']=='HIGHER_PRECISION_REPLAY_PASSED' and replay['precision_bits']==512
    assert replay['certificate_sha256']==sha(certificate_path)
    sources(certificate,ROOT)
    record=read(HERE/'NO_CONSTANT_MAP.json');sources(record,ROOT)
    original=read(ROOT/'workstreams/rate_quarter_bch/SMALLER_OUTWARD_WITNESSES.json')
    weights=[w for w,n in enumerate(g.tv.smaller_outer.spectrum()) if w and n]
    for key in ('q2_bands','adaptive_bands'):partition(original[key],weights)
    assert len(certificate['results'])==2
    for row,delta,minimum,target in zip(certificate['results'],('33/200','19/100'),(65,129),(40,30)):
        assert row['distance_target']==delta and row['target_bits']==target
        assert row['bad_weight']==int(F(delta)*(1<<22)) and row['maximum_sparse']==minimum-1
        cover=read(HERE/f'NO_CONSTANT_DENSE_COVER_{delta.replace("/","_")}.json')
        partition(cover['bands'],weights)
        assert cover['bands']==[[w for w in weights if w!=128],[128]]
        g.check_coverage(cover['selected_boxes'],32768,minimum)
        assert len(row['dense_terms'])==len(cover['selected_boxes'])
        assert len(row['adaptive_terms'])==minimum-3
        terms=[row['q1_upper']]+row['q2_terms']+row['adaptive_terms']+row['dense_terms']
        check_dyadic_union(terms,row['union_upper'],target)
    manifest=read(HERE/'BALANCED_IMPLEMENTATION.json');sources(manifest,ROOT)
    old=read(HERE/'MIXER_IMPLEMENTATION.json');sources(old,ROOT,historical_build=True)
    header,stats=generated.header(record['columns'])
    for candidate in manifest['candidates']:
        directory=HERE/'generated'/candidate['name']
        for name,digest in candidate['source_sha256'].items():assert sha(directory/name)==digest,name
        assert (directory/'BalancedMap.h').read_text()==header
        assert all(candidate[k]==v for k,v in stats.items())
    quarter=ROOT/'workstreams/rate_quarter_bch/implementation/generated'
    outer=read(quarter/'MANIFEST.json');sources(outer,ROOT,historical_build=True)
    assert outer['outer']==certificate['outer']==g.tv.smaller_outer.construction()
    for name,digest in outer['generated_sha256'].items():assert sha(quarter/name)==digest,name
    rows=[int(a,16)|(int(b,16)<<64) for a,b in re.findall(r'\{0x([0-9a-f]+)ULL,0x([0-9a-f]+)ULL\}',(quarter/'QuarterCircuit.h').read_text())]
    raw=[int(r,16) for r in certificate['outer']['generator_rows_hex']]
    assert len(rows)==32 and rows==[int(r,16) for r in outer['generator_rows_hex']]
    assert g.tv.smaller_outer.r.rank(rows)==g.tv.smaller_outer.r.rank(rows+raw)==32
    performance_path=HERE/'BALANCED_PERFORMANCE.json';performance=read(performance_path);sources(performance,ROOT)
    assert performance['status']=='SERIAL_BALANCED_MAP_PERFORMANCE'
    files=[Path(__file__),certificate_path,replay_path,performance_path,HERE/'BALANCED_IMPLEMENTATION.json',
           HERE/'NO_CONSTANT_MAP.json',quarter/'MANIFEST.json',HERE/'CMakeLists.txt',
           HERE/'test_outward.py',HERE/'BALANCED_RESULT.md',HERE/'requirements.txt']
    result=dict(status='CERTIFICATE_IMPLEMENTATION_AND_PERFORMANCE_BOUND',
                message_bits=1<<20,output_bits=1<<22,
                margin_bits=[r['margin_bits_diagnostic'] for r in certificate['results']],
                exact_dyadic_unions_checked=True,exact_outer_span_checked=True,
                exact_generated_map_checked=True,all_weight_bands_checked=True,
                source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files})
    (HERE/'BALANCED_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
