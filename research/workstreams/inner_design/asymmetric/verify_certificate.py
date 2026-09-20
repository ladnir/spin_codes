"""Check exact union, coverage, independent maps, and implementation binding."""
import hashlib
import json
from pathlib import Path
import re

import search_feedback as search
import certify
import generate_balanced as balanced
from verify_balanced_artifact import check_dyadic_union,partition

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def bind(record):
    for name,digest in record['source_sha256'].items():assert sha(ROOT/name)==digest,name


def main():
    path=HERE/'ASYMMETRIC_MARGIN_CERTIFICATE.json'
    replay_path=HERE/'ASYMMETRIC_MARGIN_CERTIFICATE_REPLAY.json'
    certificate=read(path);replay=read(replay_path)
    assert certificate['status']=='OUTWARD_ASYMMETRIC_ALL_OCCUPANCY_CERTIFICATE'
    assert certificate['precision_bits']==256 and certificate['inner']==certify.INNER
    assert (certificate['message_bits'],certificate['output_bits'])==(1<<20,1<<22)
    assert replay['status']=='HIGHER_PRECISION_REPLAY_PASSED' and replay['precision_bits']==512
    assert replay['certificate_sha256']==sha(path)
    bind(certificate)
    a_path=HERE.parent/'NO_CONSTANT_MAP.json'
    a_record=read(a_path);engine=certify.Engine(a_record)
    witness_path=search.g.tv.fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json'
    witnesses=read(witness_path)
    construction=search.g.tv.smaller_outer.construction()
    assert certificate['outer']==construction
    weights=[w for w,n in enumerate(search.g.tv.smaller_outer.spectrum()) if w and n]
    for key in ('q2_bands','adaptive_bands'):partition(witnesses[key],weights)
    assert len(certificate['results'])==2
    for row,delta,target,maximum in zip(certificate['results'],('33/200','19/100'),(40,30),(64,128)):
        assert row['distance_target']==delta and row['target_bits']==target
        assert row['maximum_sparse']==maximum and row['exact_occupation_coverage_checked']
        d=search.g.F(delta);assert row['bad_weight']==(1<<22)*d.numerator//d.denominator
        cover=read(HERE/f'DENSE_greedy3_2_{delta.replace("/","_")}.json')
        assert cover['candidate']=='greedy3_2' and cover['occupation_min']==maximum+1 and cover['occupation_max']==32768
        partition(cover['bands'],weights)
        search.g.check_coverage(cover['selected_boxes'],32768,maximum+1)
        assert len(row['dense_terms'])==len(cover['selected_boxes'])
        assert len(row['adaptive_terms'])==maximum-2
        expected=list(search.g.tv.fixed.composition.compositions(len(witnesses['q2_bands']),2))
        assert len(row['q2_terms'])==len(expected)
        check_dyadic_union([row['q1_upper']]+row['q2_terms']+row['adaptive_terms']+row['dense_terms'],row['union_upper'],target)
    implementation=read(HERE/'IMPLEMENTATION.json');bind(implementation)
    selected=next(c for c in implementation['candidates'] if c['name']=='asymmetric_greedy3_2_sparse')
    directory=HERE.parent/'generated'/selected['name']
    for name,digest in selected['source_sha256'].items():assert sha(directory/name)==digest,name
    a_header,_=balanced.header(engine.a_columns)
    expected=a_header+'\nnamespace bare_spin {\nstruct AsymmetricMap : BalancedMap {\n'
    expected+='static constexpr auto feedbackColumns=BalancedMap::columns;\n'
    expected+='static constexpr std::array<std::uint32_t,T> columns{'+','.join(map(hex,engine.columns))+'};\n'
    expected+='static constexpr auto groupedColumns=columns;\n'
    expected+='static constexpr std::array<unsigned,S> groupOrder{'+','.join(map(str,range(19)))+'};\n};\n}\n'
    assert (directory/'AsymmetricMap.h').read_text()==expected
    quarter=ROOT/'workstreams/rate_quarter_bch/implementation/generated'
    outer_manifest=read(quarter/'MANIFEST.json')
    assert outer_manifest['outer']==construction
    for name,digest in outer_manifest['generated_sha256'].items():assert sha(quarter/name)==digest,name
    rows=[int(a,16)|(int(b,16)<<64) for a,b in re.findall(r'\{0x([0-9a-f]+)ULL,0x([0-9a-f]+)ULL\}',(quarter/'QuarterCircuit.h').read_text())]
    raw=[int(r,16) for r in construction['generator_rows_hex']]
    assert len(rows)==32 and rows==[int(r,16) for r in outer_manifest['generator_rows_hex']]
    assert search.rank(rows)==search.rank(rows+raw)==32
    log=HERE/'measurements/asymmetric-correctness.log'
    assert '100% tests passed, 0 tests failed out of 4' in log.read_text()
    sources=[Path(__file__),path,replay_path,a_path,HERE/'IMPLEMENTATION.json',HERE/'AsymmetricInner.h',directory/'AsymmetricMap.h',directory/'CandidateSpin.cpp',log,
             quarter/'MANIFEST.json',quarter/'QuarterCircuit.h',quarter/'QuarterCircuit.cpp']
    result=dict(status='ASYMMETRIC_CERTIFICATE_AND_IMPLEMENTATION_BOUND',
                message_bits=1<<20,output_bits=1<<22,selected_candidate=selected['name'],
                margins_bits=[r['margin_bits_diagnostic'] for r in certificate['results']],
                exact_union_and_coverage_checked=True,exact_maps_reconstructed=True,
                exact_outer_row_space_checked=True,
                higher_precision_replay_checked=True,
                scope='Quarter-rate BCH [128,32,32], t=128,s=19,K=2^20 only; no asymptotic claim',
                source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sources})
    (HERE/'CERTIFICATE_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
