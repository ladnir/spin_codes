"""Independently reconstruct v30 bookkeeping and its positive event union."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import shutil
from flint import arb,ctx
import hd_field128_gamma_selected_index_v30 as selected

HERE=Path(__file__).resolve().parent
INDEX='hd_field128_gamma30_d40_partial_cover_v30_20261010.json'
SELECTION='hd_field128_gamma30_d40_v30_selection_20261010.json'
PINS={
 INDEX:'875c572f4ae84b5bc139a4c21be9b0b4e18aed51ab1f288c45c125a260aab88f',
 SELECTION:'1917e2990b6084348d845371ad00242552f4ac19c8ea2ff3abfa37822338f8a2',
 'hd_field128_gamma_selected_index_v30.py':'6e914264b8b42c162c93f17d174a979b1b450b99cfa107b4161907e981edc78b',
 'hd_field128_gamma_selected_index_v30_tests.py':'7d00ada745c8e541170fd2e1266fe79f775c4969dd67141d81386cf16dd12113',
 'hd_field128_gamma30_d40_v30_record_20261010.py':'440307691b144d38e085fea268c188202aaf10e1991370b5d540635f9ba16eb4'}
LIMIT,RESERVE=16<<10,20<<30

def require(test,message):
    if not test:raise ValueError(message)

def main():
    require(shutil.disk_usage(HERE).free-LIMIT>=RESERVE,'Preserve20GiB')
    original={}
    def read(name,digest):
        path=HERE/name;raw=path.read_bytes()
        require(sha256(raw).hexdigest()==digest,'Changed input '+name)
        original[path]=raw
        return json.loads(raw) if name.endswith('.json') else None
    for name,digest in PINS.items():read(name,digest)
    saved=json.loads(original[HERE/INDEX]);selection=json.loads(original[HERE/SELECTION])
    parent=read(selection['parent']['path'],selection['parent']['sha256'])
    read(selection['parent_independent_check']['path'],selection['parent_independent_check']['sha256'])
    for spec in selection['event_specs']:
        for pin in spec.values():read(pin['path'],pin['sha256'])
    rebuilt=selected.extend_v29(HERE/selection['parent']['path'],
        HERE/selection['parent_independent_check']['path'],selection['event_specs'])
    require(saved==rebuilt,'Saved object differs from rebuilt validator output')
    require(selected.canonical(saved['entries'][:260])==selected.canonical(parent['entries']),'Parent changed')
    covered=set()
    for entry in saved['entries']:
        lo,hi=entry['weight_range'];weights=set(range(lo,hi+1))
        require(1<=lo<=hi<=8192 and not covered&weights,'Bad or overlapping event')
        covered|=weights
    require(len(saved['entries'])==265 and len(covered)==8131,'Coverage count')
    pending=set(range(1,8193))-covered
    require(pending==set(range(521,548))|set(range(549,576))|set(range(917,920))|set(range(1146,1150)),'Pending domain')
    def union(entries):
        with ctx.workprec(384):
            total=arb(0);ln2=arb(2).log()
            for entry in entries:
                if entry.get('zero',False):
                    require(entry['log2_upper_rational'] is None,'Zero has logarithm')
                else:
                    f=F(entry['log2_upper_rational'])
                    total+=(ln2*arb(f.numerator)/f.denominator).exp()
            value=total.log()/ln2
            bound=F(int((value.upper()*10**12).upper().ceil().fmpz()),10**12)
            return str(bound),str(value)
    bound,enclosure=union(saved['entries'])
    added,added_enclosure=union(saved['entries'][260:])
    require(bound==saved['covered_events_union_log2_upper_rational'] and F(bound)<=-40,'Union or target')
    require(added==saved['added_events_union_log2_upper_rational'],'Added union')
    require(saved['whole_code_certificate'] is False,'Partial is not whole-code')
    selected.frozen.producer.check_source_pins(saved['aggregation_source_sha256'])
    report=dict(schema='gamma-saved-partial-index-independent-check-v30',audit_result='PASS',
        index=dict(path=INDEX,sha256=PINS[INDEX],bytes=len(original[HERE/INDEX])),
        selection=dict(path=SELECTION,sha256=PINS[SELECTION],bytes=len(original[HERE/SELECTION])),
        construction=saved['construction'],state_dimension=30,target_distance='2/5',
        full_saved_object_matches_validator_rebuild=True,parent_entry_count=260,
        parent_canonical_entry_bytes_preserved=True,parent_entries_and_order_preserved=True,
        receipt_count=265,covered_weight_count=8131,pending_weight_count=61,
        pending_weight_ranges=saved['pending_weight_ranges'],nonoverlapping_weight_ranges=True,
        independent_positive_union_precision_bits=384,
        independent_positive_union_method='Sum exp(log(2)*exact rational ceiling), then directed log2 and upward1e-12 ceiling',
        covered_events_union_log2_upper_rational=bound,covered_events_union_log2_upper_enclosure=enclosure,
        added_events_union_log2_upper_rational=added,added_events_union_log2_upper_enclosure=added_enclosure,
        covered_events_failure_bits=40,meets_covered_events_target=True,whole_code_certificate=False,
        compact_inputs_verified=[dict(path=p.relative_to(HERE).as_posix(),sha256=sha256(raw).hexdigest(),bytes=len(raw))
            for p,raw in sorted(original.items())],
        aggregation_source_pins_rechecked=len(saved['aggregation_source_sha256']),
        source_and_compact_input_bytes_rechecked_at_exit=True,
        new_mathematical_audits_inherited_from_remote=True,new_archive_verifications_inherited_from_remote=True,
        local_tree_replay=False,local_archive_revalidation=False,remote_network_access=False,
        historical_origin_verified=False,coefficient_replay=False,normalized_regional_replay=False,transfer_recomputed=False,
        auditor=dict(path=Path(__file__).name,sha256=sha256(Path(__file__).read_bytes()).hexdigest()),
        scope='Saved object, exact support coverage and positive384 event union; no coefficient or transfer replay.',
        maximum_output_bytes=LIMIT)
    raw=(json.dumps(report,indent=2,allow_nan=False)+'\n').encode()
    require(len(raw)<=LIMIT and shutil.disk_usage(HERE).free-len(raw)>=RESERVE,'Output budget')
    require(all(p.read_bytes()==data for p,data in original.items()),'Input changed at exit')
    target=HERE/'hd_field128_gamma30_d40_partial_cover_v30_20261010_independent_check.json'
    with target.open('xb') as f:f.write(raw)
    require(target.read_bytes()==raw,'Saved audit changed')
    print(json.dumps(dict(path=target.name,sha256=sha256(raw).hexdigest(),bytes=len(raw),
        audit_result='PASS',covered_weight_count=8131,pending_weight_count=61,
        covered_events_union_log2_upper_rational=bound,whole_code_certificate=False)),flush=True)
if __name__=='__main__':main()
