"""Authenticate and combine only the explicitly tested scaling occupancies."""
from __future__ import annotations
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import scaling_screen as screen

HERE=Path(__file__).resolve().parent
FILES=('lowq_v2.json','near_one_lowq_v1.json','alpha_k18_v1.json',
       'alpha_k20_v1.json','matched_k18_q64_v1.json')


def merge_results(receipts):
    merged={}
    for receipt in receipts:
        if (receipt['whole_code_certificate'] or receipt['all_occupancies_checked']
                or not receipt['source_pins_verified_at_finish']):
            raise ValueError('completed partial screen required')
        for power,item in receipt['results'].items():
            record=merged.setdefault(power,dict(geometry=item['geometry'],endpoints={},
                raw_margins={},q1_support=None))
            if record['geometry']!=item['geometry']:
                raise ValueError('geometry mismatch')
            for q,value in item['selected_endpoints'].items():
                endpoint=screen.cert.read_dyadic(value)
                record['endpoints'][q]=min(record['endpoints'].get(q,Fraction(1)),endpoint)
                record['raw_margins'][q]=max(record['raw_margins'].get(q,float('-inf')),
                    item['selected_raw_margins_display'][q])
            if 'q1_selected_support_uppers' in item:
                support=tuple(map(screen.cert.read_dyadic,item['q1_selected_support_uppers']))
                record['q1_support']=support if record['q1_support'] is None else screen.q1.min_supports(record['q1_support'],support)
    for power,record in merged.items():
        if record['q1_support'] is not None:
            value=screen.q1.combine_supports(record['q1_support'],groups=record['geometry']['groups'])
            record['endpoints']['1']=min(Fraction(1),value)
            record['raw_margins']['1']=screen.cert.bits(value)
        record.pop('q1_support')
        record['tested_q']=sorted(map(int,record['endpoints']))
        record['untested_q_count']=record['geometry']['groups']-len(record['tested_q'])
        record['endpoints']={q:screen.cert.dyadic_record(value) for q,value in record['endpoints'].items()}
        record['whole_code_certificate']=False
    return merged


def run():
    pins={}
    receipts=[]
    for name in FILES:
        path=HERE/name
        saved=json.loads(path.read_text(encoding='utf-8'))
        if not screen.cert.checked_pins(saved['source_sha256']):
            raise ArithmeticError(f'authentication failed: {name}')
        screen.cert.merge_pins(pins,saved['source_sha256'])
        pins[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
        receipts.append(saved)
    for path in (Path(__file__),HERE/'test_summarize_screen.py'):
        pins[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
    output=HERE/'summary_v1.json'
    if output.exists():raise ValueError('fresh output required')
    summary=dict(schema='packet8-wider24-scaling-SCREEN-summary-1',
        whole_code_certificate=False,all_occupancies_checked=False,
        results=merge_results(receipts),source_sha256=pins,
        source_pins_verified_at_finish=screen.cert.checked_pins(pins),
        superseded_development_receipts=['lowq_v1.json'])
    if not summary['source_pins_verified_at_finish']:raise ArithmeticError('final pin check')
    output.write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    for power,item in summary['results'].items():
        print(f'K{power}: {item["raw_margins"]}; untested={item["untested_q_count"]}',flush=True)
    print(f'{len(pins)} pins authenticated',flush=True)
    return summary


if __name__=='__main__':run()
