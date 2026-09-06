"""Compare tested states and epoch sizes at explicit occupation coverage.

These are circuit-count proxies, not measured runtime predictions. Different
outer constituents remain separate because their encoding costs differ.
"""
import csv
import hashlib
import json
from pathlib import Path

import read_grid_receipts
import run_complete_q1_grid as grid


def frontier(candidates):
    """Minimize independent XOR reductions and state-update frequency."""
    return [r for r in candidates if not any(
        a['xor_per_bit']<=r['xor_per_bit'] and a['updates_per_bit']<=r['updates_per_bit']
        and (a['xor_per_bit']<r['xor_per_bit'] or a['updates_per_bit']<r['updates_per_bit'])
        for a in candidates)]


def map_cost(path):
    payload=json.loads(path.read_text())
    t,s=payload['step_bits'],payload['state_bits']
    generators=[int(g,16) for g in payload['generator_words_hex']]
    a=sum(max(0,sum((g>>j)&1 for g in generators)-1) for j in range(t))
    b=sum(max(0,g.bit_count()-1) for g in generators)
    return dict(xor_per_bit=(a+b+t)/t,updates_per_bit=1/t,state_bits=s,
                a_xors_per_epoch=a,b_xors_per_epoch=b,output_xors_per_epoch=t)


def main():
    planned,observations,_=read_grid_receipts.snapshot()
    source=grid.HERE/'complete_grid_unions.csv'
    with source.open(newline='') as handle:unions={grid.key(r):r for r in csv.DictReader(handle)}
    native=[r for r in planned if r['native']]
    if set(unions)!={grid.key(r) for r in native}:raise ValueError('union snapshot does not cover the native grid')
    map_costs={};groups={}
    for candidate in native:
        key=grid.key(candidate);q1=observations[key];union=unions[key]
        if q1['map_tag']!=union['map_tag']:raise ValueError('map mismatch in union report')
        path=grid.pilot.ROOT/q1['map_source']
        if q1['map_tag']!='nested-'+grid.pilot.sha(path)[:16]:raise ValueError('map source changed')
        if q1['map_tag'] not in map_costs:map_costs[q1['map_tag']]=map_cost(path)
        margins={'q1':float(q1['margin_bits'])}
        for q in (4,16,32,64):
            value=union.get(f'q1_q{q}_union_margin_bits','')
            margins[f'q1_q{q}']=float(value) if value!='' else None
        value=union['full_union_margin_bits'];margins['full']=float(value) if value!='' else None
        for coverage,margin in margins.items():
            row=dict(candidate,**map_costs[q1['map_tag']],map_tag=q1['map_tag'],margin_bits=margin,
                     occupation_coverage=coverage,full_status=union['full_status'])
            groups.setdefault((candidate['series'],candidate['message_exponent'],coverage),[]).append(row)
    output=[]
    for (series,exponent,coverage),rows in groups.items():
        for threshold in (0,20,40,60):
            passing=[r for r in rows if r['margin_bits'] is not None and r['margin_bits']>=threshold]
            pareto={(r['step_bits'],r['state_bits']) for r in frontier(passing)}
            for t in grid.STEPS:
                tested=[r for r in rows if r['step_bits']==t]
                available=[r for r in tested if r['margin_bits'] is not None]
                accepted=[r for r in passing if r['step_bits']==t]
                best=min(accepted,key=lambda r:r['state_bits']) if accepted else None
                output.append(dict(series=series,message_exponent=exponent,occupation_coverage=coverage,
                    target_margin_bits=threshold,step_bits=t,evaluated_states=len(available),
                    native_states=len(tested),smallest_tested_passing_state=best['state_bits'] if best else '',
                    attained_margin_bits=best['margin_bits'] if best else '',
                    xor_per_output_bit=best['xor_per_bit'] if best else '',
                    state_updates_per_output_bit=1/t,
                    on_proxy_pareto_frontier=int((t,best['state_bits']) in pareto) if best else '',
                    status=('passing_diagnostic' if coverage=='full' else 'passing_partial_screen') if best else
                           ('coverage_missing' if not available else 'no_tested_state_passes'),
                    evidence='nearest binary64; no outward certificate; operation proxies exclude outer and field arithmetic'))
    csv_path=grid.HERE/'tested_parameter_frontiers.csv'
    with csv_path.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(output[0]));writer.writeheader();writer.writerows(output)
    summary=dict(scope='Every native finite-grid tuple, with outer constituents compared separately',
        native_parameter_count=len(native),frontier_rows=len(output),
        union_snapshot_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        csv_sha256=hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        maps=map_costs,
        full_passing_40_bit_frontiers=[r for r in output if r['occupation_coverage']=='full'
            and r['target_margin_bits']==40 and r['status']=='passing_diagnostic'],
        limitations=['A missing full result does not imply code failure.',
            'The smallest tested state is not a proof of global optimality.',
            'Circuit counts exclude GF multiplication, update addition, outer encoding, permutations, SIMD, and circuit sharing.',
            'No certificate or runtime extrapolation is inferred from partial screens.'])
    (grid.HERE/'tested_parameter_frontiers.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('maps','full_passing_40_bit_frontiers')},indent=2))


if __name__=='__main__':main()
