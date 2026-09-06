"""Verify completed batches, append them to the index, and report coverage."""
import argparse
import csv
import json

import run_complete_q1_grid as grid
from read_grid_receipts import snapshot
import spectrum_events


def common_event(rows, registry):
    events = {r.get('setup_event_id') for r in rows if r.get('setup_event_id')}
    for row in rows:
        event_id = row.get('setup_event_id')
        if not event_id: continue
        event = registry.get(event_id)
        if (event is None or event['setup_failure_bits'] != int(row['setup_failure_bits'])
                or (event['block_bits'], event['dimension']) !=
                   (int(row['block_bits']), int(row['dimension']))):
            raise ValueError('coverage row disagrees with authenticated spectrum event')
    if not events: return ''
    candidates = [a for a in events if all(spectrum_events.implies(registry[a], registry[b]) for b in events)]
    if not candidates:
        raise ValueError('coverage requires a checked common spectrum event')
    return min(candidates, key=lambda a: (-registry[a]['setup_failure_bits'], a))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--register',action='store_true')
    parser.add_argument('--require-complete',action='store_true')
    args=parser.parse_args()
    planned,observations,specs=snapshot()
    expected=sum(r['native'] for r in planned)
    if args.require_complete and len(observations)!=expected:
        raise ValueError(f'Q1 incomplete: {len(observations)}/{expected}')
    catalog_path=grid.HERE/'catalog.json'
    catalog=json.loads(catalog_path.read_text())
    event_registry=spectrum_events.load_registered_events()
    if args.register:
        known={s['path'] for s in catalog['activation_sources']}
        for spec in specs:
            if spec['path'] not in known:
                catalog['activation_sources'].append(spec)
        catalog_path.write_text(json.dumps(catalog,indent=2)+'\n')
    higher={};ranges={};composition_coverage={}
    for spec in catalog['activation_sources']:
        if spec.get('occupation')==1:
            continue
        receipt=json.loads((grid.HERE/spec['manifest']).read_text())
        path=grid.HERE/spec['path']
        if grid.pilot.sha(path)!=receipt['csv_sha256']:
            raise ValueError('changed occupation CSV')
        for dependency,digest in receipt['source_sha256'].items():
            if grid.pilot.sha(grid.pilot.ROOT/dependency)!=digest:
                raise ValueError(f'changed occupation dependency: {dependency}')
        with path.open(newline='') as handle:
            for row in csv.DictReader(handle):
                q=int(row.get('occupation') or spec.get('occupation') or 1)
                q_min=int(row.get('occupation_min') or q)
                q_max=int(row.get('occupation_max') or q)
                if q_min!=q_max:
                    if not 1<=q_min<=q_max<=int(row['outer_rows']): raise ValueError('invalid range coverage')
                    ranges.setdefault(grid.key(row),[]).append((q_min,q_max,row))
                    continue
                if q<2: continue
                if 'composition' in spec['path']:
                    composition_coverage.setdefault(grid.key(row),set()).add(q)
                group=higher.setdefault(grid.key(row),{})
                group.setdefault(q,[]).append(row)
    rows=[]
    for candidate in planned:
        observation=observations.get(grid.key(candidate))
        alternatives=higher.get(grid.key(candidate),{})
        group={q:max(bounds,key=lambda r:float(r['margin_bits'])) for q,bounds in alternatives.items()}
        interval_rows=ranges.get(grid.key(candidate),[])
        pair=group.get(2)
        if pair and (not observation or pair['map_tag']!=observation['map_tag']):
            raise ValueError('Q2 map has no matching Q1 map')
        length=(1<<candidate['message_exponent'])//candidate['dimension']
        through=int(observation is not None)
        while True:
            extended=max([through]+([through+1] if through+1 in group else [])+
                         [hi for lo,hi,r in interval_rows if lo<=through+1<=hi])
            if extended==through:break
            through=extended
        evidence_rows=[r for bounds in alternatives.values() for r in bounds]+[r for lo,hi,r in interval_rows]
        if any(not observation or r['map_tag']!=observation['map_tag'] for r in evidence_rows):
            raise ValueError('higher occupation map has no matching Q1 map')
        event=common_event(evidence_rows,event_registry)
        row=dict(candidate,outer_rows=length,
                 q1_status='complete' if observation else ('pending' if candidate['native'] else 'inapplicable'),
                 q1_margin_bits=observation['margin_bits'] if observation else '',
                 q1_witness_edge=observation['dominant_witness_at_grid_edge'] if observation else '',
                 q2_status='complete' if pair else ('pending' if candidate['native'] else 'inapplicable'),
                 q2_margin_bits=pair['margin_bits'] if pair else '',
                 evaluated_contiguous_through=through,
                 setup_event_id=event,
                 setup_failure_bits=event_registry[event]['setup_failure_bits'] if event else '',
                 required_through=length if candidate['native'] else '',
                 full_occupation_status=('evaluated' if through==length else 'pending') if candidate['native'] else 'inapplicable',
                 evidence='nearest binary64 diagnostic/reference; no current certificate')
        rows.append(row)
    with (grid.HERE/'complete_grid_coverage.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    summary=dict(candidate_count=len(planned),applicable_count=expected,
                 q1_complete=len(observations),q1_pending=expected-len(observations),
                 q2_complete=sum(2 in group for group in higher.values()),
                 through_q16_complete=sum(r['evaluated_contiguous_through']>=16 for r in rows),
                 through_q32_complete=sum(r['evaluated_contiguous_through']>=32 for r in rows),
                 through_q64_complete=sum(r['evaluated_contiguous_through']>=64 for r in rows),
                 composition_q2_q4_complete=sum({2,3,4}.issubset(qs) for qs in composition_coverage.values()),
                 full_occupation_complete=sum(r['full_occupation_status']=='evaluated' for r in rows),
                 q1_dominant_grid_edges=sum(int(r['dominant_witness_at_grid_edge']) for r in observations.values()),
                 new_completed_batches=len(specs),
                 registered_grid_batches=sum(s['path'].startswith('activation_grid_v1/') for s in catalog['activation_sources']),
                 note='Coverage is not a useful positive bound or an outward certificate; see complete_grid_unions.json.')
    (grid.HERE/'complete_grid_coverage.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
