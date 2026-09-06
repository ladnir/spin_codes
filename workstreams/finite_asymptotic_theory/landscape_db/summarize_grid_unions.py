"""Combine eligible occupation bounds without hiding gaps or setup events."""
import csv
import json
import math
import sqlite3

import numpy as np

import read_grid_receipts as receipts
import run_complete_q1_grid as grid
from build_landscape_db import infer_family
import spectrum_events


def best_cover(length,bounds,event_implications=None,event_budgets=None):
    """Return the least interval-cover union found, with one setup charge.

    Overlap is conservative: an interval's whole bound is charged even when
    only its uncovered suffix is needed. No per-occupation division occurs.
    """
    events={r.get('setup_event_id') for r in bounds if r.get('setup_event_id')}
    candidates=[]
    for event in [None]+sorted(events):
        allowed=[r for r in bounds if not r.get('setup_event_id') or r['setup_event_id']==event
                 or (event_implications or {}).get((event,r['setup_event_id']),False)]
        relevant=[r for r in allowed if r['lo']<=length and r['hi']>=1]
        nodes=sorted({1,length+1}|{min(length,r['hi'])+1 for r in relevant})
        dp={(1,False):(-math.inf,[])}
        for start in nodes:
            for used in (False,True):
                prior=dp.get((start,used))
                if prior is None or start>length:continue
                for row in relevant:
                    if not row['lo']<=start<=row['hi']:continue
                    end=min(length,row['hi'])+1
                    new_used=used or bool(row.get('setup_event_id'))
                    value=float(np.logaddexp(prior[0],row['log_bound']))
                    key=(end,new_used)
                    if key not in dp or value<dp[key][0]:
                        dp[key]=(value,prior[1]+[row['result_id']])
        for used in (False,True):
            result=dp.get((length+1,used))
            if result is None:continue
            value,ids=result
            bits=None
            if used:
                budgets=({int(event_budgets[event])} if event_budgets and event in event_budgets else
                         {int(r['setup_failure_bits']) for r in allowed if r.get('setup_event_id')==event})
                if len(budgets)!=1:raise ValueError('one setup event has inconsistent failure budgets')
                bits=budgets.pop()
                value=float(np.logaddexp(value,-bits*math.log(2)))
            candidates.append(dict(log_bound=value,result_ids=ids,setup_event_id=event if used else None,
                                   setup_failure_bits=bits))
    return min(candidates,key=lambda r:r['log_bound']) if candidates else None


def main():
    planned,observations,_=receipts.snapshot()
    events=spectrum_events.load_registered_events()
    implications=spectrum_events.containment(events)
    budgets={key:value['setup_failure_bits'] for key,value in events.items()}
    grouped={}
    db=sqlite3.connect(grid.HERE/'spin_landscape.sqlite3');db.row_factory=sqlite3.Row
    try:
        for r in db.execute('SELECT * FROM landscape WHERE comparison_eligible=1'):
            if r['setup_event_id']:
                event=events.get(r['setup_event_id'])
                if (event is None or event['setup_failure_bits']!=r['setup_failure_bits']
                        or (event['block_bits'],event['dimension'])!=(r['block_bits'],r['dimension'])):
                    raise ValueError('conditional row disagrees with its authenticated setup event')
            key=(r['outer_family'],r['block_bits'],r['dimension'],r['step_bits'],r['state_bits'],r['message_exponent'],r['map_tag'])
            grouped.setdefault(key,[]).append(dict(lo=r['occupation_min'],hi=r['occupation_max'],
                log_bound=-r['margin_bits']*math.log(2),result_id=r['result_id'],
                setup_event_id=r['setup_event_id'],setup_failure_bits=r['setup_failure_bits']))
    finally:db.close()
    rows=[];details=[]
    for candidate in planned:
        if not candidate['native']:continue
        q1=observations[grid.key(candidate)]
        key=(infer_family(candidate['series']),candidate['block_bits'],candidate['dimension'],
             candidate['step_bits'],candidate['state_bits'],candidate['message_exponent'],q1['map_tag'])
        bounds=grouped.get(key,[])
        length=(1<<candidate['message_exponent'])//candidate['dimension']
        row={k:v for k,v in candidate.items() if k!='native'}
        row['map_tag']=q1['map_tag']
        for maximum in (2,4,16,32,64):
            union=best_cover(maximum,[r for r in bounds if r['hi']<=maximum],implications,budgets)
            row[f'q1_q{maximum}_union_margin_bits']=-union['log_bound']/math.log(2) if union else ''
        full=best_cover(length,bounds,implications,budgets)
        row['full_union_margin_bits']=-full['log_bound']/math.log(2) if full else ''
        row['full_status']=('positive_diagnostic_union' if full['log_bound']<0 else 'evaluated_bound_too_loose') if full else 'incomplete_occupation_coverage'
        row['full_setup_failure_bits']=full['setup_failure_bits'] if full else ''
        row['full_setup_event_id']=full['setup_event_id'] if full else ''
        row['evidence']='nearest binary64; no outward certificate'
        rows.append(row)
        if full:details.append(dict(parameter_key=list(grid.key(candidate)),**full))
    with (grid.HERE/'complete_grid_unions.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary=dict(parameter_count=len(rows),full_evaluated=sum(r['full_status']!='incomplete_occupation_coverage' for r in rows),
                 full_positive_diagnostics=sum(r['full_status']=='positive_diagnostic_union' for r in rows),
                 sparse_q1_q4_above_40=sum(r['q1_q4_union_margin_bits']!='' and r['q1_q4_union_margin_bits']>=40 for r in rows),
                 full_certificates=0,full_union_sources=details,
                 spectrum_events_checked=len(events),
                 checked_event_implications=[dict(stronger=a,weaker=b) for (a,b),valid in implications.items() if a!=b and valid],
                 note='Full coverage, a useful positive bound, and an outward certificate are distinct outcomes.')
    (grid.HERE/'complete_grid_unions.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='full_union_sources'},indent=2))


if __name__=='__main__':main()
