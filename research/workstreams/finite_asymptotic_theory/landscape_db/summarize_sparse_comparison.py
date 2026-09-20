"""Join Q1 and Q2 for identical maps and exact outer spectra."""
import csv
import json
import math
import sqlite3
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    with sqlite3.connect(HERE/'spin_landscape.sqlite3') as db:
        db.row_factory=sqlite3.Row
        candidates=db.execute("SELECT * FROM landscape WHERE comparison_eligible=1 AND occupation_min=2 AND occupation_max=2 AND outer_model_kind='fixed_exact_spectrum' ORDER BY margin_bits DESC").fetchall()
        selected={}
        for row in candidates:
            selected.setdefault((row['inner_id'],row['outer_id'],row['message_bits'],row['bad_weight']),row)
        q2rows=list(selected.values())
        result=[]
        for row in q2rows:
            matches=db.execute("SELECT margin_bits FROM landscape WHERE comparison_eligible=1 AND occupation_min=1 AND occupation_max=1 AND inner_id=? AND outer_id=? AND message_bits=? AND bad_weight=?",
                               (row['inner_id'],row['outer_id'],row['message_bits'],row['bad_weight'])).fetchall()
            if len(matches)!=1:
                raise ValueError('Q2 needs a unique Q1 match with identical map and cutoff')
            q1margin=float(matches[0][0]); q2margin=float(row['margin_bits'])
            combined=min(q1margin,q2margin)-math.log2(1+2**(-abs(q1margin-q2margin)))
            result.append(dict(outer=row['outer_label'],message_exponent=row['message_exponent'],
                               t=row['step_bits'],s=row['state_bits'],q1_margin_bits=q1margin,q2_margin_bits=q2margin,
                               q1_q2_union_margin_bits=combined,dominant_q2_pair=[row['dominant_weight'],row['dominant_second_weight']],
                               q2_witness=row['dominant_log_surprisal'],status='Q1+Q2 diagnostic only; higher occupations tracked separately'))
    if not result: raise ValueError('no Q2 observations')
    with (HERE/'sparse_comparison.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(result[0])); writer.writeheader(); writer.writerows(result)
    (HERE/'sparse_comparison.json').write_text(json.dumps(dict(rows=result,full_distance=False),indent=2)+'\n')
    for row in result:
        if row['outer'].startswith('RM(4,9)'):
            print(row)


if __name__=='__main__': main()
