"""Replay the retained smaller-outer dense witnesses and complete union."""
import argparse
import json
import math
from pathlib import Path
from fractions import Fraction

import numpy as np
import evaluate_smaller_margins as evaluate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=evaluate.fixed.HERE/'SMALLER_MARGIN_CLOSED.json')
    path = parser.parse_args().input
    payload = json.loads(path.read_text())
    require = evaluate.fixed.outer.require
    require(payload['construction'] == evaluate.smaller_outer.construction(), 'construction mismatch')
    require(payload['inner'] == 'existing optimized t128_s19' and not payload['inner_reoptimized'], 'inner mismatch')
    for name,digest in payload['source_sha256'].items():
        require(evaluate.fixed.sha(evaluate.fixed.ROOT/name) == digest, 'changed source: '+name)
    values = evaluate.smaller_outer.spectrum()
    counts = {w:n for w,n in enumerate(values) if w and n}
    require(payload['spectrum'] == {str(w):str(n) for w,n in counts.items()}, 'spectrum mismatch')
    a,kernel,_ = evaluate.fixed.load_inner()
    for row in payload['results']:
        lo,hi = row['covered_occupations']
        require(lo==1 and len(row['occupation_margins_bits'])==hi, 'incomplete sparse range')
        require('dense' in row and row['dense']['occupation_min']==hi+1, 'missing dense range')
        require(row['output_bits']==128*row['outer_rows']==4*row['message_bits'], 'wrong rate')
        delta = Fraction(row['distance_target'])
        require(row['bad_weight']==row['output_bits']*delta.numerator//delta.denominator,'wrong distance cutoff')
        sparse = float(np.logaddexp.reduce(-np.array(row['occupation_margins_bits'])*math.log(2)))
        require(abs(-sparse/math.log(2)-row['sparse_union_margin_bits'])<1e-9,'sparse union mismatch')
        dense = evaluate.replay_dense(row['dense'],counts,a,kernel,row)
        margin = -float(np.logaddexp(sparse,dense))/math.log(2)
        require(abs(margin-row['combined_margin_bits'])<1e-9,'combined union mismatch')
        print(dict(distance=row['distance_target'],margin_bits=margin,
            exact_all_occupation_coverage=True,dense_replay=True,
            status='BINARY64_DIAGNOSTIC_NOT_OUTWARD_CERTIFICATE'),flush=True)


if __name__ == '__main__':
    main()
