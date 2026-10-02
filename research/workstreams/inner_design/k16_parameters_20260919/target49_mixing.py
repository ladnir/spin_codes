"""Separate diagnostic: add a second transvection, retaining fixed A/B maps."""
import argparse
import json
from pathlib import Path
from types import FunctionType
import math
import numpy as np
import model


def evaluate(record,rounds):
    core=model.study.core
    def epoch_logs(record,tilts,refresh=False,sharp=False):
        assert sharp and not refresh
        entries=[]
        for tilt in tilts:
            zero,one,n=core.mixing.transfers({int(w):n for w,n in record['spectrum'].items()},[v[0] for v in record['cancellation']],
                [v[1] for v in record['cancellation']],math.exp(-math.exp(float(tilt))),2.**-rounds)
            entries.append(np.array([zero,one]).reshape(2,n,n))
        with np.errstate(divide='ignore'): logs=np.log(np.array(entries))
        return logs[:,0],logs[:,1]
    evaluator=FunctionType(core.evaluate.__code__,{**core.evaluate.__globals__,'epoch_logs':epoch_logs},
                           argdefs=core.evaluate.__defaults__)
    return evaluator(record,256,16,sharp=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maps',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    source=json.loads(args.maps.read_text())
    cases=[('selected12819',model.study.record_for(256,128))]
    cases += [(f'sub64s{r["s"]}',r['inner']) for r in source['rows'] if r['found']]
    rows=[]
    for name,record in cases:
        for rounds in (1,2):
            diagnostic=evaluate(record,rounds)
            rows.append(dict(name=name,rounds=rounds,diagnostic=diagnostic))
            print(name,'rounds',rounds,diagnostic['q1_margin_bits'],flush=True)
    model.base.base.write_new(args.output,dict(status='K16_TARGET49_MIXING_BINARY64_DIAGNOSTIC',
        full_distance_proved=False,rows=rows,source_sha256={**model.base.sources(),
            Path(__file__).resolve().relative_to(model.base.ROOT).as_posix():model.base.base.sha(Path(__file__)),
            args.maps.resolve().relative_to(model.base.ROOT).as_posix():model.base.base.sha(args.maps)}))
