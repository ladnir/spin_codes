"""Retain the existing syndrome-density refinement at variable step sizes."""
import argparse
from functools import lru_cache
from fractions import Fraction as F
import json
import math
from pathlib import Path
from flint import ctx
import model
import activation_density


class Engine(model.Engine):
    @lru_cache(maxsize=24)
    def epoch(self,log_lam):
        rows=super().epoch(log_lam)
        z=(-model.number(log_lam).exp()).exp()
        output=[]
        for j,old in enumerate(rows):
            total,cap=math.comb(self.t,j),int(self.caps[j]['cap'])
            if self.m*cap<=4*total and total>self.kernel[j]:
                row=activation_density.zero_row(self.spectrum,self.kernel[j],total,cap,z**j)
                output.append(row+old[self.n:])
            else:
                output.append(old)
        return output


def run(a):
    output=a.output.resolve()
    if output.exists(): raise FileExistsError(output)
    ctx.prec=256
    rows=[]
    for s in a.s:
        e=Engine(64,s)
        if a.mode=='q1':
            row=model.study.core.evaluate(e.record,256,16,sharp=True)
            print('diagnostic Q1',s,row['q1_margin_bits'],flush=True)
        else:
            values=[]
            for q in a.q:
                value=model.study.core.grid.ladder.screen.sparse(e,q)
                values.append(value)
                print('activated sparse',s,q,value['margin_bits'],flush=True)
            row=dict(s=s,instance=e.identity(),sparse=values)
        rows.append(row)
    result=dict(status='K16_ACTIVATION_EXPLORATION',mode=a.mode,rows=rows,
        full_distance_proved=False,source_sha256=model.base.sources())
    output.parent.mkdir(parents=True,exist_ok=True)
    model.base.base.write_new(output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode',choices=['q1','sparse'],required=True)
    p.add_argument('--s',type=int,nargs='+',default=[15])
    p.add_argument('--q',type=int,nargs='+',default=[32,63])
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args())
