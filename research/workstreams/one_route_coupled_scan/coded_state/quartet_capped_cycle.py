"""FLOAT witness inside the relaxed capped operator, NOT a codeword.

Serial Peach,90s,2GiB conservative arrays,no generated files.
"""
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_mixed_counts as mixed
import quartet_fresh_count_screen as fresh
from quartet_history_screen import logchoose
from quartet_joint_moment_caps import geometry,moment_logs


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        intervals=[(a,min(512,a+63)) for a in range(0,513,64)]
        base,scale,cells,edge_ids,labels,target_map=mixed.build(intervals,only_interval=(384,447))
        selected=[i for i,c in enumerate(cells) if c[:2]==(384,447)]
        reverse={i:j for j,i in enumerate(selected)}
        totals=moment_logs([cells[i] for i in selected],0.,0.,geometry())
        transitions={}
        for l in (407,408):
            ids=np.flatnonzero(labels==l);M=np.full((len(selected),len(selected)),-np.inf)
            for e in ids:
                j=reverse[target_map[e]]
                with np.errstate(divide='ignore'):M[:,j]=scale+np.log2(base[selected,e])
            peaks=M.max(axis=1);sums=peaks+np.log2(np.exp2(M-peaks[:,None]).sum(axis=1))
            cap=np.minimum(totals[:,l],sums)
            # A single edge of this size is feasible in the relaxed cap LP.
            transitions[l]=np.minimum(M,cap[:,None])
        initial=reverse[target_map[edge_ids[408,408//mixed.BIN,0]]]
        prefix=np.full(len(selected),-np.inf);prefix[initial]=0.
        for _ in range(3):prefix=np.max(prefix[:,None]+transitions[408],axis=0)
        value=logchoose(512,408)+408*math.log2(fresh.H[1])+prefix
        value+=4048*np.diag(transitions[408])+44*np.diag(transitions[407])
        value+=mixed.terminal_logs([cells[i] for i in selected])
        value+=44*(logchoose(6527,435)-logchoose(8192,2100))
        value+=212*(logchoose(6528,436)-logchoose(8192,2100))
        value+=logchoose(8192,2100)-fresh.LOGB*165556
        best=int(np.argmax(value))
        print(dict(float_relaxed_fixed_path_log2=float(value[best]),cell=cells[selected[best]],
                   initial_count=408,transient_updates=3,self_updates_408=4048,self_updates_407=44,
                   exact_regional_payment_counts=True,actual_codeword=False,directed_certificate=False,
                   seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
