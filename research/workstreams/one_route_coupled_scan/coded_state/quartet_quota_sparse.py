"""Diagnostic: exact selected input counts and finer joint age buckets.

Not full count coverage or a distance certificate. Serial Peach,1200s,
768MiB conservative matrices plus <512MiB counter workspace,no files.
"""
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_joint_regional as joint
import quartet_mixed_counts as mixed
import quartet_fresh_count_screen as fresh
from quartet_quota_potential import normalize,fit


def main(align_old=False):
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(1200,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        if align_old:
            before=float(np.log2(fresh.H[4]/fresh.H[3]))
            fresh.H[3]=fresh.H[4]*2**(fresh.LOGB*(2*fresh.LAMBDA-fresh.TAU))
            fresh.D[2]=fresh.H[3]
            assert abs(fresh.D[3]/fresh.H[3]-fresh.D[4]/fresh.H[4])<1e-12
            print(dict(proof_potential_only=True,old_log2_h4_over_h3=before,
                       aligned_log2_h4_over_h3=float(np.log2(fresh.H[4]/fresh.H[3]))),flush=True)
        counts=(0,312,384,392,407,408,424,512);joint.BIN=16;mixed.BIN=16
        dimension=sum(len(joint.states(a)[0]) for a in counts)
        assert 6*dimension*dimension*8+4*513*dimension*8+32*(1<<20)<768*(1<<20)
        assert 2*8193*(dimension+1)*8<512*(1<<20)
        base,scale,labels,offsets=joint.build(counts)
        cells=[];edge_ids=np.full((513,33,33),-1,dtype=int)
        for i,a in enumerate(counts):
            for j,(rl,rh,sl,sh) in enumerate(joint.states(a)[0]):
                assert len(cells)==offsets[i]+j
                edge_ids[a,rl//16,sl//16]=len(cells);cells.append((a,a,rl,rh,sl,sh))
        targets=np.arange(len(cells))
        P,x,logv,growth,initial_v=normalize(base,scale,cells,edge_ids,labels,targets)
        result=fit(P,labels,targets,cells,x,logv,growth,initial_v,iterations=4)
        result.update(total_seconds=time.monotonic()-start,cells=len(cells),age_bucket_width=16,
                      aligned_old_potential=align_old)
        print(dict(restricted_best=result),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
