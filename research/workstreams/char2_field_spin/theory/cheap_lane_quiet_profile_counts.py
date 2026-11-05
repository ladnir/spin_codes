"""Singleton-deletion source counts and a fixed three-mark quiet-window screen.

For M routing pairs, weight r, and matching0 full count f, let s=r-2f.
Delete the s singleton pairs from the indexed pair cycle. Each deleted
01 or10 pair decreases the other matching's11-edge count by0 or1.
The remaining00/11 sequence has length L=M-s, weight f, and run count
c in [f-other,f-other+s]. Thus a profile count is at most
 C(M,s)*2**s*sum_c circular_run_count(L,f,c).
Selecting the singleton positions, their orientations, and the compressed
sequence recovers the original support. The compressed sequence starts at
the least remaining original position; no extra rotation factor is paid.
When L=0, use the safe2**s bound instead. The all-zero/all-one compressed
cycles have one sequence with zero runs, as in the pinned count helper.

Take the minimum of the old count and both deletion orientations. Exact
small-cycle tests compare every weight/profile with direct counts. The
numerical screen then repeats the same three complete rank-mark bounds
at size256, with unchanged source and route interfaces. Only the source
count changes. A positive-coefficient FLOAT result is not a certificate.

One fixed48-second invocation,20GiB free-space guard, estimated peak
below64MiB. No frozen source edits, generated files, or bytecode.
"""

import hashlib
import heapq
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
PIN="5d1402f5869d7cd3ffcfd97e8414fd58f983a6721b45ec723aa934b869bbeb37"


def load_marks(deadline):
    path=Path(__file__).with_name("cheap_lane_quiet_profile_marks.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen three-mark helper hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_quiet_profile_marks",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module,module.load_screen(deadline)


class DeletionCounts:
    """Exact integer upper bounds; cumulative run counts are cached in RAM."""
    def __init__(self,original):
        self.original=original
        self.records={}

    def orientation(self,pairs,weight,full,other):
        single=weight-2*full
        if single<0 or full<0 or single+full>pairs:
            return 0
        length=pairs-single
        if length==0:
            return 1<<single
        key=(pairs,weight,full)
        if key not in self.records:
            maximum=min(full,length-full)
            prefix=[0]
            for runs in range(maximum+1):
                prefix.append(prefix[-1]+self.original.circular_run_count(length,full,runs))
            self.records[key]=(math.comb(pairs,single)*(1<<single),prefix)
        factor,prefix=self.records[key]
        low=max(0,full-other)
        high=min(len(prefix)-2,full-other+single)
        if low>high:
            return 0
        return factor*(prefix[high+1]-prefix[low])

    def profile(self,pairs,weight,first,second):
        old=self.original.profile_count_upper(pairs,weight,first,second)
        if not old:
            return 0
        return min(old,self.orientation(pairs,weight,first,second),
                   self.orientation(pairs,weight,second,first))


def checks(counter,screen):
    compared=positive=zero=0
    for pairs in range(1,8):
        exact=counter.original.exact_profiles(pairs)
        for weight in range(2*pairs+1):
            for first in range(pairs+1):
                for second in range(pairs+1):
                    target=exact.get((weight,first,second),0)
                    old=counter.original.profile_count_upper(pairs,weight,first,second)
                    new=counter.profile(pairs,weight,first,second)
                    assert target<=new<=old
                    assert new==counter.profile(pairs,weight,second,first)
                    positive+=target>0; zero+=target==0; compared+=1
        screen.guard()
    assert counter.profile(4096,256,128,128)==0
    assert counter.profile(4096,256,128,127)==4096
    assert counter.profile(4096,0,0,0)==1
    assert counter.profile(4096,8192,4096,4096)==1
    return {"all_small_weight_profile_comparisons":compared,"positive_exact_profiles":positive,
            "zero_exact_profiles":zero,"large_boundary_checks":4}


def main():
    if len(sys.argv)!=1:
        raise ValueError("fixed source-count screen accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the actual source volume")
    started=time.monotonic(); marks,screen=load_marks(started+48)
    counter=DeletionCounts(screen.load_counts())
    tests=checks(counter,screen)
    tests["inherited_coefficient_mode_checks"]=screen.checks(screen.terms(.75))
    data=[]
    for name,y in marks.MARKS:
        rows=screen.terms(y); regional=[]
        for full in range(129):
            row=screen.optimize(rows,full)
            regional.append(row["regional_plain_log2"]
                            +screen.multinomial_mode_log(row["u"],row["v"])/screen.LN2)
            screen.guard()
        alpha=-math.log2(y)/8
        common=-32768*math.log2(y)+alpha*(8192*math.log2(256/255)-math.log2(255))+256*math.log2(85)
        data.append((name,common,regional))
    total=old_total=-math.inf; worst=[]; old_worst=[]
    retained=removed=positive=negative=improved=0
    usage={name:0 for name,_ in marks.MARKS}
    max_saved=0.; saving_profile=None
    selected={c:{"count":0,"minimum":math.inf,"maximum":-math.inf}
              for c in (1,4,16,64,80,128,196,256)}
    for first in range(129):
        for second in range(129):
            if first+second>255:
                continue
            old=counter.original.profile_count_upper(4096,256,first,second)
            if not old:
                continue
            base,name=min((common+128*(regional[first]+regional[second]),name)
                          for name,common,regional in data)
            old_score=base+math.log2(old)
            old_total=screen.logadd(old_total,old_score)
            marks.keep_worst(old_worst,(old_score,first,second))
            new=counter.profile(4096,256,first,second)
            if not new:
                removed+=1
                continue
            retained+=1; usage[name]+=1
            saved=math.log2(old)-math.log2(new)
            improved+=new<old
            if saved>max_saved:
                max_saved=saved; saving_profile=[first,second]
            score=base+math.log2(new); runs=256-first-second
            positive+=score>=0; negative+=score<0
            total=screen.logadd(total,score)
            marks.keep_worst(worst,(score,first,second,runs,name,saved))
            if runs in selected:
                row=selected[runs]; row["count"]+=1
                row["minimum"]=min(row["minimum"],score); row["maximum"]=max(row["maximum"],score)
        screen.guard()
    print(json.dumps({"scope":"FLOAT actual size256 quiet-window union with audited singleton-deletion counts",
        "whole_code_certificate":False,"directed":False,
        "log2_old_count_union":old_total,"log2_new_count_union":total,
        "retained_profiles":retained,"removed_impossible_profiles":removed,
        "improved_count_profiles":improved,"positive_profiles":positive,"negative_profiles":negative,
        "mark_usage":usage,"maximum_source_count_saving_bits":max_saved,
        "maximum_saving_profile":saving_profile,
        "worst_profiles":[{"log2_bound":score,"full0":f0,"full1":f1,"runs":c,
                            "mark":name,"source_count_saving_bits":saved}
                           for score,f0,f1,c,name,saved in sorted(worst,reverse=True)],
        "old_worst_profiles":[{"log2_bound":score,"full0":f0,"full1":f1,
            "new_count_log2":math.log2(counter.profile(4096,256,f0,f1)),
            "count_saving_bits":math.log2(counter.original.profile_count_upper(4096,256,f0,f1))
                                -math.log2(counter.profile(4096,256,f0,f1))}
                              for score,f0,f1 in sorted(old_worst,reverse=True)],
        "selected_run_counts":selected,"tests":tests,"marks_helper_sha256":PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "seconds":time.monotonic()-started,"generated_bytes":0,"estimated_peak_bytes":64*(1<<20)
    },indent=2))


if __name__=="__main__":
    main()
