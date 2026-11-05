"""Complete r2048 FLOAT calibrated bank using the fine rational rho cap.

Thin scoped orchestration of frozen helpers: build coarse/fine caps once,
alias the FINISHED fine129_bucket_sum cap into the bank's selected key,
then fit all31 centers and evaluate every marginal using that same cap.
Default z23/100,alpha3/32; explicit rational overrides select one witness,
not a search. Source/count/terminal and ordinary fallback behavior remain
the frozen complete-bank implementation. Fine tests augment old tests.

The private print adapter changes final provenance text only, because the
frozen bank's literal cap label describes its original coarse constructor.
No numerical result is modified. All module fields, hooks and argv are
restored. Parent owns execution; inherited180/300s,512MiB,20GiB,stdout-only
guards apply. No generated files or whole-code/certificate claim.
"""

import argparse
import builtins
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
BANK_PIN="e9ff3adda8aa363383386b8e8665b40da5f1075a8e43ddf8ffc1854ddea60776"
FINE_PIN="ca8a73442ff87ca46ce64192dd999967716cf145c65a9712553c0e81a9a69701"


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__,add_help=False)
    parser.add_argument("--z",type=F,default=F(23,100))
    parser.add_argument("--alpha",type=F,default=F(3,32))
    args,remaining=parser.parse_known_args()
    if not 0<args.z<1 or not 0<args.alpha<=1: parser.error("require0<z<1,0<alpha<=1")
    if max(abs(n).bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64bits")
    bank=load("cheap_lane_rho_calibrated_bank.py",BANK_PIN,"fine_bank_frozen_bank")
    fine=load("cheap_lane_rho_fine_buckets.py",FINE_PIN,"fine_bank_constructor")
    restorations=[]
    def replace(obj,name,value):
        exists=hasattr(obj,name)
        restorations.append((obj,name,exists,getattr(obj,name,None))); setattr(obj,name,value)
    original_load=bank.load
    def intercepted_load(name,pin,label):
        module=original_load(name,pin,label)
        if name=="cheap_lane_rho_bucket_control.py":
            original_builder,original_tests=module.bucket_pair_logs,module.exact_tests
            def builder(*values,**kwargs):
                tables,report=fine.fine_bucket_pair_logs(*values,**kwargs,coarse_builder=original_builder)
                tables["coarse_dense33_bucket_sum"]=tables["dense33_bucket_sum"]
                tables["dense33_bucket_sum"]=tables["fine129_bucket_sum"]
                report["selected_key_alias"]="dense33_bucket_sum now holds fine129_bucket_sum"
                report["all_center_fits_and_replays_use_fine_cap"]=True
                return tables,report
            def tests(pair,correction,np,guard):
                return {"coarse":original_tests(pair,correction,np,guard),
                    "fine":fine.exact_tests(pair,correction,np,guard)}
            replace(module,"bucket_pair_logs",builder); replace(module,"exact_tests",tests)
        return module
    emitted_final=False
    def emit(*values,**kwargs):
        nonlocal emitted_final
        if len(values)==1 and isinstance(values[0],str):
            row=json.loads(values[0])
            if isinstance(row,dict) and "result" in row and "pair_bound" in row:
                row["frozen_pair_bound_label_overridden"]=row["pair_bound"]
                row["pair_bound"]="fine129 disjoint-left rational buckets, intersected with coherent frozen coarse/Jensen upper cap"
                row["cap_override"]={"all_fits_marginals_and_join_use_fine":True,
                    "coarse_selected_key_is_private_alias":True,"source_count_terminal_unchanged":True}
                row["wrapper_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
                row["frozen_bank_pin"]=BANK_PIN; row["fine_constructor_pin"]=FINE_PIN
                values=(json.dumps(row,indent=2),); emitted_final=True
        builtins.print(*values,**kwargs)
    original_argv=sys.argv[:]
    try:
        replace(bank,"Z",args.z); replace(bank,"ALPHA",args.alpha)
        replace(bank,"load",intercepted_load); replace(bank,"print",emit)
        sys.argv=[sys.argv[0],*remaining]
        bank.main()
        if not any(x in ("-h","--help") for x in remaining): assert emitted_final
    finally:
        sys.argv=original_argv
        for obj,name,exists,old in reversed(restorations):
            if exists: setattr(obj,name,old)
            else: delattr(obj,name)


if __name__=="__main__":
    main()
