"""Complete r2048 FLOAT bank using the ordered two-sided rho correction.

Thin scoped wrapper around the frozen fine-bank orchestration. Build the
coherent coarse/N4/ordered tables once. Alias the FINISHED ordered cap
into the inherited fine129 selected key before all31 marker fits and all
1025 marginal extractions. Exact run masses and every ordered finite
chain remain in the frozen complete join; no selected term is added.

Explicit --z/--alpha and other options are forwarded unchanged. Defaults
remain z23/100,alpha3/32. The inherited output-only terminal, source once,
own target normalizers, FFT reliability rules and explicit ordinary
fallbacks are unchanged. No extra alpha or source count is introduced.

Only private modules/hooks are patched, then restored. The final print
adapter corrects provenance text only; no numerical field is changed.
FLOAT reconnaissance, not a certificate or other-weight claim. Parent
execution only; inherited180/300s,512MiB,20GiB,stdout-only guards apply.
No numerical work on import, generated files, or frozen source edits.
"""

import builtins
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

sys.dont_write_bytecode=True
BANK_PIN="a82044f15e0cfdf20dd85bd92818f49d42c5a24be707ad55eedec077c8f83625"
ORDERED_PIN="f816addd50969a83ad951567944d905938ece966c9260175e3b4a08b1be53255"


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def main():
    bank=load("cheap_lane_rho_fine_bank.py",BANK_PIN,"ordered_bank_frozen_bank")
    ordered=load("cheap_lane_rho_ordered_buckets.py",ORDERED_PIN,"ordered_bank_constructor")
    restorations=[]
    def replace(obj,name,value):
        restorations.append((obj,name,getattr(obj,name))); setattr(obj,name,value)
    original_load=bank.load
    def intercepted_load(name,pin,label):
        module=original_load(name,pin,label)
        if name=="cheap_lane_rho_fine_buckets.py":
            original_tests=module.exact_tests
            def builder(*values,**kwargs):
                # ordered loads its own unmodified frozen N4 module. This
                # avoids recursive invocation of this private substitution.
                tables,report=ordered.ordered_bucket_pair_logs(*values,**kwargs)
                tables["retained_N4_fine129_bucket_sum"]=tables["fine129_bucket_sum"]
                tables["fine129_bucket_sum"]=tables["ordered129_bucket_sum"]
                report["inherited_fine_key_holds_ordered_finished_cap"]=True
                report["all_center_fits_and_marginals_use_ordered_cap"]=True
                return tables,report
            def tests(pair,correction,np,guard):
                return {"frozen_N4":original_tests(pair,correction,np,guard),
                    "ordered":ordered.exact_tests(module,pair,correction,np,guard)}
            replace(module,"fine_bucket_pair_logs",builder); replace(module,"exact_tests",tests)
        return module
    emitted=False
    def emit(*values,**kwargs):
        nonlocal emitted
        if len(values)==1 and isinstance(values[0],str):
            row=json.loads(values[0])
            if isinstance(row,dict) and "result" in row and "pair_bound" in row:
                row["N4_wrapper_pair_bound_label_overridden"]=row["pair_bound"]
                row["pair_bound"]="ordered two-sided129 ratio buckets: exact off-diagonal minima, diagonal rho witness; finished upper retains N4/coarse/Jensen ceilings"
                row["cap_override"]={"all_fits_marginals_and_join_use_ordered":True,
                    "legacy_fine129_key_is_private_alias":True,
                    "source_count_terminal_unchanged":True,"selected_anchor_not_added":True}
                row["ordered_wrapper_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
                row["frozen_fine_bank_pin"]=BANK_PIN; row["ordered_constructor_pin"]=ORDERED_PIN
                values=(json.dumps(row,indent=2),); emitted=True
        builtins.print(*values,**kwargs)
    try:
        replace(bank,"load",intercepted_load)
        # This is a private module binding, not the process builtins module.
        replace(bank,"builtins",SimpleNamespace(print=emit))
        bank.main()
        if not any(x in ("-h","--help") for x in sys.argv[1:]): assert emitted
    finally:
        for obj,name,value in reversed(restorations): setattr(obj,name,value)


if __name__=="__main__":
    main()
