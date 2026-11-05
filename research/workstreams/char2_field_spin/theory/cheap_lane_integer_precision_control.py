"""Selected high-rectangle precision control for the fixed r1024 witness.

Build the pinned coherent local pair table ONCE at256-bit Arb precision.
Run the unchanged exact high rectangle(center472,F512,S158) at128 and192
integer coefficient bits. Compare its own normalized regular/terminal
entries at f471 and f512 and both ordered finite chains. Markers, z, alpha,
local caps, degree truncation and source factor are identical.

The common matrix mantissa is NOT a per-coefficient relative precision
claim. Tiny coefficients can incur larger relative upward-rounding loss.
This control measures that loss without reoptimization. It evaluates only
one regional table and two ordered profiles, not the complete weight sum.
Charging a selected chain all of run41's source mass is a labeled scalar
diagnostic; it does not bound other profiles in that run.

No frozen edits, coefficient files or numerical work on import.
One parent-run invocation:180s,512MiB,20GiB free, stdout only.
"""

from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
DRIVER_PIN="a2a2a6d784b5ae23f26b93e2c4e3caec87f7667598f8bd76d70ab3b49614e9a0"
BITS=(128,192)
FULLS=(471,512)
SECONDS,MEMORY=180,512*(1<<20)


def load_driver():
    path=Path(__file__).with_name("cheap_lane_pair_regions_integer.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=DRIVER_PIN:
        raise RuntimeError("frozen integer regional driver mismatch")
    spec=importlib.util.spec_from_file_location("precision_control_driver",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def log2(v,value):
    if not value>0: raise ArithmeticError("selected positive coefficient not resolved")
    return value.log()/v.arb(2).log()


def main():
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+SECONDS
    runner=load_driver()
    driver=runner.load("cheap_lane_pair_regions_directed.py",runner.DRIVER_PIN,"precision_frozen_schedule")
    integer=runner.load("cheap_lane_integer_regions.py",runner.INTEGER_PIN,"precision_integer_arithmetic")
    anchor=driver.load("cheap_lane_cone_cauchy_anchors.py",driver.ANCHOR_PIN,"precision_anchor")
    local=driver.load("cheap_lane_pair_local_directed.py",driver.LOCAL_PIN,"precision_local")
    join=driver.load("cheap_lane_directed_run_mass.py",driver.JOIN_PIN,"precision_join")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("180-second precision control deadline; no retry")
        if peak()>MEMORY: raise MemoryError("512MiB precision control budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"elapsed_seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    built=local.build_local_pairs(driver.Z,driver.ALPHA,anchor,guard,deadline,256,progress=progress)
    v,core=built["verifier"],built["core"]
    del built["coefficients"]
    tests={"integer_arithmetic":integer.exact_tests(guard),"run_join":join.exact_tests(guard)}
    item=next(x for x in driver.SCHEDULE if x[0]==472)
    assert item==(472,433,512,"exact") and all(item[1]<=f<=item[2] for f in FULLS)
    h,c_h,eta=anchor.source_factor(driver.WEIGHT)
    source=(v.exact(eta)*v.exact(core.BETA)**driver.WEIGHT
            *v.exact(driver.Z)**(-core.LIMIT))**v.exact(driver.ALPHA)
    runs=driver.WEIGHT-sum(FULLS)
    assert runs==41
    run_mass=join.run_counts()[runs]
    retained={}; reports=[]; original_bits=runner.BITS
    try:
        for bits in BITS:
            guard(); runner.BITS=bits
            def tagged_progress(phase,**values):
                progress(phase,coefficient_bits=bits,**values)
            regions,table_report=runner.one_table(v,built["rational_coefficients"],item,
                driver,integer,guard,tagged_progress)
            assert table_report["arithmetic"]["mantissa_bits"]==bits
            selected={f:regions[f] for f in FULLS}
            del regions
            chains=join.two_orientations(selected[FULLS[0]],selected[FULLS[1]])
            row={"coefficient_bits":bits,"table":table_report,
                "normalized_entries_log2_intervals":{str(f):{
                    "regular":[str(log2(v,x)) for x in selected[f][0]],
                    "terminal":[str(log2(v,x)) for x in selected[f][1]]} for f in FULLS},
                "ordered_profiles":[list(FULLS),list(reversed(FULLS))],
                "conditional_chain_log2_intervals":[str(log2(v,x)) for x in chains],
                "source_times_selected_chain_log2_intervals":[str(log2(v,source*x)) for x in chains],
                "full_run_mass_charged_to_selected_chain_log2_intervals":[str(log2(v,source*run_mass*x)) for x in chains],
                "run_mass_diagnostic_is_not_a_run_cover":True}
            retained[bits]=(selected,chains); reports.append(row)
            progress("selected_rectangle_precision_complete",**row)
    finally:
        runner.BITS=original_bits
    low,high=(retained[b] for b in BITS)
    entry_gains={str(f):{
        "regular":[str(log2(v,x/y)) for x,y in zip(low[0][f][0],high[0][f][0])],
        "terminal":[str(log2(v,x/y)) for x,y in zip(low[0][f][1],high[0][f][1])]} for f in FULLS}
    chain_gains=[str(log2(v,x/y)) for x,y in zip(low[1],high[1])]
    # Exponents and marker normalization come from the frozen implementation;
    # no empirical correction is applied to either output.
    guard()
    print(json.dumps({"scope":"selected high rectangle and two ordered r1024 profile controls only",
        "directed":True,"complete_weight_bound":False,"whole_code_certificate":False,
        "parameters":{"weight":driver.WEIGHT,"z":str(driver.Z),"alpha":str(driver.ALPHA),
            "local_Arb_precision":256,"coefficient_bits":BITS,"center":472,
            "full_exact_limit":512,"single_exact_limit":158,"selected_fulls":FULLS,
            "markers":[str(x) for x in driver.marker(472)]},
        "results":reports,"gain_bits_128_over_192":{"normalized_entries":entry_gains,
            "ordered_conditional_chains":chain_gains},
        "same_coherent_raw_and_local_pair_table":True,
        "coefficient_sha256":built["report"]["rational_coefficient_sha256"],
        "local_report":built["report"],"tests":tests,
        "source_log2_interval":str(log2(v,source)),
        "joint_missing":{"h":h,"c_h":str(c_h),"eta":str(eta)},
        "diagnostic_run":runs,"exact_run_mass_sha256":hashlib.sha256(str(run_mass).encode()).hexdigest(),
        "source_factor_reoptimization":False,"marker_reoptimization":False,
        "driver_pin":DRIVER_PIN,"integer_pin":runner.INTEGER_PIN,"local_pin":driver.LOCAL_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
