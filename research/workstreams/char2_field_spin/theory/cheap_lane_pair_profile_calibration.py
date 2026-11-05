"""Selected r2048 pair-marker diagnosis; no complete-size claim.

Fix z1/4, alpha119/1280, and the actual equal profile(768,768), singles512.
Build the unchanged ell1/2/3 pair caps once. Compare inherited center704,
own-f center768, and at most12 guided coordinate proposals in log markers.
The positive ordinary objective is e0 H^16383 J divided by d^256: H=M*M,
J=M*T, d=M4096(768,512)u^768v^512. Thus32767 updates and the output-only
terminal are retained. No additional alpha is applied during routing.

Positive six-moment products give analytic marker gradients/covariances.
Coordinate steps and line searches only choose witnesses. The chosen
markers are positive rationals, then ordinary/m256/m512 are replayed
independently with strict FFT controls. Ordinary stationarity is not
claimed for projected bounds. All outputs remain FLOAT diagnostics.

Every comparison charges the SAME exact indexed run512 mass, rather
than querying the exact selected profile count. This bounds the selected
profile contribution but is NOT a cover of other profiles in that run.
No source/cap/parameter search or probability-law change. No files/jobs
on import. Parent execution only:180s,512MiB,20GiB free, stdout only.

main(augment=None, augmentation_label=None) is a future wrapper seam:
augment(caps,z,alpha,guard) may replace the cached cap maps before the
pair construction. The default path never invokes an augmentation.
"""

from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from math import comb
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
SIZES_PIN="7a7024b78912abb78eaaef9e90f766eefdd7e7f807d8dfe6143f38c86c233afd"
CALIBRATION_PIN="e9f90bdda919247df80e71e5d2555f24ef7ca12323d95b1df622b22556327c56"
WEIGHT,FULL,SINGLE,RUNS=2048,768,512,512
Z,ALPHA=F(1,4),F(119,1280)
SECONDS,MEMORY,MAX_PROPOSALS=180,512*(1<<20),12
MODULI=(256,512)
LN2=math.log(2)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


class Objective:
    def __init__(self,logs,constant,calibration,np,guard):
        self.np,self.guard=np,guard; self.jets=calibration.Jets(2,np)
        self.fd=np.asarray([f for f in range(65) for s in range(65-f)],dtype=float)
        self.sd=np.asarray([s for f in range(65) for s in range(65-f)],dtype=float)
        self.local=np.asarray([logs[int(f),int(s)] for f,s in zip(self.fd,self.sd)])
        self.features=self.jets.moments(np.asarray([self.fd,self.sd]))
        self.targets=256*np.asarray([FULL,SINGLE],dtype=float); self.constant=constant

    def evaluate(self,x):
        self.guard(); np,jets=self.np,self.jets
        logs=self.local+(self.fd*x[0]+self.sd*x[1])[:,None]
        shift=float(np.max(logs)); mass=np.exp(logs-shift)
        if np.any(np.isfinite(logs)&(mass==0)):
            raise ArithmeticError("positive ordinary moment term underflow")
        values=self.features@mass
        h=jets.normalize(values[:,:4].reshape(jets.channels,2,2),shift)
        j=jets.normalize(values[:,4:].reshape(jets.channels,2,1),shift)
        result,scale=jets.multiply(jets.power(h,16383),j)
        normalized=result[:,0,0]/result[0,0,0]
        means=normalized[1:3]; covariance=np.empty((2,2))
        for channel,(a,b) in enumerate(jets.pairs,3):
            covariance[a,b]=covariance[b,a]=normalized[channel]-means[a]*means[b]
        score=scale+math.log(result[0,0,0])-float(self.targets@x)+self.constant
        if not math.isfinite(score) or not np.isfinite(covariance).all():
            raise ArithmeticError("nonfinite ordinary proposal objective")
        return {"x":x.copy(),"objective":score,"means":means,
                "gradient":means-self.targets,"covariance":covariance}


def coordinate_fit(objective,start,np,guard):
    """At most12 positive-moment evaluations, including rejected proposals."""
    history=[]; best=None
    def evaluate(x):
        nonlocal best
        if len(history)>=MAX_PROPOSALS: raise RuntimeError("proposal budget exceeded")
        point=objective.evaluate(x)
        history.append({"evaluation":len(history)+1,"log2_ordinary_FLOAT":point["objective"]/LN2,
            "markers_FLOAT":[math.exp(float(a)) for a in x],
            "tilted_region_means_FLOAT":[float(a/256) for a in point["means"]],
            "log_marker_gradient_FLOAT":[float(a/LN2) for a in point["gradient"]]})
        if best is None or point["objective"]<best["objective"]: best=point
        return point
    current=evaluate(np.asarray(start,dtype=float)); status="proposal budget exhausted"
    while len(history)<MAX_PROPOSALS:
        guard()
        if np.max(np.abs(current["gradient"]))/256<1e-5:
            status="ordinary first-order tolerance reached"; break
        diagonal=np.maximum(np.diag(current["covariance"]),1.0)
        axis=int(np.argmax(np.abs(current["gradient"])/np.sqrt(diagonal)))
        direction=float(np.clip(-current["gradient"][axis]/diagonal[axis],-2,2))
        if direction==0: status="zero proposal step"; break
        accepted=False; step=1.0
        while len(history)<MAX_PROPOSALS:
            trial=current["x"].copy(); trial[axis]=np.clip(trial[axis]+step*direction,-20,8)
            displacement=float(trial[axis]-current["x"][axis])
            if abs(displacement)<1e-12: status="coordinate box or tiny step"; break
            candidate=evaluate(trial)
            if candidate["objective"]<=current["objective"]+1e-4*current["gradient"][axis]*displacement:
                current=candidate; accepted=True; break
            step*=0.5
        if not accepted: break
    assert best is not None and len(history)<=MAX_PROPOSALS
    return best,history,status


def replay(log_coefficients,markers,modulus,source_count,broad,cells,fft,screen,np,guard):
    """One target only; no bulk marginal assignment or ordinary substitution."""
    guard(); fd=np.asarray([f for f in range(65) for s in range(65-f)])
    sd=np.asarray([s for f in range(65) for s in range(65-f)])
    logs=np.asarray([log_coefficients[f,s] for f,s in zip(fd,sd)])
    lu,lv=map(broad.logq,markers); logs=logs+(fd*lu+sd*lv)[:,None]
    channel_logs=[]
    for j in range(6):
        peak=float(np.max(logs[:,j]))
        channel_logs.append(peak+math.log(float(np.sum(np.exp(logs[:,j]-peak)))))
    def normalized(values):
        scale=max(values); entries=tuple(math.exp(x-scale) for x in values)
        if not all(math.isfinite(x) and x>0 for x in entries):
            raise ArithmeticError("ordinary selected channel underflow")
        return screen.normalized_matrix(entries,scale)
    p,t=normalized(channel_logs[:4]),normalized(channel_logs[4:])
    ordinary=(screen.matrix_power(p,64),screen.matrix_vector(screen.matrix_power(p,63),t))
    denominator=math.log(comb(4096,FULL)*comb(4096-FULL,SINGLE))+FULL*lu+SINGLE*lv
    def report(region):
        conditional=screen.paired_chains(region,region)[0]
        return {"status":"resolved","log2_selected_run_mass_bound_FLOAT":(source_count+conditional)/LN2,
                "log2_conditional_chain_FLOAT":conditional/LN2}
    result={"ordinary":report((*ordinary,denominator))}
    try:
        shift=float(np.max(logs)); mass=np.exp(logs-shift)
        if np.any(np.isfinite(logs)&(mass==0)):
            raise ArithmeticError("positive selected FFT input underflow")
        degree=np.zeros((modulus,modulus,6)); np.add.at(degree,(fd%modulus,sd%modulus),mass)
        spectra=fft.regional_spectra(degree[:,:,:4].reshape(modulus,modulus,2,2),
            degree[:,:,4:].reshape(modulus,modulus,2,1),np,guard,doublings=6)
        chosen=[]; diagnostics=[]
        for (raw,scale),expected in zip(spectra,ordinary):
            scale+=64*shift
            fft.scaled_check(tuple(float(x.real) for x in raw[0,0].ravel()),scale,expected)
            residue,_,_,diagnostic=fft.recover((raw,scale),np)
            thresholds=cells.selection_thresholds(fft,residue,raw,np)
            values,minimum=cells.select_entries(residue,(FULL%modulus,SINGLE%modulus),thresholds)
            chosen.append(screen.normalized_matrix(values,scale))
            diagnostics.append({**diagnostic,"selected_minimum_over_DC":minimum})
        result["residue"]={**report((*chosen,denominator)),"FFT_controls":diagnostics}
    except ArithmeticError as error:
        result["residue"]={"status":"unresolved","reason":str(error)}
    guard(); return result


def main(augment=None,augmentation_label=None):
    if (augment is None)!=(augmentation_label is None):
        raise ValueError("augmentation callable and explicit provenance label must be supplied together")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+SECONDS
    sizes=load("cheap_lane_pair_sizes_screen.py",SIZES_PIN,"selected_sizes")
    calibration=load("cheap_lane_broad_profile_calibration.py",CALIBRATION_PIN,"selected_moments")
    bank=sizes.load("cheap_lane_pair_bank_screen.py",sizes.BANK_PIN,"selected_bank")
    broad=sizes.load("cheap_lane_broad_peeling_screen.py",sizes.BROAD_PIN,"selected_broad")
    composition=bank.load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"selected_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"selected_triple")
    pair=correction.load_pair()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"selected_neighbor")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"selected_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("180-second selected diagnostic; no automatic retry")
        if peak()>MEMORY: raise MemoryError("512MiB selected diagnostic budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,"peak_process_bytes":peak(),**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    tests={"exact_positive_moments":calibration.exact_tests(np,screen,old),
           "old_six_moment_tests":old.proposal_tests(np),"FFT_exact_toys":fft.exact_toy(np,guard),
           "exact_scale_tests":cells.exact_scale_tests()}
    caps=composition.build_caps(Z,ALPHA,triple,anchor,guard,deadline)
    if augment is not None: caps=augment(caps,Z,ALPHA,guard)
    verifier,core=caps["verifier"],caps["core"]
    tables,single,controls=correction.corrected_pair_logs(caps["powered"],float(ALPHA),
        tuple(range(-16,17,4)),pair,old,np,guard,raw=caps["raw"])
    del caps["raw"],caps["powered"],single
    selected=tables["concavity_Jensen_minimum"]; del tables
    missing,c_h,eta=anchor.source_factor(WEIGHT)
    log_source=float((verifier.exact(ALPHA)*(verifier.exact(eta).log()+WEIGHT*verifier.exact(core.BETA).log()
        -core.LIMIT*verifier.exact(Z).log())).mid())
    numerator=8192*comb(WEIGHT-1,RUNS-1)*comb(8192-WEIGHT-1,RUNS-1)
    run_mass,remainder=divmod(numerator,RUNS); assert remainder==0
    assert RUNS==WEIGHT-2*FULL and SINGLE==WEIGHT-2*FULL
    source_count=log_source+math.log(run_mass)
    constant=source_count-256*math.log(comb(4096,FULL)*comb(4096-FULL,SINGLE))
    objective=Objective(selected,constant,calibration,np,guard)
    center=min(broad.default_centers(WEIGHT),key=lambda c:(abs(c-FULL),c)); assert center==704
    inherited=sizes.markers(WEIGHT,center); own=sizes.markers(WEIGHT,FULL)
    best,history,status=coordinate_fit(objective,np.asarray(list(map(broad.logq,own))),np,guard)
    fitted=tuple(max(F(1,1<<30),F(math.exp(float(x))).limit_denominator(1<<20)) for x in best["x"])
    progress("selected_coordinate_proposals_complete",proposals=len(history),status=status,
             fitted_markers=list(map(str,fitted)),same_caps=True)
    results=[]
    for label,markers in (("inherited_bank",inherited),("own_full_curve",own),("fitted_rational",fitted)):
        direct=objective.evaluate(np.asarray(list(map(broad.logq,markers))))
        checks=[]
        for modulus in MODULI:
            row=replay(selected,markers,modulus,source_count,broad,cells,fft,screen,np,guard)
            difference=abs(row["ordinary"]["log2_selected_run_mass_bound_FLOAT"]-direct["objective"]/LN2)
            assert difference<2e-5
            row.update({"modulus":modulus,"ordinary_two_path_control_bits":difference})
            checks.append(row)
            progress("selected_marker_residue_complete",witness=label,markers=list(map(str,markers)),**row)
        results.append({"witness":label,"markers":list(map(str,markers)),
            "ordinary_log_marker_gradient_FLOAT":[float(x/LN2) for x in direct["gradient"]],
            "tilted_region_means_FLOAT":[float(x/256) for x in direct["means"]],
            "ordinary_covariance_FLOAT":direct["covariance"].tolist(),"replays":checks})
    guard()
    print(json.dumps({"scope":"selected actual r2048 equal profile768,768 only; each term charged all run512 mass",
        "directed":False,"complete_weight_bound":False,"run_class_cover":False,"whole_code_certificate":False,
        "negative_is_not_a_certificate":True,"parameters":{"weight":WEIGHT,"fulls":[FULL,FULL],
            "singles":[SINGLE,SINGLE],"z":str(Z),"alpha":str(ALPHA),"moduli":MODULI,
            "inherited_center":center,"max_proposals":MAX_PROPOSALS,"seconds":SECONDS,"memory_budget_bytes":MEMORY},
        "fixed_source":{"log2_source_FLOAT":log_source/LN2,"log2_run_mass":math.log(run_mass)/LN2,
            "exact_run_mass_sha256":hashlib.sha256(str(run_mass).encode()).hexdigest(),
            "same_count_in_every_comparison":True,"exact_profile_queries":0,
            "joint_missing_h":missing,"c_h":str(c_h),"eta":str(eta)},
        "fit":{"method":"bounded analytic coordinate Newton with Armijo backtracking",
            "status":status,"evaluations":len(history),"history":history,
            "ordinary_stationarity_is_not_projected_stationarity":True},
        "results":results,"augmentation":augmentation_label,
        "coset_receipt":caps.get("coset_receipt"),"controls":controls,"tests":tests,
        "field_tests":caps["tests"],"block_caps":caps["block_report"],
        "sizes_pin":SIZES_PIN,"moment_pin":CALIBRATION_PIN,"composition_pin":bank.COMPOSITION_PIN,
        "correction_pin":composition.CORRECTION_PIN,"triple_pin":composition.TRIPLE_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
