"""Bounded FLOAT reconnaissance for one complete source weight, 7<=r<=4096.

For a fixed positive marker bank, extract every marginal (f,r-2f) from
positive bivariate residues. Each target keeps its own route denominator.
Local caps are powered once, before route averaging. The finite chain is
e0 (R0 R1)^127 R0 V1, including its actual output-only terminal. Ordered
profiles are streamed by run count c=r-f0-f1; opposite orientations are
not identified. All conditional chains are evaluated, without a Perron
replacement or an O(r^2) retained table.

Within each run, sort conditional chain caps. Query exact symmetric BEST
counts only while the remaining POST-source tail exceeds 2^-80/r and the
count budget remains. After any prefix, the valid envelope is
 sum(queried count*chain)+(exact run mass-queried mass)*largest remaining cap.
Count exhaustion retains this entire residual envelope. It never removes
unqueried sources. A hard overall deadline aborts without a complete report.

All numerical chain/FFT values are FLOAT, so this is NOT a certificate.
The source counts and remaining masses are exact integers. The mathematical
envelope is valid; rounding is not directed. Unresolved selected residues
abort unless --ordinary-fallback explicitly requests an ordinary positive
Cauchy numerator for those marginals. Such fallbacks are reported, never
clipped, silently substituted, or described as residue gains.

The default bank scales the successful r512 fifteen-center bank. Its curve
is u=cf*f/(4096-r+f), v=cs*(r-2f)/(4096-r+f). A zero denominator is replaced
by 1; zero numerators use 2^-30. These are declared positive witnesses,
not approximations to a degenerate distribution. No fresh calibration is
performed. Broader weights may require different supplied z/alpha/markers.
180s default, at most300s;256MiB;20GiB free; stdout only; no data files.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations_with_replacement
import json
import math
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
NEIGHBOR_PIN = "c7cce6b8763dd8145e525061c82a9c2ad51139f2d0d66d8c870de476bb374870"
PEEL_PIN = "f0de323a986f98c2d1928cfd2f6606b289f87fc39f5e3cbf993615e34f8441dd"
REFERENCE_CENTERS = (16,48,80,112,144,176,208,224,232,236,240,244,248,252,254)
ZERO_MARKER = F(1, 1 << 30)
TAIL_BITS = 80
LN2 = math.log(2)


def load(name, pin, tag):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
        raise RuntimeError("frozen dependency hash mismatch: "+name)
    spec = importlib.util.spec_from_file_location(tag, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def logq(x):
    return math.log(x.numerator)-math.log(x.denominator)


def marker_curve(r, center, cf, cs):
    empty = 4096-r+center
    assert empty >= 0 and 0 <= 2*center <= r
    denominator = max(1, empty)
    return (max(ZERO_MARKER, cf*F(center, denominator)),
            max(ZERO_MARKER, cs*F(r-2*center, denominator)))


def default_centers(r):
    # Exact nearest integer, with half-integer ties sent down; r512 replays.
    maximum = r//2
    return tuple(sorted({(2*c*maximum+255)//512 for c in REFERENCE_CENTERS}))


def run_masses(n, r):
    """Exact recurrence for n/c*C(r-1,c-1)*C(n-r-1,c-1)."""
    value = n
    for c in range(1, min(r, n-r)+1):
        yield c, value
        if c < min(r, n-r):
            value, rem = divmod(value*(r-c)*(n-r-c), c*(c+1))
            assert rem == 0


def batch_product(a, b, np):
    """Normalized, noncommuting2x2 products; one scale per profile."""
    x, lx = a; y, ly = b
    out = np.empty_like(x)
    out[:,0] = x[:,0]*y[:,0]+x[:,1]*y[:,2]
    out[:,1] = x[:,0]*y[:,1]+x[:,1]*y[:,3]
    out[:,2] = x[:,2]*y[:,0]+x[:,3]*y[:,2]
    out[:,3] = x[:,2]*y[:,1]+x[:,3]*y[:,3]
    scale = np.max(out, axis=1)
    if not np.isfinite(out).all() or np.any(scale <= 0) or np.any(out < 0):
        raise ArithmeticError("invalid finite batched matrix product")
    out /= scale[:,None]
    # All input regional matrices are strictly positive. The identity used
    # by powering has zeros, but its products with positive matrices do not.
    if np.any(out == 0):
        raise ArithmeticError("positive finite-chain entry underflowed")
    return out, lx+ly+np.log(scale)


def batch_power(base, exponent, np):
    assert exponent >= 1
    answer = None
    while exponent:
        if exponent & 1:
            answer = base if answer is None else batch_product(answer, base, np)
        exponent >>= 1
        if exponent:
            base = batch_product(base, base, np)
    return answer


def batched_chains(f0, f1, regular, regular_log, terminal, terminal_log, np, exponent=127):
    left = (regular[f0], regular_log[f0])
    right = (regular[f1], regular_log[f1])
    pair = batch_product(left, right, np)
    chain, scale = batch_product(batch_power(pair, exponent, np), left, np)
    end = terminal[f1]
    value = chain[:,0]*end[:,0]+chain[:,1]*end[:,1]
    if np.any(value <= 0) or not np.isfinite(value).all():
        raise ArithmeticError("unresolved finite output-only terminal")
    return scale+terminal_log[f1]+np.log(value)


def exact_tests(np, screen, peel):
    tests = peel.exact_tests()
    recurrence = 0
    for n in range(4, 15, 2):
        for r in range(1, n//2+1):
            rows = list(run_masses(n, r))
            assert sum(value for _, value in rows) == comb(n, r)
            for c, value in rows:
                assert value*c == n*comb(r-1,c-1)*comb(n-r-1,c-1)
                recurrence += 1
    assert default_centers(512) == REFERENCE_CENTERS
    assert marker_curve(4096,0,F(1),F(1)) == (ZERO_MARKER,F(4096))
    assert marker_curve(4096,2048,F(1),F(1)) == (F(1),ZERO_MARKER)
    matrices = ((F(1,3),F(2,5),F(3,7),F(5,11)),
                (F(7,13),F(11,17),F(13,19),F(17,23)),
                (F(19,29),F(23,31),F(29,37),F(31,41)))
    ends = ((F(2,3),F(3,5)),(F(5,7),F(7,11)),(F(11,13),F(13,17)))
    normalized = [screen.normalized_matrix(tuple(map(float,x))) for x in matrices]
    vnorm = [screen.normalized_matrix(tuple(map(float,x))) for x in ends]
    regular=np.asarray([x[0] for x in normalized]); rl=np.asarray([x[1] for x in normalized])
    terminal=np.asarray([x[0] for x in vnorm]); tl=np.asarray([x[1] for x in vnorm])
    a=np.asarray([0,1,2,1]); b=np.asarray([1,0,2,2])

    def product(x,y):
        return (x[0]*y[0]+x[1]*y[2],x[0]*y[1]+x[1]*y[3],
                x[2]*y[0]+x[3]*y[2],x[2]*y[1]+x[3]*y[3])

    checks = 0
    for exponent in (1,2,4):
        logs=batched_chains(a,b,regular,rl,terminal,tl,np,exponent)
        for k,(i,j) in enumerate(zip(a,b)):
            pair=product(matrices[i],matrices[j]); value=(F(1),F(0),F(0),F(1))
            for _ in range(exponent): value=product(value,pair)
            value=product(value,matrices[i])
            truth=value[0]*ends[j][0]+value[1]*ends[j][1]
            assert abs(logs[k]-logq(truth)) < 2e-12
            checks += 1
    logs=batched_chains(a,b,regular,rl,terminal,tl,np)
    for k,(i,j) in enumerate(zip(a,b)):
        reference=screen.paired_chains((normalized[i],vnorm[i],0),
                                       (normalized[j],vnorm[j],0))[0]
        assert abs(logs[k]-reference) < 2e-10
        checks += 1
    tests.update(exact_run_recurrence_checks=recurrence,
                 noncommuting_batched_terminal_checks=checks,
                 reference_bank_and_boundary_markers=True)
    return tests


def marginals(coefficients, r, centers, cf, cs, modulus, ordinary_fallback,
              cells, fft, screen, np, guard, progress):
    assignment=[min(centers,key=lambda c:(abs(f-c),c)) for f in range(r//2+1)]
    indices=[(f,s) for f in range(33) for s in range(33-f)]
    fd=np.asarray([f for f,_ in indices]); sd=np.asarray([s for _,s in indices])
    local=np.asarray([coefficients[f,s] for f,s in indices])
    if np.any(local < 0) or not np.isfinite(local).all():
        raise ArithmeticError("invalid powered local coefficient")
    positive=local>0
    loglocal=np.full(local.shape,-math.inf)
    loglocal[positive]=np.log(local[positive])
    result=[None]*(r//2+1); reports=[]; fallbacks=[]
    stats={"max_imaginary_over_DC":0.0,"max_negative_over_DC":0.0,
           "min_selected_over_DC":1.0}
    for center in centers:
        guard()
        targets=[f for f,c in enumerate(assignment) if c==center]
        if not targets: continue
        u,v=marker_curve(r,center,cf,cs); lu,lv=logq(u),logq(v)
        logs=loglocal+(fd*lu+sd*lv)[:,None]
        shift=float(np.max(logs)); mass=np.exp(logs-shift)
        # A common scalar factor is restored128 times for either regional
        # object: the terminal also contains128 packet factors.
        ordinary_logs=[]
        for j in range(6):
            maximum=float(np.max(logs[:,j]))
            ordinary_logs.append(maximum+math.log(float(np.sum(np.exp(logs[:,j]-maximum)))))

        def from_logs(values):
            scale=max(values); entries=tuple(math.exp(x-scale) for x in values)
            if not all(x>0 for x in entries):
                raise ArithmeticError("ordinary packet channel underflowed")
            return screen.normalized_matrix(entries,scale)

        op=from_logs(ordinary_logs[:4]); ot=from_logs(ordinary_logs[4:])
        ordinary=(screen.matrix_power(op,128),
                  screen.matrix_vector(screen.matrix_power(op,127),ot))
        tables=[]; cell_error=None
        try:
            if np.any(positive & (mass==0)):
                raise ArithmeticError("positive local route term underflowed in FFT input")
            degree=np.zeros((modulus,modulus,6))
            np.add.at(degree,(fd%modulus,sd%modulus),mass)
            p=degree[:,:,:4].reshape(modulus,modulus,2,2)
            t=degree[:,:,4:].reshape(modulus,modulus,2,1)
            spectra=fft.regional_spectra(p,t,np,guard)
            for kind,(raw,scale),expected in zip(("regular","terminal"),spectra,ordinary):
                scale += 128*shift
                fft.scaled_check(tuple(float(x.real) for x in raw[0,0].ravel()),scale,expected)
                residue,_,_,diagnostic=fft.recover((raw,scale),np)
                thresholds=cells.selection_thresholds(fft,residue,raw,np)
                tables.append((residue,scale,thresholds))
                for key in ("max_imaginary_over_DC","max_negative_over_DC"):
                    stats[key]=max(stats[key],diagnostic[key])
        except ArithmeticError as error:
            if not ordinary_fallback: raise ArithmeticError(f"center{center}: {error}") from error
            cell_error=str(error)
        for f in targets:
            guard(); s=r-2*f; error=cell_error
            projected=[]
            if error is None:
                try:
                    for residue,scale,thresholds in tables:
                        values,minimum=cells.select_entries(residue,(f%modulus,s%modulus),thresholds)
                        projected.append(screen.normalized_matrix(values,scale))
                        stats["min_selected_over_DC"]=min(stats["min_selected_over_DC"],minimum)
                except ArithmeticError as failure:
                    if not ordinary_fallback:
                        raise ArithmeticError(f"weight{r},center{center},f{f},s{s}: {failure}") from failure
                    error=str(failure)
            if error is not None:
                projected=ordinary
                fallbacks.append({"full":f,"single":s,"center":center,"reason":error})
            denominator=math.log(comb(4096,f)*comb(4096-f,s))+f*lu+s*lv
            result[f]=(projected[0],projected[1],denominator,u,v)
        report={"center":center,"u":str(u),"v":str(v),"first_full":targets[0],
                "last_full":targets[-1],"marginals":len(targets),
                "ordinary_fallbacks":sum(x["center"]==center for x in fallbacks)}
        reports.append(report); progress("marker_cell_complete",**report)
        # Do not retain any previous cell's coefficient or frequency arrays.
        tables.clear()
        if "spectra" in locals(): del spectra
        if "degree" in locals(): del degree,p,t
        if "residue" in locals(): del residue,raw
    assert all(x is not None for x in result)
    return result,reports,stats,fallbacks


def join(r,regions,counts,screen,np,log_source,query_limit,count_seconds,
         deadline,guard,progress):
    regular=np.asarray([x[0][0] for x in regions])
    terminal=np.asarray([x[1][0] for x in regions])
    rl=np.asarray([x[0][1]-x[2] for x in regions])
    tl=np.asarray([x[1][1]-x[2] for x in regions])
    target=-TAIL_BITS*LN2-math.log(r)
    count_deadline=min(deadline-25,time.monotonic()+count_seconds)
    queries=attempts=peeled=hits=ordered=unmet=0; count_elapsed=0.0; stop_reason=None
    total=all_tails=coarse=-math.inf; source_total=0
    largest=[]; coarse_largest=[]; known_profiles=[]
    maximum=r//2
    for c,number in run_masses(8192,r):
        guard()
        first=np.arange(max(0,r-c-maximum),min(maximum,r-c)+1,dtype=int)
        second=r-c-first
        assert len(first)>0 and np.all((second>=0)&(second<=maximum))
        chain=batched_chains(first,second,regular,rl,terminal,tl,np)
        if not np.isfinite(chain).all(): raise ArithmeticError("nonfinite conditional chain")
        order=np.lexsort((first,-chain))
        ordered += len(first); source_total += number
        coarse_run=log_source+math.log(number)+float(chain[order[0]])
        coarse=screen.logadd(coarse,coarse_run)
        remaining=number; k=0; prefix=-math.inf; cached={}
        while remaining:
            assert k<len(order)
            position=int(order[k]); a,b=int(first[position]),int(second[position])
            current=float(chain[position])
            tail=log_source+math.log(remaining)+current
            if tail<=target: break
            key=(min(a,b),max(a,b))
            if key in cached:
                value=cached[key]; hits+=1
            else:
                if stop_reason is None:
                    if queries>=query_limit: stop_reason="exact-query limit"
                    elif time.monotonic()>=count_deadline: stop_reason="exact-count time reserve"
                if stop_reason is not None: break
                counts.DEADLINE=count_deadline
                before=time.monotonic()
                attempts+=1
                try:
                    value=counts.profile_count(4096,r,*key)
                except RuntimeError as error:
                    if "deadline exhausted" not in str(error) or time.monotonic()<count_deadline:
                        raise
                    stop_reason="exact-count deadline during a query; partial count discarded"
                    break
                finally:
                    count_elapsed+=time.monotonic()-before
                    counts.DEADLINE=deadline
                queries+=1; cached[key]=value
            assert 0<=value<=remaining
            if value:
                weighted=log_source+math.log(value)+current
                known_profiles.append({"full":[a,b],"single":[r-2*a,r-2*b],
                    "runs":c,"log2_weighted_profile_FLOAT":weighted/LN2,
                    "log2_exact_count":math.log(value)/LN2,
                    "log2_conditional_chain_FLOAT":current/LN2})
                known_profiles=sorted(known_profiles,key=lambda x:x["log2_weighted_profile_FLOAT"],reverse=True)[:6]
                prefix=screen.logadd(prefix,math.log(value)+current)
            remaining-=value; k+=1; peeled+=1
        if remaining:
            tail=log_source+math.log(remaining)+float(chain[int(order[k])])
        else:
            tail=-math.inf
        upper=screen.logadd(log_source+prefix,tail)
        assert upper<=coarse_run+1e-7
        met=tail<=target; unmet+=not met
        total=screen.logadd(total,upper); all_tails=screen.logadd(all_tails,tail)
        row={"runs":c,"log2_bound_FLOAT":upper/LN2,
             "post_source_tail_log2_FLOAT":tail/LN2 if math.isfinite(tail) else None,
             "tail_budget_met_FLOAT":met,"peeled_ordered_profiles":k,
             "ordered_profiles":len(first),"remaining_count_log2":math.log(remaining)/LN2 if remaining else None,
             "residual_maximizing_full":([int(first[order[k]]),int(second[order[k]])] if remaining else None)}
        largest.append(row); largest=sorted(largest,key=lambda x:x["log2_bound_FLOAT"],reverse=True)[:6]
        coarse_largest.append({"runs":c,"log2_run_envelope_FLOAT":coarse_run/LN2,
            "maximizing_full":[int(first[order[0]]),int(second[order[0]])]})
        coarse_largest=sorted(coarse_largest,key=lambda x:x["log2_run_envelope_FLOAT"],reverse=True)[:6]
        if c%max(32,r//16)==0 or c==r:
            progress("run_prefix_complete",through_runs=c,exact_queries=queries,
                     ordered_profiles=ordered,query_stop_reason=stop_reason)
    assert source_total==comb(8192,r)
    assert ordered==(r//2+1)**2-(r%2==0)
    assert total<=coarse+1e-7
    if unmet==0: assert all_tails<=-TAIL_BITS*LN2+1e-7
    return {"weight":r,"complete_source_weight":True,"ordered_profiles":ordered,
        "log2_peeling_bound_FLOAT":total/LN2,"log2_unpeeled_run_envelope_FLOAT":coarse/LN2,
        "post_source_tail_log2_FLOAT":all_tails/LN2 if math.isfinite(all_tails) else None,
        "tail_budget_bits":TAIL_BITS,"per_run_tail_log2":target/LN2,
        "tail_budget_met_for_all_runs_FLOAT":unmet==0,"runs_missing_tail_budget":unmet,
        "completed_unique_exact_queries":queries,"exact_query_attempts":attempts,
        "symmetric_count_cache_hits":hits,
        "peeled_ordered_profiles":peeled,"count_seconds":count_elapsed,
        "query_stop_reason":stop_reason,"largest_run_bounds":largest,
        "top_exact_count_weighted_profiles":known_profiles,
        "largest_unpeeled_run_envelopes":coarse_largest,
        "exact_count_checks":{"run_recurrence_sum_equals_binomial":True,
            "remaining_masses_nonnegative":True,"all_ordered_chains_retained":True,
            "both_terminal_orientations":True,"updates":32767},
        "powered_source_log2_FLOAT":log_source/LN2}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight",type=int,default=512)
    parser.add_argument("--modulus",type=int,choices=(128,256),default=128)
    parser.add_argument("--centers")
    parser.add_argument("--z",default="17447/25000")
    parser.add_argument("--alpha",default="119/1280")
    parser.add_argument("--curve-full",default="1562175/990029")
    parser.add_argument("--curve-single",default="1140210/880327")
    parser.add_argument("--seconds",type=float,default=180)
    parser.add_argument("--count-seconds",type=float,default=30)
    parser.add_argument("--max-queries",type=int,default=10000)
    parser.add_argument("--ordinary-fallback",action="store_true")
    args=parser.parse_args(); r=args.weight
    try:
        z,alpha,cf,cs=map(F,(args.z,args.alpha,args.curve_full,args.curve_single))
        centers=default_centers(r) if args.centers is None else tuple(sorted(map(int,args.centers.split(","))))
    except (ValueError,ZeroDivisionError): parser.error("invalid rational or center list")
    if not 7<=r<=4096: parser.error("weight must be in7..4096")
    if not 1<=len(centers)<=16 or len(set(centers))!=len(centers) or not 0<=centers[0]<=centers[-1]<=r//2:
        parser.error("choose1..16 distinct centers in0..floor(r/2)")
    if not 0<z<1 or not 0<alpha<=1 or not all(F(1,256)<=x<=256 for x in (cf,cs)):
        parser.error("require0<z<1,0<alpha<=1,curve multipliers in[1/256,256]")
    if max(abs(n).bit_length() for x in (z,alpha,cf,cs) for n in (x.numerator,x.denominator))>64:
        parser.error("rational parameters limited to64bits")
    if not 30<=args.seconds<=300 or not 0<=args.count_seconds<=min(120,args.seconds-25) or not 0<=args.max_queries<=100000:
        parser.error("seconds30..300,count-seconds0..min(120,seconds-25),queries0..100000")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    neighbor=load("cheap_lane_neighbor_screen.py",NEIGHBOR_PIN,"broad_neighbor")
    peel=load("cheap_lane_profile_peeling.py",PEEL_PIN,"broad_peeling_tests")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    calibration=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"broad_calibration")
    gradient=calibration.load_gradient(); anchor=gradient.load_anchor()
    peak_memory=anchor.peak_memory_reader()

    def guard():
        if time.monotonic()>deadline: raise TimeoutError("overall screen deadline; no complete bound or automatic retry")
        if peak_memory()>256*(1<<20): raise MemoryError("256MiB screen budget exceeded")

    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)

    import numpy as np
    from flint import ctx
    ctx.prec=256
    counts=screen.load("cheap_lane_fast_profile_counts.py",screen.COUNT_PIN,"broad_counts")
    counts.DEADLINE=deadline
    tests=exact_tests(np,screen,peel)
    tests["exact_count_small_checks"]=counts.checks()
    tests["exact_FFT_toy_coefficients"]=fft.exact_toy(np,guard)
    tests["exact_scale_checks"]=cells.exact_scale_tests()
    adapter=anchor.load_adapter(); core=adapter.load_verifier()
    core.Z,core.ALPHA=z,alpha

    class BroadVerifier(core.Verifier):
        def guard(self): guard()

    verifier=BroadVerifier(2,deadline)
    caps=verifier.local_caps(list(combinations_with_replacement(range(17),4)))
    exponent=verifier.exact(alpha)
    powered={occ:tuple(verifier.positive(x)**exponent for x in row) for occ,row in caps.items()}
    del caps
    for row in powered.values():
        for value in row:
            approximate=float(value.mid())
            if not math.isfinite(approximate) or approximate<0 or (value>0 and approximate==0):
                raise ArithmeticError("positive powered cap is not representable in FLOAT")
    coefficients=calibration.float_coefficients(powered,np,guard); del powered
    regions,bank,stats,fallbacks=marginals(coefficients,r,centers,cf,cs,args.modulus,
        args.ordinary_fallback,cells,fft,screen,np,guard,progress)
    del coefficients
    h,c_h,eta=anchor.source_factor(r)
    source_log=exponent*(verifier.exact(eta).log()+r*verifier.exact(core.BETA).log()
                         -core.LIMIT*verifier.exact(z).log())
    log_source=float(source_log.mid())
    report=join(r,regions,counts,screen,np,log_source,args.max_queries,args.count_seconds,
                deadline,guard,progress)
    replay=(r==512 and args.modulus==128 and centers==REFERENCE_CENTERS and not fallbacks
        and z==F(17447,25000) and alpha==F(119,1280) and cf==F(1562175,990029)
        and cs==F(1140210,880327) and report["tail_budget_met_for_all_runs_FLOAT"])
    if replay:
        assert abs(report["log2_peeling_bound_FLOAT"]-(-59.24900640651474))<2e-4
    guard()
    print(json.dumps({"scope":"complete upper-envelope FLOAT screen for this source weight only",
        "directed":False,"negative_is_not_a_certificate":True,"whole_code_certificate":False,
        "other_weights_covered":False,"fresh_marker_or_z_calibration":False,
        "parameters":{"weight":r,"modulus":args.modulus,"centers":centers,
            "z":str(z),"alpha":str(alpha),"curve_full":str(cf),"curve_single":str(cs),
            "marker_reference_weight":r,"boundary_zero_marker":str(ZERO_MARKER),
            "boundary_denominator_floor":1,"seconds":args.seconds,
            "count_seconds":args.count_seconds,"max_queries":args.max_queries},
        "joint_missing":{"h":h,"c_h_exact":str(c_h),"eta_exact":str(eta)},
        "result":report,"bank":bank,"FFT_controls":stats,
        "ordinary_fallback":{"explicitly_enabled":args.ordinary_fallback,
            "marginals_used":len(fallbacks),"first16":fallbacks[:16],
            "all_fallback_full_indices":[x["full"] for x in fallbacks],
            "method":"ordinary positive Cauchy numerator; no residue or coefficient clipping"},
        "r512_peeling_control_replayed":replay,"tests":tests,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak_memory(),"generated_bytes":0,
        "neighbor_pin":NEIGHBOR_PIN,"peeling_tests_pin":PEEL_PIN,
        "FFT_pin":cells.FFT_PIN,"screen_pin":fft.SCREEN_PIN,"count_pin":screen.COUNT_PIN,
        "calibration_pin":screen.CALIBRATION_PIN,"core_pin":anchor.CORE_PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))


if __name__=="__main__":
    main()
