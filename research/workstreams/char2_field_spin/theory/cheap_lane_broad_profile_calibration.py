"""Selected-profile FLOAT marker differential; no complete-weight claim.

Default r1024 profiles are(384,500) and(384,384), at fixed rational
z304397809/625000000 and alpha119/1280. Cache the frozen powered local
coefficients once. For each actual ordered profile compare three witnesses:
the inherited marker-bank centers, its own-f curve, and at most12 guided
ordinary-Cauchy marker proposals starting from the own-f curve.

For unequal profiles use four log markers and the finite numerator
 H=e0(P0^128 P1^128)^127 P0^128 P1^127 T1.
The denominator is(d0*d1)^128, d_i=M4096(f_i,s_i)u_i^f_i v_i^s_i.
Positive mass/first/second moments give the four-coordinate gradient and
Hessian. Their targets are128(f0,s0,f1,s1). Equal profiles enforce shared
markers and use P^32767 T, with two targets256(f,s). These are distinct
optimizations; no averaging of unequal profiles is used.

Newton steps are bounded proposals only. Round the chosen markers to
positive rationals, then compare ordinary Cauchy and residue128/256 at
the SAME exact source count, z,alpha,caps and finite ordered terminal.
Residues use the frozen FFT noise controls. An unresolved projection is
reported, not clipped or replaced by an ordinary bound. Ordinary-Cauchy
stationarity is not claimed for a residue projection.

All reported bounds are FLOAT, not certificates. There is no fresh code,
refresh law, source-density substitution, complete-profile envelope, or
automatic parameter sweep.180s/256MiB/20GiB; stdout only, no generated data.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations_with_replacement, product
import json
import math
from math import comb
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
BROAD_PIN = "df955f7f8458b0e22314c82311b2b713c2142cd9c24b1d188441edf1a208c33f"
MAX_PROPOSALS, SECONDS, MEMORY_LIMIT = 12, 180, 256*(1 << 20)
ZERO_MARKER, LN2 = F(1, 1 << 30), math.log(2)
MODULI = (128, 256)


def load_broad():
    path = Path(__file__).with_name("cheap_lane_broad_peeling_screen.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != BROAD_PIN:
        raise RuntimeError("frozen broad-screen hash mismatch")
    spec = importlib.util.spec_from_file_location("selected_broad_screen", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def logq(value):
    return math.log(value.numerator)-math.log(value.denominator)


class Jets:
    """Positive count moments; only the final covariance subtracts means."""
    def __init__(self, axes, np):
        self.axes, self.np = axes, np
        self.pairs = [(i,j) for i in range(axes) for j in range(i,axes)]
        self.channels = 1+axes+len(self.pairs)

    def normalize(self, value, scale=0.0):
        largest = float(self.np.max(value[0]))
        if not largest>0 or not self.np.isfinite(value).all() or self.np.any(value<0):
            raise ArithmeticError("invalid positive moment matrix")
        return value/largest, scale+math.log(largest)

    def multiply(self, left, right):
        a,la=left; b,lb=right
        out=self.np.empty((self.channels,a.shape[1],b.shape[2]))
        out[0]=a[0]@b[0]
        for i in range(self.axes):
            out[1+i]=a[1+i]@b[0]+a[0]@b[1+i]
        for k,(i,j) in enumerate(self.pairs,1+self.axes):
            out[k]=(a[k]@b[0]+a[0]@b[k]
                    +a[1+i]@b[1+j]+a[1+j]@b[1+i])
        return self.normalize(out,la+lb)

    def power(self, value, exponent):
        identity=self.np.zeros((self.channels,2,2)); identity[0]=self.np.eye(2)
        answer=(identity,0.0)
        while exponent:
            if exponent&1: answer=self.multiply(answer,value)
            exponent//=2
            if exponent: value=self.multiply(value,value)
        return answer

    def moments(self, counts):
        np=self.np
        return np.asarray([np.ones(counts.shape[1])]+[counts[i] for i in range(self.axes)]
                          +[counts[i]*counts[j] for i,j in self.pairs])


def exact_tests(np, screen, calibration):
    """Noncommuting four-packet paths versus exact rational count moments."""
    p0=(((0,0,F(1,3)),(1,0,F(2,5))), ((0,1,F(3,7)),),
        ((1,0,F(5,11)),), ((0,0,F(7,13)),(0,1,F(11,17))))
    p1=(((0,1,F(2,7)),), ((1,0,F(3,11)),(0,0,F(5,13))),
        ((0,0,F(7,17)),(1,0,F(11,19))), ((0,1,F(13,23)),))
    terminal=(((1,0,F(17,29)),), ((0,0,F(19,31)),(0,1,F(23,37))))
    checks=0
    for axes in (2,4):
        jets=Jets(axes,np)
        marker=(F(2,7),F(3,11)) if axes==2 else (F(2,7),F(3,11),F(5,13),F(7,17))

        def lift(terms,side):
            offset=0 if axes==2 else 2*side
            result=[F(0)]*jets.channels
            for f,s,coefficient in terms:
                count=[0]*axes; count[offset:offset+2]=[f,s]
                mass=coefficient*marker[offset]**f*marker[offset+1]**s
                for k,weight in enumerate([1]+count+[count[i]*count[j] for i,j in jets.pairs]):
                    result[k]+=mass*weight
            return result

        a=np.asarray([list(map(float,lift(x,0))) for x in p0]).T.reshape(jets.channels,2,2)
        b=np.asarray([list(map(float,lift(x,1))) for x in p1]).T.reshape(jets.channels,2,2)
        t=np.asarray([list(map(float,lift(x,1))) for x in terminal]).T.reshape(jets.channels,2,1)
        aj,bj,tj=jets.normalize(a),jets.normalize(b),jets.normalize(t)
        actual,scale=jets.multiply(jets.multiply(jets.multiply(aj,bj),aj),tj)
        for initial in (0,1):
            truth=[F(0)]*jets.channels
            for states in product((0,1),repeat=3):
                route=(initial,)+states
                choices=(p0[2*route[0]+route[1]],p1[2*route[1]+route[2]],
                         p0[2*route[2]+route[3]],terminal[route[3]])
                for selected in product(*choices):
                    count=[0]*axes; coefficient=F(1)
                    for side,(f,s,weight) in zip((0,1,0,1),selected):
                        offset=0 if axes==2 else 2*side
                        count[offset]+=f; count[offset+1]+=s; coefficient*=weight
                    mass=coefficient
                    for x,k in zip(marker,count): mass*=x**k
                    for k,weight in enumerate([1]+count+[count[i]*count[j] for i,j in jets.pairs]):
                        truth[k]+=weight*mass
            for k,value in enumerate(truth):
                got=actual[k,initial,0]*math.exp(scale)
                assert abs(got-float(value))<=3e-12*max(1,float(value))
                checks+=1
        # The mass channel also agrees with the frozen ordinary paired chain.
        ra=screen.matrix_power(screen.normalized_matrix(tuple(a[0].ravel())),1)
        rb=screen.matrix_power(screen.normalized_matrix(tuple(b[0].ravel())),1)
        tv=screen.normalized_matrix(tuple(t[0].ravel()))
        plain=screen.matrix_vector(screen.matrix_product(screen.matrix_product(ra,rb),ra),tv)
        assert abs(screen.first_log(plain)-(scale+math.log(actual[0,0,0])))<2e-12
        checks+=1
        if axes==2:
            old=calibration.normalized_product(calibration.normalized_product(
                calibration.normalized_product(aj,bj,np),aj,np),tj,np)
            assert abs(old[1]-scale)<1e-12
            assert np.allclose(old[0],actual,rtol=2e-12,atol=1e-14)
            checks+=1
    return checks


class OrdinaryObjective:
    def __init__(self, coefficients, profiles, constant, np, guard):
        self.np,self.guard=np,guard
        self.profiles=profiles; self.equal=profiles[0]==profiles[1]
        self.jets=Jets(2 if self.equal else 4,np)
        self.indices=[(f,s) for f in range(33) for s in range(33-f)]
        self.fd=np.asarray([f for f,s in self.indices],dtype=float)
        self.sd=np.asarray([s for f,s in self.indices],dtype=float)
        self.local=np.asarray([coefficients[f,s] for f,s in self.indices])
        self.constant=constant
        factor=256 if self.equal else 128
        self.targets=factor*np.asarray(profiles[0] if self.equal else profiles[0]+profiles[1],dtype=float)
        self.features=[]
        for side in range(1 if self.equal else 2):
            counts=np.zeros((self.jets.axes,len(self.indices)))
            counts[2*side]=self.fd; counts[2*side+1]=self.sd
            self.features.append(self.jets.moments(counts))

    def packet(self,x,side):
        np=self.np; offset=2*side
        monomial=self.fd*x[offset]+self.sd*x[offset+1]
        shift=float(np.max(monomial))
        weighted=self.features[side]*np.exp(monomial-shift)[None,:]
        values=weighted@self.local
        p=values[:,:4].reshape(self.jets.channels,2,2)
        t=values[:,4:].reshape(self.jets.channels,2,1)
        return self.jets.normalize(p,shift),self.jets.normalize(t,shift)

    def evaluate(self,x):
        self.guard(); np,jets=self.np,self.jets
        p0,t0=self.packet(x,0)
        if self.equal:
            result,scale=jets.multiply(jets.power(p0,32767),t0)
        else:
            p1,t1=self.packet(x,1)
            r0,r1=jets.power(p0,128),jets.power(p1,128)
            v1=jets.multiply(jets.power(p1,127),t1)
            result,scale=jets.multiply(jets.multiply(jets.power(jets.multiply(r0,r1),127),r0),v1)
        normalized=result[:,0,0]/result[0,0,0]
        means=normalized[1:1+jets.axes]
        covariance=np.empty((jets.axes,jets.axes))
        for k,(i,j) in enumerate(jets.pairs,1+jets.axes):
            covariance[i,j]=covariance[j,i]=normalized[k]-means[i]*means[j]
        value=scale+math.log(result[0,0,0])-float(self.targets@x)+self.constant
        return {"x":x.copy(),"objective":value,"means":means,
                "gradient":means-self.targets,"covariance":covariance}


def calibration(objective,start,np,guard):
    history=[]; best=None
    active=objective.targets>0

    def record(x):
        nonlocal best
        assert len(history)<MAX_PROPOSALS
        point=objective.evaluate(x)
        denominator=256 if objective.equal else 128
        history.append({"evaluation":len(history)+1,"log2_ordinary_FLOAT":point["objective"]/LN2,
            "markers_FLOAT":[math.exp(float(a)) for a in x],
            "tilted_region_means_FLOAT":[float(a/denominator) for a in point["means"]]})
        if best is None or point["objective"]<best["objective"]: best=point
        return point

    current=record(np.asarray(start,dtype=float))
    status="proposal budget exhausted"
    while len(history)<MAX_PROPOSALS:
        guard()
        if np.max(np.abs(current["gradient"][active]))/(256 if objective.equal else 128)<1e-5:
            status="ordinary-Cauchy FLOAT first-order tolerance reached"; break
        gradient=current["gradient"][active]
        covariance=current["covariance"][np.ix_(active,active)]
        if np.isfinite(covariance).all() and np.linalg.eigvalsh(covariance)[0]>1e-8:
            partial=-np.linalg.solve(covariance,gradient)
        else: partial=-gradient/np.maximum(objective.targets[active],1)
        if float(partial@gradient)>=0: partial=-gradient/np.maximum(objective.targets[active],1)
        partial/=max(1,float(np.max(np.abs(partial)))/2)
        direction=np.zeros_like(start); direction[active]=partial
        accepted=False; step=1.0
        while len(history)<MAX_PROPOSALS:
            trial=current["x"].copy()
            trial[active]=np.clip(trial[active]+step*direction[active],-20,4)
            displacement=trial-current["x"]
            if np.max(np.abs(displacement))<1e-12:
                status="proposal box or negligible step reached"; break
            candidate=record(trial)
            if candidate["objective"]<=current["objective"]+1e-4*float(current["gradient"]@displacement):
                current=candidate; accepted=True; break
            step*=0.5
        if not accepted: break
    return best,history,status


class Replay:
    """Cache selected2x2/2x1 values only; never retain past FFT tables."""
    def __init__(self,coefficients,r,cells,fft,screen,np,guard):
        self.r,self.cells,self.fft,self.screen,self.np,self.guard=r,cells,fft,screen,np,guard
        indices=[(f,s) for f in range(33) for s in range(33-f)]
        self.fd=np.asarray([f for f,s in indices]); self.sd=np.asarray([s for f,s in indices])
        local=np.asarray([coefficients[f,s] for f,s in indices])
        self.positive=local>0; self.loglocal=np.full(local.shape,-math.inf)
        self.loglocal[self.positive]=np.log(local[self.positive]); self.cache={}

    def marginal(self,full,markers,modulus):
        key=(full,markers,modulus)
        if key in self.cache: return self.cache[key]
        self.guard(); np,screen,fft=self.np,self.screen,self.fft
        u,v=markers; lu,lv=logq(u),logq(v); single=self.r-2*full
        logs=self.loglocal+(self.fd*lu+self.sd*lv)[:,None]
        ordinary_logs=[]
        for j in range(6):
            largest=float(np.max(logs[:,j]))
            ordinary_logs.append(largest+math.log(float(np.sum(np.exp(logs[:,j]-largest)))))

        def from_logs(values):
            scale=max(values); entries=tuple(math.exp(a-scale) for a in values)
            if not all(a>0 for a in entries): raise ArithmeticError("ordinary channel underflow")
            return screen.normalized_matrix(entries,scale)

        p0,t0=from_logs(ordinary_logs[:4]),from_logs(ordinary_logs[4:])
        ordinary=(screen.matrix_power(p0,128),screen.matrix_vector(screen.matrix_power(p0,127),t0))
        denominator=math.log(comb(4096,full)*comb(4096-full,single))+full*lu+single*lv
        ordinary_region=(*ordinary,denominator)
        self.cache[full,markers,0]={"status":"resolved","region":ordinary_region}
        if modulus==0: return self.cache[key]
        try:
            shift=float(np.max(logs)); mass=np.exp(logs-shift)
            if np.any(self.positive & (mass==0)):
                raise ArithmeticError("positive local FFT input underflow")
            degree=np.zeros((modulus,modulus,6))
            np.add.at(degree,(self.fd%modulus,self.sd%modulus),mass)
            spectra=fft.regional_spectra(degree[:,:,:4].reshape(modulus,modulus,2,2),
                degree[:,:,4:].reshape(modulus,modulus,2,1),np,self.guard)
            selected=[]; diagnostics=[]
            for (raw,scale),expected in zip(spectra,ordinary):
                scale+=128*shift
                fft.scaled_check(tuple(float(x.real) for x in raw[0,0].ravel()),scale,expected)
                residue,_,_,diagnostic=fft.recover((raw,scale),np)
                thresholds=self.cells.selection_thresholds(fft,residue,raw,np)
                values,minimum=self.cells.select_entries(residue,(full%modulus,single%modulus),thresholds)
                selected.append(screen.normalized_matrix(values,scale))
                diagnostics.append({**diagnostic,"selected_minimum_over_DC":minimum})
            result={"status":"resolved","region":(*selected,denominator),"FFT_controls":diagnostics}
        except ArithmeticError as error:
            result={"status":"unresolved","reason":str(error)}
        self.cache[key]=result
        return result

    def profile(self,fulls,markers,log_source_count):
        result={}
        for modulus in (0,)+MODULI:
            entries=[self.marginal(f,m,modulus) for f,m in zip(fulls,markers)]
            label="ordinary" if modulus==0 else "residue"+str(modulus)
            if any(a["status"]!="resolved" for a in entries):
                result[label]={"status":"unresolved","reasons":[a.get("reason") for a in entries]}
                continue
            conditional=self.screen.paired_chains(entries[0]["region"],entries[1]["region"])[0]
            result[label]={"status":"resolved","log2_bound_FLOAT":(log_source_count+conditional)/LN2,
                "log2_conditional_chain_FLOAT":conditional/LN2,
                "FFT_controls":[a.get("FFT_controls") for a in entries] if modulus else None}
        return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight",type=int,default=1024)
    parser.add_argument("--profiles",default="384:500,384:384")
    parser.add_argument("--z",type=F,default=F(304397809,625000000))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    args=parser.parse_args(); r=args.weight
    try: profiles=[tuple(map(int,item.split(":"))) for item in args.profiles.split(",")]
    except ValueError: parser.error("profiles must be comma-separated full0:full1 pairs")
    if not 7<=r<=4096 or not 1<=len(profiles)<=2 or len(set(profiles))!=len(profiles):
        parser.error("require r7..4096 and one or two distinct profiles")
    if any(len(pair)!=2 or not all(0<=f<=r//2 for f in pair) for pair in profiles):
        parser.error("each full count must lie in0..floor(r/2)")
    if not 0<args.z<1 or not 0<args.alpha<=1:
        parser.error("require0<z<1,0<alpha<=1")
    if max(x.bit_length() for a in (args.z,args.alpha) for x in (a.numerator,a.denominator))>64:
        parser.error("rational witnesses limited to64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+SECONDS
    broad=load_broad()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"selected_neighbor")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"selected_old_calibration")
    gradient=old.load_gradient(); anchor=gradient.load_anchor(); peak=anchor.peak_memory_reader()

    def guard():
        if time.monotonic()>deadline: raise TimeoutError("180-second calibration deadline; no automatic retry")
        if peak()>MEMORY_LIMIT: raise MemoryError("256MiB calibration budget exceeded")

    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)

    import numpy as np
    from flint import ctx
    ctx.prec=256
    counts=screen.load("cheap_lane_fast_profile_counts.py",screen.COUNT_PIN,"selected_broad_counts")
    counts.DEADLINE=deadline
    tests=exact_tests(np,screen,old)+old.proposal_tests(np)+gradient.exact_tests()
    tests+=fft.exact_toy(np,guard)+cells.exact_scale_tests()
    adapter=anchor.load_adapter(); core=adapter.load_verifier()
    core.Z,core.ALPHA=args.z,args.alpha

    class LocalVerifier(core.Verifier):
        def guard(self): guard()

    verifier=LocalVerifier(2,deadline)
    caps=verifier.local_caps(list(combinations_with_replacement(range(17),4)))
    exponent=verifier.exact(args.alpha)
    powered={occ:tuple(verifier.positive(value)**exponent for value in row) for occ,row in caps.items()}
    del caps
    for row in powered.values():
        for value in row:
            approximate=float(value.mid())
            if not math.isfinite(approximate) or approximate<0 or (value>0 and approximate==0):
                raise ArithmeticError("positive powered local cap not representable in FLOAT")
    coefficients=old.float_coefficients(powered,np,guard); del powered
    replay=Replay(coefficients,r,cells,fft,screen,np,guard)
    h,c_h,eta=anchor.source_factor(r)
    source_log=exponent*(verifier.exact(eta).log()+r*verifier.exact(core.BETA).log()
                         -core.LIMIT*verifier.exact(args.z).log())
    logsource=float(source_log.mid())
    centers=broad.default_centers(r); cf,cs=F(1562175,990029),F(1140210,880327)
    reports=[]
    progress("cached_local_caps_ready",tests=tests)
    for fulls in profiles:
        guard(); number=counts.profile_count(4096,r,*fulls)
        if not number:
            progress("empty_actual_profile",fulls=fulls); continue
        targets=tuple((f,r-2*f) for f in fulls)
        assignments=tuple(min(centers,key=lambda c:(abs(c-f),c)) for f in fulls)
        inherited=tuple(broad.marker_curve(r,c,cf,cs) for c in assignments)
        own=tuple(broad.marker_curve(r,f,cf,cs) for f in fulls)
        log_source_count=logsource+math.log(number)
        constant=log_source_count-128*sum(math.log(comb(4096,f)*comb(4096-f,s)) for f,s in targets)
        objective=OrdinaryObjective(coefficients,targets,constant,np,guard)

        def vector(markers):
            flat=markers[0] if objective.equal else markers[0]+markers[1]
            return np.asarray([logq(x) for x in flat])

        initial=objective.evaluate(vector(inherited))
        best,history,status=calibration(objective,vector(own),np,guard)
        chosen=tuple(max(ZERO_MARKER,F(math.exp(float(x))).limit_denominator(1<<20)) for x in best["x"])
        if objective.equal: fitted=(chosen,chosen)
        else: fitted=(chosen[:2],chosen[2:])
        rounded=objective.evaluate(vector(fitted))
        progress("marker_proposals_complete",fulls=fulls,axes=objective.jets.axes,
                 proposals=len(history),initial_ordinary_log2=initial["objective"]/LN2,
                 rounded_fitted_ordinary_log2=rounded["objective"]/LN2)
        comparison={}
        for label,markers in (("inherited_bank",inherited),("own_f_curve",own),("fitted",fitted)):
            values=replay.profile(fulls,markers,log_source_count)
            direct=objective.evaluate(vector(markers))
            assert abs(values["ordinary"]["log2_bound_FLOAT"]-direct["objective"]/LN2)<2e-6
            comparison[label]={"markers":[[str(x) for x in pair] for pair in markers],
                "ordinary_log_marker_gradient_FLOAT":[float(x/LN2) for x in direct["gradient"]],
                "ordinary_tilted_region_means_FLOAT":[float(x/(256 if objective.equal else 128)) for x in direct["means"]],
                "bounds":values}
            progress("selected_witness_replayed",fulls=fulls,witness=label,
                     results={name:{k:v for k,v in value.items() if k!="FFT_controls"} for name,value in values.items()})
        report={"fulls":fulls,"profiles":targets,"inherited_centers":assignments,
            "exact_source_count_log2_FLOAT":math.log(number)/LN2,
            "exact_source_count_sha256":hashlib.sha256(str(number).encode()).hexdigest(),
            "count_policy":"same exact indexed count for every comparison; one ordered profile",
            "calibration":{"axes":objective.jets.axes,"targets":list(map(int,objective.targets)),
                "maximum_proposals":MAX_PROPOSALS,"proposals":len(history),"status":status,
                "history":history,"log_marker_box":[-20,4],
                "zero_target_axes_fixed":[bool(x==0) for x in objective.targets],
                "residue_stationarity_claimed":False},"comparison":comparison}
        reports.append(report)
        progress("selected_profile_complete",**report)
    guard()
    print(json.dumps({"scope":"only the explicitly listed actual ordered profiles; FLOAT differential",
        "directed":False,"whole_code_certificate":False,"complete_weight_bound":False,
        "weight":r,"z":str(args.z),"alpha":str(args.alpha),"moduli":MODULI,
        "profiles_requested":profiles,"profiles_nonempty":len(reports),
        "source_factor":{"h":h,"c_h_exact":str(c_h),"eta_exact":str(eta),
            "powered_source_log2_FLOAT":logsource/LN2},
        "tests":{"exact_and_frozen_moment_checks":tests,"finite_updates":32767,
            "terminal_packets":1,"normalizer_power_per_matching":128,
            "final_ordinary_replay_matches_dual_objective":True},
        "broad_pin":BROAD_PIN,"calibration_pin":screen.CALIBRATION_PIN,
        "core_pin":anchor.CORE_PIN,"count_pin":screen.COUNT_PIN,"FFT_pin":cells.FFT_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
