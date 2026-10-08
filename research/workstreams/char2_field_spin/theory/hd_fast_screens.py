"""Fixed-witness replay of the first fast/high-distance geometry screen.

All coefficients are uniform nonzero GF(2^16), with fixed GF256 maps and
one independent uniform permutation per outer-coordinate region. This
does not optimize, certify intervals, time native code, or exhibit bad words.
The expected margins detect changes in the existing floating proof helper;
they are not acceptance thresholds for a construction.
"""
from fractions import Fraction as F
import argparse
import json
import math

from high_distance_expansion_screen import bound, profile
from high_distance_window_obstruction import lower_bound


# name, b, k, width, state, chronological retained-output pattern,
# sparse?, exact active count/input weight, fixed natural-log marks, margin.
CASES = (
    ("4K-small64-A64", 64,16,4,4,(2,),True,64,
     (4.917050486910787,2.6432434987824065),-380.5156131033579),
    ("4K-small64-density03", 64,16,4,4,(2,),False,125829,
     (-2.623715287439001,-.15295983596328971,2.5571011153983227,2.6176716924846035),-7783.748505450742),
    ("4K-small128-A64", 128,32,4,4,(2,),True,64,
     (5.611271889772042,2.6451252887047016),38.79758547508949),
    ("4K-small128-density03", 128,32,4,4,(2,),False,125829,
     (-2.6239781201875494,-.1529507763291939,2.554982895472363,2.647828752372469),-1121.589538374168),
    ("2K-state4-A16", 256,128,4,4,(4,),True,16,
     (6.217618013797049,2.0114849606107517),-2515.6614914215847),
    ("2K-state4-density03", 256,128,4,4,(4,),False,62915,
     (-2.513585167041394,-.07889607306903107,2.805601752915663,.27462567642365154),-220206.82306474278),
    ("2K-state8-A16", 256,128,4,8,(4,),True,16,
     (6.444132925444082,2.4164826347157238),-955.1715565355782),
    ("2K-state8-density03", 256,128,4,8,(4,),False,62915,
     (-1.8225666076884834,-.1427611827510447,4.044367472128175,1.3013421820759459),-160191.44388690064),
    ("4K-state2-A16", 256,64,4,2,(2,),True,16,
     (6.3079058631574645,3.1414555891152744),-264.3656441838304),
    ("4K-state2-density03", 256,64,4,2,(2,),False,125829,
     (-2.6307148933876627,-.14186565272376805,3.316721144839583,2.5479201859288354),-23583.597608661992),
    ("4K-state3-A16", 256,64,4,3,(2,),True,16,
     (6.353718467529401,3.0940715496286715),580.4840869423573),
    ("4K-state3-A64", 256,64,4,3,(2,),True,64,
     (6.300209131144886,2.7327987552028477),-98.3767105892588),
    ("4K-state3-density03", 256,64,4,3,(2,),False,125829,
     (-2.6311786519983986,-.1476451619275239,3.0103515692682326,2.597927078265783),-5969.707050757615),
    ("3K-width4-state8-A16", 192,64,4,8,(3,3,2),True,16,
     (6.311359063564223,3.47929681093817),1349.3460440732622),
    ("3K-width4-state8-density03", 192,64,4,8,(3,3,2),False,95109,
     (-2.679521626797551,-.1069214104689734,2.5784790472522805,2.408772269175698),-37379.576357691985),
    ("3K-width8-state8-A16", 192,64,8,8,(5,5,6),True,16,
     (6.311613850170314,3.4798332290242793),1350.5761367060995),
    ("3K-width8-state8-density03", 192,64,8,8,(5,5,6),False,95109,
     (-2.711526434219499,-.10432601378856415,2.420742938471738,2.425425584587043),-35164.815373868354),
    ("3K-outer255-A64", 255,85,8,8,(5,5,6),True,64,
     (6.262262416527541,2.3272573193266783),-2009.5679232528503),
    ("3K-outer255-density03", 255,85,8,8,(5,5,6),False,95472,
     (-2.6503852130108103,-.11558364108098056,2.220075994123187,2.4614879827882166),-30259.434628370956),
)


def replay(case, check=False):
    name,b,k,width,state,pattern,sparse,target,marks,expected = case
    options = dict(bits=16,length=b,dimension=k,pattern=pattern,width=width,
                   state=state,distance=F(2,5),k_log=20)
    g = profile(b,k,pattern,width,state,F(2,5),20)
    value = -bound(marks,target,sparse,**options)
    if check and not math.isclose(value,expected,rel_tol=1e-8,abs_tol=2e-5):
        raise AssertionError((name,value,expected))
    return dict(name=name,**g,bits=16,mode="active" if sparse else "weight",
                target=target,log_parameters=marks,margin_bits=value,
                rigorous=False,whole_code=False,
                scope="one floating fixed-witness bound; a negative margin is not a counterexample")


def window_gate(b,k,width,state,pattern,exclude_rows=0):
    """Exact equal-window necessary gate, also valid with periodic puncture.

    Choose an integral number of complete retention periods in every region.
    The retained output fraction then equals the chosen raw-input fraction.
    A nonzero source in the exit-equation kernel is a bad source even if its
    punctured output is zero: in that case the requested dimension is lost.
    Zero lower bound means inconclusive, not that all windows are safe.
    """
    if not 1 <= k < b or width < 1 or state < 0 or not pattern or any(v < 1 for v in pattern):
        raise ValueError("invalid window geometry")
    unit = math.lcm(64,width*len(pattern))
    rows = ((2**20+k*unit-1)//(k*unit))*unit
    if not 0 <= exclude_rows < rows:
        raise ValueError("invalid excluded row count")
    g = dict(rows=rows)
    period_positions = width*len(pattern)
    window = ((2*g['rows']-1)//5//period_positions)*period_positions
    mean,_ = lower_bound(b,k,g['rows'],window,state)
    mean *= F(rows-exclude_rows,rows)
    capacity = b*state
    maximum = k*window
    deterministic = max(0,b*window-(b-k)*rows-k*exclude_rows)
    probability = (F(1) if deterministic > capacity else
                   max(F(0),(mean-capacity)/(maximum-capacity)) if maximum > capacity else F(0))
    return dict(outer=[b,k],width=width,state=state,pattern=pattern,
                rows=g['rows'],window=window,
                excluded_full_rows=exclude_rows,
                expected_shortened_dimension=float(mean),exit_capacity=capacity,
                deterministic_shortened_dimension_lower_bound=deterministic,
                deterministic_kernel_dimension_lower_bound=max(0,deterministic-capacity),
                failure_probability_lower_bound=str(probability),
                failure_probability_lower_bound_float=float(probability),
                exact=True,scope="one fixed equal-window family only")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check",action="store_true")
    p.add_argument("--name",help="replay only cases whose names contain this text")
    p.add_argument("--windows",action="store_true")
    a = p.parse_args()
    selected = [c for c in CASES if a.name is None or a.name in c[0]]
    if not selected:
        p.error("no matching fixed witnesses")
    for case in selected:
        print(json.dumps(replay(case,a.check)),flush=True)
    if a.windows:
        geometries={(c[1],c[2],c[3],c[4],c[5]) for c in selected}
        # Known failing small-outer control checks that the gate is informative.
        geometries.add((64,32,4,4,(4,)))
        geometries.add((192,128,4,4,(5,5,6)))
        geometries.add((256,192,4,4,(6,)))
        for args in sorted(geometries):
            print(json.dumps(window_gate(*args)),flush=True)
        print(json.dumps(window_gate(256,144,8,8,(9,),18)),flush=True)


if __name__ == "__main__":
    main()
