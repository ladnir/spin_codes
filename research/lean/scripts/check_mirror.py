"""Compare the Python `Fix` mirror against Lean, numerator by numerator.

Run `lake env lean scripts/mirror_cases.lean` and pipe its output here, or let
this script run it.  Any disagreement is a port bug and must be fixed before
the generator is trusted to choose a partition.
"""
import subprocess, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixmirror as F

S = F.SCALE
lg = F.directLog(10)

A1 = (32000000000000000000000000000, 35000000000000000000000000000)
B1 = (190000000000000000000000000000, 195000000000000000000000000000)
W1 = (99000000000000000000000000000, 102000000000000000000000000000)
A2 = (300000000000000000000000000000, 400000000000000000000000000000)
C2 = (400000000000000000000000000000, 500000000000000000000000000000)
A3 = (0, 200000000000000000000000000000)
A4 = (400000000000000000000000000000, 400000000000000000000000000000)
C4 = (200000000000000000000000000000, 210000000000000000000000000000)
U1 = 330000000000000000000000000000
THR = -76800000000000000000000


def fmt(name, v):
    if v is None:
        return f"{name} NONE"
    return f"{name} {v[0]} {v[1]}"


expected = [
    fmt("log2", F.LOG2),
    fmt("flogA_9_10", F.flogA(9, 10, 9)),
    fmt("flogA_5_4", F.flogA(5, 4, 9)),
    fmt("flogA_3_2", F.flogA(3, 2, 9)),
    fmt("flogQ_1_3", F.flogQ(1, 3, 10)),
    fmt("flogQ_big", F.flogQ(4096000000000000000000000000000000, S, 10)),
    fmt("golay1", F.fGolayG(F.ofInt(1))),
    fmt("fxlogx_01_02", F.fxlogxL(lg, (100000000000000000000000000000,
                                       200000000000000000000000000000))),
    fmt("fxlogx_0_02", F.fxlogxL(lg, A3)),
    fmt("fpiEval_A2C2", F.fpiEvalL(lg, A2, C2)),
    fmt("fpiCentered_A2C2", F.fpiCenteredG(lg, A2, C2)),
    fmt("fpiBest_A2C2", F.fpiBestL(lg, A2, C2)),
    fmt("fpiBest_A4C4", F.fpiBestL(lg, A4, C4)),
    fmt("fDpa_A2C2", F.fDpaL(lg, A2, C2)),
    fmt("fDpc_A2C2", F.fDpcL(lg, A2, C2)),
    fmt("fgObj_U1A1", F.fgObjL(lg, F.sc(U1), A1)),
    fmt("fpiBest_A1B1", F.fpiBestL(lg, A1, B1)),
    fmt("fpiBest_B1W1", F.fpiBestL(lg, B1, W1)),
    fmt("fpiClamp_A2C2", F.fpiEvalClampL(lg, A2, C2)),
    fmt("fpiBestClamp_straddle", F.fpiBestClampL(lg,
        (8000000000000000000000000000, 500000000000000000000000000000),
        (0, 500000000000000000000000000000))),
    f"infeas_A1B1W1 {str(F.infeasible(A1, B1, W1)).lower()}",
    f"leafOK_A1B1W1 {str(F.leafOK(lg, THR, U1, A1, B1, W1)).lower()}",
]

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
out = subprocess.run(["lake", "env", "lean", os.path.join("scripts", "mirror_cases.lean")],
                     cwd=root, capture_output=True, text=True, shell=True)
lines = [l.strip() for l in out.stdout.splitlines() if l.strip()]
if out.returncode != 0 and not lines:
    print("lean failed:\n", out.stdout[-3000:], out.stderr[-3000:])
    raise SystemExit(1)

got = {l.split()[0]: l for l in lines}
bad = 0
for e in expected:
    key = e.split()[0]
    g = got.get(key)
    if g is None:
        print(f"MISSING in Lean output: {key}")
        bad += 1
    elif g != e:
        print(f"MISMATCH {key}\n  lean   {g}\n  python {e}")
        bad += 1
print(f"compared {len(expected)} cases, {bad} mismatches")
raise SystemExit(1 if bad else 0)
