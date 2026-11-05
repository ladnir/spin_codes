"""FLOAT input-weight diagnostic on the frozen four-packet tile operator.

The actual exact-T outer shadow has input weight at least33r. Introduce
v>=1 before the scalar relaxation. Each ordinary packet cap is multiplied
by v^ell, where ell is its allowed input-coordinate count. Quiet packets
instead use their exact scalar kernel weight enumerator. For one block,

    quiet_ell(v,z) = 256^-ell W_{ker C_I}(v*z).

Every four columns of C are independent. Thus the number a_w of kernel
vectors nonzero on every coordinate of a fixed w-subset is determined by
sum_{j<=w} C(w,j)a_j=256^max(w-4,0). The polynomial is
W_ell(t)=sum_w C(ell,w)a_w t^w. All its coefficients are nonnegative.
For a quiet emission-only terminal, the exact moment is
[(1+255*v*z)/256]^ell. Neither replacement includes an outgoing mask count.

The frozen tile4 procedure averages raw route coefficients before alpha.
The input-weight threshold contributes -33*r*alpha*log2(v) to the final
log2 bound. Both v=1 and1.01 use the same exact quiet replacements.
This is a FLOAT artificial repeated-paired-marginal comparison, not a
support-class bound or certificate. It does not change the construction.

One invocation runs the two fixed marks serially, with a48-second shared
deadline,20GiB free-space preflight, and estimated peak memory below256MiB.
Run python -B. Stdout only; no generated files or frozen helper edits.
"""

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
TILES_SHA256 = "a96bcfc95c927501264f9787692e509b03f64a8f968570bc0ac631f581d8f34e"


def load_tiles():
    path = Path(__file__).with_name("cheap_lane_mask_tiles.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != TILES_SHA256:
        raise RuntimeError("frozen mask-tiles helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_mask_tiles", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def kernel_coefficients():
    full = []
    rows = []
    for ell in range(17):
        value = 256**max(ell-4, 0)-sum(math.comb(ell, j)*full[j] for j in range(ell))
        assert isinstance(value, int) and value >= 0
        full.append(value)
        row = [math.comb(ell, w)*full[w] for w in range(ell+1)]
        assert sum(row) == 256**max(ell-4, 0)
        assert all(isinstance(x, int) and x >= 0 for x in row)
        rows.append(row)
    assert full[0] == 1 and full[1:5] == [0]*4 and full[5] == 255
    return rows


def marked_tables(m, t, degrees, v, z, rows, np):
    # raw_tables orders the153 unordered half counts as(a,a,b,b).
    marker = np.power(v, 2*degrees)
    mm = m*marker[:, None, None]
    tt = t*marker[:, None, None]
    moment = [sum(coefficient*(v*z)**w for w, coefficient in enumerate(row))/256**ell
              for ell, row in enumerate(rows)]
    terminal = (1+255*v*z)/256
    row = 0
    for a in range(17):
        for b in range(a, 17):
            assert degrees[row] == a+b
            mm[row, 0, 0] = moment[a]**2*moment[b]**2
            tt[row, 0, 0] = terminal**(2*a+2*b)
            row += 1
    assert row ==153
    assert np.isfinite(mm).all() and np.isfinite(tt).all()
    return mm, tt


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed diagnostic accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 <<30):
        raise RuntimeError("requires20GiB free on actual source volume")
    start = time.monotonic()
    tiles = load_tiles(); grouped = tiles.load_grouped()
    singleton = grouped.load_singleton()
    paired = singleton.load_block().load_paired(start+48)
    tests = {**tiles.toy_test(grouped, paired), **grouped.toy_test(paired)}
    rows = kernel_coefficients()
    tests["nonnegative_integer_MDS_coefficients"] = sum(len(row) for row in rows)
    tests["MDS_kernel_cardinality_checks"] =17
    args = SimpleNamespace(r=256, z=.85, alpha=.523, activity_mark=1.0)
    filtered, middle = paired.imports()
    m, t, degrees, weights = grouped.raw_tables(args, singleton, paired, filtered)
    np = paired.np
    mm, tt = marked_tables(m, t, degrees, 1.0, args.z, rows, np)
    ratios = mm[:, 0, 0]/m[:, 0, 0]
    terminal_ratios = tt[:, 0, 0]/t[:, 0, 0]
    assert np.all(ratios <=1+2e-12) and np.all(terminal_ratios <=1+2e-12)
    assert np.array_equal(mm[:, 1:], m[:, 1:])
    assert np.array_equal(tt[:, 1:], t[:, 1:])
    tests["v1_max_quiet_regular_ratio"] = float(ratios.max())
    tests["v1_max_quiet_terminal_ratio"] = float(terminal_ratios.max())
    results = []
    for v in (1.0, 1.01):
        paired.guard()
        stage = time.monotonic()
        mm, tt = marked_tables(m, t, degrees, v, args.z, rows, np)
        a, b = tiles.raw_polynomials(mm, tt, degrees, weights, np, middle)
        a, b = tiles.double_pair(a, b, 2, 128, paired)
        a = tiles.jensen(a,128,[math.comb(16,j) for j in range(17)],args.alpha,np)
        b = tiles.jensen(b,128,[1],args.alpha,np)
        result = grouped.evaluate(args,a,b,5,paired,filtered,middle)
        correction = -33*args.r*args.alpha*math.log2(v)
        chain = result["artificial_chain"]
        results.append({"input_mark":v,"regional_log2":result["regional_log2"],
            "raw_nonquiet_log2":chain["nonquiet_log2_bound"],
            "global_input_threshold_log2_correction":correction,
            "adjusted_nonquiet_log2":chain["nonquiet_log2_bound"]+correction,
            "adjusted_quiet_log2":chain["quiet_log2_bound"]+correction,
            "adjusted_restart_floor_log2":result["artificial_restart_floor_log2"]+correction,
            "adjusted_source_log2_prefactor":chain["source_log2_prefactor"]+correction,
            "largest_first_active_column":chain["largest_first_active_columns"][0]["column_zero_based"],
            "seconds":time.monotonic()-stage})
    print(json.dumps({"scope":"FLOAT artificial repeated-marginal input-weight diagnostic",
        "whole_code_certificate":False,"actual_support_class_bound":False,
        "parameters":vars(args),"tile_packets":4,"results":results,
        "adjusted_score_change":results[1]["adjusted_nonquiet_log2"]-results[0]["adjusted_nonquiet_log2"],
        "tests":tests,"estimated_peak_bytes":64*(1 <<20),
        "seconds":time.monotonic()-start,"generated_bytes":0,
        "tiles_source_sha256":TILES_SHA256},indent=2))


if __name__ == "__main__":
    main()
