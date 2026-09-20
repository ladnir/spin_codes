"""Dense-tail discovery with a fixed interval of all-one row counts.

For each ordinary-row count d, the moment envelope is uniform over h in
the interval. The count C(8192,d)*2^(8192-d) includes every constant/zero
assignment, so it safely bounds the restricted h interval as well.
"""
import argparse
import math
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize_scalar
import bridge as base
import screen_larger_state_v2 as cache
import occupation_three as groups
import screen_exponential_modes as cap_loader
import scaled_adaptive as scaled
import positive_line_hull as hull
import general_occupancy as general
from screen_constant_row_split import moment,log_power


def box_envelope(region,lo,hi,maximum):
    assert 0 <= lo <= hi <= 8192 and maximum <= 8192-lo
    padded = np.pad(region,((0,hi-lo),(0,0),(0,0)),constant_values=-np.inf)
    windows = np.lib.stride_tricks.sliding_window_view(padded[lo:hi+maximum+1],hi-lo+1,axis=0)
    return windows.max(axis=-1)


def screen(qmin,boxes,tilt,tag):
    assert 1 <= qmin <= 8192 and tag.isidentifier()
    output = base.HERE/'generated'/f'larger_constant_boxes_{tag}_screen.json'
    assert not output.exists()
    region = cache.region_logs('t64_s20',tilt)
    caps,sources = cap_loader.latest_caps();bands = groups.BANDS[:-1]
    correction = base.CUTOFF*math.exp(tilt/10)
    floors = []
    for floor in (512,1024,1536,2048,3072,4096):
        # Every codeword with h>=floor has at most8192-floor ordinary rows.
        # sum_d C(L,d)2^(L-d) a^d <= 3^L a^(L-floor), for a>=1.
        value = log_power(region[floor:].max(axis=0))+correction
        value += 8192*math.log(3)+(8192-floor)*math.log((1<<128)-2)
        floors.append(dict(all_one_rows_minimum=floor,margin_bits_diagnostic=-value/math.log(2)))
    print('Direct floor screens',floors,flush=True)
    rows = []
    for lo,hi in boxes:
        dmin,dmax = max(0,qmin-hi),8192-lo
        assert dmin <= dmax
        envelope = box_envelope(region,lo,hi,dmax)
        anchor = (dmin+dmax)//2
        ps,gammas = [],[]
        for band in bands:
            w = np.array(band)
            v = np.array([math.log(caps[x])-math.log(math.comb(256,x)) for x in band])
            def objective(theta):
                p = 1/(1+math.exp(-theta))
                gamma = float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p)))
                return moment(envelope[:anchor+1],anchor,p)+anchor*gamma
            opt = minimize_scalar(objective,bounds=(-6,14),method='bounded',options={'xatol':1e-6})
            p = 1/(1+math.exp(-float(opt.x)))
            ps.append(p);gammas.append(float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))))
        ps = np.array(ps);roots = np.exp(np.array(gammas)/256)
        left,right = roots*(1-ps),roots*ps;keep = hull.indices(left,right)
        exponents = np.floor(envelope.max(axis=(1,2))/math.log(2)).astype(np.int64)
        mantissas = np.exp(envelope-exponents[:,None,None]*math.log(2))
        values = []
        if dmin == 0:
            values.append((0,log_power(envelope[0])+correction+8192*math.log(2)))
        for d,matrix,exponent in scaled.matrices(mantissas,exponents,left[keep],right[keep],outward=False):
            if d < dmin:continue
            value = general.log_power(matrix)+256*exponent*math.log(2)+correction
            value += math.log(math.comb(8192,d))+(8192-d)*math.log(2)+d*math.log(len(bands))
            values.append((d,value))
        worst = max(values,key=lambda pair:pair[1])
        row = dict(all_one_interval=[lo,hi],ordinary_interval=[dmin,dmax],anchor=anchor,
            worst_ordinary_count=worst[0],minimum_margin_bits_diagnostic=-worst[1]/math.log(2),
            union_margin_bits_diagnostic=-(float(np.logaddexp.reduce([v for d,v in values])))/math.log(2),
            p=[base.encode(F.from_float(float(p))) for p in ps])
        rows.append(row);print('Box',lo,hi,'worst',row['minimum_margin_bits_diagnostic'],'at d',worst[0],flush=True)
    files = [Path(__file__),Path(cache.__file__),Path(groups.__file__),Path(cap_loader.__file__),Path(scaled.__file__),
        Path(hull.__file__),Path(general.__file__),base.HERE/'screen_constant_row_split.py',
        base.HERE/f'generated/larger_v2_t64_s20_region_{tilt}.json']+sources
    base.write_new(output,dict(status='LARGER_CONSTANT_BOX_SCREEN_ONLY',configuration='t64_s20',
        total_occupancy_minimum=qmin,tilt=tilt,rows=rows,floors=floors,bands=bands,
        used_caps={str(w):caps[w] for w in base.WEIGHTS},cap_receipts=[str(p.relative_to(base.HERE)) for p in sources],
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in files}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--qmin',type=int,required=True);parser.add_argument('--boxes',nargs='+',required=True)
    parser.add_argument('--tilt',type=int,default=8);parser.add_argument('--tag',required=True)
    args = parser.parse_args();screen(args.qmin,[tuple(map(int,b.split(':'))) for b in args.boxes],args.tilt,args.tag)
