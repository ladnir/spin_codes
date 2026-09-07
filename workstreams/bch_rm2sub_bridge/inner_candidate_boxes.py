"""Bounded complete dense-type diagnostics using the incoming four-state engine.

Local wrapper only: no shared producer or frozen receipt is changed.
Replay checks fixed witnesses and exact coverage, not outward arithmetic.
"""
import argparse
import itertools
import math
from pathlib import Path
import time

import numpy as np

import inner_candidate_screen as screen

base = screen.base
engine = screen.refresh


def plain(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k:plain(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [plain(v) for v in value]
    return value


def count_box(lower, upper, total):
    remaining = total-sum(lower)
    dimensions = len(lower)
    if any(a>b for a,b in zip(lower,upper)): return 0
    result = 0
    for mask in range(1 << dimensions):
        shifted = remaining-sum(upper[i]-lower[i]+1 for i in range(dimensions) if mask>>i&1)
        if shifted >= 0:
            result += (-1)**mask.bit_count()*math.comb(shifted+dimensions-1,dimensions-1)
    return result


def check_cover(boxes, length, minimum):
    dimensions = len(boxes[0]['lower'])
    total = 0
    for box in boxes:
        lo,hi = box['lower'],box['upper']
        assert len(lo)==len(hi)==dimensions and all(0<=a<=b<=length for a,b in zip(lo,hi))
        assert hi[0] <= length-minimum
        total += count_box(lo,hi,length)
    for a,b in itertools.combinations(boxes,2):
        lo = [max(x,y) for x,y in zip(a['lower'],b['lower'])]
        hi = [min(x,y) for x,y in zip(a['upper'],b['upper'])]
        assert count_box(lo,hi,length)==0, 'Overlapping integer type boxes'
    expected = math.comb(length+dimensions-1,dimensions-1)-math.comb(minimum+dimensions-2,dimensions-1)
    assert total==expected, 'Incomplete integer type cover'
    return total


def replay(saved):
    for path,digest in saved['source_sha256'].items():
        assert base.sha(base.ROOT/path)==digest, path
    name = saved['configuration']; t,s,spectrum,kernel = screen.load(name)
    length = saved['rows']; result = saved['result']
    covered = check_cover(result['selected_boxes'],length,result['occupation_min'])
    epochs = engine.Epochs(t,s,spectrum,[kernel.get(j,0) for j in range(t+1)])
    logs = []
    for box in result['selected_boxes']:
        witness = box['witness']; lam = math.exp(witness['log_surprisal'])
        corners = engine.typed.vertices(box['lower'],box['upper'],length)
        proposal = np.array(witness['proposal']); probabilities = np.array(witness['probabilities'])
        moment = engine.terminal_logs(engine.epoch_mixture(epochs.at(lam),t,[proposal@probabilities]),256*length//t)[0]
        value = float(max(engine.typed.point_logs(corners,length,256,witness['log_density_costs'],
            proposal,moment,256*length//10,lam)))+engine.typed.lattice_log_count(box['lower'],box['upper'])
        assert abs(value-box['own_log_bound']) < 1e-7
        logs.append(value)
    union = float(np.logaddexp.reduce(logs))
    assert abs(union-result['log_union_upper']) < 1e-7
    return dict(status='DENSE_TYPE_FLOAT_WITNESSES_AND_EXACT_COVER_REPLAYED',
                margin_bits=-union/math.log(2),integer_types=str(covered),boxes=len(logs),
                outward_certificate=False)


def run(name, exponent, minimum, nodes, output):
    if output.exists(): raise FileExistsError('Use a fresh result path')
    t,s,spectrum,kernel = screen.load(name)
    if exponent < 7: raise ValueError('Complete outer rows required')
    length = 1 << (exponent-7)
    counts = screen.caps_module.caps()
    # Four nonzero categories; keep the constant all-one row separate.
    bands = [[w for w in counts if 38<=w<64], [w for w in counts if 64<=w<=192],
             [w for w in counts if 192<w<256], [256]]
    started = time.monotonic()
    model = engine.DenseModel(counts,256,t,s,spectrum,[kernel.get(j,0) for j in range(t+1)],
        length,[-4.,-3.,-2.,-1.,0.,.5,1.],bands=bands)
    result = plain(model.search(minimum=minimum,maximum_nodes=nodes,target_bits=70))
    covered = check_cover(result['selected_boxes'],length,minimum)
    hashes = screen.source_hashes()
    for path in (Path(__file__),Path(engine.typed.__file__),Path(engine.composition.__file__),
                 Path(engine.boxes.__file__),Path(engine.balanced.__file__)):
        hashes[path.relative_to(base.ROOT).as_posix()] = base.sha(path)
    data = dict(status='BCH256_INNER_COMPLETE_DENSE_TYPE_DIAGNOSTIC',configuration=name,
        message_exponent=exponent,rows=length,result=result,integer_types=str(covered),
        margin_bits=-result['log_union_upper']/math.log(2),seconds=time.monotonic()-started,
        shell_caps={str(w):n for w,n in counts.items()},source_sha256=hashes,
        full_distance_proved=False,limitations=['Complete dense type cover, nearest binary64 only.',
        'Occupancies below the stated minimum are not included.',
        'Weak upper bounds do not establish an obstruction or encoder failure.'])
    base.write_new(output,data)
    print({k:v for k,v in data.items() if k in ('status','configuration','margin_bits','seconds')},flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--configuration',choices=(*screen.maps.NAMES,*base.CONFIGS))
    parser.add_argument('--m',type=int,default=20)
    parser.add_argument('--minimum',type=int,default=1024)
    parser.add_argument('--nodes',type=int,default=63)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    if args.verify: print(replay(base.read(args.output)),flush=True)
    else:
        if args.configuration is None: parser.error('--configuration required for discovery')
        run(args.configuration,args.m,args.minimum,args.nodes,args.output)
