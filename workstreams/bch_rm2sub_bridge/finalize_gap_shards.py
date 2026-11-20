"""512-bit replay and exact extraction of a complete interval from partial shards."""
import argparse
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from flint import arb,ctx
import bridge as base
import larger_state_maps as maps
import close_larger_state_gap as gap
import screen_exponential_modes as cap_loader
import occupation_three as groups
import polynomial_regions as poly


def run(directory_name,lower,upper,tag):
    assert Path(directory_name).name == directory_name and tag.isidentifier()
    assert 2 <= lower <= upper <= 8192
    directory = base.HERE/'generated'/directory_name
    files = sorted(directory.glob('shard*.json'));assert files
    output = base.HERE/'generated'/f'larger_range_{tag}_outward.json'
    assert not output.exists()
    best = {};hashes = {}
    for path in files:
        saved = base.read(path)
        for source,digest in saved['local_sha256'].items():
            assert base.sha(base.HERE/source) == digest
            assert source not in hashes or hashes[source] == digest
            hashes[source] = digest
        lo,hi = saved['occupancy_range']
        assert [row['occupation'] for row in saved['rows']] == list(range(lo,hi+1))
        for row in saved['rows']:
            q,power = row['occupation'],row['upper_power_of_two']
            if lower <= q <= upper and power <= -70:best[q] = min(best.get(q,power),power)
    assert sorted(best) == list(range(lower,upper+1))
    sources = files+[Path(__file__)]
    hashes.update({str(p.relative_to(base.HERE)):base.sha(p) for p in sources})
    total = sum((gap.as_fraction(p) for p in best.values()),F(0))
    base.write_new(output,dict(status='LARGER_STATE_OUTWARD_RANGE',configuration='t64_s20',occupancy_range=[lower,upper],
        parameters=dict(message_bits=1<<20,output_bits=1<<21,step_bits=64,state_bits=20,distance_cutoff=base.CUTOFF),
        rows=[dict(occupation=q,upper=base.encode(gap.as_fraction(best[q]))) for q in sorted(best)],
        range_upper=base.encode(total),range_below_2_to_minus_60=total<F(1,1<<60),
        all_occupations_certified=False,local_sha256=hashes))
    ctx.prec=512;t,s,spectrum,kernel=maps.load('t64_s20');caps,_=cap_loader.latest_caps()
    maximum=max(base.read(path)['occupancy_range'][1] for path in files)
    @lru_cache(maxsize=8)
    def region(tilt):
        lam=(arb(tilt)/10).exp()
        return poly.regions(t,s,spectrum,kernel,(-lam).exp(),maximum)
    for index,path in enumerate(files):
        saved=base.read(path);lo,hi=saved['occupancy_range']
        actual=gap.evaluate(region(saved['witness']['tilt']),saved['witness'],lo,hi,groups.BANDS,caps)
        assert len(actual)==len(saved['rows'])
        for (q,bound,_),row in zip(actual,saved['rows']):
            assert row['occupation']==q and bound<=gap.as_fraction(row['upper_power_of_two'])
        print('Extracted gap replay shard',index,'passed',lo,hi,flush=True)
    base.write_new(output.with_name(f'larger_range_{tag}_replay.json'),dict(
        status='LARGER_STATE_512_BIT_RANGE_REPLAY_PASSED',rows_checked=len(best),
        producer_sha256=base.sha(output),verifier_sha256=base.sha(Path(__file__))))
    print('Extracted complete interval replayed',lower,upper,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',required=True)
    parser.add_argument('--lower',type=int,required=True);parser.add_argument('--upper',type=int,required=True)
    parser.add_argument('--tag',required=True);args=parser.parse_args()
    run(args.directory,args.lower,args.upper,args.tag)
