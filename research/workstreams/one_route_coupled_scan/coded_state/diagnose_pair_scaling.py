"""Selected larger-budget tests of the fixed pair proof, not a certificate.

Serial standard-library diagnostic, bounded to 45s, no output files.
"""
import time
import argparse
from diagnose_pair_budgeted_reset import selected
from diagnose_pair_defects import profile, volume_bound
from diagnose_pair_region_rank import rank_cap
from diagnose_pair_episodes import volume_bits


def block_comparison(h,a,c,T,slots,block,cap,dimension):
    """Hypothetical independent uniform superregular block maps, not SL2.

    Cyclic disjoint-bit partitions cover nine state-index bits in ceil(9/k)
    phases for block=2^k. Only supports closed under the current partition
    are clean. The same decreasing-support argument bounds proper states.
    """
    assert block in (2,4,8)
    phase_count=(9+(block.bit_length()-2))//(block.bit_length()-1)
    e=h-c+a
    proper=min(T,phase_count*(h+a)-c)
    full=T-proper
    # No proper activation claim when the existing cap permits512 inputs.
    runs=e if cap<512 else h+1
    bits=26*(h+2)+(2*(h+2)-1).bit_length()+9+80
    regions=min(256,(full+30*runs)//16)
    rank=512*full if regions<=128 else rank_cap(full,runs-1,bits)['rank_lower']
    # Union-bound Schwartz-Zippel rejection for every square minor of A and
    # A+delta I. Uniform matrices before conditioning give uniform Au at u!=0.
    import math
    rejected=2*sum(k*math.comb(block,k)**2 for k in range(1,block+1))
    maps=(1<<21)//block
    # -log2(1-x) <= (3/2)*x/(1-x), valid since log2(e)<3/2.
    density=(3*maps*rejected+2*((1<<32)-rejected)-1)//(2*((1<<32)-rejected))
    # For the actual SL2 control use its already established global factor2,
    # not a claim that its setup is the uniform-matrix comparison distribution.
    if block==2:density=1
    placement=26*(h+1)+(2*(h+1)-1).bit_length()
    W=460*(1<<21)//1000
    inverse=density+volume_bits(slots,W)-32*rank+placement
    random_min=512*full+block*proper
    forward=density+volume_bits(random_min,W)+32*(dimension-random_min)+placement
    growth=0; size=cap
    assert cap>0
    while size<512:
        size=block*size+cap;growth+=1
    if growth>e:
        reset=-(1<<50)
    else:
        window=phase_count*(e-growth+1)
        reset=density+32*dimension+12+(window**512-1).bit_length()-32*512+placement
    upper=min(inverse,forward,reset)
    return dict(block=block,phases=phase_count,proper=proper,full_runs=runs,
                density_bits=density,inverse=inverse,forward=forward,reset=reset,upper=upper)


def run():
    began=time.monotonic()
    # Separate many activations from many internal defects. Fixed T avoids
    # reporting a favorable length search as general history coverage.
    for h,c in ((8,8),(12,12),(16,16),(12,6),(16,8),(16,4)):
        a=1; T=1967
        base=profile(h,a,c,T,460)
        result=selected(h,a,c,T,False)
        bits=26*(h+2)+(2*(h+2)-1).bit_length()+9+80
        # Diagnostic ideal: each nonzero incoming state would supply all512
        # random coordinates. This is not an assumption about the pair map.
        rank=rank_cap(T,c-1,bits)['rank_lower']
        placement=26*(h+1)+(2*(h+1)-1).bit_length()
        ideal=1+volume_bound(base['slots'],460*(1<<21)//1000)-32*rank+placement
        alternatives=[block_comparison(h,a,c,T,base['slots'],b,result['batch_cap'],result['dimension'])
                      for b in (4,8)]
        pair=block_comparison(h,a,c,T,base['slots'],2,result['batch_cap'],result['dimension'])
        print(dict(h=h,c=c,T=T,e=h-c+a,S=result['batch_cap'],
                   dimension=result['dimension'],growth=result['growth_required'],
                   proper=result.get('proper'),full_regions=result.get('full_region_cap'),
                   upper=min(result['upper'],pair['upper']),ideal_dense_upper=ideal,
                   uniform_block_comparisons=alternatives,
                   comparison_only=True,seconds=time.monotonic()-began),flush=True)
        if time.monotonic()-began>45:raise TimeoutError('scaling diagnostic 45s cap')


def test():
    began=time.monotonic()
    # Backward closure of any surviving coordinate over every cyclic window
    # reaches all512 coordinates. This checks the actual bit-partition schedule.
    schedules={2:tuple((i,) for i in range(9)),
               4:((0,1),(2,3),(4,5),(6,7),(8,0)),
               8:((0,1,2),(3,4,5),(6,7,8))}
    cases=0
    for block,phases in schedules.items():
        for start in range(len(phases)):
            for coordinate in range(512):
                support={coordinate}
                for offset in reversed(range(len(phases))):
                    for bit in phases[(start+offset)%len(phases)]:
                        support.update(tuple(x^(1<<bit) for x in support))
                assert len(support)==512
                cases+=1
    # The pair clean-run length eight is attainable as support geometry:
    # the lower half is closed under bits0..7 but not under bit8.
    half=set(range(256))
    assert all({x^(1<<bit) for x in half}==half for bit in range(8))
    assert {x^256 for x in half}!=half
    print(dict(backward_closure_cases=cases,pair_eight_phase_half_support=True),flush=True)
    # Strengthen the old count by retaining the minimum mass of each clean
    # proper run. On width4 the loss is <=2 per clean run and2 per defect.
    from diagnose_pair_defect_envelope import partner
    states={(mask,1,0,0):0 for mask in range(1,16)}; checks=0
    for t in range(7):
        nxt_states={}
        for (mask,c,e,T),mass in states.items():
            live=mask|partner(mask,t%2,4)
            for nxt in range(16):
                cc=c+int(mask==0 and nxt!=0)
                ee=e+int(mask!=0 and (live!=mask or bool(nxt&~live)))
                tt=T+int(mask!=0);mm=mass+live.bit_count()
                if cc>3 or ee>3:continue
                assert mm>=4*tt-2*cc-4*ee
                key=(nxt,cc,ee,tt)
                nxt_states[key]=min(nxt_states.get(key,1<<30),mm)
                checks+=1
        states=nxt_states
    assert 512*1967-3586*12-4096==959976<460*(1<<21)//1000
    print(dict(clean_run_mass_checks=checks,h12_c12_T1967_random_floor=959976,
               target_weight=460*(1<<21)//1000),flush=True)
    # A cheap two-layer butterfly is not a uniform four-coordinate map.
    # If the first active pair collapses to one coordinate, the crossing
    # superregular second layer produces exactly two nonzero coordinates.
    from diagnose_uniform_pairs import field
    from fractions import Fraction
    for q,poly in ((8,11),(16,19)):
        mul=field(q,poly);delta=2
        good=[(a,b,c,d) for a in range(1,q) for b in range(1,q)
              for c in range(1,q) for d in range(1,q)
              if mul[a][d]^mul[b][c]==1 and a!=delta and d!=delta
              and mul[a^delta][d^delta]^mul[b][c]!=0]
        receipts=[]
        for x,y in ((1,0),(1,1)):
            collapsed=sum(((mul[a][x]^mul[b][y])==0)
                          or ((mul[c][x]^mul[d][y])==0) for a,b,c,d in good)
            if y==0:assert collapsed==0
            else:assert collapsed>0
            receipts.append(dict(input=(x,y,0,0),two_nonzero_outputs=str(Fraction(collapsed,len(good)))))
        print(dict(q=q,accepted_pairs=len(good),two_layer_support_test=receipts),flush=True)
        if time.monotonic()-began>20:raise TimeoutError('scaling tests 20s cap')


def window():
    began=time.monotonic()
    for h in (12,16):
        for T in (1936,1952,1967,1984,2000):
            c=h;a=1
            base=profile(h,a,c,T,460)
            result=selected(h,a,c,T,False)
            others=[block_comparison(h,a,c,T,base['slots'],b,result['batch_cap'],
                                     result['dimension']) for b in (4,8)]
            pair=block_comparison(h,a,c,T,base['slots'],2,result['batch_cap'],result['dimension'])
            print(dict(h=h,c=c,T=T,cap=result['batch_cap'],pair=min(result['upper'],pair['upper']),
                       block4=others[0]['upper'],block8=others[1]['upper'],
                       selected_profiles_only=True),flush=True)
            if time.monotonic()-began>45:raise TimeoutError('window diagnostic 45s cap')
    print(dict(seconds=time.monotonic()-began,whole_code_certificate=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--test',action='store_true')
    parser.add_argument('--window',action='store_true')
    args=parser.parse_args()
    if args.test:test()
    elif args.window:window()
    else:run()
