"""Batch-local eligible rows for the no-full-reset frontier, fixed SL2 model.

Preserve the older episode counts. This is an alternative support-growth
exclusion, not an assertion of whole-code or extension-field closure.
"""
import argparse
import time
from fractions import Fraction as F

from diagnose_pair_defects import profile
from diagnose_pair_region_rank import early_joint,joint_cell
from diagnose_pair_full_reset_frontier import localized_cap

NEG=-(1<<50)


def cover(h=7):
    began=time.monotonic(); checked=old_residual=rescued=0; worst=NEG; witness=None
    groups={}; examples=[]
    placement=26*(h+1)+(2*(h+1)-1).bit_length()
    old_bits=placement+9+80
    local_bits=26*(h+2)+(2*(h+2)-1).bit_length()+12+80
    for initial in (0,1):
        for episodes in range(1,h+initial+1):
            for live in range(episodes,4096):
                base=profile(h,initial,episodes,live,460,True)
                checked+=1
                if base is None:continue
                value=base['log2_upper']
                if value+placement>-80:
                    local=NEG
                    for early,used,regions in early_joint(h,initial,episodes):
                        cell=joint_cell(h,initial,episodes,live,early,used,regions,old_bits)
                        if cell is not None:
                            local=max(local,value if 'unresolved' in cell else min(value,cell['upper']))
                        if time.monotonic()-began>45:raise TimeoutError('localized cover 45s limit')
                    value=local
                if value+placement>-80:
                    old_residual+=1
                    # Full incoming support must occur, forcing the terminal
                    # interval to end at Q_(L-1) in this no-full-reset branch.
                    if live>base['partial_cap']:
                        extra=h-episodes+initial
                        cap=localized_cap(F(live+episodes,16),
                                          2*episodes-initial-1,local_bits)
                        state_cap=((1<<(extra+1))-1)*cap['cap']
                        if state_cap<512:
                            rescued+=1;value=NEG
                        else:
                            key=(initial,episodes)
                            item=groups.setdefault(key,dict(count=0,first_T=live,last_T=live,
                                                            worst=NEG,largest_batch_cap=0))
                            item['count']+=1;item['last_T']=live
                            item['worst']=max(item['worst'],value+placement)
                            item['largest_batch_cap']=max(item['largest_batch_cap'],cap['cap'])
                            if len(examples)<5:
                                examples.append(dict(a=initial,c=episodes,T=live,extra=extra,
                                                     batch_cap=cap['cap'],state_cap=state_cap,
                                                     upper=value+placement))
                    else:
                        raise AssertionError('residual without guaranteed full state')
                if value>worst:worst,witness=value,(initial,episodes,live)
                if time.monotonic()-began>45:raise TimeoutError('localized cover 45s limit')
    upper=max(worst+placement,-80)+2
    print(dict(kind='localized_no_full_reset',h=h,checked=checked,
               old_residual=old_residual,rescued=rescued,residual=old_residual-rescued,
               residual_groups=groups,examples=examples,good_route_upper=worst+placement,
               worst=witness,route_exceptions=(-80,-80),family_upper=upper,
               closed=upper<=-40,whole_code_certificate=False,
               seconds=time.monotonic()-began),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--defects',type=int,default=7)
    args=p.parse_args()
    assert 0<=args.defects<=16
    cover(args.defects)
