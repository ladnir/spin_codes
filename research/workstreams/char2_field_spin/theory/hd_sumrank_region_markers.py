"""FLOAT entry-specific regional coefficient witnesses; no certification.

Each matrix entry gets its own positive coefficient marker. This preserves
the exact active count in each region without explicitly convolving a long
matrix polynomial. A selected marker always gives an upper bound; numerical
optimization need not find its minimum. Directed replay is still required.
"""

from fractions import Fraction as F
from math import comb,log2

from scipy.optimize import minimize_scalar

from hd_sumrank_transfer import fractional_packets
from hd_randomized_outer_state_profile import power,total


def screen(rows,mark,alpha,*,state=4,mixer='feedback_refresh',input_mark=F(1),tiles=8):
    width=8*tiles
    packets=fractional_packets(F(mark),F(alpha),mixer=mixer,input_mark=F(input_mark),state=state,tiles=tiles)
    size=len(packets[0])
    choose=log2(comb(8192,rows))
    cache={}

    def at(odds):
        if odds not in cache:
            local=tuple(tuple(total(log2(comb(width,j))+j*odds+packets[j][r][s]
                                    for j in range(width+1)) for s in range(size)) for r in range(size))
            block=power(local,8192//width)
            cache[odds]=tuple(tuple(x-rows*odds-choose for x in row) for row in block)
        return cache[odds]

    entries,markers=[],[]
    for r in range(size):
        row,marks=[],[]
        for s in range(size):
            found=minimize_scalar(lambda x:at(x)[r][s],bounds=(-32,32),method='bounded',options={'xatol':1e-7})
            row.append(float(found.fun))
            marks.append(float(found.x))
        entries.append(tuple(row))
        markers.append(tuple(marks))
    bound=(choose+float(alpha)*(rows*(4096-128*log2(65535))-log2(65535)
                              -838860*log2(mark)-129*rows*log2(input_mark))
           +total(power(tuple(entries),256)[0]))
    return dict(active_rows=rows,mark=str(mark),alpha=str(alpha),state=state,mixer=mixer,
                input_mark=str(input_mark),tiles=tiles,regional_log2_entries=entries,
                marker_log_odds=markers,log2_upper_screen=bound,directed=False,whole_code_certificate=False)


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=int,nargs='+',default=(16,128,512,1024,2048))
    parser.add_argument('--state',type=int,default=4)
    parser.add_argument('--tiles',type=int,default=8)
    parser.add_argument('--mixer',choices=('rank_one','uniform_pre','feedback_refresh'),default='feedback_refresh')
    parser.add_argument('--marks',type=F,nargs='+',default=(F(63,64),F(7,8),F(1,2),F(1,4),F(1,8),F(1,16)))
    parser.add_argument('--alphas',type=F,nargs='+',default=(F(1,32),F(1,16),F(1,12),F(1,8),F(1,4)))
    args=parser.parse_args()
    best={r:None for r in args.rows}
    for z in args.marks:
        for alpha in args.alphas:
            for r in args.rows:
                result=screen(r,z,alpha,state=args.state,mixer=args.mixer,tiles=args.tiles)
                if best[r] is None or result['log2_upper_screen']<best[r]['log2_upper_screen']:
                    best[r]=result
        print(json.dumps(dict(completed_mark=str(z),best=list(best.values()))),flush=True)
