"""FLOAT proposals for less frequent eight-check feedback.

Keep the certified global RS emission and rank-one shared refresh. Repeat
the packet schedule (4,4,4,8), where each number is the local feedback rank.
Each pair retains the existing occupancy-dependent matrix product before
fractional powers. Compose the two resulting positive matrix polynomials.
This is a valid, potentially loose proposal bound, not a certificate.
"""
import argparse
import json
import numpy as np
from hd_bucket_feedback_screen import paired_coefficients
from hd_systematic_mds_screen import optimize


def compose(left,right):
    """Multiply two positive 2-by-2 matrix polynomials in logarithmic form."""
    peak_left=float(np.max(left));peak_right=float(np.max(right))
    a=np.exp2(left-peak_left);b=np.exp2(right-peak_right)
    out=np.zeros((len(a)+len(b)-1,2,2))
    for r in range(2):
        for s in range(2):
            for h in range(2):out[:,r,s]+=np.convolve(a[:,r,h],b[:,h,s])
    with np.errstate(divide='ignore'):return np.log2(out)+peak_left+peak_right


def coefficients(z,alpha,state,_width,strong_pairs=1):
    if strong_pairs not in (1,3):raise ValueError('Use two or four pairs per period')
    left=paired_coefficients(z,alpha,state,64,cutoff=64,local_rank=4)
    right=paired_coefficients(z,alpha,state,64,cutoff=64,local_rank=4,right_rank=8)
    for _ in range(strong_pairs):left=compose(left,right)
    return left


def screen(rows,iterations=100,regional=False,strong_pairs=1):
    # The generic evaluator's width denotes half of a composition unit.
    # Four actual 64-symbol packets form this 256-symbol unit.
    period=128*(strong_pairs+1)
    result=optimize(rows,16,iterations,regional,width=period//2,
                    pair_source=lambda z,a,s,w:coefficients(z,a,s,w,strong_pairs))
    result.update(construction='global-emission-sparser-eight-check-feedback',
                  packet_width=64,composition_width=period,
                  feedback_schedule=[4,4]+[4,8]*strong_pairs,whole_code_certificate=False)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=int,nargs='+',default=[256,512,2048,3072])
    parser.add_argument('--iterations',type=int,default=100)
    parser.add_argument('--regional',action='store_true')
    parser.add_argument('--strong-pairs',type=int,choices=(1,3),default=1)
    args=parser.parse_args()
    for rows in args.rows:print(json.dumps(screen(rows,args.iterations,args.regional,args.strong_pairs)),flush=True)
