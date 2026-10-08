"""FLOAT paired-region screen for the systematic candidate, never a certificate."""
from fractions import Fraction
import argparse,json
from hd_identity_feedthrough_screen import packet_logs
from hd_sumrank_profile_pairs import coefficients,exact_region_screen
from hd_sumrank_certificate_witnesses import PAIRED


def screen(rows,state=16):
    record=min(PAIRED,key=lambda x:abs(rows-(x[0]+x[1])/2))
    z=float(Fraction(record[2]))**(rows/((record[0]+record[1])/2));base=float(Fraction(record[3]))
    alphas=sorted({min(1,base*x) for x in (.75,1,1.25)})
    best=None
    for scale in (.9,1,1.1):
        mark=z**scale
        cs=coefficients(mark,alphas,state=state,cutoff=record[5],packet_source=packet_logs)
        for alpha,coeff in zip(alphas,cs):
            result=exact_region_screen(rows,mark,alpha,coeff)
            result.update(state=state,feedthrough='identity',feedback='ideal-transitive',cutoff=record[5])
            if best is None or result['log2_upper_screen']<best['log2_upper_screen']:best=result
    return best


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows',type=int,nargs='+',default=(16,128,512))
    p.add_argument('--state',type=int,default=16)
    args=p.parse_args()
    for r in args.rows:print(json.dumps(screen(r,args.state)),flush=True)
