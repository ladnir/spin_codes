"""Compare the local overlap bound with the triangle bound, not a certificate."""
import argparse

from mean import build,zero_refine
from spectral_feedback import build as feedback_build
from full_feedback_refinement import coefficients
from group_moment import maps
from occupancy_memory import C,Z
from flint import arb_mat,ctx


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilts',nargs='+',default=['.068','.072','.076'])
    parser.add_argument('--penalty',default='.9')
    args=parser.parse_args();ctx.prec=192
    overlap=build(10)
    feedback=feedback_build(10)
    spectrum=maps()[2]
    for tilt in args.tilts:
        triangle=coefficients(feedback,spectrum,tilt,args.penalty)
        base=[arb_mat([[1]*11 for _ in range(11)]) for _ in range(33)]
        for j in range(5,11):base[j][C,Z]=triangle[j,C,Z]
        changed=zero_refine(base,overlap,feedback,tilt,args.penalty)
        print('LOCAL ONLY shape-maximized C-to-Z; tilt',tilt,flush=True)
        for j in range(5,11):
            print(j,'triangle',base[j][C,Z],'chord',changed[j][C,Z],
                  'ratio upper',(changed[j][C,Z]/base[j][C,Z]).upper(),flush=True)
    print('No region or support cover was computed; not a distance certificate.',flush=True)


if __name__=='__main__':
    main()
