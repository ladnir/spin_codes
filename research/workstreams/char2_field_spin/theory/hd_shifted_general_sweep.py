"""Small deliberately varied profile sweep; never a coverage certificate."""
import argparse
import json
from hd_shifted_general_screen import proposal,bound

PROFILES = (
    (1,0,0,0), (8,8,0,0), (64,64,1,1), (256,256,16,16),
    (512,512,64,64), (512,512,400,400), (512,480,450,420),
    (768,256,48,48), (896,128,100,100), (1024,0,0,0),
    (1024,1024,256,256), (2048,2048,1024,1024),
)

# Fixed FLOAT witnesses make the discriminating cases replayable without
# rerunning an optimizer. Rigorous receipts live in the certificate helper.
CHECKPOINTS = {
    'singleton': ((1,0,0,0), .9991587655850359, .999999515263104,
        ((-11.99964773652837,-40.,-40.),(-40.,-11.99964773652837,-40.))),
    'mixed': ((512,512,400,400), .4826158925972776, .1026151436509744,
        ((-4.499283928159036,-4.49927835554462,-2.30991334083137),)*2),
    'asymmetric': ((512,480,450,420), .49367570922234283, .09802545440976157,
        ((-5.408942950886019,-6.4540781593556,-2.2088270710712097),
         (-5.44207015077779,-4.8278390848501544,-2.296457984764048))),
    'paired64': ((32,32,32,28), .9540638902942163, .10184711565334122,
        ((-40.,-40.,-6.9886846867721655),
         (-9.544465570910006,-9.544528368948043,-6.3944959528860315))),
    'near256': ((128,128,124,116), .8350122129577269, .08918303443460507,
        ((-10.080831746527448,-9.977415804200383,-4.3144662039236135),
         (-8.118270772195903,-8.290861583210848,-4.4003996747508065))),
    'near512': ((256,256,248,240), .6957946023787829, .08652665164588119,
        ((-8.901459997925475,-8.9838430652593,-3.2616285318132805),
         (-7.7737116499633565,-7.903324695309117,-3.306504793208383))),
}


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first',type=int,default=0)
    parser.add_argument('--last',type=int,default=len(PROFILES))
    parser.add_argument('--iterations',type=int,default=130)
    parser.add_argument('--replay',choices=tuple(CHECKPOINTS),nargs='+')
    parser.add_argument('--grid',type=int,default=128)
    parser.add_argument('--paired',action='store_true',help='Use selective two-step composition in fixed replays')
    args = parser.parse_args()
    if args.replay:
        for name in args.replay:
            profile,z,alpha,odds = CHECKPOINTS[name]
            if args.paired:
                from hd_shifted_good_pair_screen import bound as paired_bound
                result = paired_bound(profile,z,alpha,odds,grid=args.grid,exact_count=True)
            else:
                result = bound(profile,z,alpha,odds,method='alias',
                               grid=args.grid,exact_count=True)
            result['checkpoint'] = name
            print(json.dumps(result),flush=True)
        raise SystemExit(0)
    for profile in PROFILES[args.first:args.last]:
        result = proposal(profile,args.iterations,method='hybrid')
        if result['log2_upper_screen']>0:
            result = proposal(profile,args.iterations,method='hybrid',start=result)
        print(json.dumps(result),flush=True)
        if result['log2_upper_screen']>-80:
            try:
                print(json.dumps(bound(profile,result['mark'],result['alpha'],result['marker_log_odds'],method='alias')),flush=True)
            except ArithmeticError as error:
                print(json.dumps({'profile':profile,'alias_rejected':str(error),'whole_code_certificate':False}),flush=True)
