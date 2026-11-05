"""Bounded selected-family map; NOT a whole-domain certificate or distance test.

Preserves every existing verifier. Uses exact integer endpoint bounds with the
retained bank. Failed bounds do not establish feasible histories or bad words.
Run serially; outputs only compact summaries. Peak arrays are below 32 MiB.
"""
import argparse
from collections import Counter
import hashlib
from pathlib import Path
from functools import lru_cache
import time
import numpy as np
import cover_layered_location_regions as c
from layered_mixed_gate import occupation_markers, occupation_floor, no_empty_point, choose_log2

@lru_cache(maxsize=1)
def source_binomials():
    # Compute the large integer binomials once, not once per proof marker.
    return np.array([c.inverse_binomial(s) for s in range(1,c.M+1)],dtype=np.int64)


def checkpoint():
    root=Path(__file__).parent
    files=sorted(p for p in root.glob('*.py') if not p.name.startswith('diagnose_'))
    h=hashlib.sha256()
    for p in files:
        h.update(p.name.encode()+b'\0'+p.read_bytes())
    return dict(files=len(files),sha256=h.hexdigest())


class PointMap:
    def __init__(self,R,all_sources=False):
        self.R=R; self.P=(R*c.M-c.b)//2
        self.H,self.E,self.u=c.region_env(R,'factorial')
        self.coeff=[c.location_coeff(R,i) for i in range(len(c.BANK))]
        self.all_sources=all_sources
        self.s=np.arange(1 if all_sources else 5901,c.M+1,dtype=np.int64)
        slopes=np.array([slope for _,_,_,slope in self.coeff],dtype=np.int64)
        self.source=slopes[:,None]*self.s[None,:]+(R-1)*source_binomials()[self.s-1][None,:]
        self.occ=occupation_markers(R)

    def necessary(self,x):
        z,r,f,e,B=x; l=self.P-z-r-f; V=self.P-z+f
        if min(z,r,l,f,e,B)<0 or e>f or V==0 or B>(400*c.N//1000+c.b)//2:return False
        # Count zero-coordinate episodes, now including both zeros of 00 pairs.
        # H=P+z-f; its head/tail counts are z+r and z+l. The same temporal
        # min-plus lower bound applies with effective full count f-z.
        return min(z+r,z+l)>=occupation_floor(self.R,f-z,B,self.occ)

    def score(self,x):
        z,r,f,e,B=x; P=self.P; R=self.R; D=P-z+r+f-2*e
        degrees=np.minimum.reduce((np.full(c.M+1,D,dtype=np.int64),
            P-z-e+self.E,self.H+r+f-e,R*self.u))
        # Twice the same half-integral cap used by joint_dimension.
        d2=min(D,2*(R-128)*c.M,2*max(0,int((degrees-128*self.u).max())))
        constants=np.array([const+sum(a*v for a,v in zip(aa,x))+rate*d2//2
                            for const,aa,rate,_ in self.coeff],dtype=np.int64)
        values=self.source+constants[:,None]
        if self.all_sources:
            # Valid source-specific sharpening: d_S <= min(dbar,(R-128)*s).
            correction=np.minimum(0,2*(R-128)*self.s-d2)
            rates=np.array([rate for _,_,rate,_ in self.coeff],dtype=np.int64)
            values+=rates[:,None]*correction[None,:]//2
        per_s=values.min(axis=0)
        # The small-source bound is affine in the same dimension cap.
        initial,internal,terminal,mark=c.plain_logs(); m=28
        const=c.POWER*(1+m*c.b+(m-32)*P+(32-m)*(R-128)*5900+c.M)
        const+=(c.b//4)*(initial+terminal)+(P-c.b//2)//16*internal+(400*c.N//1000)*mark[5]
        aa=[c.POWER*v+vlog for v,vlog in zip((32-m,-32,m-32,32,0),mark[:5])]
        small=-((-(const+sum(a*v for a,v in zip(aa,x))+m*c.POWER*d2//2))//c.POWER)
        if self.all_sources:per_s[self.s<=5900]=np.minimum(per_s[self.s<=5900],small*c.POWER)
        pos=int(per_s.argmax());bank=int(values[:,pos].argmin())
        large=-((-int(per_s[pos]))//c.POWER)
        location=(large if self.all_sources else max(small,large))+13
        # Independent coarse complete-profile bounds, not marker extrapolation.
        # Output supports cost <=2^u, exceptional masks <=2^P; each fixed
        # inverse history costs <=2*q^(d-D-e). Forward: union exceptional
        # positions over <=2^P masks, paying one fresh root per exception.
        u=R*c.M; W=400*c.N//1000
        inverse=u+P+32*(W-D-e)+16*d2+1
        forward=P+32*((R-128)*c.M-e)+1
        bounds={'location':location,'output':inverse,'exception_forward':forward}
        l=P-z-r-f; V=P-z+f
        pattern=choose_log2(P,z)+choose_log2(P-z,r)+choose_log2(P-z-r,f)+choose_log2(f,e)
        positions=u-2*c.b; changes=min(positions,2*B+c.b)
        change_bits=positions if 2*changes>=positions else (changes+1).bit_length()+choose_log2(positions,changes)
        pattern=min(pattern,c.b+change_bits+choose_log2(f,e))
        # All four pair types: retain b, not b/2, in boundary losses.
        # Each one-sided pair has one initially-zero free coordinate. Apart
        # from births and the terminal batch, its next-state zero pays a root.
        roots=max(e,l+r+e-B-c.b//2)
        bounds['support_forward']=pattern+32*((R-128)*c.M-roots)+1
        # Independent zero-output equations lie on temporal paths. Dropping
        # exact edge counts, <=2^V choices and a<=V+b-C are conservative.
        C=max(0,V-W,V+B-c.b-W)
        bounds['support_inverse']=pattern+V+32*(c.b+e-r-C)+16*d2+1
        if z==0:
            old=no_empty_point(R,f,r,e,(self.H,self.E,self.u),births=B,dimension=d2//2)
            bounds['old_no_empty']=old['best']
        # DIFFERENT ensemble: rejected-uniform SL2 pair maps. Nonzero pairs
        # have two random inverse coordinates, with total density loss <2.
        # No exceptional label exists. Project the class onto (z,r,f,B).
        idealD=2*(P-z)
        ideald=min(idealD//2,(R-128)*c.M,max(0,int((np.minimum(idealD,2*self.E)-128*self.u).max())))
        idealpattern=pattern-choose_log2(f,e)
        idealoutput=u+32*(W-idealD+ideald)+1
        idealstate=idealpattern+V+32*(V+c.b-C+ideald-idealD)+1
        idealroots=max(0,idealD-V-c.b)
        idealforward=idealpattern+32*((R-128)*c.M-idealroots)+1
        ideal=min(idealoutput,idealstate,idealforward)
        method=min(bounds,key=bounds.get)
        return dict(bound=bounds[method],location=location,method=method,small=small,large=large,
                    s=int(self.s[pos]),bank=c.BANK[bank][0],dimension_twice=d2,
                    ideal_pair_bound=ideal,ideal_pair_dimension=ideald)


def scan(seconds,all_sources=False):
    start=time.monotonic(); before=checkpoint(); stats=[]; examples=[]
    R_values=(129,144,160,176,192,200,208,224,240,256)
    for R in R_values:
        env=PointMap(R,all_sources); P=env.P; seen=set(); counts=Counter(); failed_by=Counter(); worst=[]
        for znum in (0,1,4,8):
            z=P*znum//16
            for fnum in (0,4,8,12,16):
                f=(P-z)*fnum//16
                for rnum in (0,1,2):
                    r=(P-z-f)*rnum//2
                    for enum in (0,1,2,3,4):
                        e=f*enum//4
                        for B in sorted({0,P//32,P//8,min(P//2,(400*c.N//1000+c.b)//2)}):
                            if time.monotonic()-start>seconds:raise TimeoutError('Incomplete diagnostic screen; no coverage claim')
                            x=(z,r,f,e,B)
                            if x in seen:continue
                            seen.add(x);counts['selected']+=1
                            if not env.necessary(x):counts['excluded']+=1;continue
                            value=env.score(x);bad=value['bound']>=-100
                            counts['failed' if bad else 'passed']+=1
                            if bad:counts['failed_also_ideal' if value['ideal_pair_bound']>=-100 else 'failed_resolved_by_ideal_pair']+=1
                            if bad:
                                failed_by[znum]+=1
                                worst.append((value['bound'],x,value))
                                worst=sorted(worst,reverse=True)[:3]
        item=dict(R=R,counts=dict(counts),failed_by_empty_fraction=[(key,val) for key,val in sorted(failed_by.items())],worst=worst)
        stats.append(item)
        print(dict(R=R,counts=dict(counts),worst=worst[:1]),flush=True)
        examples.extend((R,)+v for v in worst)
    after=checkpoint();assert before==after
    print(dict(summary=stats,all_sources=all_sources,elapsed_seconds=time.monotonic()-start,checkpoint=before,
        diagnostic_only=True,whole_domain_covered=False,whole_code_certificate=False),flush=True)


def slices():
    center=c.hardest()[1:6];out=[]
    for R in (192,200,224,256):
        env=PointMap(R)
        for axis,values in ((0,(0,16384,49136,65536,98304)),
                            (3,(0,184260,294816,368520,442224,589632)),
                            (4,(0,49136,98304,196544))):
            row=[]
            for val in values:
                x=list(center);x[axis]=val;x=tuple(x)
                row.append((val,env.necessary(x),env.score(x)))
            out.append(dict(R=R,axis=axis,values=row))
    print(dict(slices=out,diagnostic_only=True),flush=True)


def tests():
    for R in (192,194,208):
        env=PointMap(R)
        for x in (c.hardest()[1:6],(49137,24568,589632,368520,49136)):
            got=env.score(x);want=c.adaptive_box(R,x,x,step=1,refined='factorial')
            assert got['location']==want['per_class'],(R,x,got,want)
    print(dict(exact_point_replays=6))


def probes():
    from refine_layered_location_weights import refine
    cases=((200,c.hardest()[1:6]),
           (192,(393088,0,196544,98272,24568)),
           (160,(327552,163776,0,0,81888)))
    for R,x in cases:
        env=PointMap(R);old=env.score(x);s=old['s']
        profile=(R,)+x+(old['dimension_twice']//2,0)
        result=refine(profile,s,25,seconds=20,free_m=True,refined='factorial')
        print(dict(R=R,x=x,s=s,old=old,optimized_slice_upper=result['exact_upper'],
                   m=result['m'],seconds=result['seconds'],stop=result['stop'],
                   whole_profile_closed=False),flush=True)


def ideal_examples():
    for R,x in ((192,c.hardest()[1:6]),(200,c.hardest()[1:6]),
                (192,(393088,0,294816,221112,98272)),
                (160,(327552,0,327552,81888,81888)),
                (256,(524160,0,524160,524160,0))):
        env=PointMap(R)
        print(dict(R=R,x=x,score=env.score(x),different_ensemble=True),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seconds',type=int,default=120)
    ap.add_argument('--test',action='store_true');ap.add_argument('--slices',action='store_true')
    ap.add_argument('--probes',action='store_true')
    ap.add_argument('--retuned-map',action='store_true')
    ap.add_argument('--ideal-examples',action='store_true')
    args=ap.parse_args();assert 1<=args.seconds<=180
    if args.test:tests()
    elif args.slices:slices()
    elif args.probes:probes()
    elif args.ideal_examples:ideal_examples()
    elif args.retuned_map:probes();scan(args.seconds,all_sources=True)
    else:scan(args.seconds)
