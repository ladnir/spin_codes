"""Small exact FULL/EMPTY projections of the frozen positive pair majorant.

The input has six exact nonnegative rational channels H00,H01,H10,H11,
J0,J1 at every pair degree (f,s), f+s<=64. Relabel e=64-f-s without
changing channels or temporal order. Positive rectangular truncation in
f,e preserves every retained coefficient under multiplication. Six
updates (H,J)<-(H*H,H*J) give H^64 and H^63 J, hence 32767 updates in
the complete 256-region chain with one output-only terminal.

At actual source weight r the regional target is (f,e=4096-r+f), with
s=r-2f. Divide by C(4096,f)C(4096-f,s), including no marker, orientation
factor, source/count factor or new alpha. "Exact projection" refers to
the selected pair majorant; Arb evaluates an enclosure of that projection.
It does not mean an exact encoder enumerator or a complete-size bound.

API edge_regions(v,degrees,max_full,guard,progress=None) returns the
actual4096 map f->(flat R4,flat V2), plus a compact report. The generic
edge_regions_many accepts at most17 actual weights and returns a map for
each. Only explicitly covered targets are returned. The verifier v must
provide exact(Fraction/int) and positive(Arb); caller owns precision,
deadline, memory, >=20GiB disk preflight and source/dependency pins.

No numerical work, file writes, imports of proof providers or standalone
numerical main occur on import. Parent executes all tests and jobs.
"""

from fractions import Fraction as F
from math import comb

LOCAL_PIN="1ed1d59d9beed36b9a1bd1a2bc65af043efe342e95e2afd5047f90c696c04c47"
PAIR_SLOTS,REGION_SLOTS,DOUBLINGS=64,4096,6
MAX_FULL,MAX_EMPTY,MAX_WEIGHTS=8,16,17


def _limits(full,empty):
    if type(full) is not int or type(empty) is not int:
        raise ValueError("integer rectangular limits required")
    if not 0<=full<=MAX_FULL or not 0<=empty<=MAX_EMPTY:
        raise ValueError("require0<=full<=8 and0<=empty<=16")


def _seed(v,degrees,full,empty,slots,guard):
    """Dense rectangular six-channel table; absent input is never zero-filled."""
    expected=(slots+1)*(slots+2)//2
    if len(degrees)!=expected:
        raise ValueError("complete local triangular table required")
    for (f,s),row in degrees.items():
        if type(f) is not int or type(s) is not int or min(f,s)<0 or f+s>slots:
            raise ValueError("invalid local pair degree")
        if len(row)!=6 or any(type(x) not in (int,F) or x<0 for x in row):
            raise ValueError("six exact nonnegative rational channels required")
    zero=v.exact(0); result=[]
    for f in range(full+1):
        guard()
        for e in range(empty+1):
            result.append(tuple(v.exact(x) for x in degrees[f,slots-f-e])
                          if f+e<=slots else (zero,)*6)
    return result


def _double(v,table,full,empty,guard):
    """One positive MM/MT convolution; explicit fixed two-state hot kernel."""
    stride=empty+1; zero=v.exact(0); result=[]
    for f in range(full+1):
        for e in range(stride):
            guard()
            h00=h01=h10=h11=j0=j1=zero
            for a in range(f+1):
                left=a*stride; right=(f-a)*stride
                for b in range(e+1):
                    x=table[left+b]; y=table[right+e-b]
                    h00+=x[0]*y[0]+x[1]*y[2]
                    h01+=x[0]*y[1]+x[1]*y[3]
                    h10+=x[2]*y[0]+x[3]*y[2]
                    h11+=x[2]*y[1]+x[3]*y[3]
                    j0+=x[0]*y[4]+x[1]*y[5]
                    j1+=x[2]*y[4]+x[3]*y[5]
            result.append(tuple(v.positive(x) for x in (h00,h01,h10,h11,j0,j1)))
    return result


def _power(v,table,full,empty,steps,guard,progress=None):
    for step in range(steps):
        table=_double(v,table,full,empty,guard)
        if progress is not None:
            progress("edge_projection_doubling",doubling=step+1,
                     regular_pairs=1<<(step+1),terminal_regular_pairs=(1<<(step+1))-1,
                     full_limit=full,empty_limit=empty,
                     retained_coefficients_per_channel=(full+1)*(empty+1))
    return table


def edge_regions_many(v,degrees,weights,guard,full_limit=8,empty_limit=16,progress=None):
    """Return only covered actual-r marginals, already own-normalized."""
    _limits(full_limit,empty_limit)
    weights=tuple(weights)
    if not 1<=len(weights)<=MAX_WEIGHTS or len(set(weights))!=len(weights):
        raise ValueError("one to17 distinct actual source weights required")
    if any(type(r) is not int or not 1<=r<8192 for r in weights):
        raise ValueError("nonconstant actual weights1..8191 required")
    targets={}
    for r in weights:
        targets[r]=tuple((f,REGION_SLOTS-r+f,r-2*f) for f in range(full_limit+1)
            if 0<=REGION_SLOTS-r+f<=empty_limit and 0<=r-2*f<=REGION_SLOTS-f)
        if not targets[r]:
            raise ValueError("requested weight has no target in the explicit rectangle")
    guard(); before=getattr(v,"intersections",0)
    table=_seed(v,degrees,full_limit,empty_limit,PAIR_SLOTS,guard)
    table=_power(v,table,full_limit,empty_limit,DOUBLINGS,guard,progress)
    regions={}; coverage={}; normalizers=0
    for r in weights:
        own={}; details=[]
        for f,e,s in targets[r]:
            guard()
            if f+e+s!=REGION_SLOTS or 2*f+s!=r:
                raise ArithmeticError("inconsistent actual-source target")
            denominator=comb(REGION_SLOTS,f)*comb(REGION_SLOTS-f,s)
            if denominator!=comb(REGION_SLOTS,e)*comb(REGION_SLOTS-e,f):
                raise ArithmeticError("full/empty multinomial disagreement")
            divisor=v.exact(denominator)
            row=table[f*(empty_limit+1)+e]
            normalized=tuple(v.positive(x/divisor) for x in row)
            own[f]=(normalized[:4],normalized[4:])
            details.append({"full":f,"empty":e,"single":s})
            normalizers+=1
        regions[r]=own; coverage[r]=details
    guard()
    return regions,{"method":"positive directed FULL/EMPTY rectangular projection",
        "actual_weights":weights,"coverage":coverage,"complete_size_coverage":False,
        "full_limit":full_limit,"empty_limit":empty_limit,
        "retained_coefficients_per_channel":(full_limit+1)*(empty_limit+1),
        "doublings":DOUBLINGS,"regular":"H^64","terminal":"H^63 J",
        "local_capacity":PAIR_SLOTS,"regional_capacity":REGION_SLOTS,
        "local_degree_relabeling":"(f,s)->(f,64-f-s)",
        "target_empty":"4096-actual_weight+full",
        "normalization":"C(4096,f) C(4096-f,s); own integer target",
        "normalizer_checks":normalizers,"marker_applied":False,"new_alpha_applied":False,
        "source_or_count_applied":False,"extra_orientation_factor":False,
        "no_residue_aliases":True,"projection_of_frozen_majorant_not_encoder":True,
        "nonnegative_intersections":getattr(v,"intersections",0)-before,
        "expected_local_source_pin":LOCAL_PIN,"generated_bytes":0}


def edge_regions(v,degrees,max_full,guard,progress=None):
    """Convenience API: actual4096, exact f=e=0..max_full projections."""
    result,report=edge_regions_many(v,degrees,(4096,),guard,
        full_limit=max_full,empty_limit=max_full,progress=progress)
    own=result[4096]
    if tuple(own)!=tuple(range(max_full+1)):
        raise ArithmeticError("incomplete requested4096 endpoint strip")
    return own,report


def exact_tests(guard=lambda:None):
    """Independent rational products, terminal order and actual-r targets."""
    class Rational:
        @staticmethod
        def exact(x): return F(x)
        @staticmethod
        def positive(x):
            if x<0: raise ArithmeticError("negative exact toy coefficient")
            return x
    v=Rational(); checks=0
    def mm(a,b):
        return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                           for j in range(len(b[0]))) for i in range(len(a)))
    def add(*matrices):
        return tuple(tuple(sum(m[i][j] for m in matrices) for j in range(len(matrices[0][0])))
                     for i in range(len(matrices[0])))
    def split(row): return ((row[0],row[1]),(row[2],row[3])),((row[4],),(row[5],))
    def poly(a,b):
        result={}
        for (f,e),left in a.items():
            for (g,h),right in b.items():
                guard(); value=mm(left,right); key=f+g,e+h
                result[key]=add(result[key],value) if key in result else value
        return result
    def flat(matrix): return tuple(x for row in matrix for x in row)
    for slots in (1,2):
        degrees={}
        for f in range(slots+1):
            for s in range(slots-f+1):
                weight=comb(slots,f)*comb(slots-f,s)
                degrees[f,s]=tuple(F(weight*(3+2*f+3*s+j),7+j) for j in range(6))
        # The independent product retains every coefficient, without truncation.
        hp={}; jp={}
        for (f,s),row in degrees.items(): hp[f,slots-f-s],jp[f,slots-f-s]=split(row)
        hh,jj=hp,jp; table=_seed(v,degrees,2,2,slots,guard)
        for step in range(3):
            if step:
                jj=poly(hh,jj); hh=poly(hh,hh)
                table=_double(v,table,2,2,guard)
            for f in range(3):
                for e in range(3):
                    expected=flat(hh.get((f,e),((F(0),F(0)),(F(0),F(0)))))
                    expected+=flat(jj.get((f,e),((F(0),),(F(0),))))
                    assert table[3*f+e]==expected; checks+=6
            # Zero full/empty degree is precisely the frozen top-singleton path.
            a,b=hp[0,0],jp[0,0]
            er=((F(1),F(0)),(F(0),F(1)))
            for _ in range((1<<step)-1): er=mm(er,a)
            assert table[0]==flat(mm(er,a))+flat(mm(er,b)); checks+=6
        if slots==2:
            one=_double(v,_seed(v,degrees,1,1,slots,guard),1,1,guard)
            a,b,c,d=(hp[key] for key in ((0,0),(1,0),(0,1),(1,1)))
            t,u,w,x=(jp[key] for key in ((0,0),(1,0),(0,1),(1,1)))
            expected=flat(add(mm(a,d),mm(d,a),mm(b,c),mm(c,b)))
            expected+=flat(add(mm(a,x),mm(d,t),mm(b,w),mm(c,u)))
            assert one[3]==expected; checks+=6
            assert mm(a,b)!=mm(b,a); checks+=1
            assert mm(a,t)!=((sum(t[k][0]*a[k][0] for k in range(2)),),
                              (sum(t[k][0]*a[k][1] for k in range(2)),)); checks+=1
    # Exercise the production six-doubling wrapper on a tiny four-atom table.
    degrees={(f,s):(F(0),)*6 for f in range(65) for s in range(65-f)}
    for f,e in ((0,0),(1,0),(0,1),(1,1)):
        degrees[f,64-f-e]=tuple(F(2+f+e+j,3+j) for j in range(6))
    own,report=edge_regions(v,degrees,1,guard)
    neighbors,nreport=edge_regions_many(v,degrees,(4095,4096,4097),guard,1,2)
    assert tuple(own)==(0,1) and tuple(neighbors[4097])==(1,); checks+=2
    assert own==neighbors[4096]; checks+=1
    assert report["doublings"]==6 and report["normalizer_checks"]==2; checks+=2
    assert nreport["normalizer_checks"]==5; checks+=1
    assert comb(4096,1)*comb(4095,4094)==4096*4095; checks+=1
    assert nreport["coverage"][4095][0]=={"full":0,"empty":1,"single":4095}; checks+=1
    assert nreport["coverage"][4097][0]=={"full":1,"empty":0,"single":4095}; checks+=1
    for arguments in ((9,1),(1,17),(-1,0)):
        try: _limits(*arguments)
        except ValueError: checks+=1
        else: raise AssertionError("invalid rectangle accepted")
    return {"exact_fraction_checks":checks,"independent_untruncated_products":True,
        "full_empty_mapping":True,"f0_top_singleton_identity":True,
        "f1_same_pair_and_two_ordered_split_placements":True,"noncommuting_terminal":True,
        "six_doublings_and_own_normalizers":True,"neighbor_actual4097_rejects_full0":True}
