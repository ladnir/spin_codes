"""Full support covers using weighted CDF caps inside each support interval.

Avoid multiplying the largest count in an interval by its worst binomial
denominator when those two extrema occur at different supports.
"""
import argparse
from fractions import Fraction
from heapq import heappush,heappop
from itertools import product
from math import comb,factorial,log
from collections import Counter
import numpy as np
from scipy.optimize import minimize,minimize_scalar
from scipy.special import gammaln,logsumexp
from flint import arb,arb_mat,arb_poly,ctx

from bch_joint_support import authenticated_caps,support_caps
from basis_lattice import improve_caps
from shortened_bound import dimension_caps
from occupancy_memory import prepare,build,epoch_operators,rounded,C,TERMINAL
from occupancy_model import local_data,placement
from occupancy_adaptive import volume,split,geometry_test
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_power_moment,log_binomial_mass
from occupancy_memory_verify import DENOMINATOR
from group_rank_one_verify import up


def bernoulli_masses(numerators):
    """Outward coefficients of the Bernoulli-count polynomial.

    Equal exact probabilities are grouped before Arb polynomial powering.
    This evaluates the same product as the scalar recurrence, including
    the endpoint probabilities zero and one.
    """
    if any(type(n) is not int or not 0<=n<=DENOMINATOR for n in numerators):
        raise ValueError('integer probability numerators in the fixed denominator required')
    polynomial=arb_poly([1])
    for numerator,count in Counter(numerators).items():
        p=arb(numerator)/DENOMINATOR
        polynomial*=arb_poly([1-p,p])**count
    return [polynomial[i] for i in range(len(numerators)+1)]


def increments(counts,lo,hi):
    return [counts[lo]]+[counts[u]-counts[u-1] for u in range(lo+1,hi+1)]


def fold_log(counts,lo,hi,p):
    weights=np.array([-log_binomial_mass(256,u,p) for u in range(lo,hi+1)])
    majorant=np.maximum.accumulate(weights[::-1])[::-1]
    return float(logsumexp([log(c)+w for c,w in zip(increments(counts,lo,hi),majorant) if c]))


def fold_function(counts,lo,hi):
    """Precompute a floating proposal objective; outward replay is unchanged."""
    supports=np.arange(lo,hi+1)
    choose=gammaln(257)-gammaln(supports+1)-gammaln(257-supports)
    logs=np.array([log(c) if c else -np.inf for c in increments(counts,lo,hi)])
    def evaluate(p):
        weights=-choose-supports*np.log(p)-(256-supports)*np.log1p(-p)
        majorant=np.maximum.accumulate(weights[::-1])[::-1]
        return float(logsumexp(logs+majorant))
    return evaluate


def fold_arb(counts,lo,hi,p):
    weights=[]
    for u in range(lo,hi+1):
        mass=arb(comb(256,u))*p**u*(1-p)**(256-u)
        lower=arb(mass.lower())
        assert lower>0
        weights.append(up(1/lower))
    for i in range(len(weights)-2,-1,-1):
        weights[i]=max(weights[i],weights[i+1])
    return up(sum((c*w for c,w in zip(increments(counts,lo,hi),weights)),arb(0)))


def self_test():
    checks=0
    for shells in product(range(3),repeat=4):
        caps=[sum(shells[:j+1])+j+1 for j in range(4)]
        for lo in range(4):
            for hi in range(lo,4):
                for raw in ((1,4,2,3),(4,3,2,1),(1,2,3,4)):
                    weights=[Fraction(x,7) for x in raw[lo:hi+1]]
                    envelope=[max(weights[j:]) for j in range(len(weights))]
                    bound=sum(c*w for c,w in zip(increments(caps,lo,hi),envelope))
                    actual=sum(shells[j]*weights[j-lo] for j in range(lo,hi+1))
                    assert actual<=bound<=caps[hi]*max(weights)
                    checks+=1
    print('CDF folding:',checks,'exact weighted interval inequalities',flush=True)
    counts=[u*u+1 for u in range(257)]
    for lo,hi in ((38,256),(66,92),(128,128),(200,256)):
        fast=fold_function(counts,lo,hi)
        for p in (.0004,.1,.5,.9,.999999):
            assert abs(fast(p)-fold_log(counts,lo,hi,p))<1e-10
    print('Vectorized witness objectives match scalar reference',flush=True)


class RetainedCover:
    """A node may certify its entire box or use all of its children."""
    def __init__(self):
        self.nodes={}

    def add(self,item,parent=None):
        identifier=item[1]
        self.nodes[identifier]=dict(item=item,parent=parent,children=[],best=-item[0],use_children=False)
        if parent is not None:
            self.nodes[parent]['children'].append(identifier)

    def update(self,identifier):
        while identifier is not None:
            node=self.nodes[identifier]
            own=-node['item'][0]
            children=node['children']
            other=float(logsumexp([self.nodes[c]['best'] for c in children])) if children else float('inf')
            best=min(own,other)
            assert best<=node['best']+1e-10
            node['best']=best
            node['use_children']=other<own
            identifier=node['parent']

    def selected(self,identifier=1):
        node=self.nodes[identifier]
        if not node['use_children']:
            return [node['item']]
        return [item for child in node['children'] for item in self.selected(child)]


def retained_test():
    tree=RetainedCover()
    for identifier,bound,parent in ((1,4,None),(2,2,1),(3,3,1)):
        tree.add((-log(bound),identifier),parent)
    tree.update(1)
    assert [i[1] for i in tree.selected()]==[1]
    tree.add((-log(.2),4),2);tree.add((-log(.3),5),2)
    tree.update(2)
    assert [i[1] for i in tree.selected()]==[4,5,3]
    tree.add((-log(2),6),3);tree.add((-log(2),7),3)
    tree.update(3)
    assert [i[1] for i in tree.selected()]==[4,5,3]
    assert abs(tree.nodes[1]['best']-log(3.5))<1e-14
    print('Retained-cover parent/child and monotonicity checks passed',flush=True)


def placement_prefix_test():
    """Higher retained degrees cannot affect a lower occupancy coefficient."""
    from flint import fmpq_mat
    local=[fmpq_mat([[1,1],[0,1]]),fmpq_mat([[1,0],[1,1]]),
           fmpq_mat([[2,1],[0,1]])]
    full=placement(local,epochs=3,windows=2,matrix=fmpq_mat,
                   rounding=lambda x:x,maximum_groups=6)
    checks=0
    for degree in range(7):
        short=placement(local,epochs=3,windows=2,matrix=fmpq_mat,
                        rounding=lambda x:x,maximum_groups=degree)
        assert short==full[:degree+1]
        checks+=len(short)
    print('Placement prefix reuse:',checks,'exact noncommuting coefficient checks',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=7,choices=range(1,2049),metavar='Q')
    parser.add_argument('--occupancies',type=int,nargs='+',choices=range(1,2049),metavar='Q',
                        help='Replay each listed occupancy using shared region operators; overrides --groups')
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--target-bits',type=int,default=55)
    parser.add_argument('--max-splits',type=int,default=4000)
    parser.add_argument('--screen-only',action='store_true')
    parser.add_argument('--fresh',action='store_true')
    parser.add_argument('--window-average',action='store_true')
    parser.add_argument('--multi-average',action='store_true')
    parser.add_argument('--prefix-flags',action='store_true')
    parser.add_argument('--shortening-moments',action='store_true',help='Exact higher-rank count bounds from shortened-code containment moments')
    parser.add_argument('--positive-shortening',action='store_true',help='Additional exact positive-polynomial shortening dimension bounds')
    parser.add_argument('--dual-shortening',action='store_true',help='Replay BCH dual-distance premise and bound complementary shortened dimensions')
    parser.add_argument('--dual-moments',action='store_true',help='Positive dual-complement identities for averaged shortening counts')
    parser.add_argument('--dual-sandwich',action='store_true',help='Exact repaired BCH-sandwich witnesses for selected dual shells; enables dual moments')
    parser.add_argument('--window-histogram',type=int,choices=range(1,13),help='Exact distinct-window output moments through this occupancy, using expansion histograms')
    parser.add_argument('--joint-cancellation',action='store_true',help='Joint feedback/output moments through three active windows')
    parser.add_argument('--fresh-collision',action='store_true')
    parser.add_argument('--zero-moment',action='store_true')
    parser.add_argument('--full-feedback',type=int,choices=range(2,11),help='Exact feedback distributions through this local occupancy')
    parser.add_argument('--mature-tail',type=int,choices=(56,64))
    parser.add_argument('--pair-tail',action='store_true')
    parser.add_argument('--class-tail',action='store_true')
    parser.add_argument('--mixing-rounds',type=int,choices=range(1,9),default=1)
    parser.add_argument('--weight-tilt',default='1',help='Total input-weight tilt >=1; jointly coupled to all-one penalties')
    parser.add_argument('--retain-parents',action='store_true')
    parser.add_argument('--joint-witness',action='store_true')
    parser.add_argument('--joint-top',type=int,default=2)
    parser.add_argument('--probe-supports',type=int,nargs='+',default=[])
    parser.add_argument('--probe-vector',type=int,nargs='+',action='append',default=[])
    parser.add_argument('--penalties',nargs='+',default=['1'])
    parser.add_argument('--tilts',nargs='+',default=['.0016','.0025','.0032','.004','.005','.0064','.008'])
    args=parser.parse_args()
    occupancies=sorted(set(args.occupancies or [args.groups]))
    args.groups=max(occupancies)
    assert arb(args.weight_tilt)>=1
    self_test();geometry_test()
    retained_test();placement_prefix_test()
    if args.mixing_rounds!=1:
        import mixing_rounds
        mixing_rounds.self_test()
        print('ALTERNATIVE INNER:',args.mixing_rounds,'independent transvections per step; certificate scope excludes performance.',flush=True)
    if args.fresh_collision:
        import fresh_collision
        fresh_collision.self_test()
    spectrum=authenticated_caps()
    dimensions=dimension_caps()
    if args.positive_shortening:
        import shortening_polynomial
        shortening_polynomial.self_test()
        dimensions=shortening_polynomial.improve_dimensions(dimensions)
    if args.dual_shortening:
        import dual_shortening
        dual_shortening.self_test()
        dimensions=dual_shortening.improve_dimensions(dimensions)
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    if args.shortening_moments:
        import shortening_moments
        shortening_moments.self_test()
        caps,_=shortening_moments.improve(caps,dimensions)
        print('Outer counts refined by exact shortened-code containment moments',flush=True)
    if args.dual_moments or args.dual_sandwich:
        import dual_moments
        dual_moments.self_test()
        dual_spectrum=None
        if args.dual_sandwich:
            import dual_sandwich
            dual_sandwich.self_test()
            dual_spectrum=dual_sandwich.refined_spectrum()
        caps=dual_moments.refine_bch(caps,dimensions,dual_spectrum=dual_spectrum)
    counts=[sum(row[u] for row in caps) for u in range(257)]
    assert all(a<=b for a,b in zip(counts,counts[1:]))
    prepared=prepare(local_data(4))
    terminal=TERMINAL
    pair_tail_probabilities=None
    assert not args.pair_tail or args.mature_tail
    assert not args.class_tail or args.mature_tail
    if args.mature_tail:
        import mature_tail
        tail_inputs=mature_tail.prepare_inputs()
        tail_probabilities=mature_tail.census(tail_inputs,prepared)
        terminal=mature_tail.TAIL_TERMINAL
        if args.pair_tail:
            import pair_tail
            pair_tail_probabilities=pair_tail.census(tail_inputs,prepared)
    size=len(terminal)
    if args.zero_moment:
        import zero_moment
        zeros=zero_moment.census()
    if args.full_feedback:
        import full_feedback_refinement
        feedback=full_feedback_refinement.census(args.full_feedback)
        full_feedback_refinement.check_pairs(feedback)
    if args.joint_cancellation:
        import cancellation_joint
        cancellation_data=cancellation_joint.census()
        cancellation_joint.check_fresh(cancellation_data,prepared)
        if args.full_feedback and args.full_feedback>=3:
            cancellation_joint.check_feedback(cancellation_data,feedback)
    from occupancy_allones import weighted_cdf
    count_sets={p:weighted_cdf(spectrum,caps,dimensions,p,prefix_flags=args.prefix_flags) for p in args.penalties}
    if args.weight_tilt!='1':
        import total_weight_caps
        ctx.prec=args.precision
        total_weight_caps.self_test()
        import joint_weight_caps
        joint_weight_caps.self_test()
        count_sets={p:(total_weight_caps.caps(spectrum,counts,args.weight_tilt) if p=='1' else
                       joint_weight_caps.caps(spectrum,caps,dimensions,args.weight_tilt,p,args.prefix_flags))
                    for p in args.penalties}
    if args.fresh or args.multi_average:
        import occupancy_fresh_moment as fresh_moment
        fresh=fresh_moment.fresh_census(prepared)
    if args.multi_average:
        import occupancy_multi_average as multi_average
        multi_average.self_test()
    if args.window_average:
        import occupancy_window_average as window_average
        windows=window_average.prepare_inputs()
    if args.window_histogram:
        import window_histogram
        window_histogram.self_test()
        histogram_data=window_histogram.census(prepared)
    ctx.prec=args.precision
    # The detailed census includes four local occupancies. Extra degrees
    # are harmless when the requested region coefficient has degree <4.
    local_groups=max(4,min(args.groups,32))
    operators={}
    for tilt in args.tilts:
        if args.window_histogram:
            histogram_moments=window_histogram.moments(histogram_data,tilt,args.window_histogram)
        if args.window_average:
            values=window_average.averages(windows,tilt)
        for penalty in args.penalties:
            if args.fresh:
                ops=fresh_moment.refine(prepared,fresh,tilt,local_groups,penalty,args.weight_tilt)
            else:
                ops=epoch_operators(prepared,tilt,maximum_groups=local_groups,pair_conditioned=True,full_penalty=penalty,input_penalty=args.weight_tilt)
            if args.window_average:
                ops=window_average.refine(prepared,values,tilt,local_groups,penalty,ops,args.weight_tilt)
            if args.multi_average:
                ops=multi_average.refine(ops,fresh,tilt,penalty,args.weight_tilt)
            if args.fresh_collision:
                ops=fresh_collision.refine(ops,prepared,tilt,penalty,args.weight_tilt)
            if args.zero_moment:
                ops=zero_moment.refine(ops,zeros,tilt,penalty,args.weight_tilt)
            if args.mature_tail:
                ops=mature_tail.lift(ops,prepared[0][0][0],tail_probabilities,tilt,penalty,args.mature_tail,pair_tail_probabilities,args.class_tail,args.weight_tilt)
                mature_tail.self_test(tail_inputs,ops,tilt,penalty,input_penalty=args.weight_tilt)
            if args.mixing_rounds!=1:
                ops=mixing_rounds.transform(ops,prepared[0][0][0],tilt,args.mixing_rounds)
                if args.mature_tail:
                    mature_tail.self_test(tail_inputs,ops,tilt,penalty,args.mixing_rounds,args.weight_tilt)
            if args.window_histogram:
                ops=window_histogram.refine(ops,histogram_data,histogram_moments,penalty,
                                            args.mixing_rounds,args.weight_tilt)
            if args.full_feedback:
                ops=full_feedback_refinement.refine(ops,feedback,prepared[0][0][0],tilt,penalty,
                                                    args.mixing_rounds,args.weight_tilt)
            if args.joint_cancellation:
                ops=cancellation_joint.refine(ops,cancellation_data,tilt,penalty,
                                              args.mixing_rounds,args.weight_tilt)
            exact=placement(ops,rounding=rounded,maximum_groups=args.groups)
            arrays=[np.array([[float(t[i,j]) for j in range(size)] for i in range(size)]) for t in exact]
            operators[tilt,penalty]=exact,arrays
    print('Built CDF-cover operators',list(operators),flush=True)
    verified=[]
    for groups in occupancies:
        args.groups=groups
        print('Starting occupancy cover',groups,flush=True)
        verified.append(cover(args,operators,count_sets,terminal))
    if args.occupancies:
        if any(value is None for value in verified):
            print('Batch incomplete: at least one requested occupancy has no outward certificate.',flush=True)
        else:
            total=arb(0)
            for value in verified:total=up(total+value)
            assert total>0
            print('VERIFIED OCCUPANCY BATCH',occupancies,'precision',args.precision,
                  'upper',total,'margin',-total.log()/arb(2).log(),flush=True)
            if occupancies!=list(range(1,2049)):
                print('Unlisted occupancies remain; not a full-code certificate.',flush=True)


def cover(args,operators,count_sets,terminal,*,cutoff=209715):
    """Independent support cover and outward replay for args.groups."""
    if type(cutoff) is not int or not 0<=cutoff<(1<<21):
        raise ValueError('integer output-weight cutoff in [0,2^21) required')
    size=len(terminal)
    assert all(len(exact)>args.groups and len(region)>args.groups
               for exact,region in operators.values())
    cache={}
    def proposal(choice,interval):
        key=choice,interval
        if key not in cache:
            lo,hi=interval
            region=operators[choice][1]
            counts=count_sets[choice[1]]
            fold=fold_function(counts,lo,hi)
            def objective(z):
                p=1/(1+np.exp(-z))
                return log_power_moment(matrix_for_probabilities(region,[p]*args.groups),256,terminal)+args.groups*fold(p)
            fit=minimize_scalar(objective,bounds=(-8.,16.),method='bounded',options={'xatol':1e-7})
            p=1/(1+np.exp(-fit.x))
            numerator=max(1,min(DENOMINATOR-1,round(p*DENOMINATOR)))
            cache[key]=numerator,fold(numerator/DENOMINATOR)
        return cache[key]
    def witness(box):
        best=None
        candidates=[]
        for choice,(_,region) in operators.items():
            tilt,penalty=choice
            items=[proposal(choice,interval) for interval in box]
            nums=[n for n,_ in items]
            value=log_power_moment(matrix_for_probabilities(region,[n/DENOMINATOR for n in nums]),256,terminal)
            value+=float(tilt)*cutoff+sum(w for _,w in items)
            if best is None or value<best[0]:
                best=value,choice,nums
            candidates.append((value,choice,nums))
        if args.joint_witness:
            intervals=sorted(set(box))
            indices=[intervals.index(interval) for interval in box]
            for _,choice,nums in sorted(candidates)[:args.joint_top]:
                region=operators[choice][1]
                counts=count_sets[choice[1]]
                functions=[fold_function(counts,lo,hi) for lo,hi in intervals]
                start=np.array([nums[box.index(interval)]/DENOMINATOR for interval in intervals])
                def objective(logits):
                    ps=1/(1+np.exp(-logits))
                    value=log_power_moment(matrix_for_probabilities(region,[ps[i] for i in indices]),256,terminal)
                    folds=[f(p) for f,p in zip(functions,ps)]
                    return value+float(choice[0])*cutoff+sum(folds[i] for i in indices)
                if getattr(args,'analytic_gradient',False):
                    from gf16_packets.cdf_gradient import JointObjective
                    analytic=JointObjective(region,box,counts,terminal,float(choice[0])*cutoff)
                    fit=minimize(analytic,np.log(start/(1-start)),jac=True,method='L-BFGS-B',
                                 bounds=[(-8.,16.)]*len(start),options={'maxiter':60,'ftol':1e-11})
                else:
                    fit=minimize(objective,np.log(start/(1-start)),method='L-BFGS-B',bounds=[(-8.,16.)]*len(start),
                                 options={'maxiter':60,'ftol':1e-11})
                ps=1/(1+np.exp(-fit.x))
                numbers=[max(1,min(DENOMINATOR-1,round(p*DENOMINATOR))) for p in ps]
                quantized=np.array(numbers)/DENOMINATOR
                value=objective(np.log(quantized/(1-quantized)))
                if value<best[0]:
                    best=value,choice,[numbers[i] for i in indices]
        return best
    if args.probe_supports or args.probe_vector:
        vectors=[(u,)*args.groups for u in args.probe_supports]+[tuple(v) for v in args.probe_vector]
        for vector in vectors:
            assert len(vector)==args.groups and all(38<=u<=256 for u in vector)
            value,choice,nums=witness(tuple((u,u) for u in sorted(vector)))
            multiplicity=factorial(args.groups)
            for n in Counter(vector).values():
                multiplicity//=factorial(n)
            score=(value+log(comb(2048,args.groups)*multiplicity))/log(2)
            print('POINT support vector',vector,'log2 upper',score,'choice',choice,'numerators',nums,flush=True)
        print('Selected-point binary64 diagnostic with the cover CDF; no complete coverage or outward certificate.')
        return
    heap=[];serial=0
    tree=RetainedCover() if args.retain_parents else None
    locations=comb(2048,args.groups)
    def push(box,mult,parent=None):
        nonlocal serial
        score,tilt,ps=witness(box)
        serial+=1
        item=(-score-log(locations*mult),serial,box,mult,tilt,ps)
        heappush(heap,item)
        if tree is not None:
            tree.add(item,parent)
    push(((38,256),)*args.groups,1)
    # This buffer only decides when to attempt the outward replay. It is
    # not part of the proof bound: the final Arb result must still be
    # strictly below 2**(-target_bits). Preserve the existing default.
    proposal_buffer=getattr(args,'proposal_buffer_bits',2)
    if not 0<=proposal_buffer<=64:
        raise ValueError('proposal buffer must be between zero and 64 bits')
    threshold=-(args.target_bits+proposal_buffer)*log(2)
    check_interval=getattr(args,'check_interval',25)
    if type(check_interval) is not int or check_interval<1:
        raise ValueError('positive integer cover check interval required')
    score=-heap[0][0]
    steps=0
    while heap and steps<args.max_splits and score>threshold:
        item=heappop(heap)
        _,_,box,mult,_,_=item
        if all(lo==hi for lo,hi in box):
            if tree is not None:
                # Freeze this leaf; other branches can still improve.
                continue
            heappush(heap,item)
            print('Dominating singleton',box,flush=True)
            break
        for child,child_mult in split(box,mult):
            push(child,child_mult,item[1] if tree is not None else None)
        if tree is not None:
            tree.update(item[1])
        steps+=1
        if steps%check_interval==0:
            score=tree.nodes[1]['best'] if tree is not None else float(logsumexp([-x[0] for x in heap]))
            if steps%250==0:
                print('CDF-cover splits',steps,'leaves',len(heap),'log2 union',score/log(2),flush=True)
    if tree is not None:
        heap=tree.selected()
    score=float(logsumexp([-x[0] for x in heap]))
    assert sum(volume(x[2],x[3]) for x in heap)==219**args.groups
    print('Full-domain CDF cover:',steps,'splits;',len(heap),'leaves; log2 union',score/log(2),flush=True)
    print('Largest leaves:',[(-x[0]/log(2),x[2],x[4]) for x in sorted(heap)[:3]],flush=True)
    if score>threshold or args.screen_only:
        print('No outward certificate from this run.',flush=True)
        return
    total=arb(0)
    # Precision and count_sets stay fixed throughout this replay. Many
    # groups in a box have the same interval and exact probability; avoid
    # reconstructing their identical outward folds once per group.
    outward_folds={}
    for index,(_,_,box,mult,choice,nums) in enumerate(heap,1):
        tilt,penalty=choice
        ps=[arb(n)/DENOMINATOR for n in nums]
        masses=bernoulli_masses(nums)
        matrix=sum((m*r for m,r in zip(masses,operators[choice][0])),arb_mat(size,size))**256
        term=sum((matrix[0,j] for j in range(size) if terminal[j]),arb(0))*(arb(tilt)*cutoff).exp()*locations*mult
        for (lo,hi),numerator,p in zip(box,nums,ps):
            key=penalty,lo,hi,numerator
            if key not in outward_folds:
                outward_folds[key]=fold_arb(count_sets[penalty],lo,hi,p)
            term=up(term*outward_folds[key])
        total=up(total+term)
        if index%1000==0:
            print('Outward CDF-cover leaves',index,'/',len(heap),flush=True)
    assert 0<total<arb(2)**-args.target_bits
    print('VERIFIED all supports at occupancy',args.groups,'cutoff',cutoff,'precision',args.precision,'upper',total,
          'margin',-total.log()/arb(2).log(),flush=True)
    print('Other occupancies remain; not a full-code certificate.')
    return total


if __name__=='__main__':
    main()
