"""Common-state-potential bounds summed over one matching marginal.

Fix an active-group weight r and a positive two-state vector v.  If the
regional matrix R(r,f) satisfies R(r,f)v <= c(r,f)v for every full-pair
count f, positivity bounds the ordered 256-region product by the product
of its c factors, with initial factor v[0]/min(v).  Cauchy--Schwarz over
indexed binary supports gives

 sum_U c(r,a(U))**128 c(r,b(U))**128
 <= sum_f one_pairing_count(M,r,f) * c(r,f)**256.

Both alternating matchings have this same marginal law; their counts are
not independent.  The vector may depend on r, but never on f within this
sum.  This reduces the complete (a,b) domain to a one-dimensional sum.
The scalar-code support convention and original projective normalization
are unchanged. The saved directed checkpoints use cyclic-run convolution
from hd_shifted_run_potential, which retains more joint geometry than this
unweighted Cauchy bound. FLOAT screens are not directed certificates.
"""
from fractions import Fraction as F
from math import log2

import numpy as np
from flint import arb, ctx
from scipy.optimize import minimize, minimize_scalar

from hd_shifted_alias_certificate import rational
from hd_shifted_coverage import rational_ceiling
from hd_shifted_pairs_count import one_pairing_count
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION
from hd_sumrank_paired_certificate import scalar_power


def full_counts(pairs, weight):
    if (not isinstance(pairs, int) or not isinstance(weight, int)
            or pairs < 1 or not 0 <= weight <= 2*pairs):
        raise ValueError('Expected a valid integer pair count and weight')
    return range(max(0, weight-pairs), weight//2+1)


def positive_vector(vector):
    if len(vector) != 2 or any(isinstance(v, float) for v in vector):
        raise ValueError('Two positive exact rational potential entries required')
    vector = tuple(map(F, vector))
    if min(vector) <= 0:
        raise ValueError('Potential entries must be positive')
    return vector


def contraction_arb(entries, vector):
    """Directed c satisfying Av <= cv, from four nonnegative entries."""
    vector = positive_vector(vector)
    entries = tuple(rational(x) if isinstance(x, F) else arb(x) for x in entries)
    if len(entries) != 4 or any(not x >= 0 for x in entries):
        raise ValueError('A nonnegative flattened two-state matrix is required')
    v0, v1 = map(rational, vector)
    a, b, c, d = entries
    rows = ((a+b*v1/v0).upper(), (c*v0/v1+d).upper())
    return max(rows)


def marginal_sum_arb(pairs, weight, regions, vector, *, repetitions=256):
    """Bound the full support sum using one potential for all histograms.

    regions[f] must upper-enclose the SAME regional operator for each f.
    The default exponent is 256 chronological regions, 128 per matching.
    Smaller even exponents support independent exact test fixtures.
    """
    vector = positive_vector(vector)
    if not isinstance(repetitions, int) or repetitions < 2 or repetitions % 2:
        raise ValueError('A positive even number of regional repetitions is required')
    value = arb(0)
    for full in full_counts(pairs, weight):
        contraction = contraction_arb(regions[full], vector)
        value += arb(one_pairing_count(pairs, weight, full))*contraction**repetitions
    return (rational(vector[0]/min(vector))*value).upper()


def weighted_transfer_arb(pairs, run_tilt, singleton_marker, full_marker):
    """trace(K**pairs), marking t**(even11+odd11) on indexed cycles.

    K[l,v] sums the two choices of a pair's first bit u.  Its factor is
    S**[u!=v] F**[u=v=1] t**(l*u+u*v), where l is the preceding last bit.
    The trace closes the labelled cycle; it does not quotient rotations.
    """
    if not isinstance(pairs,int) or pairs<1:
        raise ValueError('Positive integer pair count required')
    if any(isinstance(x,float) for x in (run_tilt,singleton_marker,full_marker)):
        raise ValueError('Exact rational count markers required')
    t,s,f = map(F,(run_tilt,singleton_marker,full_marker))
    if min(t,s,f)<=0:
        raise ValueError('Positive count markers required')
    t,s,f = map(rational,(t,s,f))
    power = scalar_power((1+s,s+f*t,1+s*t,s+f*t*t),pairs)
    return power[0]+power[3]


def weighted_count_upper_arb(pairs, weight, full, run_tilt, count_markers):
    """Upper bound N_f(t)=sum_{U:a(U)=f,|U|=r} t**(a(U)+b(U))."""
    if full not in full_counts(pairs,weight) or len(count_markers)!=2:
        raise ValueError('Valid histogram and two count markers required')
    single = weight-2*full
    value = weighted_transfer_arb(pairs,run_tilt,*count_markers)
    s,f = map(rational,map(F,count_markers))
    return (value/(s**single*f**full)).upper()


def weighted_marginal_sum_arb(pairs, weight, regions, vector, run_tilt,
                              count_markers, *, repetitions=256):
    """Cauchy under the SAME positive cycle weight t**(a+b) for every f.

    Each N_f(t) has its own positive coefficient markers, used only to
    upper-bound that count. The tilt t and the state potential v are common
    across f. No independence of the two matching statistics is assumed.
    """
    vector = positive_vector(vector)
    if not isinstance(repetitions,int) or repetitions<2 or repetitions%2:
        raise ValueError('A positive even repetition count is required')
    if isinstance(run_tilt,float) or F(run_tilt)<=0:
        raise ValueError('A positive exact rational cycle tilt is required')
    tilt,value = rational(F(run_tilt)),arb(0)
    for full in full_counts(pairs,weight):
        count = weighted_count_upper_arb(pairs,weight,full,run_tilt,count_markers[full])
        contraction = contraction_arb(regions[full],vector)
        value += count*contraction**repetitions/tilt**(2*full)
    return (rational(vector[0]/min(vector))*value).upper()


def marginal_log_screen(pairs, weight, log_regions, log_ratio, *, repetitions=256):
    """FLOAT log2 marginal sum, including the initial-state potential cost."""
    values = []
    for full in full_counts(pairs, weight):
        a, b, c, d = np.asarray(log_regions[full], dtype=float).ravel()
        contraction = max(np.logaddexp2(a, b+log_ratio),
                          np.logaddexp2(c-log_ratio, d))
        values.append(log2(one_pairing_count(pairs, weight, full))
                      +repetitions*contraction)
    return float(np.logaddexp2.reduce(values)+max(0., -log_ratio))


def optimize_potential_screen(pairs, weight, log_regions, *, bounds=(-256., 256.),
                              repetitions=256):
    """Optimize only one state ratio; regional matrices and witness stay fixed."""
    objective = lambda ratio: marginal_log_screen(pairs, weight, log_regions, ratio,
                                                  repetitions=repetitions)
    result = minimize_scalar(objective, bounds=bounds, method='bounded',
                             options={'xatol': 1e-9, 'maxiter': 100})
    ratio = float(result.x)
    return dict(log2_potential_ratio=ratio, log2_marginal_screen=objective(ratio),
                evaluations=int(result.nfev), directed=False)


def atlas_bank_screen(weight, z, alpha, centers, *, grid=128):
    """Bounded FLOAT helper: fit markers at explicitly supplied f centers.

    Each atlas bounds every f.  Normalize first, then take entrywise minima.
    The caller chooses the finite centers; this does not launch a search
    over weights, z, alpha, or support profiles.  Returned values are logs.
    """
    from hd_shifted_bulk_search import LOG_FACTORIAL
    from hd_shifted_projected_screen import alias_table
    from hd_shifted_swap_screen import positive_pair_matrix
    feasible = tuple(full_counts(4096, weight))
    centers = tuple(sorted(set(centers)))
    if not centers or len(centers)>100 or any(f not in feasible for f in centers):
        raise ValueError('At most100 explicitly supplied feasible marker centers required')
    full = np.array(feasible)
    single, empty = weight-2*full, 4096-weight+full
    log_mult = (LOG_FACTORIAL[4096]-LOG_FACTORIAL[empty]
                -LOG_FACTORIAL[single]-LOG_FACTORIAL[full])
    best = np.full((4,len(full)), np.inf)
    witnesses = []
    for target in centers:
        singleton, vacant = weight-2*target, 4096-weight+target
        if vacant == 0:
            raise ValueError('This medium-weight marker proposer requires an empty pair')
        initial = np.array((log2(max(1,singleton)/vacant)+.5,
                            log2(max(1,target)/vacant)+.8))
        bounds = [(-30.,5.), (-30.,5.)]
        if singleton == 0:
            initial[0],bounds[0] = -64.,(-64.,-64.)
        if target == 0:
            initial[1],bounds[1] = -64.,(-64.,-64.)
        def objective(odds):
            positive = positive_pair_matrix(z,alpha,tuple(odds))
            scale = float(np.max(positive))
            a,b,c,d = (positive/scale).ravel()
            rho = (a+d+np.sqrt((a-d)**2+4*b*c))/2
            return 64*(log2(scale)+log2(rho))-singleton*odds[0]-target*odds[1]
        result = minimize(objective, initial, method='L-BFGS-B', bounds=bounds,
                          options={'maxiter':60, 'ftol':1e-11})
        odds = tuple(map(float,result.x))
        table = alias_table(z,alpha,odds,grid)
        candidate = (table[:,single%grid,full%grid]-log_mult
                     -single*odds[0]-full*odds[1])
        np.minimum(best,candidate,out=best)
        witnesses.append(dict(full_pairs=int(target),odds=odds,
                              objective_evaluations=int(result.nfev)))
    return {int(f): best[:,i].reshape(2,2) for i,f in enumerate(full)}, witnesses


class SwapPotentialBound:
    """Reuse directed positive/alias witnesses for complete weight classes.

    Each marker/parity atlas bounds the same swap-averaged regional matrix.
    Taking their entrywise minimum AFTER normalization remains an upper
    bound and can improve the coefficient enclosure.  No parity-dependent
    potential is introduced by this minimum.
    """
    def __init__(self, witness):
        if witness.construction != SWAP_CONSTRUCTION:
            raise ValueError('This helper requires the independently swapped construction')
        self.witness = witness
        self.cache = {}

    def regional(self, weight, full):
        if full not in full_counts(4096, weight):
            raise ValueError('Invalid regional histogram')
        key = weight, full
        if key not in self.cache:
            w = self.witness
            with ctx.workprec(w.precision):
                hist = 4096-weight+full, weight-2*full, full
                log_mult = w.log_factorial(4096)-sum((w.log_factorial(n) for n in hist), arb(0))
                rows = []
                for parity in (0, 1):
                    numerator = (w.regions[parity] if w.grid is None or w.uniform_aliases
                                 else w.regions[parity].numerator(((hist[1], hist[1]), (full, full))))
                    normalizer = log_mult+sum((rational(n)*mark for n, mark in
                                               zip(hist[1:], w.log_markers[parity])), arb(0))
                    scale = arb(2)**(-normalizer)
                    rows.append(tuple((entry*scale).upper() for entry in numerator))
                self.cache[key] = tuple(min(a, b) for a, b in zip(*rows))
        return self.cache[key]

    def receipt(self, weight, vector):
        if not 1 <= weight <= 8192:
            raise ValueError('Nonzero scalar support weight required')
        vector = positive_vector(vector)
        w = self.witness
        with ctx.workprec(w.precision):
            regions = {f: self.regional(weight, f) for f in full_counts(4096, weight)}
            value = marginal_sum_arb(4096, weight, regions, vector)
            bound = value.log()/arb(2).log()+w.slope*weight+w.constant
            ceiling = rational_ceiling(bound)
        return dict(active_groups=weight, potential=list(map(str, vector)),
                    log2_upper_rational=str(ceiling), log2_upper_enclosure=str(bound),
                    directed=True, covers_all_matching_profiles=True,
                    covers_all_even_odd_splits=True, whole_code_certificate=False,
                    construction=w.construction, field_order=256, forbidden_weight=734003,
                    support_count_method='common-potential-cauchy-marginal',
                    regional_histograms=len(regions), witness=w.witness())


def directed_bank_regions(weight,z,alpha,markers,*,grid=128,precision=256,bits=60,progress=None):
    """Normalized regional enclosures from a finite exact marker bank.

    All atlases use the same z, alpha, coefficient rounding, and swap law.
    Construct each atlas once and discard it after taking the normalized
    entrywise minimum. This minimum upper-bounds the same exact operator.
    """
    from hd_shifted_swap_certificate import swap_unmarked
    from hd_shifted_unoriented_certificate import unsigned_components
    from hd_shifted_good_pair_certificate import component_power
    from hd_shifted_coverage import RegionalAliases
    feasible = tuple(full_counts(4096,weight))
    if (weight==0 or isinstance(z,float) or isinstance(alpha,float)
            or not 0<F(z)<1 or not 0<F(alpha)<=1 or precision<128 or bits<16
            or not isinstance(grid,int) or grid<64 or grid&(grid-1)):
        raise ValueError('Invalid exact directed atlas-bank witness')
    rows = []
    for row in markers:
        if len(row)!=2 or any(isinstance(x,float) for x in row):
            raise ValueError('Exact two-entry regional markers required')
        row = tuple(map(F,row))
        if min(row)<=0:
            raise ValueError('Positive regional markers required')
        if row not in rows:
            rows.append(row)
    if not rows:
        raise ValueError('A nonempty marker bank is required')
    z,alpha = F(z),F(alpha)
    with ctx.workprec(precision):
        single,paired = swap_unmarked(z,alpha,precision,bits)
        log_two = arb(2).log()
        factorial = {n:(arb(n)+1).lgamma()/log_two for n in
                     {4096}|{v for f in feasible for v in (4096-weight+f,weight-2*f,f)}}
        log_mult = {f:factorial[4096]-sum((factorial[n] for n in
                    (4096-weight+f,weight-2*f,f)),arb(0)) for f in feasible}
        best = {}
        for index,row in enumerate(rows):
            pieces = unsigned_components(single,paired,row)
            powered,vector,rho = component_power(*pieces,grid,dimensions=2)
            atlas = RegionalAliases.from_power(powered,vector,rho)
            logs = tuple(rational(m).log()/log_two for m in row)
            for full in feasible:
                singleton = weight-2*full
                numerator = atlas.numerator(((singleton,singleton),(full,full)))
                scale = arb(2)**(-log_mult[full]-singleton*logs[0]-full*logs[1])
                values = tuple((entry*scale).upper() for entry in numerator)
                best[full] = values if full not in best else tuple(min(a,b) for a,b in zip(best[full],values))
            if progress:
                progress(dict(stage='directed-regional-atlas',weight=weight,precision=precision,
                              completed=index+1,total=len(rows)))
    return best


def run_bank_receipt(weight,z,alpha,markers,vector,*,grid=128,precision=256,bits=60,progress=None):
    """Directed complete-weight receipt using root's positive run convolution."""
    from hd_grouped_packet_screen import LIMIT
    from hd_shifted_cover_backend import profile_count
    from hd_shifted_run_potential import run_sum_arb
    vector = positive_vector(vector)
    regions = directed_bank_regions(weight,z,alpha,markers,grid=grid,precision=precision,
                                    bits=bits,progress=progress)
    z,alpha = F(z),F(alpha)
    with ctx.workprec(precision):
        value = run_sum_arb(4096,weight,regions,vector)
        log = lambda x:x.log()/arb(2).log()
        prefix = rational(alpha)*(weight*(2048-128*log(arb(255)))-log(arb(255))-LIMIT*log(rational(z)))
        bound = log(value)+prefix
        ceiling = rational_ceiling(bound)
    ranges = ((weight,weight),(0,4096),(0,4096))
    witness = dict(z=str(z),alpha=str(alpha),markers=[list(map(str,map(F,row))) for row in markers],
                   potential=list(map(str,vector)),grid=grid,precision=precision,
                   packet_mantissa_bits=bits,regional_minimum='after-histogram-normalization')
    return dict(active_groups=weight,ranges=ranges,integer_profile_triples=profile_count(ranges),
                log2_upper_rational=str(ceiling),log2_upper_enclosure=str(bound),witness=witness,
                directed=True,covers_all_supports_at_weight=True,covers_all_even_odd_splits=True,
                whole_code_certificate=False,construction=SWAP_CONSTRUCTION,pair_swaps=True,
                field_order=256,outer_groups=8192,source_dimension=1048576,output_length=2097152,
                forbidden_weight=LIMIT,support_count_method='common-potential-run-convolution',
                coefficient_method='minimum-normalized-directed-alias-bank',regional_histograms=len(regions))


# Exact selected witnesses, not saved operator values or probability bounds.
# Replays below regenerate every packet coefficient and every alias atlas.
SAVED_RUN_WITNESSES = {
    128: dict(z='9106085449535871/10000000000000000',
              alpha='1878066517809951/20000000000000000',
              potential='4808642444/418075895493',markers=(
        '12249948778/282683253799 1/18446744073709551616',
        '15650021973/483988410832 2774711160/405666644623',
        '23110270223/1076963431914 14303224219/1051487276839',
        '7143748786/668317609815 18880210059/930406732954',
        '454707733/85210235729 13724145229/581217416188',
        '2265348929/849642639497 22160035071/877016045686',
        '1092562161/819842133386 16812314195/644306604816',
        '398939823/598823159957 5591052683/210932160619',
        '1/18446744073709551616 24454259239/908440550132')),
    512: dict(z='6918720502163723/10000000000000000',
              alpha='9447819854595409/100000000000000000',
              potential='4168407278/596457008347',markers=(
        '59417882768/309709900899 1/18446744073709551616',
        '201812822104/1056356808593 401880885/851558688881',
        '164772873407/866120123277 584689399/619633773277',
        '174958093384/927486939273 1966549621/1042636041934',
        '155100255415/836432258001 3354808987/890364586774',
        '186687493595/1042710209119 721316973/95937455659',
        '171291200345/1029725637773 11970335517/799646890463',
        '54859639713/388279749724 17881652805/602509416316',
        '101481765713/869985580217 17526769967/396978194341',
        '67651814751/732171270166 58613641044/1003472856257',
        '7446726490/108630500969 72689031967/1002942927763',
        '9094255955/201288233666 90129833450/1044097212759',
        '17004746361/759424938500 46691828118/467486427829',
        '1234101024/110330911127 35814077351/336181576340',
        '5827373941/1040866341558 95944707418/873543548959',
        '2720525993/970916752699 122260557893/1096696112758',
        '125438627/89479922874 39161420042/348710419895',
        '367281305/523819444216 102659048938/910785775861',
        '1/18446744073709551616 6819521031/60282525847')),
    1024: dict(z='4824969675687487/10000000000000000',
               alpha='10257997854894543/100000000000000000',
               potential='934811780/346161200607',markers=(
        '397541866217/871886858341 1/18446744073709551616',
        '432896428609/1023320767978 1491808961/82572791389',
        '383878181278/982274301931 25438293476/711610592835',
        '390280247001/1086378389006 49064392019/924683608692',
        '42023075553/127988845241 76359638149/1090504739735',
        '295095174142/990134539531 59423423580/685813105781',
        '195307149079/727862155760 101985644227/990658113470',
        '236148269998/987267682485 130744019707/1099265020979',
        '163757730466/777569506381 140804871569/1045849571264',
        '139454374161/764028535961 129008976559/859782117622',
        '36126204046/233182496101 166844738399/1009886641623',
        '121359655627/949845865351 82181714623/456175470581',
        '100708471738/996850821893 214029527499/1098082472122',
        '76206111611/1019446425378 152265400499/726796791407',
        '34522123431/702757244614 234663653615/1048217219651',
        '23382775852/960728856337 207002680507/870224388926',
        '10029106906/822633969921 116474435268/475970228915',
        '5107059833/834752916510 264146556976/1064698598719',
        '1668057375/543838403764 256997234626/1028882617459',
        '679927807/442667993773 155598043876/620837728593',
        '372390307/484479532628 65696530415/261689799672',
        '1/18446744073709551616 72804904735/289518741617')),
}


def checkpoint_receipt(weight,*,precision=256,progress=None):
    if weight not in SAVED_RUN_WITNESSES:
        raise ValueError('No saved exact run-potential witness at this weight')
    witness = SAVED_RUN_WITNESSES[weight]
    markers = tuple(tuple(map(F,row.split())) for row in witness['markers'])
    return run_bank_receipt(weight,F(witness['z']),F(witness['alpha']),markers,
                            (F(1),F(witness['potential'])),precision=precision,progress=progress)


if __name__ == '__main__':
    import argparse,json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weights',type=int,nargs='+',default=sorted(SAVED_RUN_WITNESSES))
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--progress',action='store_true')
    args = parser.parse_args()
    progress = (lambda row:print(json.dumps(row),flush=True)) if args.progress else None
    for weight in args.weights:
        print(json.dumps(checkpoint_receipt(weight,precision=args.precision,progress=progress)),flush=True)
