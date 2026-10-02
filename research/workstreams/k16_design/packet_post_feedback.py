"""Isolated post-feedback mixing candidates; no retained certificate is reused.

Fix A:F_2^s -> F_2^t and C=A^T with CA=0. For each physical step, setup
samples an independent uniform M in GL(s,2). The first candidate emits
y=x+A*a, then sets a'=M*(a+C*x). State starts at zero and is not flushed.
This moves the existing uniform update; it does not add a second GL map.

Fix the complete input sequence independently of all M. The current output
is determined before using the current M. Conditional on a+C*x being nonzero,
M*(a+C*x) is uniform nonzero, independently of the weighted emitted history.
Consequently the tilted state measure is exactly a combination of zero mass
and uniform-nonzero mass after every step. This averages over setup matrices;
it does not condition on their realized values or assert seedwise uniformity.

Put m=2^s-1. For fixed x, let c=C*x and
    H_x(z) = sum_{a!=0} z^wt(x+A*a) / m,
    G_x(z) = 1_{c!=0} z^wt((I+A*C)*x) / m.
The exact two-state transfer has rows
    [1_{c=0} z^wt(x), 1_{c!=0} z^wt(x)], [G_x, H_x-G_x].
The return event and emitted weight are correlated. In particular, the
original recurrence's H/m return bound cannot be substituted for G here.

Physical occupancy j means a uniform j-subset of four-bit packet positions,
with independent uniform nonzero labels. We compute G_j exactly through a
requested small occupancy. Higher occupancies use valid lower/upper bounds.
The return map I+AC is an involution, not generally weight preserving.
For nonzero Cx, wt((I+AC)x)>=max(1,d_A-wt(x)); its parity equals wt(x).
If the constant word lies in im(A), Cx=0 additionally forces even input parity.
These facts, total return probability, and Cauchy give conservative bounds.
An upper-envelope row uses [G_upper,H_upper-G_lower], never H-G_upper.

The second candidate retains a'=M*a+C*x, then applies independent sampled
transvections after feedback. Each fixes zero and acts on a nonzero vector
as half identity plus half uniform refresh. Positive composition preserves
the old arbitrary/class/density envelope. It is not an exact two-state model.

All preparation reuses caller-supplied fresh map profiles. Small return
histograms are regenerated locally; no full state census runs here. Tail
comparisons cover only explicit q>=3 occupancies and cannot certify a code.
"""
from __future__ import annotations

from dataclasses import asdict
from fractions import Fraction as Q
import hashlib
from itertools import combinations, product
import json
from math import comb

from flint import arb, arb_mat, arb_poly, ctx
import packet_q1 as q1
import packet_rs_state_sparse as sparse
from packet_regional_power import placement_power
from packet_uniform_tail import regional_uniform
from rs_uniform_envelope import UniformInputEnvelope

aq, up = q1.kernel_t64.aq, q1.kernel_t64.up


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _map_structure(physical):
    images = q1.kernel_t64.kernel_maps.authenticate(physical)
    bits, width = physical['bits'], 4*physical['windows']
    if not 1 <= bits <= 22 or not 1 <= physical['windows'] <= 16:
        raise ValueError('state dimensions 1..22 and at most sixteen physical packets required')
    rows = tuple(int(images[1 << j]) for j in range(bits))
    columns = tuple(map(int, physical['columns']))
    actual = tuple(sum(((row >> i) & 1) << j for j, row in enumerate(rows)) for i in range(width))
    if columns != actual or any((a & b).bit_count() & 1 for a in rows for b in rows):
        raise ValueError('post-feedback candidate requires C=A^T and CA=0')
    rank = q1.kernel_t64.s16_maps.binary_rank
    if rank(rows) != bits:
        raise ValueError('injective expansion required')
    if (sum(map(int, physical['histogram_multiplicities'])) != (1 << bits)-1
            or sum(map(int, physical['multiplicities'])) != 1 << bits
            or any(sum(map(int, h)) != width//4 for h in physical['histograms'])
            or any(sum(map(int, h)) != width//4 for h in physical['records'])):
        raise ValueError('complete fresh image and character profiles required')
    distance = min(sum(i*int(n) for i, n in enumerate(h)) for h in physical['histograms'])
    if distance <= 0 or distance % 2:
        raise ValueError('injective even-weight expansion required')
    constant = rank([*rows, (1 << width)-1]) == bits
    return rows, columns, distance, constant


def return_histograms(physical, *, exact_through=3):
    """Count return-output weights without enumerating any GL matrices.

    Each j contributes C(W,j)*15^j inputs. The two histograms partition them
    according to Cx=0. On that event (I+AC)x=x, so its histogram also gives
    exact tilted zero-state self-loops. No floating arithmetic enters counts.
    """
    _map_structure(physical)
    W, width = physical['windows'], 4*physical['windows']
    if type(exact_through) is not int or not 0 <= exact_through <= min(3, W):
        raise ValueError('exact return cutoff must be between zero and min(3,W)')
    columns, images = physical['columns'], physical['map_images']
    packet = []
    for position in range(W):
        choices = []
        for label in range(1, 16):
            feedback = 0
            for bit in range(4):
                if label >> bit & 1:
                    feedback ^= int(columns[4*position+bit])
            x = label << (4*position)
            choices.append((feedback, x ^ int(images[feedback])))
        packet.append(choices)
    returned, zero, denominators = [], [], []
    for j in range(exact_through+1):
        rh, zh = [0]*(width+1), [0]*(width+1)
        for positions in combinations(range(W), j):
            for selected in product(*(packet[p] for p in positions)):
                feedback, transformed = 0, 0
                for c, y in selected:
                    feedback ^= c
                    transformed ^= y
                (rh if feedback else zh)[transformed.bit_count()] += 1
        denominator = comb(W, j)*15**j
        if (sum(rh)+sum(zh) != denominator or Q(sum(zh), denominator)
                != Q(physical['zero_probabilities'][j]) or rh[0]):
            raise ArithmeticError('return census disagrees with exact feedback probability or invertibility')
        returned.append(tuple(rh)); zero.append(tuple(zh)); denominators.append(denominator)
    return dict(exact_through=exact_through, return_counts=tuple(returned),
                zero_counts=tuple(zero), denominators=tuple(denominators))


def _profile_polynomial(histogram, inactive, active):
    result = arb_poly([1])
    for w, n in enumerate(histogram):
        result *= arb_poly([inactive[w], active[w]])**int(n)
    return result


def _moments(physical, z, *, zero_feedback=True):
    """Arb balls for all H_j and Z_j; no endpoint subtraction is hidden here."""
    W, S = physical['windows'], 1 << physical['bits']
    powers = [z**w for w in range(5)]
    active = [sum((comb(4, k)-int(k == w))*powers[k] for k in range(5))/15
              for w in range(5)]
    image = arb_poly([0])
    for histogram, multiplicity in zip(physical['histograms'], physical['histogram_multiplicities']):
        image += int(multiplicity)*_profile_polynomial(histogram, powers, active)
    means = [image[j]/((S-1)*comb(W, j)) for j in range(W+1)]
    if not zero_feedback:
        return means, None
    # Walsh inversion at syndrome zero. Signed intermediate coefficients are
    # enclosed as balls; only the final nonnegative probability is clipped.
    fourier = [((1+z)**(4-r)*(1-z)**r-1)/15 for r in range(5)]
    zero = arb_poly([0])
    for profile, multiplicity in zip(physical['records'], physical['multiplicities']):
        zero += int(multiplicity)*_profile_polynomial(profile, [arb(1)]*5, fourier)
    return means, [zero[j]/(S*comb(W, j)) for j in range(W+1)]


def _bounds(value):
    if not value.is_finite():
        raise ArithmeticError('finite interval required; increase proof precision')
    lo, hi = max(arb(0), value.lower()), value.upper()
    if hi < 0 or lo > hi:
        raise ArithmeticError('nonnegative mass has an inconsistent interval')
    return lo, max(arb(0), hi)


class PostFeedbackModel:
    """New one-GL recurrence, sharing maps but not old transition claims."""
    def __init__(self, physical, *, exact_through=3, source_map_record=None):
        rows, columns, self.distance, self.constant = _map_structure(physical)
        self.physical, self.bits, self.width = physical, physical['bits'], 4*physical['windows']
        self.counts = return_histograms(physical, exact_through=exact_through)
        self.source_map_record = source_map_record
        self.metadata = dict(schema='post-feedback-uniform-gl-local-model-1',
            recurrence='y=x+A*a; next_a=M*(a+C*x)', distribution='uniform_gl',
            physical_t=self.width, state_bits=self.bits, matrix_coordinates=['zero', 'uniform_nonzero'],
            physical_updates_independent=True, zero_initial_state=True, final_flush=False,
            maps_shared_with_original_only=True, source_map_record=source_map_record,
            map_sha256=physical['map_sha256'], expansion_rows_hex=list(map(hex, rows)),
            feedback_columns=list(columns), feedback_times_expansion_zero=True,
            contains_constant_expansion_word=self.constant, minimum_expansion_weight=self.distance,
            exact_return_through=exact_through, return_histogram_sha256=_digest(self.counts),
            return_map='I+A*C; an involution, not a weight isometry',
            higher_return_bounds='distance/parity, total return probability, and Cauchy; lower bounds retained',
            whole_code_certificate=False)

    def physical_operators_at_z(self, z, *, diagnostics=False):
        z = arb(z)
        if not z.is_finite() or not 0 < z <= 1:
            raise ValueError('finite output-weight variable in (0,1] required')
        W, m = self.physical['windows'], (1 << self.bits)-1
        means, zeros = _moments(self.physical, z)
        means2, _ = _moments(self.physical, z*z, zero_feedback=False)
        powers = [z**i for i in range(self.width+1)]
        packet = sum(comb(4, i)*powers[i] for i in range(1, 5))/15
        matrices, details = [], []
        for j in range(W+1):
            total = packet**j
            if j <= self.counts['exact_through']:
                denominator = self.counts['denominators'][j]
                zero = sum((n*powers[w] for w, n in enumerate(self.counts['zero_counts'][j])), arb(0))/denominator
                returned = sum((n*powers[w] for w, n in enumerate(self.counts['return_counts'][j])), arb(0))/(m*denominator)
                glo, ghi = _bounds(returned)
            else:
                zero = zeros[j]
                nonzero = 1-Q(self.physical['zero_probabilities'][j])
                probability = aq(nonzero)/m
                hlo, hhi = _bounds(means[j])
                _, h2hi = _bounds(means2[j])
                minimum = max(1, self.distance-4*j)
                glo = _bounds(probability*powers[self.width])[0]
                caps = [hhi, up(probability*powers[minimum]), up((h2hi*probability).sqrt())]
                if self.constant:
                    odd = (1-(-Q(1, 15))**j)/2
                    even_nonzero = nonzero-odd
                    if even_nonzero < 0:
                        raise ArithmeticError('constant feedback row must detect every odd input')
                    odd_min = minimum if minimum & 1 else minimum+1
                    even_min = minimum+1 if minimum & 1 else minimum
                    caps.append(up((aq(odd)*powers[odd_min]+aq(even_nonzero)*powers[even_min])/m))
                    glo = max(glo, _bounds((aq(odd)*powers[self.width-1]
                                           +aq(even_nonzero)*powers[self.width])/m)[0])
                ghi = min(caps)
            zlo, zhi = _bounds(zero)
            _, thi = _bounds(total)
            _, hhi = _bounds(means[j])
            zhi, ghi = min(zhi, thi), min(ghi, hhi)
            if glo > ghi or zlo > zhi:
                raise ArithmeticError('inconsistent return or zero-feedback bounds')
            zn = up(total-zlo)
            un = up(means[j]-glo)
            if zn < 0 or un < 0:
                raise ArithmeticError('negative complementary mass; increase precision')
            matrix = arb_mat([[zhi, zn], [ghi, un]])
            if any(not matrix[a, b].is_finite() or not matrix[a, b] >= 0
                   for a in range(2) for b in range(2)):
                raise ArithmeticError('finite nonnegative two-state envelope required')
            matrices.append(matrix)
            details.append(dict(occupancy=j, exact_return=j <= self.counts['exact_through'],
                                G_lower=glo, G_upper=ghi, H_upper=hhi, Z_lower=zlo, Z_upper=zhi))
        return (matrices, details) if diagnostics else matrices

    def physical_operators(self, tilt, *, diagnostics=False):
        tilt = Q(tilt)
        if tilt <= 0:
            raise ValueError('positive rational tilt required')
        return self.physical_operators_at_z((-aq(tilt)).exp(), diagnostics=diagnostics)

    def local_operators(self, tilt):
        return q1.kernel_t64.convolve(self.physical_operators(tilt))


class PostTransvectionModel:
    """Original uniform GL update followed by independent cheap transvections."""
    def __init__(self, data, map_record, *, updates=1):
        bits = sparse.validated(data, map_record)
        if type(updates) is not int or not 1 <= updates <= 8:
            raise ValueError('one through eight independent post-feedback transvections required')
        self.data, self.physical = data, q1.kernel_t64.authenticate(data)
        self.source_map_record, self.bits, self.updates = map_record, bits, updates
        self.width = 4*self.physical['windows']
        self.metadata = dict(schema='gl-plus-post-transvection-local-model-1',
            recurrence='y=x+A*a; next_a=R_post*(M*a+C*x)',
            distribution='uniform_gl_then_independent_transvections', physical_t=self.width,
            state_bits=bits, post_transvections=updates, lazy_weight=str(Q(1, 2**updates)),
            source_map_record=map_record, map_sha256=data['map_sha256'],
            zero_initial_state=True, final_flush=False, physical_updates_independent=True,
            sampler='u uniform nonzero; v uniform in u-perp including zero; R=I+u*v^T',
            matrix_coordinates='unchanged zero/arbitrary/uniform-density/expansion-weight classes',
            whole_code_certificate=False)

    def physical_operators(self, tilt):
        tilt = Q(tilt)
        if tilt <= 0:
            raise ValueError('positive rational tilt required')
        z = (-aq(tilt)).exp()
        old = q1.kernel_t64.sparse_kernel.outward_at_z(self.physical, z)
        if self.physical.get('birth_density', 'classes') == 'capped':
            old = q1.kernel_t64.kernel_birth_density.refine_local(self.physical, old, z, Q(1, 2))
        size = old[0].nrows()
        alpha, beta = aq(Q(1, 2**self.updates)), aq(1-Q(1, 2**self.updates))
        post = arb_mat(size, size)
        post[0, 0], post[2, 2] = 1, 1
        for index in (1, *range(3, size)):
            post[index, index], post[index, 2] = alpha, beta
        return [q1.rounded(operator*post) for operator in old]

    def local_operators(self, tilt):
        return q1.kernel_t64.convolve(self.physical_operators(tilt))


def prepare(data, map_record, *, exact_through=3):
    """Reuse a fresh authenticated t64 map preparation for the new recurrence."""
    sparse.validated(data, map_record)
    return PostFeedbackModel(q1.kernel_t64.authenticate(data), exact_through=exact_through,
                             source_map_record=map_record)


def _authenticate_model(model):
    """Bind the actual prepared maps to the receipt, including direct constructors."""
    if model.source_map_record is None:
        raise ValueError('authenticated source-map record required for numerical receipts')
    if isinstance(model, PostFeedbackModel):
        wrapped = q1.kernel_t64.wrap(model.physical)
        if _digest(model.counts) != model.metadata['return_histogram_sha256']:
            raise ValueError('prepared return histograms changed')
    elif isinstance(model, PostTransvectionModel):
        wrapped = model.data
        if wrapped.get('physical_data') is not model.physical:
            raise ValueError('post-transvection model no longer refers to its prepared physical maps')
    else:
        raise ValueError('prepared post-feedback candidate required')
    bits = sparse.validated(wrapped, model.source_map_record)
    if (bits != model.bits or wrapped['map_sha256'] != model.metadata['map_sha256']
            or model.width != wrapped['physical_step_bits']):
        raise ValueError('candidate metadata no longer matches the prepared maps')


def run_tail(model, *, K, occupancies, tilts, precision=192):
    """Fresh q>=3 comparison; returns a scoped receipt without writing files."""
    if not isinstance(model, (PostFeedbackModel, PostTransvectionModel)) or model.width != 64:
        raise ValueError('prepared physical t64 candidate required')
    if type(K) is not int or K <= 0 or K % 4096:
        raise ValueError('K must be a positive multiple of 4096')
    geometry = q1.Geometry(K//128, 64, 128)
    qs = tuple(occupancies)
    if (not qs or any(type(q) is not int or not 3 <= q <= geometry.group_count for q in qs)
            or tuple(sorted(set(qs))) != qs or type(precision) is not int or precision < 192
            or model.source_map_record is None):
        raise ValueError('ordered occupancies in3..L, precision>=192, and authenticated source maps required')
    tilts = tuple(map(Q, tilts))
    if not tilts or min(tilts) <= 0 or len(set(tilts)) != len(tilts):
        raise ValueError('distinct positive rational tilts required')
    _authenticate_model(model)
    sources = sparse.source_snapshot(model.source_map_record)
    beta = UniformInputEnvelope(16, 8, 4, 4).beta
    best, choices = dict.fromkeys(qs), dict.fromkeys(qs)
    previous = ctx.prec
    try:
        ctx.prec = precision
        for tilt in tilts:
            local = model.local_operators(tilt)
            regional = placement_power(local, epochs=geometry.macros_per_region,
                                        windows=32, maximum_groups=max(qs))
            factor = (aq(tilt)*(geometry.N//10)).exp()
            for q in qs:
                matrix = regional_uniform(regional, q)**64
                moment = sum((matrix[0, j] for j in range(matrix.ncols())), arb(0))
                value = up(comb(geometry.group_count, q)*aq(beta)**q*factor*moment)
                if not value.is_finite() or not value > 0:
                    raise ArithmeticError('finite positive occupancy endpoint required')
                if best[q] is None or value < best[q]:
                    best[q], choices[q] = value, str(tilt)
            if sources != sparse.source_snapshot(model.source_map_record):
                raise RuntimeError('loaded mathematical sources changed during comparison')
            _authenticate_model(model)
        total = up(sum(best.values(), arb(0)))
        return dict(schema='post-feedback-rs16-tail-comparison-1', K=K, N=geometry.N,
            geometry=asdict(geometry), threshold=geometry.N//10, distance='1/10',
            model=model.metadata, state_bits=model.bits, occupancy_covered=list(qs),
            physical_steps=geometry.N//64, macros_per_region=geometry.macros_per_region,
            placement_backend='binary', zero_initial_state=True, final_flush=False,
            terminal='all state-envelope mass; no reset between physical steps or regions',
            precision=precision, beta=str(beta), tilts=list(map(str, tilts)),
            source_sha256=sources, occupancy_choices=choices,
            occupancy_uppers={str(q): q1.endpoint(v) for q, v in best.items()},
            occupancy_margin_bits={str(q): str(-v.log()/arb(2).log()) for q, v in best.items()},
            union_upper=q1.endpoint(total), margin_bits=str(-total.log()/arb(2).log()),
            fresh_computation=True, whole_code_certificate=False,
            scope='New inner/setup distribution with the same RS16 outer maps and routing. '
                  'Only listed q>=3 first-moment contributions are covered; no q1/q2, '
                  'whole-code, performance, or seeded guarantee follows.')
    finally:
        ctx.prec = previous
