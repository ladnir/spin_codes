"""Map-bound regional refinements without legacy 19-state entry points.

The retained regional proposal/replay bodies are rebound to this module's
local-operator factory. No imported module or process-wide global is patched.
Fresh map preparation supplies exact profiles. Joint-return integer counts
are rebuilt per model; numerical caps are rebuilt at the current precision.
"""
from fractions import Fraction as Q
from types import FunctionType

from flint import arb
import kernel_maps
import regional_count as retained
import occupancy_birth_classes
import return_moment
import lazy_density
import single_packet
import fiber_density

aq = kernel_maps.aq


def local_operators(model, witness, *, _proposal_scratch=None):
    data = model.data
    images = kernel_maps.authenticate(data)
    lam = Q(witness['parameters'][0]); z = (-aq(lam)).exp()
    if lam <= 0:
        raise ValueError('positive output tilt required')
    feedback = retained.feedback_class_interval(witness, data['windows'])
    exact_zero = witness.get('regional_exact_zero', False)
    if type(exact_zero) is not bool:
        raise ValueError('boolean exact-zero option required')
    birth_density = data.get('birth_density', 'classes')
    activity = Q(witness.get('s16_birth_density_activity', '1/2'))
    if not 0 <= activity <= 1:
        raise ValueError('valid birth-row selection activity required')
    key = (data['map_sha256'], data['distribution'], data.get('updates'), lam, exact_zero,
           birth_density, activity, data.get('profile_partition', False))
    saved = None if _proposal_scratch is None else _proposal_scratch.get('s16_local_base')
    if data['distribution'] == 'uniform_gl':
        import sparse_kernel
        # With alpha=0 the joint/lazy/fiber refinements have no contribution.
        # Do not invoke the legacy helpers, which require an updates field.
        if saved is None or saved[0] != key:
            if data.get('profile_partition', False):
                import kernel_profiles
                local = kernel_profiles.outward_local_at_z(data, z, activity)
            else:
                local = sparse_kernel.outward_at_z(data, z)
                if birth_density == 'capped':
                    import kernel_birth_density
                    local = kernel_birth_density.refine_local(data, local, z, activity)
            saved = (key, local)
            if _proposal_scratch is not None:
                _proposal_scratch['s16_local_base'] = saved
        return saved[1]
    if saved is not None and saved[0] == key:
        local = saved[1]
    else:
        local = occupancy_birth_classes.outward_at_z(data, z)
        if exact_zero:
            local = occupancy_birth_classes.refine_zero(data, local, z)
        if _proposal_scratch is not None:
            _proposal_scratch['s16_local_base'] = (key, local)
    if 'regional_joint_return_through' in witness:
        through = witness['regional_joint_return_through']
        if type(through) is not int or not 0 <= through <= min(4, data['windows']):
            raise ValueError('joint-return occupancy must be in 0..min(4,W)')
        key = (data['map_sha256'], through)
        cache = getattr(model, 's16_joint_return_cache', None)
        if cache is None:
            cache = model.s16_joint_return_cache = {}
        if key not in cache:
            cache[key] = return_moment.census(images, data['columns'], data['bits'], through)
        local = return_moment.refine_class_returns(data, local, cache[key], z)
    feedback_source = local
    if 'regional_lazy_density_through' in witness:
        through = witness['regional_lazy_density_through']
        if type(through) is not int or not 0 <= through <= data['windows']:
            raise ValueError('lazy-density occupancy must be in 0..W')
        key = (data['map_sha256'], lam, through)
        saved = None if _proposal_scratch is None else _proposal_scratch.get('s16_lazy_caps')
        if saved is None or saved[0] != key:
            single = single_packet.census(images, data['columns'], z) if through else None
            saved = (key, lazy_density.density_caps(data, z, single))
            if _proposal_scratch is not None:
                _proposal_scratch['s16_lazy_caps'] = saved
        caps = saved[1]
        local = lazy_density.candidate(data, local, z, caps, through)
    if feedback is not None:
        attached = getattr(model, 's16_fiber_density_data', None)
        if (attached is None or attached['map_sha256'] != data['map_sha256']
                or attached.get('updates') != data.get('updates')):
            attached = model.s16_fiber_density_data = fiber_density.attach(data)
        uniform = witness.get('regional_feedback_uniform_classes', False)
        replace = witness.get('regional_feedback_uniform_replace', False)
        options = {'include_uniform': True} if uniform and not replace else {}
        local = fiber_density.candidate(attached, local, z, *feedback, allocation='classes', **options)
        if replace:
            local = fiber_density.replace_uniform_classes(attached, local, feedback_source, z, *feedback)
    return local


# Retain all count-MGF, interval validation and outward arithmetic unchanged.
# Only the local-kernel dependency changes, in this private globals dictionary.
_bindings = dict(retained.__dict__, local_operators=local_operators)
_propose = FunctionType(retained.propose.__code__, _bindings, 'propose', retained.propose.__defaults__)
outward = FunctionType(retained.outward.__code__, _bindings, 'outward', retained.outward.__defaults__)
prepare_witness = retained.prepare_witness


def propose(model, cell, witness):
    # The potential only selects among valid complete rows. Record its input
    # so outward replay uses the same proposal activity without trusting caps.
    if model.data.get('birth_density', 'classes') == 'capped':
        weights, _ = model.weights(cell, model.tilt)
        witness = dict(witness, s16_birth_density_activity=str(weights[1]/sum(weights)))
    return _propose(model, cell, witness)
