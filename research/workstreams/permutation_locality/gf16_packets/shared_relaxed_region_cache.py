"""Worker-scoped exact regional-polynomial cache for search only.

No saved numerical record is imported. A cache miss freshly constructs the
outward polynomial from the actual inner. A hit reuses it only for the same
retained data object, precision, output tilt, and local refinements. Cell
counts and outer bounds still use the original outward implementation.
Final certificate replay does not install this context.
"""
from collections import OrderedDict
from contextlib import contextmanager
from fractions import Fraction as Q
from types import FunctionType

from flint import arb_mat, ctx
import regional_count as regional


LOCAL_SELECTORS = (
    'regional_exact_zero', 'regional_joint_return_through', 'regional_lazy_density_through',
    'regional_feedback_classes_from', 'regional_feedback_classes_through',
    'regional_feedback_uniform_classes', 'regional_feedback_uniform_replace',
)
COUNT_SELECTORS = (
    'regional_count_parts', 'regional_direct_counts', 'regional_tilted_atom',
    'regional_fine_tilts', 'regional_tilted_variance',
)


def local_key(model, witness):
    """Fail closed if a future regional selector is not classified here."""
    if not isinstance(witness, dict) or not isinstance(model.data, dict):
        raise ValueError('actual data object and rational witness required')
    if any(key.startswith('regional_') and key not in LOCAL_SELECTORS+COUNT_SELECTORS
           for key in witness):
        raise ValueError('unknown regional selector: cache classification required')
    windows = model.data.get('windows')
    if type(windows) is not int or windows < 1:
        raise ValueError('positive actual packet geometry required')
    parameters = witness.get('parameters')
    if not isinstance(parameters, list) or len(parameters) != 3 or Q(parameters[0]) <= 0:
        raise ValueError('positive output tilt and outer duals required')
    for key in ('regional_exact_zero', 'regional_feedback_uniform_classes',
                'regional_feedback_uniform_replace'):
        if type(witness.get(key, False)) is not bool:
            raise ValueError('boolean local refinement selector required')
    for key, limit in (('regional_joint_return_through', min(4, windows)),
                       ('regional_lazy_density_through', windows)):
        if key in witness and (type(witness[key]) is not int or not 0 <= witness[key] <= limit):
            raise ValueError('valid local occupancy cutoff required')
    regional.feedback_class_interval(witness, windows)
    selectors = tuple((key, key in witness, witness.get(key)) for key in LOCAL_SELECTORS)
    return id(model.data), ctx.prec, regional.sc.G, Q(parameters[0]), selectors


def matrix_copies(polynomial):
    # Downstream code currently multiplies without mutation. Copies also
    # make cache correctness robust to later in-place matrix operations.
    return [arb_mat(matrix) for matrix in polynomial]


class RegionalCache:
    def __init__(self, maximum_entries=4):
        if type(maximum_entries) is not int or not 1 <= maximum_entries <= 32:
            raise ValueError('cache capacity in 1..32 required')
        self.maximum_entries = maximum_entries
        self.entries = OrderedDict()
        self.hits = 0
        self.misses = 0

    def polynomial(self, model, witness):
        key = local_key(model, witness)
        saved = self.entries.get(key)
        if saved is not None:
            if saved[0] is not model.data:
                raise ArithmeticError('cached data identity mismatch')
            self.entries.move_to_end(key)
            self.hits += 1
            return matrix_copies(saved[1])
        precision = ctx.prec
        result = regional.placement(regional.local_operators(model, witness))
        if ctx.prec != precision or local_key(model, witness) != key:
            raise ArithmeticError('local precision, data, or selectors changed during construction')
        if len(result) != regional.sc.G+1 or not result:
            raise ArithmeticError('complete regional polynomial required')
        size = result[0].nrows()
        if size < 1 or any(matrix.nrows() != size or matrix.ncols() != size for matrix in result):
            raise ArithmeticError('matching square coefficient matrices required')
        if any(not matrix[i,j].is_finite() or matrix[i,j] < 0
               for matrix in result for i in range(size) for j in range(size)):
            raise ArithmeticError('finite nonnegative outward polynomial required')
        # Retain the object itself, so its id cannot be recycled while cached.
        self.entries[key] = model.data, tuple(result)
        while len(self.entries) > self.maximum_entries:
            self.entries.popitem(last=False)
        self.misses += 1
        return matrix_copies(result)

    def evaluate(self, original, model, cell, witness):
        """Run the unchanged outward body with a private polynomial provider.

        Cloning the function's globals avoids process-global monkeypatches
        of placement/local_operators and avoids duplicating proof formulas.
        Only those two providers are substituted for this invocation.
        """
        if (not isinstance(original, FunctionType)
                or not {'local_operators', 'placement'} <= set(original.__code__.co_names)):
            raise ValueError('expected original regional outward function required')
        polynomial = self.polynomial(model, witness)
        marker = object()
        calls = [0, 0]

        def local_provider(actual_model, actual_witness):
            if actual_model is not model or actual_witness is not witness or calls[0]:
                raise ArithmeticError('unexpected local-operator call in cached outward body')
            calls[0] += 1
            return marker

        def placement_provider(local):
            if local is not marker or calls[0] != 1 or calls[1]:
                raise ArithmeticError('unexpected placement call in cached outward body')
            calls[1] += 1
            return polynomial

        namespace = dict(original.__globals__, local_operators=local_provider,
                         placement=placement_provider)
        body = FunctionType(original.__code__, namespace, original.__name__,
                            original.__defaults__, original.__closure__)
        body.__kwdefaults__ = original.__kwdefaults__
        result = body(model, cell, witness)
        if calls != [1, 1]:
            raise ArithmeticError('outward body did not consume its complete cached polynomial')
        return result


@contextmanager
def install(maximum_entries=4):
    """Temporarily enable a cache in one search worker; always restore it."""
    original = regional.outward
    if getattr(original, '_shared_search_region_cache', False):
        raise ValueError('nested regional cache contexts are not supported')
    cache = RegionalCache(maximum_entries)

    def cached(model, cell, witness):
        return cache.evaluate(original, model, cell, witness)

    cached._shared_search_region_cache = True
    regional.outward = cached
    try:
        yield cache
    finally:
        regional.outward = original
