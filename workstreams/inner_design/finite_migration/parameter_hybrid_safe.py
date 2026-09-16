"""Reject rounded endpoint proposals during hybrid IMT witness discovery."""
import math

import parameter_hybrid_dense as prior


class Mixed(prior.maps.dense.Dense):
    def moment(self,theta,tilt):
        if not 0 < theta < 1:
            return math.inf
        return super().moment(theta,tilt)


class Dense(prior.Dense):
    def __init__(self,*args):
        super().__init__(*args)
        self.mixed = Mixed(*args)

    def moment(self,theta,tilt):
        if not 0 < theta < 1:
            return math.inf
        return super().moment(theta,tilt)

    def evaluate(self,lower,upper):
        result = super().evaluate(lower,upper)
        if result is not None and not math.isfinite(result['own_log_bound']):
            raise ArithmeticError('No finite interior-probability witness')
        return result
