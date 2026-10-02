"""Refresh envelope requiring nonzero feedback for a return to zero."""
import refresh_kernel as base

actual=base.actual
prepare=base.prepare


def floating(data,probabilities,tilt):
    return base.floating(data,probabilities,tilt,feedback_aware=True)


def outward_at_z(data,probabilities,z):
    return base.outward_at_z(data,probabilities,z,feedback_aware=True)


def outward(data,probabilities,tilt):
    return base.outward(data,probabilities,tilt,feedback_aware=True)
