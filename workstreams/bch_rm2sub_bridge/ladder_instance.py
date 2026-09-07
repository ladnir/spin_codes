"""Extend instance identities through m24 without modifying frozen m<=20 code."""
import hashlib
import json

import certificate_search_core as core


def instance(m):
    core.require(type(m) is int and 16 <= m <= 24, 'Ladder supports message exponents 16..24')
    original = core.instance('t128_s19', min(m, 20))
    if m <= 20:
        return original
    value = dict(original)
    del value['fingerprint']
    rows = 1 << (m-7)
    core.require(rows % value['step_bits'] == 0, 'Incomplete inner epochs')
    value.update(message_exponent=m, rows=rows, cutoff=256*rows//10)
    value['fingerprint'] = hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()
    return value
