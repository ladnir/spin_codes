"""Frozen FLOAT-proposed literals and a deterministic219-class assignment.

Scope: one dense group plus sixteen two-row groups, ordered source(8,0),
mu16 only. Fresh directed arithmetic is required before certification.
The mask partition is selected before any setup randomness is observed.
"""
from fractions import Fraction as F

from hd_xor_paired_gate import CONSTRUCTION


PROPOSAL_SHA256 = '0abded252784870277a1b285d93e522b9bdb717288a42448ce49469097fc1fc5'
SOURCE_MATCHING_COUNTS = (8, 0)
MU_RANGE = (16, 16)
SPARSE_SUPPORT_FLOOR = 544
DENSE_COLUMN_FLOOR = 17

# label, rational z, rational log2(t), assigned fixed dense row masks
SAVED_CLASSES = (
    ('a8_fit', '8970988817957597/9007199254740992', '-2266954815/33554432', (
        25,26,28,41,42,44,49,50,52,56,57,58,60,61,62,63,73,74,76,81,82,84,
        88,89,90,92,93,94,95,105,106,108,109,110,111,113,114,116,117,118,119,
        120,121,122,123,124,125,126,127,137,138,140,145,146,148,152,153,154,
        156,157,158,159,169,170,172,173,174,175,177,178,180,181,182,183,184,
        185,186,187,188,189,190,191,207,215,219,221,222,223,231,235,237,238,
        239,243,245,246,247,249,250,251,252,253,254,255)),
    ('a3_mask7_fit', '1121881298265941/1125899906842624', '-66751297/1048576', (
        7,11,13,14,15,19,21,22,23,27,35,37,38,39,43,51,67,69,70,71,75,83,
        99,131,133,134,135,139,147,163,195)),
    ('a4_mask240_fit', '8973939342378351/9007199254740992', '-2386334719/33554432', (
        29,30,31,45,46,47,53,54,55,59,77,78,79,85,86,87,91,97,98,100,101,
        102,103,104,107,112,115,141,142,143,149,150,151,155,161,162,164,165,
        166,167,168,171,176,179,193,194,196,197,198,199,200,201,202,203,204,
        205,206,208,209,210,211,212,213,214,216,217,218,220,224,225,226,227,
        228,229,230,232,233,234,236,240,241,242,244,248)),
)


def checked_classes():
    seen, result = set(), []
    for label, mark, potential, masks in SAVED_CLASSES:
        z, logp = F(mark), F(potential)
        if not 0 < z < 1 or logp > 0 or set(masks) & seen or len(set(masks)) != len(masks):
            raise ValueError('Invalid literal or overlapping fixed dense classes')
        if any(type(mask) is not int or not 0 < mask < 256 or mask.bit_count() < 3 for mask in masks):
            raise ValueError('Only arity3..8 dense row sets are in scope')
        seen.update(masks)
        result.append((label, z, logp, masks))
    if seen != {mask for mask in range(256) if mask.bit_count() >= 3}:
        raise ValueError('The frozen assignment must cover exactly219 dense row sets')
    return tuple(result)
