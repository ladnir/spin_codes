"""Reuse the exact exhaustive state-law tests on the tighter envelope."""
import unittest
import test_general_occupancy as tests
import tightened_occupancy as tightened

if __name__=='__main__':
    tests.model=tightened
    unittest.main(module=tests)
