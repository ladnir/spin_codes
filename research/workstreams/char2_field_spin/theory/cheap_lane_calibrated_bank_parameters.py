"""Explicit rational proof parameters for the frozen complete r2048 bank.

Set private Z/ALPHA before rebuilding all local caps, powers, corrections,
source factors and marker witnesses. Every actual profile and terminal is
unchanged. Delegate other arguments to the bounded frozen bank. FLOAT only;
stdout only, no generated files or automatic parameter search.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
BANK_PIN = "e9ff3adda8aa363383386b8e8665b40da5f1075a8e43ddf8ffc1854ddea60776"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z", type=F, required=True)
    parser.add_argument("--alpha", type=F, required=True)
    args, remaining = parser.parse_known_args()
    if not F(1,10) <= args.z <= F(9,10) or not F(1,64) <= args.alpha <= F(1,4):
        parser.error("require 1/10<=z<=9/10 and 1/64<=alpha<=1/4")
    if max(n.bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator)) > 64:
        parser.error("rational witnesses limited to 64 bits")
    path = Path(__file__).with_name("cheap_lane_rho_calibrated_bank.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != BANK_PIN:
        raise RuntimeError("frozen calibrated-bank hash mismatch")
    spec = importlib.util.spec_from_file_location("parameter_calibrated_bank", path)
    bank = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bank)
    bank.Z, bank.ALPHA = args.z, args.alpha
    print(json.dumps({"phase":"calibrated_bank_parameter_provenance",
        "z":str(args.z),"alpha":str(args.alpha),"bank_pin":BANK_PIN,
        "wrapper_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}), flush=True)
    previous = sys.argv
    try:
        sys.argv = [previous[0], *remaining]
        bank.main()
    finally:
        sys.argv = previous


if __name__ == "__main__":
    main()
