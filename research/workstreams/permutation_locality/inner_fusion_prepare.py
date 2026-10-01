"""Extract the authenticated prior winner for isolated implementation A/B tests."""
import argparse
import hashlib
from pathlib import Path

PRIOR_SHA256 = "1f796772fc380519fecec1dd45418fcbaeab701f56e9fcd71d3bdd9faaaeb733"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    raw = Path(__file__).with_name("inner_packet_probe.cpp").read_bytes()
    if hashlib.sha256(raw).hexdigest() != PRIOR_SHA256:
        raise ArithmeticError("prior winner source changed")
    source = raw.decode().replace("\r\n", "\n")
    marker = "\ntemplate<int Mode> static void experiment("
    if source.count(marker) != 1:
        raise ArithmeticError("ambiguous reference extraction")
    args.output.write_text("// Authenticated packet winner: " + PRIOR_SHA256 +
                           "\n#pragma once\n" + source.split(marker)[0] + "\n}\n",
                           encoding="utf-8", newline="\n")
    print("Retained packet source SHA256", PRIOR_SHA256)


if __name__ == "__main__":
    main()
