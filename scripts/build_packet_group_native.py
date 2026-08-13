#!/usr/bin/env python3
"""Build the optional native packet-group diagnostic backend."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import platform
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--portable", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    directory = Path(__file__).resolve().parent
    source = directory / "packet_group_shared_drive_native.cpp"
    system = platform.system()
    if args.output is not None:
        output = args.output.resolve()
    elif system == "Windows":
        output = directory / "packet_group_native.dll"
    elif system == "Darwin":
        output = directory / "libpacket_group_native.dylib"
    else:
        output = directory / "libpacket_group_native.so"

    compiler_name = Path(args.compiler).name.lower()
    if compiler_name in {"cl", "cl.exe"}:
        command = [
            args.compiler,
            "/nologo",
            "/O2",
            "/DNDEBUG",
            "/std:c++17",
            "/LD",
            str(source),
            f"/Fe:{output}",
        ]
    else:
        command = [
            args.compiler,
            "-O3",
            "-DNDEBUG",
            "-std=c++17",
            "-fPIC",
            "-shared",
        ]
        if not args.portable:
            command.append("-march=native")
        command.extend([str(source), "-o", str(output)])

    output.parent.mkdir(parents=True, exist_ok=True)
    print(" ".join(command), flush=True)
    subprocess.run(command, check=True)
    print(output)


if __name__ == "__main__":
    main()
