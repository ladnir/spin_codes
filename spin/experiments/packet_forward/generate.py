"""Compatibility entry point for the production generator."""
import runpy
from pathlib import Path
runpy.run_path(str(Path(__file__).resolve().parents[2]/"tools/generate_packet_forward.py"),run_name="__main__")
