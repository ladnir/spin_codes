"""Emit one candidate leaf for the compact Lean checker resource pilot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from emit_scalable_dense_geometry_pilot import qinput


def emit(candidate: Path, destination: Path, box_id: int) -> dict:
    payload = json.loads(candidate.read_text(encoding="utf-8"))
    assert payload["status"] == "candidate_complete"
    box = payload["boxes"][box_id]
    witness = payload["witnesses"][str(box["witness"])]
    slope, intercept = box["outer_line"]
    a0, a1 = box["alpha"]
    x0, x1 = box["row_density"]
    vector = ", ".join(qinput(value) for value in witness["vector"])
    source = f"""import SpinCodes.UniformMixer.ScalableDensePackedLeaf

namespace Spin.UniformMixer.ScalableDensePackedPilot

open Spin.Structured.DenseOccupationFixed

set_option maxHeartbeats 1000000

def line : DenseLine := ⟨{qinput(slope)}, {qinput(intercept)}⟩
def witness : DenseWitness := ⟨{qinput(witness['p'])},
  {qinput(witness['y'])}, {qinput(witness['rho'])},
  {qinput(witness['z'])}, ![{vector}]⟩
def box : DenseBox := ⟨{qinput(a0)}, {qinput(a1)},
  {qinput(x0)}, {qinput(x1)}⟩

theorem checked : densePackedLeafCheck line witness box = true := by
  decide +kernel

theorem claim :
    DenseClaim (DenseTarget (1 / 10000000) 10000000) box.rect :=
  densePackedLeafCheck_sound line witness box checked

#print axioms claim

end Spin.UniformMixer.ScalableDensePackedPilot
"""
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(source.encode("utf-8"))
    return {"box_id": box_id, "witness": box["witness"],
            "source_bytes": len(source.encode("utf-8")),
            "destination": str(destination)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--box-id", type=int, default=0)
    args = parser.parse_args()
    print(json.dumps(emit(args.candidate, args.destination, args.box_id)))


if __name__ == "__main__":
    main()
