"""Optional ctypes loader for packet-group diagnostic hot kernels."""

from __future__ import annotations

import ctypes
import os
from pathlib import Path

import numpy as np


_DOUBLE_ARRAY = np.ctypeslib.ndpointer(
    dtype=np.float64, ndim=1, flags=("C_CONTIGUOUS", "ALIGNED")
)
_INT32_ARRAY = np.ctypeslib.ndpointer(
    dtype=np.int32, ndim=1, flags=("C_CONTIGUOUS", "ALIGNED")
)
_LOADED_SHARED_DRIVE_APPLY = None
_LOADED_POINT_CAPS = None
_LOAD_ATTEMPTED = False


def _candidate_paths() -> list[Path]:
    configured = os.environ.get("PACKET_GROUP_NATIVE_LIB")
    directory = Path(__file__).resolve().parent
    candidates = []
    if configured:
        candidates.append(Path(configured))
    candidates.extend(
        [
            directory / "libpacket_group_native.so",
            directory / "packet_group_native.dll",
            directory / "libpacket_group_native.dylib",
        ]
    )
    return candidates


def load_shared_drive_apply():
    """Return the native function when present, otherwise ``None``."""

    global _LOADED_SHARED_DRIVE_APPLY, _LOADED_POINT_CAPS, _LOAD_ATTEMPTED
    if os.environ.get("PACKET_GROUP_DISABLE_NATIVE") == "1":
        return None
    if _LOAD_ATTEMPTED:
        return _LOADED_SHARED_DRIVE_APPLY
    _LOAD_ATTEMPTED = True
    for path in _candidate_paths():
        if not path.is_file():
            continue
        library = ctypes.CDLL(str(path))
        function = library.packet_group_shared_drive_apply
        function.argtypes = [
            _DOUBLE_ARRAY,
            _DOUBLE_ARRAY,
            _INT32_ARRAY,
            _DOUBLE_ARRAY,
            _DOUBLE_ARRAY,
            _DOUBLE_ARRAY,
            _INT32_ARRAY,
            _DOUBLE_ARRAY,
        ]
        function.restype = ctypes.c_int
        # Keep the library alive through the function object.
        function._packet_group_library = library
        function._packet_group_library_path = str(path)
        _LOADED_SHARED_DRIVE_APPLY = function
        try:
            point_caps = library.packet_group_point_caps
        except AttributeError:
            # An older optional library may contain only the apply kernel.
            point_caps = None
        if point_caps is not None:
            point_caps.argtypes = [
                ctypes.c_int,
                _DOUBLE_ARRAY,
                _DOUBLE_ARRAY,
            ]
            point_caps.restype = ctypes.c_int
            point_caps._packet_group_library = library
            point_caps._packet_group_library_path = str(path)
            _LOADED_POINT_CAPS = point_caps
        return _LOADED_SHARED_DRIVE_APPLY
    return None


def load_point_caps():
    """Return the native point-cap builder when present, otherwise ``None``."""

    global _LOADED_POINT_CAPS
    load_shared_drive_apply()
    return _LOADED_POINT_CAPS
