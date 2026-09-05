#!/usr/bin/env python3
"""Load the local NumPy supermodel pickle through an explicit type allowlist."""

from __future__ import annotations

import json
import math
import os
import pickle
import stat
import sys
from pathlib import Path
from typing import Any

import numpy as np


MAX_PICKLE_BYTES = 256 * 1024 * 1024
MAX_ARRAY_ELEMENTS = 100_000
MAX_DEPTH = 100
MAX_CONVERTED_ITEMS = 6_000_000
MAX_CONVERTED_JSON_BYTES = 64 * 1024 * 1024

ALLOWED_GLOBALS = {
    ("numpy", "dtype"),
    ("numpy", "ndarray"),
    ("numpy.core.multiarray", "_reconstruct"),
    ("numpy.core.multiarray", "scalar"),
    ("numpy._core.multiarray", "_reconstruct"),
    ("numpy._core.multiarray", "scalar"),
}


class RestrictedUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str) -> Any:
        if (module, name) not in ALLOWED_GLOBALS:
            raise pickle.UnpicklingError(f"disallowed pickle global: {module}.{name}")
        return super().find_class(module, name)


class ConversionBudget:
    def __init__(self) -> None:
        self.items = 0
        self.json_bytes = 0
        self.active_containers: set[int] = set()
        self.seen_containers: dict[int, Any] = {}

    def charge(self, *, items: int = 0, json_bytes: int = 0) -> None:
        self.items += items
        self.json_bytes += json_bytes
        if self.items > MAX_CONVERTED_ITEMS:
            raise ValueError("converted data exceeds the item limit")
        if self.json_bytes > MAX_CONVERTED_JSON_BYTES:
            raise ValueError("converted data exceeds the JSON-size limit")

    def enter_container(self, value: Any, *, remember: bool = True) -> int:
        identity = id(value)
        if identity in self.active_containers:
            raise ValueError("cyclic container references are not allowed")
        if remember and identity in self.seen_containers:
            raise ValueError("aliased container references are not allowed")
        self.active_containers.add(identity)
        if remember:
            # Retain the reference so Python cannot reuse its id during conversion.
            self.seen_containers[identity] = value
        self.charge(items=1, json_bytes=2)
        return identity

    def leave_container(self, identity: int) -> None:
        self.active_containers.remove(identity)

    def charge_scalar(self, value: Any) -> None:
        encoded = json.dumps(value, allow_nan=False, separators=(",", ":"))
        self.charge(items=1, json_bytes=len(encoded.encode("utf-8")))


def convert(
    value: Any,
    budget: ConversionBudget,
    depth: int = 0,
    *,
    track_aliases: bool = True,
) -> Any:
    if depth > MAX_DEPTH:
        raise ValueError("pickle nesting exceeds the conversion limit")

    if value is None or isinstance(value, (str, bool, int)):
        budget.charge_scalar(value)
        return value
    if isinstance(value, float):
        converted = value if math.isfinite(value) else None
        budget.charge_scalar(converted)
        return converted
    if isinstance(value, complex):
        return convert({"real": value.real, "imag": value.imag}, budget, depth + 1)
    if isinstance(value, np.generic):
        return convert(value.item(), budget, depth + 1)
    if isinstance(value, np.ndarray):
        if value.dtype.hasobject:
            raise ValueError("object-dtype NumPy arrays are not allowed")
        remaining_items = MAX_CONVERTED_ITEMS - budget.items
        if value.size > MAX_ARRAY_ELEMENTS or value.size > remaining_items:
            raise ValueError("NumPy array exceeds the element limit")
        identity = budget.enter_container(value)
        try:
            return convert(
                value.tolist(),
                budget,
                depth + 1,
                track_aliases=False,
            )
        finally:
            budget.leave_container(identity)
    if isinstance(value, dict):
        identity = budget.enter_container(value, remember=track_aliases)
        try:
            converted: dict[str, Any] = {}
            for key, item in value.items():
                converted_key = str(key)
                budget.charge_scalar(converted_key)
                budget.charge(json_bytes=1)
                converted[converted_key] = convert(
                    item,
                    budget,
                    depth + 1,
                    track_aliases=track_aliases,
                )
            return converted
        finally:
            budget.leave_container(identity)
    if isinstance(value, (list, tuple, set)):
        identity = budget.enter_container(value, remember=track_aliases)
        try:
            converted_items = []
            for item in value:
                budget.charge(json_bytes=1)
                converted_items.append(
                    convert(
                        item,
                        budget,
                        depth + 1,
                        track_aliases=track_aliases,
                    )
                )
            return converted_items
        finally:
            budget.leave_container(identity)
    raise ValueError(f"unsupported value type: {type(value).__module__}.{type(value).__name__}")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: load-supermodel-json.py <supermodel_data.pkl>")

    source = Path(sys.argv[1])
    source_lstat = source.lstat()
    if stat.S_ISLNK(source_lstat.st_mode) or not stat.S_ISREG(source_lstat.st_mode):
        raise ValueError("pickle source must be a regular, non-symlink file")
    if source_lstat.st_size > MAX_PICKLE_BYTES:
        raise ValueError("pickle source exceeds the 256 MiB limit")

    descriptor = os.open(source, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        with os.fdopen(descriptor, "rb") as handle:
            loaded = RestrictedUnpickler(handle).load()
    except Exception:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise

    json.dump(
        convert(loaded, ConversionBudget()),
        sys.stdout,
        allow_nan=False,
        separators=(",", ":"),
    )


if __name__ == "__main__":
    main()
