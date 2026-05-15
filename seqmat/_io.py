"""File I/O helpers: pickle and JSON."""
from __future__ import annotations

import json
import pickle
from pathlib import Path
from typing import Any, Union


def dump_pickle(path: Union[str, Path], data: Any) -> None:
    """Save data to a pickle file."""
    with open(path, "wb") as f:
        pickle.dump(data, f)


def unload_pickle(path: Union[str, Path]) -> Any:
    """Load data from a pickle file."""
    with open(path, "rb") as f:
        return pickle.load(f)


def dump_json(path: Union[str, Path], data: Any) -> None:
    """Save data to a JSON file (indent=2)."""
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def unload_json(path: Union[str, Path]) -> Any:
    """Load data from a JSON file."""
    with open(path, "r") as f:
        return json.load(f)


__all__ = ["dump_pickle", "unload_pickle", "dump_json", "unload_json"]
