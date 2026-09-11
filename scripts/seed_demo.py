#!/usr/bin/env python3
"""Install the bundled example into a fresh workspace; never replace user data."""
import os
import shutil
from pathlib import Path

root = Path(__file__).resolve().parents[1]
data = Path(os.environ.get("DYNAMOL_DATA_DIR", str(root / "data")))
for name in ("demo", "ubiquitin-start"):
    source = root / "examples" / name
    target = data / "datasets" / name
    if source.is_dir() and not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, target)
        print(f"Installed bundled example: {name}")
