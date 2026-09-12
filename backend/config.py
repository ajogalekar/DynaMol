"""Local storage and deliberately bounded v0.1 resource limits."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = Path(os.environ.get("DYNAMOL_DATA_DIR", str(ROOT / "data"))).resolve()
DATASETS_DIR = DATA_ROOT / "datasets"
JOBS_DIR = DATA_ROOT / "jobs"
MAX_UPLOAD_BYTES = 250 * 1024 * 1024
MAX_ATOMS = int(os.environ.get("DYNAMOL_MAX_ATOMS", "100000"))
if not 1 <= MAX_ATOMS <= 200_000:
    raise ValueError("DYNAMOL_MAX_ATOMS must be between 1 and 200000; the default is 100000.")
MAX_FRAMES = 10_000
MAX_COORD_BYTES = 256 * 1024 * 1024
CPU_THREADS = max(1, min(4, int(os.environ.get("DYNAMOL_CPU_THREADS", "2"))))
for directory in (DATASETS_DIR, JOBS_DIR):
    directory.mkdir(parents=True, exist_ok=True)
