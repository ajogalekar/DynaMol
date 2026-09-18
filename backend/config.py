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
MAX_DURATION_PS = float(os.environ.get("DYNAMOL_MAX_DURATION_PS", "10000"))
if not 0 < MAX_DURATION_PS <= 1_000_000:
    raise ValueError("DYNAMOL_MAX_DURATION_PS must be between 0 and 1000000; the default is 10000 (10 ns).")
MAX_STEPS = int(os.environ.get("DYNAMOL_MAX_STEPS", "10000000"))
if not 1 <= MAX_STEPS <= 2_000_000_000:
    raise ValueError("DYNAMOL_MAX_STEPS must be between 1 and 2000000000; the default is 10000000.")
MAX_FRAMES = 10_000
MAX_COORD_BYTES = 256 * 1024 * 1024
CPU_THREADS = max(1, min(4, int(os.environ.get("DYNAMOL_CPU_THREADS", "2"))))
# Above this prepared-atom count, GBn2 implicit solvent (NoCutoff, O(N^2)) is slow
# on CPU; preparation recommends explicit TIP3P/PME (O(N), also more accurate) for
# such protein-only systems. Advisory only: implicit stays selectable.
RECOMMEND_EXPLICIT_ATOMS = int(os.environ.get("DYNAMOL_RECOMMEND_EXPLICIT_ATOMS", "4000"))
for directory in (DATASETS_DIR, JOBS_DIR):
    directory.mkdir(parents=True, exist_ok=True)


def is_app_bundle() -> bool:
    """True when running inside the packaged macOS .app rather than a source checkout."""
    return any(part.endswith(".app") for part in ROOT.parts)


def setup_guidance(component: str, source_command: str, detail: str = "") -> str:
    """User-facing recovery text for a missing bundled engine/runtime.

    Installed-app users can never act on source-checkout instructions, so tailor the
    message: reinstall for the .app (first launch re-unpacks the engines), the install
    script for a source checkout. Keeps the diagnostic detail only where it is useful.
    """
    if is_app_bundle():
        return (f"{component} isn't set up in this copy of DynaMol. Quit DynaMol, download "
                "the latest version, and drag it into Applications to replace this copy, then "
                "reopen it. The first launch finishes setting up the bundled engines, which "
                "can take a minute or two.")
    prefix = detail.strip().rstrip(".") + ". " if detail.strip() else ""
    return f"{prefix}Run {source_command} from your DynaMol source checkout, then restart the app."
