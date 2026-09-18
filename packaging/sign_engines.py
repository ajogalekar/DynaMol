#!/usr/bin/env python3
"""Sign every Mach-O inside the bundled engine archives so the app can be notarized.

Each engine archive (Contents/Resources/engines/<name>.tar.gz) is extracted, every Mach-O
inside is Developer ID signed (executables get the hardened-runtime entitlements; libraries
get hardened runtime only), nested .app/.framework bundles are sealed deepest-first, the
archive is repackaged, and its sha256 in manifest.json is updated to match.

Env: IDENTITY, ENT (entitlements plist), APP (path to DynaMol.app), WORK (scratch dir).
"""
import hashlib, json, os, shutil, subprocess, sys
from pathlib import Path

IDENTITY = os.environ["IDENTITY"]
ENT = os.environ["ENT"]
APP = Path(os.environ["APP"])
WORK = Path(os.environ["WORK"])
MANIFEST = APP / "Contents/Resources/manifest.json"
ENGDIR = APP / "Contents/Resources/engines"

MACHO_MAGIC = {b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf",
               b"\xfe\xed\xfa\xce", b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca"}

def is_macho(p: Path) -> bool:
    try:
        with open(p, "rb") as f:
            return f.read(4) in MACHO_MAGIC
    except OSError:
        return False

def codesign(path: Path, entitlements: bool) -> subprocess.CompletedProcess:
    cmd = ["/usr/bin/codesign", "--force", "--timestamp", "--options", "runtime", "--sign", IDENTITY]
    if entitlements:
        cmd += ["--entitlements", ENT]
    cmd.append(str(path))
    return subprocess.run(cmd, capture_output=True, text=True)

def sha256(path: Path) -> str:
    with open(path, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()

def sign_tree(root: Path) -> int:
    machos = []
    for dp, dn, fn in os.walk(root):
        for name in fn:
            p = Path(dp) / name
            if not p.is_symlink() and is_macho(p):
                machos.append(p)
    libs = exes = fail = 0
    for p in machos:
        ft = subprocess.run(["/usr/bin/file", "-b", str(p)], capture_output=True, text=True).stdout
        ent = "executable" in ft
        r = codesign(p, ent)
        if r.returncode == 0:
            exes += ent; libs += (not ent)
        else:
            fail += 1
            if fail <= 15:
                print(f"   FAIL {p}: {r.stderr.strip()[:90]}", flush=True)
    bundles = [Path(dp) / d for dp, dn, _ in os.walk(root) for d in dn if d.endswith((".app", ".framework"))]
    bundles.sort(key=lambda x: len(x.parts), reverse=True)  # deepest first
    bfail = 0
    for b in bundles:
        if codesign(b, entitlements=b.suffix == ".app").returncode:
            bfail += 1
    print(f"   signed libs={libs} exes={exes} bundles={len(bundles)} (bundle_fail={bfail}) macho_fail={fail}", flush=True)
    return fail + bfail

def main() -> int:
    manifest = json.loads(MANIFEST.read_text())
    total_fail = 0
    for eng, info in manifest["engines"].items():
        arc = info["archive"]
        src = ENGDIR / arc
        print(f"=== {eng}: {arc} ({src.stat().st_size/1e6:.0f} MB) ===", flush=True)
        ex = WORK / eng
        if ex.exists():
            shutil.rmtree(ex)
        ex.mkdir(parents=True)
        print("  extracting...", flush=True)
        subprocess.run(["/usr/bin/tar", "-xzf", str(src), "-C", str(ex)], check=True)
        subprocess.run(["/usr/bin/xattr", "-cr", str(ex)], capture_output=True)
        print("  signing Mach-O...", flush=True)
        total_fail += sign_tree(ex)
        print("  repackaging...", flush=True)
        newarc = WORK / arc
        subprocess.run(["/usr/bin/tar", "-czf", str(newarc), "-C", str(ex), "."], check=True)
        digest = sha256(newarc)
        shutil.move(str(newarc), str(src))          # replace archive in the app
        info["sha256"] = digest                      # update manifest hash
        print(f"  new sha256={digest}", flush=True)
        shutil.rmtree(ex)                            # free space
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"=== engine signing complete; total_fail={total_fail} ===", flush=True)
    return 1 if total_fail else 0

if __name__ == "__main__":
    sys.exit(main())
