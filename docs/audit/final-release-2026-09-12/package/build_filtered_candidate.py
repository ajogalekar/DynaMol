from pathlib import Path
import importlib.util
root=Path.cwd()
spec=importlib.util.spec_from_file_location("dynamol_package_candidate",root/"scripts/package_macos.py")
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
original=module.ignore_copy
audit=root/"docs/audit/final-release-2026-09-12"
def filtered(directory,names):
    ignored=original(directory,names)
    path=Path(directory).resolve()
    if path==audit and "package" in names:ignored.append("package")
    if path.is_relative_to(audit):
        ignored.extend(name for name in names if Path(name).suffix.lower() in {".zip",".mp4",".bin"})
    return sorted(set(ignored))
module.ignore_copy=filtered
module.build(root/"build/releases/final-audit-2026-09-12/filtered",app_only=True)
