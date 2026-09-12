"""Exercise the local DynaMol API; keep every job and stage for diagnosis."""
import argparse
import json
from pathlib import Path
import time
import urllib.error
import urllib.request


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--base-url", default="http://127.0.0.1:8792")
    parser.add_argument("--actions", default="{}")
    parser.add_argument("--build-loops", action="store_true")
    parser.add_argument("--prepared", action="store_true")
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent
    stamp = time.strftime("%Y%m%d-%H%M%S")
    record_path = folder / f"{args.case}-run-{stamp}.json"
    record = {"case": args.case, "source_dataset_id": args.dataset, "seed": args.seed,
              "purpose": "Short explicit-water software stability test; no equilibration or biological claim", "jobs": []}

    def save():
        record_path.write_text(json.dumps(record, indent=2) + "\n")

    def request(path, data=None):
        req = urllib.request.Request(args.base_url + path,
            data=json.dumps(data).encode() if data is not None else None,
            headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=180) as response:
                return json.load(response)
        except urllib.error.HTTPError as error:
            raise RuntimeError(error.read().decode()) from error

    def await_job(job):
        record["jobs"].append(job)
        save()
        started, last = time.monotonic(), None
        while True:
            job = request("/api/jobs/" + job["id"])
            record["jobs"][-1] = job
            save()
            state = (job["status"], job["stage"], int(job["progress"] // 10))
            if state != last:
                print(json.dumps({k: job.get(k) for k in ("id", "name", "status", "stage", "progress", "elapsed_seconds", "error")}), flush=True)
                last = state
            if job["status"] in {"completed", "failed", "cancelled", "interrupted"}:
                if job["status"] != "completed":
                    raise RuntimeError(job.get("error", job["status"]))
                return job.get("output_dataset_id") or job["dataset_id"]
            if time.monotonic() - started > 7200:
                request("/api/jobs/" + job["id"] + "/cancel", {})
                raise RuntimeError("Two-hour stage deadline exceeded; requested cancellation, artifacts retained")
            time.sleep(5)

    try:
        prepared = args.dataset
        if not args.prepared:
            settings = {"dataset_id": args.dataset, "name": f"Pressure test · {args.case} · prepared",
                        "ph": 7, "build_missing_residues": args.build_loops, "optimize_sidechains": True,
                        "ligand_actions": json.loads(args.actions), "seed": args.seed}
            record["preparation_settings"] = settings
            readiness = request(f"/api/datasets/{args.dataset}/readiness", {"mode": "preparation", "settings": settings})
            record["preparation_readiness"] = readiness
            save()
            if not readiness["ready"]:
                raise RuntimeError(str(readiness["blockers"]))
            prepared = await_job(request("/api/preparations", settings))
        record["prepared_dataset_id"] = prepared
        solvated = await_job(request(f"/api/datasets/{prepared}/solvate", {"padding_nm": 1, "ph": 7, "seed": args.seed}))
        record["solvated_dataset_id"] = solvated
        metadata = request(f"/api/datasets/{solvated}")
        atoms = metadata["atoms"]
        ca = [atom for atom in atoms if atom["category"] == "protein" and atom["name"] == "CA"]
        center = ca[len(ca) // 3]
        group = [atom for atom in atoms if atom["chain"] == center["chain"] and atom["resid"] == center["resid"]]
        by_name = {atom["name"]: atom["index"] for atom in group}
        # Use bonded backbone atoms for angle/dihedral, independent of dataset atom numbering.
        bonds = metadata["bonds"]
        next_n = next((b if a == by_name["C"] else a for a, b in bonds
                       if (a == by_name["C"] or b == by_name["C"])
                       and atoms[b if a == by_name["C"] else a]["name"] == "N"), None)
        if next_n is None:
            raise RuntimeError("No connected backbone N found for the live dihedral")
        measurements = [
            {"id": "backbone-distance", "kind": "distance", "atoms": [ca[len(ca)//4]["index"], ca[len(ca)//2]["index"]], "label": "Backbone separation", "trackDuringRun": True},
            {"id": "backbone-angle", "kind": "angle", "atoms": [by_name[n] for n in ("N", "CA", "C")], "label": "Backbone N–CA–C", "trackDuringRun": True},
            {"id": "backbone-dihedral", "kind": "dihedral", "atoms": [by_name[n] for n in ("N", "CA", "C")] + [next_n], "label": "Backbone psi", "trackDuringRun": True},
        ]
        settings = {"dataset_id": solvated, "name": f"Pressure test · {args.case} · 100 ps", "engine": "openmm",
                    "duration_ps": 100, "timestep_fs": 2, "report_interval": 500, "temperature_k": 300,
                    "seed": args.seed, "solvent": "explicit", "padding_nm": 1, "minimize": True,
                    "equilibration_steps": 5000, "measurements": measurements}
        record["simulation_settings"] = settings
        record["simulation_readiness"] = request(f"/api/datasets/{solvated}/readiness", {"mode": "simulation", "settings": settings})
        save()
        if not record["simulation_readiness"]["ready"]:
            raise RuntimeError(str(record["simulation_readiness"]["blockers"]))
        record["output_dataset_id"] = await_job(request("/api/jobs", settings))
        record["status"] = "completed"
        save()
        print(json.dumps({"completed": args.case, "record": str(record_path), "dataset": record["output_dataset_id"]}), flush=True)
    except Exception as error:
        record["status"] = "failed"
        record["error"] = str(error)
        save()
        raise


if __name__ == "__main__":
    main()
