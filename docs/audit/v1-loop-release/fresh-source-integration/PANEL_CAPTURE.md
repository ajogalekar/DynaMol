# Private panel capture helper

`panel_capture.py` captures one fresh common preparation immediately at the
`_refine_loops` call boundary. The complete loop candidate validator must run
separately. A snapshot is identity/serialization evidence; it is not structural
acceptance, simulation readiness, or app publication.

Run with the stable focused Python, using an explicit JSON configuration:

```sh
/Users/ashujo/.cache/dynamol-runtimes/covalent-v1-focused/bin/python -B \
  /Users/ashujo/Documents/Science/DynaMol/docs/audit/v1-loop-release/fresh-source-integration/panel_capture.py \
  --config /absolute/path/capture.json
```

The configuration requires every field shown below. `source_repo` is optional
and defaults to `/Users/ashujo/Documents/Science/DynaMol`.

```json
{
  "case": "CASE-ID",
  "root": "/absolute/path/to/new-attempt-root",
  "dataset_id": "raw_dataset_id",
  "expected_raw_atoms": 1234,
  "seed": 2026
}
```

An equivalent environment-only invocation sets `DYNAMOL_CAPTURE_CASE`,
`DYNAMOL_CAPTURE_ROOT`, `DYNAMOL_CAPTURE_DATASET`, `DYNAMOL_CAPTURE_RAW_ATOMS`,
and `DYNAMOL_CAPTURE_SEED`. `DYNAMOL_CAPTURE_CONFIG` can instead name the JSON.
The child reads the configuration copy frozen inside the attempt root.

The new attempt root must already contain:

- `workspace/datasets/<dataset_id>/`: exactly one unprepared, single-frame dataset,
  including its `topology.pdb`, `metadata.json`, ordinary raw dataset sidecars,
  all available `original-source/` files, and the sequence evidence the backend
  will actually select. Evidence must be local inside this dataset.
- `overlay/`: the previously qualified pure-Python `dimorphite_dl` and `loguru`
  dependencies, including their data/metadata. Native dependencies come from the
  focused runtime.
- `raw-source-provenance.json` and `overlay-provenance.json`: each has a `files`
  array. Every row has absolute `source` and `copy` paths plus the SHA256 string
  `sha256`. The raw manifest must enumerate every file in the raw dataset, and
  each `copy` must name its actual location in this new root. Extra provenance
  fields are preserved. The source paths must remain readable for post-run
  verification. Overlay rows must point to their verified new copies.

Minimal manifest shape:

```json
{
  "files": [
    {
      "source": "/absolute/original/raw-file",
      "copy": "/absolute/new-attempt-root/workspace/datasets/raw_dataset_id/raw-file",
      "sha256": "64 lowercase hexadecimal characters"
    }
  ]
}
```

Do not copy prepared coordinates, previous job folders, or ligand parameter
outputs into the new root. The helper freezes backend Python, bundled residue
data and project metadata before entry. `DYNAMOL_DATA_DIR` selects the private
workspace; explicit installed ProMod3 and AmberTools prefixes keep native
resolution independent of the copied source root. No repository `.venv` or
`.tools` runtime is selected. The same qualified AmberTools file pins are
rechecked; fresh Python and ProMod3 import/startup evidence is recorded inside
the supervised child before preparation.

The inherited supervisor lifecycle and `OwnedTree` class are unchanged: one
600-second attempt, 4 GiB aggregate descendant RSS limit sampled once per
second, two threads, disk reserves, and cleanup of owned descendants. Fresh
AM1-BCC/SQM is permitted only as ordinary ligand preparation. MD and complete
loop refinement are not performed. Unsupported chemistry, missing evidence,
partial sequence coverage, altered atom identities, or inability to reach the
capture boundary fail with explicit evidence. This helper does not change
backend support or acceptance thresholds.

Outputs include `result.json`, `capture-result.json`, `snapshot/snapshot.json`,
raw-to-common movement/mapping reports, all parameter sidecars, recursive base
force-field files, the actual selected sequence evidence, every raw dataset
file, frozen implementation sources, native logs, resource samples, cleanup
status, and an artifact manifest. `early-failure.json` records a failure before
normal supervision when a valid existing root is available. Use a new attempt
root after any attempt; existing plans or frozen source directories are refused.

The pure-Python helper checks use `test_panel_capture.py`; they execute no
preparation or native tools. Native replay of each case remains separate work.
