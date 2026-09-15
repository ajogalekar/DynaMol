# Isolated quantum worker

`backend/qm_worker.py` runs in the private `.tools/qm` Python environment. The main application does not import PySCF. Initial supported methods are closed-shell RHF and RKS with explicitly specified molecular charge, spin `0`, basis, and (for RKS) functional. `B3LYPG` is the explicit Gaussian-VWN B3LYP variant; the ambiguous `B3LYP` alias is rejected.

```sh
.tools/qm/bin/python backend/qm_worker.py input.json fresh-output-directory
```

The output directory must not exist. Failures leave `result.json` with `accepted:false`, the input, native log and traceback; no `arrays.npz` is retained as an accepted result. Numerical acceptance is not force-field validation or evidence that the chosen chemical state is correct.

## Input

```json
{
  "schema_version": 1,
  "atom_ids": ["water/O", "water/H1", "water/H2"],
  "elements": ["O", "H", "H"],
  "coords_angstrom": [[0,0,0], [0.75,0,0.59], [-0.76,0,0.60]],
  "charge": 0,
  "spin": 0,
  "method": "RHF",
  "basis": "sto-3g",
  "operations": ["gradient", "hessian", "esp"],
  "esp_points_bohr": [[3,2,4], [-4,3,2]],
  "threads": 2,
  "max_memory_mb": 8000
}
```

Supply exactly one of `coords_angstrom` or `coords_bohr`. Atom IDs must be stable and unique. `basis` and optional `ecp` accept a named string or element-to-name dictionary; absent-element keys and unresolved requested ECPs are rejected. For RKS, add `functional:"B3LYPG"` and `grid_level:3`. `density_fit` defaults to `false`; if enabled, optional `auxbasis` and the resolved auxiliary basis are recorded.

`operations` is a nonempty list from `energy`, `gradient`, `hessian`, `esp`. ESP requires explicit `esp_points_bohr`; this worker does not generate a fitting surface. A second job can evaluate ESP at a grid constructed around the first job’s optimized output coordinates.

Optional optimization uses `optimization:{"enabled":true,"max_steps":100,"freeze_atom_ids":[...]}`. Constraints freeze the named Cartesian coordinates at their input positions; they do not constrain fitted charges. The five geomeTRIC tolerances can be supplied explicitly:

| Field | Default | Units |
| --- | ---: | --- |
| `convergence_energy` | 1e-6 | Hartree |
| `convergence_grms` | 3e-4 | Hartree/bohr |
| `convergence_gmax` | 4.5e-4 | Hartree/bohr |
| `convergence_drms` | 0.0012 | Å |
| `convergence_dmax` | 0.0018 | Å |

Optimization must report convergence, every evaluated SCF must converge, and frozen atoms must remain within 1e-5 bohr. Requested final properties are recomputed at the actual final geometry. Limits include `max_scf_cycles` (100), `max_cphf_cycles` (100), and `max_wall_seconds` (7200, maximum 172800). Threads are restricted to at most two and PySCF memory to at most 8000 MB. The memory setting is an allocation hint, not an OS memory limit. Callers must impose an outer process timeout as well; a Python signal cannot necessarily interrupt a long native-library call immediately.

Explicit fixed dihedrals are supported for bounded constrained optimizations:

```json
{
  "optimization": {
    "enabled": true,
    "max_steps": 100,
    "freeze_atom_ids": ["cap/C", "cap/N"],
    "dihedral_constraints": [
      {"atom_ids": ["CB", "SG", "ligand/C25", "ligand/C24"], "target_degrees": 60}
    ],
    "dihedral_tolerance_degrees": 0.1
  }
}
```

All IDs must exist in the supplied atom inventory; the example names are illustrative. Each ordered A–B–C–D quartet must have four distinct atoms. Up to 32 explicit targets are accepted, in the range −360° to 360°. The original target is retained; its periodic equivalent in [−180°,180°] is written to geomeTRIC’s documented `$set`/`dihedral` constraint with one-based indices. Cartesian freezes retain their original behavior. Repeated/reversed duplicate quartets, conflicting targets on the same quartet, all-four-atoms-frozen redundancies, and initially undefined or nearly linear dihedrals are rejected. Other infeasible combinations must still converge and satisfy the final checks.

Initial, evaluated-step and final actual angles are recorded with the atom-ID/index map and constraint-file hash. The worker independently checks the final periodic error against `dihedral_tolerance_degrees` (default 0.1°, maximum 0.5°) and frozen-coordinate drift against 1e-5 bohr. A native convergence flag alone is insufficient. Dihedrals that become singular during optimization fail. These are constrained QM geometries/energies for caller-defined profile work, not an automatic torsion fit or proof that a molecular-mechanics trajectory is physically accurate.

An optional `initial_checkpoint` reuses the density of a completed accepted worker calculation at exactly the same initial geometry:

```json
{
  "initial_checkpoint": {
    "result_path": "/absolute/previous-output/result.json",
    "result_sha256": "64 lowercase hexadecimal characters",
    "checkpoint_path": "/absolute/previous-output/scf.chk",
    "checkpoint_sha256": "64 lowercase hexadecimal characters"
  }
}
```

The hashes must be actual SHA-256 digests, not the descriptive strings above. The source result, arrays, expanded-method record and checkpoint are copied into the new output’s `initial-checkpoint/` directory. The source files are never overwritten. Checks require matching stable atom order, elements, charge, spin, electron/AO counts, PySCF version, spherical basis, resolved basis/ECP, and (for RKS) functional and grid level. Both source arrays and the checkpoint molecule must match the requested initial geometry within 1e-10 bohr. No AO projection or density transfer to a different geometry is allowed.

Checkpoint orbitals must be finite, orthonormal in the requested basis and have valid closed-shell occupancies; their density must be symmetric and have the correct electron count. The checkpoint energy must match the pinned accepted result within 1e-9 Hartree. That provenance check does not recompute the source energy. A DF calculation may seed an exact calculation, but the requested SCF must independently converge under its own settings. Older accepted worker reports without a checkpoint digest require the separately pinned input digest and all the same state/density/energy checks. New result records include `checkpoint_file` and `checkpoint_sha256`. For optimization, the requested method first converges from the verified density, then geomeTRIC makes explicit subsequent geometry changes.

The implementation uses the public [PySCF checkpoint loader](https://pyscf.org/_modules/pyscf/scf/chkfile.html) and the explicit [`SCF.kernel(dm0=...)` initial-density argument](https://pyscf.org/_modules/pyscf/scf/hf.html).

## Output

`result.json` records acceptance, energy, atom IDs, elements, electron count, SCF/optimization/response completion, input/worker/array hashes, versions and explicit units. `resolved-method.json` preserves actual expanded basis/ECP data and hashes, functional resolution and density-fitting settings. `progress.json` records the current stage; native output goes to `native.log`.

`requested_threads` and `actual_pyscf_threads` distinguish the configured request from PySCF’s reported OpenMP thread count. The original macOS wheel reports one even when two are requested; the rebuilt parallel runtime reports two. This does not measure instantaneous CPU utilization or every numerical library’s threading. If `sys.prefix/dynamol-qm-runtime.json` exists, the worker copies it to `runtime-manifest.json` and records its SHA-256 as declared build provenance. It does not claim to verify every loaded binary against that manifest.

`arrays.npz` contains:

| Key | Shape | Units |
| --- | --- | --- |
| `coords_angstrom` | N×3 | Å |
| `coords_bohr` | N×3 | bohr |
| `atomic_numbers` | N | nuclear atomic number |
| `effective_nuclear_charges` | N | e, after ECP core removal |
| `gradient_hartree_per_bohr` | N×3 | Hartree/bohr |
| `hessian_hartree_per_bohr2` | 3N×3N | Hartree/bohr²; atom-major xyz |
| `esp_points_bohr` | M×3 | bohr |
| `esp_hartree_per_e` | M | Hartree/e |

Optional arrays appear only for requested operations. The native Hessian axes `(atom_i,atom_j,xyz_i,xyz_j)` are explicitly transposed to `(atom_i,xyz_i,atom_j,xyz_j)` before flattening. The exact `pyscf.lib.param.BOHR` conversion constant is recorded.

ESP uses the public [PySCF 031-MEP example](https://github.com/pyscf/pyscf/blob/master/examples/1-advanced/031-MEP.py): positive effective nuclear potential minus electronic potential from batched `int3c2e` integrals. Deterministic sample points are checked independently with `int1e_rinv` to 1e-8 Hartree/e. ESP points coincident with nuclei are rejected. ECP calculations represent effective point cores plus explicit valence electrons; they do not reconstruct an all-electron potential close to a removed core.

## Validation

[The original checkpoint](qm-worker-checks/worker-checkpoint.json) and [original native test log](qm-worker-checks/native-final-tests.log) retain the initial 24 passing tests. The [restart checkpoint](qm-worker-checks/restart-worker-checkpoint.json) and [restart test log](qm-worker-checks/native-restart-tests.log) record **31 passing tests**, including reruns of all original numerical checks and the new restart controls. They cover all nine water RHF/STO-3G Cartesian finite differences, Hessian ordering, independent ESP integrals, positive-ion far-field sign, Zn LANL2DZ effective-core charges, B3LYPG RKS gradient/Hessian/ESP, constrained optimization and failure cases. New numerical artifacts are preserved in `qm-worker-checks/native-suite-restart-1`.

The initial full water derivative check had maximum gradient error 1.98e-9 Hartree/bohr and Hessian error 7.54e-8 Hartree/bohr² at a 1e-4-bohr displacement. The restart suite also compares DF-initialized exact RHF with a fresh exact RHF calculation (energy within 1e-10 Hartree, gradient within 1e-8 Hartree/bohr, ESP within 1e-9 Hartree/e), verifies warm constrained optimization, and rejects changed identities, geometry, charge, basis, functional/grid, artifact hashes, checkpoint energy and invalid orbitals/occupancies. These are numerical implementation tests on tiny molecules, not validation of a large metal-site model or agreement with Gaussian on that model.

The [dihedral checkpoint](qm-worker-checks/dihedral-worker-checkpoint.json) records 40 passing full-suite checks plus a separate [negative-convergence control](qm-worker-checks/dihedral-false-convergence-test.log), for 41 passing tests. Actual RHF/STO-3G H–O–O–H constrained optimizations reached 59.999999973° for a 60° target and −60.000007715° for the periodic 300° target; both frozen oxygen coordinates remained within 1.6e-7 bohr. The separate control makes an optimizer falsely report convergence at the unchanged 90° geometry and verifies that the worker rejects it without accepted arrays. The [parallel-runtime smoke result](qm-worker-checks/parallel-dihedral-smoke/result.json) also verifies the constraint, two reported PySCF threads, and copied runtime-manifest provenance. No large covalent or metal-site scan was run for this worker test.

Optimization behavior follows the public [PySCF geometry-optimization interface](https://pyscf.org/user/geomopt.html) and [geomeTRIC constraint specification](https://geometric.readthedocs.io/en/latest/constraints.html). Geometries, charges, spin states, surface grids and fitting constraints remain explicit responsibilities of the calling scientific workflow.
