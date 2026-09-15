# Private pyDiSGro 0.1.0 evaluation — 2026-09-14

**Verdict: fast coordinate generation, but no admissible repaired complex from this bounded test.** Keep this as an isolated research candidate generator. Do not vendor or enable it in DynaMol. All 40 retained proposals need substantial anchor correction; none satisfies the current 1 Å observed-backbone displacement cap before refinement. Redistribution rights for inherited code/data also remain unresolved.

## Frozen experiment and results

Two sequential, single-thread CPU workers used 256 trials, 32 distance states, 20 retained conformations, no sidechain sampling, and seeds 91401/91402. Each worker had a 180 s deadline and a supervisor enforcing a 4 GiB resident-memory ceiling (RSS sampled every 0.2 s). No MD or QM was performed. The original package source was not edited.

| Case | Missing residues | Wall time | Peak sampled RSS | Claimed closed / stored | Retained | Coherent gross geometry | Coherent CCTBX | Within 1 Å anchor cap | Geometry after restoring original anchors |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1UA2 | A44–55, KLGHRSEAKDGI | 5.92 s | 123.1 MiB | 100 / 73 | 20 | 20/20 | 1/20 | 0/20 | 0/20 |
| 8K5R | A177–179, AKN | 4.03 s | 122.2 MiB | 196 / 196 | 20 | 20/20 | 3/20 | 0/20 | 0/20 |

The maximum displacement of an observed anchor atom ranges from **1.868–2.657 Å** for 1UA2 and **1.621–2.631 Å** for 8K5R. A sampler `closed` flag is therefore not evidence of a loop fitting the original protein. For example, restoring 1UA2 candidate 000's original anchor creates a 3.142 Å C55–N56 gap.

Ramachandran scores use the existing native **CCTBX 2025.11** tables, residue classes, and thresholds. Each JSON includes the native-library path/hash and all scored residues, including the two anchors. Gross geometry reuses the existing independent backbone validator unchanged. The contact diagnostic checks generated backbone atoms against **every retained source heavy atom**, excluding local backbone 1–2/1–3 neighbors, and records overlaps below 1.8 Å. This limited overlap screen is not full force-field, sidechain, chirality, rotamer, or complete-complex validation. No whole-preparation pass is claimed.

These are two saved normalized repair contexts, not a survey of all deposited molecules or an assertion that every original crystallographic water was present in the context. Source context/config files were copied byte-for-byte and rechecked. Other gaps elsewhere in either protein were not rebuilt in this experiment.

## Complete-complex identity and parser limitations

The package `Structure.readPdb` keeps only `ATOM` records and canonical residue templates. In these actual inputs it omits:

| Context | Source atoms | Omitted by sampler |
|---|---:|---|
| 1UA2 | 2,323 | TPO A170: 11 atoms; ATP B1: 31 atoms; terminal ASN A311 OXT: 1 atom |
| 8K5R | 2,502 | TPO A186: 11 atoms; VQE B1: 21 atoms; terminal MET A327 OXT: 1 atom |

The sampler search and empirical ranking **do not account for those omitted components**. A later full-scaffold check does not retroactively make the search aware of the complete chemical environment. The isolated adapter retains all source records in the immutable scaffold, extracts only identity-mapped missing backbone coordinates, and keeps changes proposed for observed anchors separately. It never uses the raw sampler PDB as a replacement complex.

Independent parser probes reproduce silent omission of an insertion-coded atom and acceptance of a one-residue sequence for a twelve-residue gap. The adapter requires an exact sequence/range match and unambiguous source identities before running. The empty parser result still reports a default `numRes=1`, so actual populated atoms, not that count, are audited.

The raw writer removes chain identifiers and `CONECT` records and renumbers residues. Both raw candidate-000 outputs contain **zero heterogen atoms**. They are retained with a `RAW-LOSSY-PARSER-OUTPUT` filename for evidence only. No zero-coordinate heavy atoms were observed in those two inspected raw outputs; potential behavior on other incomplete templates remains untested.

## Provenance and parity

The inspected PyPI source archive has SHA256:

`766da438d5fcbef5bfeb65543a72562350b22acbc1235dded2d8b64a46e19680`

Its explicit MIT notice names the 2026 port author. All **eight bundled parameter/distribution files plus their README** are byte-identical to the original C++ repository pinned at commit `6c25da054ff391e93b7d41c221cbc0086d665371`. The full, nontruncated repository tree contains no LICENSE/COPYING/copyright-named file; the inspected README and source header supply no redistribution grant. That is an unresolved rights chain, not a conclusion that the port author's MIT notice clears the inherited material. The upstream README describes the original command-line workflow, while its data README identifies the empirical tables. [Pinned original README](https://github.com/aa14k/Disgro/blob/6c25da054ff391e93b7d41c221cbc0086d665371/Readme.txt), [pinned data README](https://github.com/aa14k/Disgro/blob/6c25da054ff391e93b7d41c221cbc0086d665371/data/readme.txt).

Every downloaded upstream file has its source URL, Git blob ID, SHA256, and comparison in `upstream-provenance/comparison.json`; all extracted package files have a checked source manifest. No author was contacted.

This is **not numerical parity validation** against the published C++ implementation. The port explicitly defaults away from the original sidechain-energy accumulation behavior, uses NumPy's random generator, and leaves some original CLI features unsupported. Neither sidechain behavior nor C++/Python candidate equivalence was tested. The original method is described by Tang, Zhang, and Liang; its published performance does not establish the correctness or speed of this new port. [Original paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003539).

## Artifacts and handoff

Project-owned scripts and compact evidence:

- `evaluate_private.py`: bounded worker/supervisor and identity-preserving proposal adapter.
- `inspect_saved.py`: source hash verification, synthetic parser probes, and saved-output audit; does not resample.
- `result-summary.json`: machine-readable combined results, provenance, and candidate paths.

Private source/provenance and failed first adapter attempt:

`/Users/ashujo/.cache/dynamol-research/v1-loop-release/disgro-evaluation-v1/`

Successful frozen evaluation:

`/Users/ashujo/.cache/dynamol-research/v1-loop-release/disgro-evaluation-v2/`

The first attempt failed before sampling because the adapter imported `SMC` from the package root instead of `pydisgro.smc`. Its original script, protocol, logs, and failures are retained in v1. Only the private adapter was corrected. The v2 directory contains its protocol, all 512 per-trial events, all stored proposals in NPZ files, all 40 JSON candidates, raw lossy outputs, native CCTBX/geometry details, and resource/exit reports. A first inspection mistakenly used the empty parser's default residue count; that inspection is retained as `inspection-v1.json`, and the corrected probe counts populated atoms. No failed sampler result was discarded or rerun with a looser gate.

Each `candidate-NNN.json` uses:

```json
{
  "schema_version": 1,
  "immutable_scaffold": "/absolute/path/original-context.pdb",
  "source_sha256": "...",
  "loop": {"chain_id": "A", "start_residue": 177, "end_residue": 179, "sequence": "AKN", "left_anchor": 176, "right_anchor": 180},
  "generated_atoms": [{"identity": ["A", "177", "", "ALA", "N"], "xyz_angstrom": [0.0, 0.0, 0.0], "element": "N"}],
  "research_context_atoms": [],
  "app_ready": false,
  "whole_preparation_pass": false
}
```

The illustrative zero coordinate above is schema notation, not a real candidate. Actual candidates require finite, non-placeholder coordinates. `generated_atoms` contains N/CA/C/O only for missing residues; `research_context_atoms` contains those backbone atoms for the two observed anchors. Restoring those anchors is a separate, audited operation. Sidechains must be rebuilt by a later validated step.

**Exact next bounded step:** if testing native restrained refinement of an independent seed family is useful, take **8K5R candidate 005** from the v2 directory. Its coherent CCTBX score passes and it has no detected source-scaffold gross overlaps, but its anchor shift is 2.041 Å. Rebuild sidechains and refine within the original complete-complex parameter bundle, then apply the original observed-atom cap, peptide geometry, all-atom contact/chirality, and CCTBX gates. It is an input to that test, not an accepted repair. Candidate 009 is a second such seed with a larger 2.631 Å anchor shift. The only CCTBX-passing 1UA2 seed, candidate 003, also has three gross overlaps and a 2.570 Å shift, making it a less favorable immediate probe. Do not expand the sampler budget or integrate it into the release based on these results.
