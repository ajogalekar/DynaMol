# 1UA2 source-supported versus reconstructed boundary torsions

**ARG57 is a pre-existing source outlier whose phi changes during the mobile-flank replays. ASN56 is a newly defined boundary conformation because its preceding ILE55 carbonyl was missing.** Neither distinction makes the current complete candidates acceptable. No threshold or acceptance-policy change is recommended from this single case.

This read-only audit compared the original selected 1UA2 source, the archived refinement-entry coordinates, and all six completed 1UA2 candidate outputs in `hybrid-complete-validation-v1`. Exact atom identities determined the availability and movement of each defining atom. Native CCTBX classifications and double-precision dihedrals were used. Observed N/CA/C/O atoms at residues 42, 43, 56, 57 and 58 are identical between the original source and the archived input (maximum difference 0 Å).

## Boundary interpretation

The missing sequence is A:44–55. Values below are in degrees. The representative final values are from `fragment_mc-2027`.

| Residue | Original experimental support | Final phi / psi | Interpretation |
|---|---|---:|---|
| LYS42 | Both torsions defined; source −105.450284 / 127.897670 | −105.450284 / 131.766054 | Source-supported psi changes because observed ILE43 N moves; remains favored. |
| ILE43 | Phi defined (−123.942760); psi undefined because LYS44 N is missing | −115.367892 / 118.868231 | Newly defined boundary pair; favored in this output. The archived seed's psi was modeled, not observed. |
| ASN56 | Phi undefined because ILE55 C is missing; psi defined (62.367211) | 80.037638 / 46.559098 | Newly defined boundary pair; remains an outlier in this output. Its archived seed classification cannot be treated as an unchanged experimental outlier. |
| ARG57 | Both torsions defined; source 67.719090 / −83.366294 | 64.672912 / −83.366294 | Pre-existing source outlier with changed phi, because ASN56 C moves. ARG57's own atoms and psi remain fixed. |

The original ARG57 general-class CCTBX score is **0.00029607199579877354**, classified as an outlier. Its phi uses ASN56 C–ARG57 N–ARG57 CA–ARG57 C. Its psi uses ARG57 N–ARG57 CA–ARG57 C–THR58 N. Thus an unmoved ARG57 residue does not imply that its complete backbone conformation is unchanged: a preceding residue contributes to phi.

For ASN56, the archived entry reports phi −35.382918 and psi 62.367211. The former already depends on the constructed ILE55 C; the complete pair has no fully observed original counterpart. In the representative result, ASN56 N, CA and C move 0.504755, 0.038479 and 0.060575 Å respectively. This changes its source-supported psi as well as the newly defined phi.

## ARG57 across the six replays

| Attempt | ASN56 C displacement (Å) | Final ARG57 phi | ARG57 in selected reference? |
|---|---:|---:|---|
| archive-fragment | 0 | 67.719090 | No; unchanged outside the selected modeled region |
| torsion_mc-2026 | 0.044593 | 65.503513 | Yes; defining atom changed |
| torsion_mc-2027 | 0.010876 | 67.166103 | Yes; defining atom changed |
| fragment_mc-2026 | 0.022005 | 67.920812 | Yes; defining atom changed |
| fragment_mc-2027 | 0.060575 | 64.672912 | Yes; defining atom changed |
| disgro-representative | 0.091913 | 62.951650 | Yes; defining atom changed |

ARG57 remains an outlier in all six coordinate sets; its psi stays exactly −83.366294°. The repair did not create its original outlier status, but the mobile-flank outputs do alter a defining coordinate. Reporting it as a newly created outlier would be inaccurate; reporting it as wholly unchanged would also be inaccurate.

Inspection of every selected reference row found **no fully unchanged, source-supported pre-existing outlier unnecessarily included in these six actual selections**. The fixed-context archive attempt excludes ARG57. In all five mobile-flank alternatives, both outside neighbors selected by the validator are actually affected: LYS42 psi depends on moved ILE43 N, and ARG57 phi depends on moved ASN56 C. This does not prove that the validator's one-neighbor expansion is minimal for every future candidate; it establishes the distinction for these saved outputs without loosening any gate.

## Provenance

Repository dataset root: `/Users/ashujo/Documents/Science/DynaMol/docs/audit/loop-fallback/runs/root-final/1UA2/workspace/datasets`.

| Input | SHA256 |
|---|---|
| `3a661cfc670048c5/topology.pdb` — original selected complete source | `c58d4cce211bab37bb61ede35799dce5534408f27cbceac7fa8b5208263b8ded` |
| `4dcb5c25f5c14559/loop-refinement-topology.cif` | `79784fdf00348d2153f25622899d521a1d95bc28466ccc7364cd84eaa0bdfa40` |
| `4dcb5c25f5c14559/loop-refinement-input.npz` | `8bd9acd39fe9bc22dbfb78efc8e25807d7784343f42a9622b4e7b1ddfe672a0a` |
| Native `mmtbx_validation_ramachandran_ext.so` | `b7a531377948dc6c1534ffe2f4041cfe28c38eac6817e12e0808e828da28935c` |

Candidate root: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/hybrid-complete-validation-v1/1UA2/search`. Each row below names `<attempt>/validation/refined.npz`; the adjacent `reference.json` and `result.json` identify the selected residue set and recorded checks.

| Attempt | Refined-coordinate SHA256 |
|---|---|
| archive-fragment | `f51467e077bf5d58b7f8289dad14841302741654c3ded8482a2f6efeac4334b3` |
| torsion_mc-2026 | `e365519da0027857ecb99e7d10645f2ea2b654ca6fc69403b5af38154665ea83` |
| torsion_mc-2027 | `c995c191c89ae7617b3edbe45717f28d573c7cd1d853cf9520643cfbbd169862` |
| fragment_mc-2026 | `f01da986bb2a57188dbd64d024b3d7f792d85b64378a822c70ee25a093b72bf7` |
| fragment_mc-2027 | `c90dfd41900804483fd437a8039fdf7601c9676d4113b45295c50cd004965ac4` |
| disgro-representative | `3dd6a192670d29ff3dc3a6c94bbf5e5a97db0f9d9e30b4d745823f1867eb2bc8` |

The audit ran with the isolated research Python and native CCTBX target. No coordinates, existing evidence, validator code or acceptance thresholds were modified; no sampling/refinement/MD/QM was performed for this interpretation.
