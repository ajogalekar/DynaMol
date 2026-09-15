# Shift the native closure anchor past the observed carbonyl

**The installed API supports this geometry. An intact control closed a suffix using a temporary first-modeled-residue anchor while preserving the preceding residue's N/CA/C/O and first peptide join exactly.** This demonstrates API feasibility, not that 1UA2 or 8K5R can be repaired under the existing gates. No target-protein search was performed.

## Construction contract

For 1UA2, start with one coherent candidate covering source ILE43 through ASN56. Align its ILE43 N/CA/C frame to the source. Match its N–CA–C–O orientation to source ILE43 using `SetPsiTorsion(0, adjusted_psi, True)`: **the final `True` is necessary** to rotate O43 and the entire downstream chain while leaving N43/CA43/C43 fixed. The default `False` instead rotates N43, which would invalidate the intended locked frame. Rotation preserves the incoming candidate's local lengths/angles; it cannot make a mismatched source triangle or carbonyl length/angle exact. Check all four original stem atoms before closure.

Create the temporary anchor from the coherent prefix and close the suffix:

```python
temporary = locked_full.Extract(0, 2).ToEntity()
virtual_first = temporary.residues[1]
suffix = locked_full.Extract(1, len(locked_full))
closer = modelling.CCDCloser(
    virtual_first, original_c_stem, suffix.GetSequence(), torsions, seed
)
closed = suffix.Copy()
success = closer.Close(suffix, closed)
combined = locked_full.Copy()
if success:
    combined.ReplaceFragment(closed, 1, False)
```

`Extract(start, end)` is zero-based and excludes `end`; its result is an independent backbone copy. `ToEntity()` produces chain A, residues 1/2, connected as a peptide, with N/CA/C/O and optional CB. These temporary identities are **not source residue numbers**: record the explicit map 1→source ILE43 and 2→source LYS44. `BackboneList.InsertInto(chain, start_resnum)` is also available when native numbering is useful. Do not export this temporary backbone entity as the prepared complete complex.

The n-stem and c-stem handles may come from different entities. The temporary LYS44's preceding residue is coherent ILE43, allowing CCD to calculate the incoming phi44 and Ile neighbor class. Supply the **original source/native-model ASN56 handle** as c-stem, not the candidate's floating ASN56. Its following original ARG57 provides psi56 and the far-end O reconstruction direction. Verify native `InSequence` for both external links, exact sequence, source maps, units and handles before closing.

Recombine with `ReplaceFragment(..., superpose_stems=False)` to copy the already closed suffix without another transform. It preserves the locked prefix. Score and validate the recombined full candidate with the unchanged retained environment and full source identity map. Never overwrite O43 alone afterward to manufacture a passing frame or peptide plane.

## Degrees of freedom and implementation limits

Relative to full-span CCD, this removes **psi43 and phi44** from closure's adjustable torsions. Their values now belong to the coherent proposal/virtual-anchor construction. CCD can adjust psi44, remaining modeled phi/psi, and phi56. It does not change omega through those rotations. Every observed-only omega basin and both newly constructed joins still require independent checks.

A virtual anchor fixed for an entire native suffix-MC attempt is possible with native sampler/closer objects. It also freezes phi44 for that entire attempt. Rebuilding the virtual anchor for every full-span proposal is a different protocol: the compiled `LoopCandidates.FillFromMonteCarloSampler` expects native sampler/closer pointers and cannot simply accept a Python closure object. That version needs a controlled explicit driver or a native wrapper with tested proposal, rejection and output-state handling. The installed Python `SampleMonteCarlo` helper has unvalidated state reassignment issues noted in [BOUNDARY-API.md](/Users/ashujo/Documents/Science/DynaMol/docs/audit/v1-loop-release/hybrid-fallback/BOUNDARY-API.md); it should not be substituted without controls.

For the secondary 8K5R case, the shifted span is source A177–180, with source LEU176 as the temporary preceding residue. It leaves only two internal residues, so native KIC's three-internal-pivot requirement prevents its use as this suffix closer. Native CCD can be constructed, but closure feasibility remains untested here. No context expansion or relaxation of the 1 Å cap follows from this design.

Native CCD still fits only the c-stem N/CA/C target, then reconstructs its O from the following observed N. That can move observed ASN56 atoms and alter ARG57 phi. The existing cap and affected-boundary reference checks remain necessary. Locking source O43 also constrains the direction of the missing peptide N44; it may conflict with allowed boundary geometry. The control does not establish that these constraints are jointly satisfiable for 1UA2.

## One intact control

The preserved six-residue fixture was used only to test transformations, slicing, temporary peptide connectivity, cross-entity target handles and recombination. A prescribed +0.1-radian first-psi rotation was undone by the source-carbonyl orientation adjustment, then exactly one shifted-anchor CCD closure was performed. There were no random conformer draws or search retries.

| Check | Result |
|---|---:|
| Coherent orientation roundtrip, maximum N/CA/C/O difference | 0.0000103 Å |
| Closure return | True |
| Locked prefix N/CA/C/O displacement during suffix closure | 0 Å |
| Virtual first-modeled N/CA/C displacement | 0 Å |
| First join C–N length, O–C–N angle, omega and incoming phi | Exactly unchanged |
| Input suffix modified by `Close(input, output)` | No |
| C-stem N/CA/C maximum displacement from original | 0.00000815 Å |
| C-stem O displacement from original | 0.175217 Å |
| C-stem O difference from native next-N reconstruction formula | 0 Å |

Evidence: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/carbonyl-lock-design-v1/intact_control.py` and `intact-control.json`, which retain original, orientation-locked and recombined coordinates. Fixture SHA256: `039779d649c37c7f01f56fceebcc5bab28ec0d21ca68df877f875d04d14af1a1`. Runtime: installed ProMod3 3.6.0 / OST 2.11.1, two-thread environment; 0.0134 seconds native activity, 122,880,000 bytes peak RSS. Native API/source provenance is recorded in [BOUNDARY-API.md](/Users/ashujo/Documents/Science/DynaMol/docs/audit/v1-loop-release/hybrid-fallback/BOUNDARY-API.md).

No app code, existing coordinates, force-field parameters, source maps or acceptance gates were modified. This is a bounded construction hypothesis for a subsequent frozen comparison, with complete-complex validation still required.
