# Bounded 8K5R coherent fragment-pool comparison

**No candidate passed all local backbone checks.** This comparison expands the smallest failing 8K5R loop (source A:177–179, AKN) from one coherent fragment to a frozen pool of 40 distinct native fragments. All raw fragments have valid gross internal backbone geometry and preserve the observed-only peptide cis/trans basins, but none fits every relevant observed N/CA/C/O atom within the original 1 Å cap.

## Frozen search and observed endpoint diagnostics

The raw native pool and exact closure IDs were saved before closure or independent CCTBX outcomes. The rule retained the first eight native database candidates per eligible context, with a global cap of 40 and at most two observed internal context residues. No candidates were substituted after a failed reference outcome. The selected subset consists of two per context with the lowest raw complete observed N/CA/C/O squared displacement. The ranking and cap both include carbonyl oxygens.

| Source span | Native database matches | Retained raw fragments | Selected for KIC |
|---|---:|---:|---:|
| 176–180 | 0 | 0 | 0 |
| 175–180 | 10 | 8 | 2 |
| 176–181 | 96 | 8 | 2 |
| 174–180 | 39 | 8 | 2 |
| 175–181 | 427 | 8 | 2 |
| 176–182 | 97 | 8 | 2 |

The pool consists of coherent native `BackboneList` copies before CCD, KIC or observed-coordinate restoration. Every retained fragment has a coordinate hash, source-identity mapping, native context/ordinal, and per-observed-atom displacement record. Experimental fragment accession metadata is not exposed by the API used here; no held-out database or independent-structure claim is made. Raw pool diagnostics show a placement mismatch already present before KIC, rather than malformed transplanted peptide joins.

## Outcomes

| Stage | Coordinates screened | Gross backbone geometry pass | Original observed omega basins preserved | Local CCTBX pass | Complete observed 1 Å cap pass | Combined pass |
|---|---:|---:|---:|---:|---:|---:|
| Raw coherent pool | 40 | 40 | 40 | 11 | 0 | 0 |
| KIC outputs from 10 selected fragments | 40 | 40 | 40 | 0 | 0 | 0 |
| Intact control KIC outputs | 12 | 12 | 12 | 4 | 4 | 4 |

Three of the ten selected raw fragments pass the local reference before closure; none of their returned KIC alternatives retains a passing complete local reference. The 40 KIC outputs are raw returned solutions, not necessarily 40 independent conformations. Each selected fragment used at most four internal pivot triples; each triple was tested unperturbed and with one independent non-pivot conformer draw. The native coil sampler seed was 2026. The observed n-stem carbonyl orientation was represented before closure; source observed-only cis/trans state and the complete N/CA/C/O cap were checked afterward. Bond lengths, bond angles and force constants were not fitted.

The closest closed result is `fragment-06-2618ffce2178`, trial 1, solution 0: LEU176 C moves **1.892 Å**, LEU176 O 1.721 Å, and SER180 N 1.348 Å. It also has reference outliers at LEU176, ALA177, LYS178 and GLN181. It therefore does not qualify for promotion to full-atom preparation.

The native batch completed in **0.84 s**, with peak sampled RSS **937,836,544 bytes** (0.87 GiB), under the 180 s / 4 GiB supervisor limits. No adapter errors, cancellations or resource-limit events occurred. All source/helper hashes were unchanged across the native run; both independent screens also verified their inputs. The new scripts passed AST and targeted whitespace checks.

## Limits and next implication

This bounded native-order pool is not an exhaustive search of the returned database candidates or protein conformations. The original rigid-stem context returned no native database matches; the extended contexts contained more matches than the retained quota. The results do not establish that the gap is impossible to model and do not justify relaxing the displacement cap. They show that simply selecting among these coherent fragments and closing them with the tested pivot/sampling plan does not yet satisfy the observed geometry and local reference together.

Only backbone coordinates were assessed. The independent screen uses exact original observed coordinates outside each candidate span and includes torsions affected by rebuilt or moved atoms. A 0.0001 Å tolerance distinguishes native coordinate roundoff when selecting affected residues; it does not change the 1 Å quality cap or native CCTBX thresholds. Unchanged distant experimental outliers are excluded from the local gate. Sidechain chirality/rotamers, complete retained-environment contacts, full force-field refinement and saved full-complex admission remain untested. No app, MD or QM run was initiated and no complete candidate was accepted.

[Compact results](fragment-pool-summary.json), [native pool/closure implementation](compare_native_fragment_pool.py), [raw reference screen](audit_raw_fragment_pool.py), and [post-closure reference screen](screen_kic_backbones.py). Full inputs, coordinates, frozen plans, implementation snapshots and checks: `/Users/ashujo/.cache/dynamol-research/v1-loop-release/8k5r-fragment-pool-v1`.
