# Public DMP charge-artifact recovery

This separate follow-up recovered four **candidate DMP charge/connectivity sets linked to Panteva's cited source**, but did not recover an unambiguously identified Panteva simulation topology. No parameters were assigned, no simulation was run, and no application or shared checkpoint was changed.

## Bounded retrieval and provenance

Two of the three authorized public database lookups were used. The observed, unauthenticated R.E.DD.B. form returned no match for `dimethyl phosphate`; `dimethyl` returned the four whole-molecule projects below. A third query and guanosine search were unnecessary. Eight observed summary/archive links were downloaded with normal TLS verification. No downloaded code was executed.

All four project pages cite **Dupradeau et al., PCCP 2010, 12, 7821–7839**, the article cited as reference 51 by [Panteva et al.](https://theory.rutgers.edu/resources/pdfs/Panteva_JChemPhysB_2015_v119_p15460.pdf). This is a direct bibliographic connection; neither Panteva's inspected paper/SI nor these project records identifies one of these project codes as the actual fitted-system input.

| Public project | Uploaded, Paris time | Last updated, Paris time | Program/orientation family | O4/O5 charge, each | MOL2 charge sum |
|---|---|---|---|---|---|
| [W-11](https://upjv.q4md-forcefieldtools.org/REDDB/projects/W-11/) | 7 Nov 2005 15:39 | 27 Sep 2011 17:31 | GAMESS-US | −0.7956 e | −1.0000 e |
| [W-12](https://upjv.q4md-forcefieldtools.org/REDDB/projects/W-12/) | 7 Nov 2005 16:04 | 27 Sep 2011 17:36 | GAMESS-US | −0.7956 e | −1.0000 e |
| [W-13](https://upjv.q4md-forcefieldtools.org/REDDB/projects/W-13/) | 7 Nov 2005 16:16 | 15 Mar 2013 10:34 | Gaussian 1998 | −0.7952 e | −1.0002 e |
| [W-14](https://upjv.q4md-forcefieldtools.org/REDDB/projects/W-14/) | 7 Nov 2005 16:18 | 27 Sep 2011 17:38 | Gaussian 1998 | −0.7952 e | −1.0002 e |

W-11, W-12 and W-14 name Francois-Yves Dupradeau, UFR de pharmacie, UPJV, as the submitter. W-13's live submitter fields are blank; its publication authors are still present. Publication authorship has not been substituted for the missing submitter identity. Dates above are displayed database metadata, not inferred from archive timestamps.

## What was recovered

Each archive provides a 13-atom, 12-bond DMP coordinate/charge graph and two RESP inputs. The pages describe a gauche/gauche conformer, HF/6-31G*, Connolly surface, two-stage RESP and one program-selected orientation. The first-stage input explicitly equates atoms 7 and 8. Graph inspection independently identifies them as O4/O5, each connected only to phosphorus; O2/O6 each connect phosphorus to a carbon. Atom count, graph references and finite numerical values were checked.

The intended RESP total charge is −1. W-13/14's four-decimal output sums to −1.0002 e; this small residual is compatible with rounding. It is preserved, not silently normalized. W-11/12 exchange charge values between the two methyl/bridging-oxygen arms; W-13/14 show the analogous exchange. These variants are not a unique fitting input.

The MOL2 type column contains only element labels (`C/H/O/P`), and every stored bond order is single. The listed AMBER converter and additional force-field files are placeholders, not usable parameter files. W-13's archive also omits the listed AMBER converter file. Consequently these archives do **not** establish the original Amber atom typing, bonded terms, Lennard-Jones parameters or exact fitting topology. A missing extra parameter file does not itself mean standard force-field terms were impossible; those terms simply are not identified here.

## Consequence for corrected 6OIM

The previously frozen native GDP model is Meagher–Redman–Carlson GDP(−3). Its Mg-coordinating O2B carries −0.9552 e, versus the recovered DMP candidates' −0.7956/−0.7952 e. The native GDP alpha-phosphate oxygens are −0.9474 e. This is a concrete charge-model difference, **not proof that transfer fails**, and does not validate transfer of Panteva's fitted correction. The exact DMP variant, complete topology and applicability to the preserved GDP/protein environment remain unresolved. This DMP-only pass supplies no new GAFF2 `c6` mapping.

Public download works; an explicit redistribution grant was not recovered from these archives. The site combines a free-project notice with a rights-reserved footer. Nothing has been selected for bundling.

## Evidence and next step

`artifact-inventory.json` records every request URL, UTC retrieval, HTTP status, content hash, archive member and charge graph. `evidence.json` records source identities, checks, placeholders, native GDP comparison and the remaining gaps. Raw source downloads are retained in `/Users/ashujo/.cache/dynamol-research/complete-complex-track/fit-artifact-recovery-red-db-v1/`. The original recovery report remains unchanged.

The next bounded source step is to inspect the Dupradeau 2010 manuscript/SI for these orientation variants and any exact DMP library identifiers. Only an explicit primary-source selection can identify the Panteva input. If none is available, retain these as candidates and require separate complete-system research validation before claiming support; do not choose a charge set from its filename or matching citation.
