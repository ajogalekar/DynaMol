# Loop-building alternatives and automatic fallback

Research review, 14 September 2026. No new generator was installed, executed
or admitted to the application by this review. The 12-residue-per-gap release
scope, frozen panel and existing acceptance checks remain unchanged.

## Recommended architecture

Use a fast fragment-based attempt followed by bounded alternative sampling
when no candidate passes the complete preparation checks. Longer gaps can
receive larger search budgets. A provisional 1–6 versus 7–12 grouping is
useful for measuring performance, but is not an established algorithmic
boundary. The existing three-residue 8K5R failure is a reason to permit
fallback at any supported length.

Keep one default "Automatic loop repair" control. All generators must retain
the exact sequence and atom mapping, respect observed context, and pass the
same join, chirality, backbone-conformation, retained-environment and saved-file
checks after refinement. A generator's own score or closed chain is not an
acceptance result. Longer loops may admit several plausible conformations;
passing geometry does not establish the experimentally missing conformation.

## Existing ProMod3 routes

`backend/promod_loop_worker.py::construct_loops` already tries database
fragments and a native torsion Monte Carlo fallback. However, the stage loop
stops when `model.gaps` becomes empty. Later full-atom refinement/quality
failure does not cause a return to candidate generation. The priority
integration change is therefore **fallback after a failed complete candidate**,
with bounded, recorded attempts and coherent observed context preserved.

The unused native fragment-insertion Monte Carlo mode changes proposals from
individual torsions to structural fragment replacements. It can reuse the
existing runtime but shares ProMod3's closure/scoring infrastructure, so it
is not an independent implementation. Its installed default fragment length
is nine; it cannot be assigned indiscriminately to every 7–12-residue gap.
`GenerateDeNovoTrajectories` is free-chain generation, not anchored repair.
The full homology-modeling pipeline is also not a drop-in complete-complex
repair path. [Official pipeline documentation](https://openstructure.org/promod3/3.5/modelling/pipeline/),
[primary ProMod3 paper](https://doi.org/10.1371/journal.pcbi.1008667).

## Independent generators

**DiSGro / pyDisgro is the first external candidate to benchmark.** Its
distance-guided sequential chain-growth sampling provides a different search
strategy. The original paper reports results for longer loops including
10–17 residues; its best-sampled and energy-selected accuracies differ, so
candidate generation and selection both need evaluation. Historical CPU
timings do not predict DynaMol's desktop latency.
[Original study](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003539).

The [pyDisgro 0.1.0 distribution](https://pypi.org/project/pydisgro/) is a
recent NumPy-only Python port. Its source archive was read without executing
package code and checked against the published SHA256
`766da438d5fcbef5bfeb65543a72562350b22acbc1235dded2d8b64a46e19680`.
It contains an MIT license notice, although its PyPI license metadata is
empty. This is not yet a redistribution clearance: the linked
[original DiSGro repository](https://github.com/aa14k/Disgro), inspected at
`6c25da054ff391e93b7d41c221cbc0086d665371`, has no license grant in its
license-file inventory or examined code/data READMEs. The port's notice
does not establish rights to inherited implementation and empirical tables.
Resolve provenance before bundling. Native algorithm parity, ranking,
complete-complex handling, cancellation and actual CPU timing remain untested.

**RCD/RCD+ is a scientific comparator with distribution gaps.** The available
sampling-core download carries BSD-style terms but supplies Linux64 binaries;
the official repository contains a README rather than implementation source.
It omits the full RCD+ pipeline, whose published all-atom results include
Rosetta refinement. Thus those results cannot be assigned to the standalone
core. No verified Mac source/build is available from this inspection.
[Official repository](https://github.com/chaconlab/RCD),
[primary RCD+ paper](https://academic.oup.com/nar/article/44/W1/W395/2499369).

**LoopTK** provides closure/sampling building blocks rather than a turnkey
ranked predictor. Its official distribution lists older source, but access
and current code/data license could not be verified in this review; no Mac
build was tested. [Distribution](https://simtk.org/frs/?group_id=214),
[documentation](https://robotics.stanford.edu/looptk/doc/index.html).

Rosetta NGK/KIC-with-fragments is a relevant separately licensed comparator.
Rosetta's own FAQ explicitly distinguishes source availability from open
source and prohibits ordinary licensee redistribution. It cannot silently
become part of the downloadable bundle.
[Rosetta licensing FAQ](https://rosettacommons.org/software/licensing-faq/).

KarmaLoop/AutoLoop are possible future learned-model candidates. This review
did not establish a license-cleared, supported Mac repair path. The examined
KarmaLoop repository includes a CUDA-oriented environment and no visible
top-level license; the article's license is not a software/weights license.
[Author repository](https://github.com/karma211225/KarmaLoop).

## Next bounded comparison

First correct candidate routing in an isolated prototype and compare native
torsion versus fragment proposals where applicable. In parallel, review the
pyDisgro source/data provenance and implementation before a private bounded
benchmark. Evaluate the same known failures and frozen panel; add masked
intact-loop controls with the removed coordinates withheld from generation
and ranking. Stratify 1–6 and 7–12 residues, while allowing fallback in either
group. Record candidate-selection accuracy separately from best-in-ensemble
accuracy and complete-preparation acceptance. No new method is selected for
shipping yet.
