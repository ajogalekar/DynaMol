# Independent review: two-anchor additive repair

No remaining scoped blocker after the separate v2 hydrogen-parent admission correction. This verdict concerns the standalone additive geometry candidate, not a prepared 1CLL system, a validated calcium model, or the newer ethanol force-field package.

I verified all 15 original artifacts and all seven v2 artifacts, the three captured original-function AST hashes, and the exact source delta: only `admit` plus the new CCD-column validator changed. The initially identified gap was that a disconnected declared H or H–H pair could be ignored by the heavy-only graph reader. v2 requires every CCD H to have exactly one heavy parent and validates complete required columns. The author’s 41 admission-only checks reproduce both original gaps, reject them after correction, and show unchanged admission for the frozen EOH specimen. I inspected those tests and results without rerunning RDKit, embedding, relaxation, fitting, QM, or MD.

Independent stdlib-only coordinate checks confirm both deposited carbon coordinates are bitwise unchanged throughout the 24 saved alternatives. The saved coordinates match the declared rigid rotations to 3.67e-15 Å, and all pair distances are preserved to 9.33e-15 Å. H21 remains an explicitly deposited hydrogen record; it was not relabeled as oxygen. The selected candidate remains the explicitly arbitrary first surviving orientation, with physical acceptance and simulation readiness false.

The full hydrogen name mapping, complete ligand mechanics, full-system environment/materialization, and physical calcium-site validation remain separate obligations. No inference of a unique experimental missing-oxygen pose follows from the coarse clash screen. The newer mapped ethanol package was not reviewed here.

An initial read-only harness compared AST digests using system Python 3.14.6 and failed because the recorded digests use Python 3.12 AST formatting. That failure is preserved. The final stdlib-only audit ran with the qualified Python 3.12.13 interpreter and passed; no scientific source or original digest changed.

Evidence: `/Users/ashujo/.cache/dynamol-research/reviews/1cll-two-anchor-independent-v1/review.json`.
