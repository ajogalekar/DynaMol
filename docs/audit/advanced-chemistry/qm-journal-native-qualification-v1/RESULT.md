# Tiny native evaluation-journal result

The isolated reviewed candidate passed this small implementation qualification in **9.56 seconds**, with 23 no-QM harness checks and 15 native child processes. Candidate/source/helper files and the installed runtime remained unchanged.

Both DF RHF/STO-3G and DF RKS/B3LYPG/STO-3G/grid2 water optimizations converged in five completed evaluations with oxygen fixed. Compared with identical-input original-worker runs, final energies, all Cartesian gradient components, and all Cartesian coordinates were exactly identical. Both journal manifests and all ten evaluation records retained explicit UNCONVERGED/unaccepted/checkpoint-not-authorized states. Every record was checked against its callback and independently replayed by the original worker at exactly the recorded Bohr coordinates. Maximum replay differences were 4.2633e-14 Hartree in energy and 1.1671e-9 Hartree/Bohr per gradient component, below the predeclared 1e-8 and 1e-7 bounds.

The separate native persistence-failure run computed a genuine first gradient, then deliberately raised ENOSPC during record publication. It correctly failed with no accepted arrays, no completed journal record and no successful callback reference. The temporary write and traceback remain available.

All 15 owned child processes/groups were absent after completion. Root-authorized existing worker 37089/controller 37084 was only read for admission; it was not signalled or changed. The tiny-only 10 GiB entry reserve retained 8 GiB running reserve, with a 256 MiB total audit budget. Final run artifacts used about 5.23 MB before final reports. The failed 35 GiB-preflight v1 is retained and launched no native work.

RSS was sampled once per second by the unchanged reviewed supervisor. Since each tiny child finished in under one second, the retained samples cover startup and **do not measure true peak memory**. This result therefore makes no full-parent or peak-memory qualification claim. It also does not validate force-field parameters, metal chemistry, a Hessian, checkpoint reuse or molecular accuracy.

Full evidence: `/Users/ashujo/.cache/dynamol-research/qm-journal-native-qualification-v2/`; see SUMMARY.json, run/result.json, run/cleanup-review.json and artifact-manifest.json. Native calculations were not repeated to generate this report.
