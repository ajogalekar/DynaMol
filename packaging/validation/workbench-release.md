# Workbench package refresh · September 11, 2026

The local self-contained Apple Silicon package **9413166b6b1729cb** passed the [application check](workbench-app.json) and [archive/source verification](workbench-archive.json). The ZIP is `build/releases/DynaMol-0.1.0-macos-arm64.zip` (536,728,341 bytes); its SHA-256 is recorded beside the archive and in the verification record.

With a fresh private runtime/data directory whose path contains spaces and a system-only PATH, the bundled interpreter launched both native workers. OpenMM and GROMACS completed real short ubiquitin simulations, each produced multiple trajectory frames, and both recorded native temperature observations. Readiness and 101-frame structural RMSD also passed through the packaged API.

The service restarted on a different allocated port. The saved workspace frame, named atom selection, structural-analysis settings and named project were identical afterward. Service reuse, progress page, authenticated shutdown and retained molecular datasets/jobs were checked as well. These calls did not use the source-development server or its data.

ZIP CRC checks passed; all 615 recorded internal files matched their manifest and current source files. Only the public bundled examples were included. These acceptance records were generated after the archive was built and are retained beside the source; they are not self-referential files inside that archive.

This remains a local prototype for Apple Silicon/macOS 14+, with an ad-hoc launcher signature. The check ran on the development Mac, not an independent clean machine. Apple Developer signing/notarization, native Quit-menu interaction on an unlocked desktop, and complete public redistribution materials remain separate release work. No molecular convergence or accuracy claim follows from these short native checks.
