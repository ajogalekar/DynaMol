# Stable local RESP runtime

The copied, previously qualified RESP executable passed relocation regression:
all printed charges in both stages of the 113-atom reference and constrained
synthetic 94-atom fixture match their retained original outputs exactly. All 18
fixed canonical charges and the neutral caps are preserved. These are execution
and relocation checks, not newly calculated physical charges.

The isolated prefix is `/Users/ashujo/.cache/dynamol-runtimes/resp-heap-local-v1`.
Only the copied executable's library-search header and ad-hoc signature changed.
The four files' file-backed Mach-O section payloads are unchanged; the three
Fortran runtime libraries remain byte-identical. Actual dyld logs show all three
dependencies resolving inside the local prefix, without Documents paths.
The original runtime and active quantum worker remain untouched.

Four successful fitting stages took 1.739 seconds total, with a largest sampled
RSS of 2,076,114,944 bytes under a 3 GB sampled guard. Every owned process exited.
The fixed workspace declares four 8000-by-8000 double arrays (2.048 GB of array
storage); sampling is not an OS allocation cap or a guarantee of the true peak.

Two earlier tests remain retained. The initial 1 GB test budget stopped during
initialization; its exact peak was not persisted. A separate long-absolute-path
test exposed native filename truncation and exit0 with no fit output. The final
test stages short names in its local working directory and requires real parsed
outputs. The original qualified DynaMol fitter already uses this convention;
no change to its numerical fitting logic was needed.

The preliminary `-h` loader probe printed an unknown-flag usage message with
exit0. It establishes dependency loading only; the four actual fitting stages
provide functional evidence.

The immutable admission JSON is `local-qualification/qualification.json` inside
the prefix, SHA `dfe338698f80029eb15668cb44480b695e875ac0b85f37794b0e422448889851`.
Its four-file inventory SHA is `80a80b1e6716aade0a42150957ee9056e56a49a16fa398e93d93212f167fb0c7`;
the 91-artifact evidence manifest SHA is `f519f2266c1de508994292303192b14c612b9610e079f72d894879fcae7cf189`.
Independent review is recorded separately; the frozen admission artifact is not
rewritten when that review completes. Actual ESP provenance, convergence,
geometry, charge-fit and model-acceptance gates remain required downstream.
