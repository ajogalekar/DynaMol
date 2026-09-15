# Standalone Mac application

Current Mac builds display DynaMol in a dedicated, resizable application window
with a Dock icon. The interface runs in Apple's system WebKit inside the app;
it does not require an open Chrome or Safari window. OpenMM, GROMACS, Python,
AmberTools and ProMod3 remain bundled, with the same supported chemistry and
resource limits as the browser version.

This candidate is separate from the published browser-opening v0.1.1 installer.
The explicit [published release page](https://github.com/ajogalekar/DynaMol/releases/tag/v0.1.1)
continues to identify the previous download until a new candidate is published.
GitHub's `latest/download` shortcut does not resolve a prerelease asset.

## Everyday use

- Open DynaMol from Applications to show its workspace. First launch displays
  progress while unpacking the included engines; it uses no package manager.
- Use Open Files for a native Mac file chooser. Exports use a native Save dialog.
- The green window button and View menu provide fullscreen. The molecule
  viewer retains its own zoom, representations, playback and measurement tools.
- Closing the window keeps DynaMol in the Dock. Click its Dock icon to reopen
  the same workspace. Use **DynaMol → Quit DynaMol** or **Command-Q** to exit.
- If work is running, Quit offers a choice to keep jobs running or stop them
  before exiting. Reopen the app to inspect retained jobs and saved results.
- External documentation links open in your browser. Molecular work stays in
  the local application service.

## Existing installations

The app keeps the same bundle identity and data directory,
`~/Library/Application Support/DynaMol`. Replacing the application preserves
structures, projects, simulation results and private runtime caches. Quit the
old app before replacing it; closing a browser tab alone does not quit its
menu-bar launcher.

Apple Silicon and macOS 14 or newer remain the target. This change does not add
Developer ID signing or notarization; the existing per-app Gatekeeper opening
instructions still apply. See [installation details](PACKAGING.md).

## Building and verification

Follow the [Mac packaging instructions](../packaging/README.md). The launcher is
Swift/AppKit with WKWebView, and uses the system WebKit framework rather than
shipping an additional browser engine. The manifest identifies the interface as
`native-macos-webkit`. The app uses the existing private service startup,
runtime integrity checks, background workers and shutdown endpoint.

Native-window acceptance is separate from the earlier Chrome audit. It must
exercise WebGL rendering, file import and download dialogs, fullscreen, window
reopening and quitting, as well as verify the exact signed bundle and container.
