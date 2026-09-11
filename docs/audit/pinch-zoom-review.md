# Pinch zoom regression check

All four focused browser tests passed on the rebuilt local application at `http://127.0.0.1:8765/` using macOS Chrome (11.6 seconds; no failures, skips, or retries). Tests measure the rendered molecule's pixel extent and centroid rather than accessing NGL or React internals. Raw results and geometry measurements are in [pinch-zoom-results.json](pinch-zoom-results.json).

- Native Playwright Ctrl+wheel increased the molecule from 228 × 254 to 292 × 327 pixels and restored it to 228 × 254. Browser viewport scale stayed at 1 and page dimensions/scroll position were unchanged.
- Synthetic WebKit `gesturestart/change/end` sequences with cumulative scales 1.1, 1.2, 1.3 and interleaved duplicate Ctrl+wheel events produced 296 × 330 pixels (approximately 1.3×), then restored the original dimensions. This tests WebKit's DOM event protocol in Chrome; it does not substitute for a physical trackpad test in Safari or the embedded macOS WebKit browser.
- Real Chrome DevTools Protocol two-touch input, increasing separation by only 1.5 pixels per event while drifting the center, grew the molecule to 330 × 369 pixels and moved its centroid 27.9 pixels right and 13.0 pixels down. The inverse gesture restored dimensions and returned the centroid within 0.7 pixels. Page scale remained unchanged.
- Ordinary wheel input still zoomed the molecule in both directions (228 → 245 → 229 pixels wide).

Run from `frontend/`:

```sh
DYNAMOL_BASE_URL=http://127.0.0.1:8765 npm run test:e2e -- tests/pinch-zoom.spec.ts
```

The old production build was replaced before this focused suite ran; no execution baseline against the old build is claimed. Initial test harness failures arose from DOM overlays entering the screenshot bounds and an overly large assumption about one ordinary wheel event's zoom amount; the final harness isolates the molecular scene and checks observable direction and expected pinch scale.
