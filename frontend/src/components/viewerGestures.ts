interface GestureControls {
  zoom(factor: number): void;
  pan(dx: number, dy: number): void;
}

type PinchEvent = Event & { scale: number };

function touchPair(touches: TouchList) {
  if (touches.length !== 2) return null;
  const [a, b] = [touches[0], touches[1]];
  return {
    distance: Math.hypot(a.clientX - b.clientX, a.clientY - b.clientY),
    x: (a.clientX + b.clientX) / 2,
    y: (a.clientY + b.clientY) / 2,
  };
}

/** Own pinch input before NGL's document listeners interpret Ctrl+wheel as clipping. */
export function bindViewerGestures(host: HTMLElement, controls: GestureControls) {
  let gestureScale: number | null = null;
  let previousTouch: ReturnType<typeof touchPair> = null;
  const options = { capture: true, passive: false };
  const consume = (event: Event) => {
    event.preventDefault();
    event.stopPropagation();
  };
  const zoom = (factor: number) => {
    if (Number.isFinite(factor) && factor > 0 && factor !== 1) controls.zoom(factor);
  };

  const wheel = (event: WheelEvent) => {
    // Chromium and Firefox expose trackpad pinch as a Ctrl-modified wheel event.
    // Leave ordinary scrolling and the other mouse controls with NGL.
    if (!event.ctrlKey) return;
    consume(event);
    if (gestureScale !== null || previousTouch) return;
    const unit = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? host.clientHeight : 1;
    zoom(Math.exp(Math.max(-1, Math.min(1, -event.deltaY * unit * 0.01))));
  };
  const gestureStart = (event: Event) => {
    consume(event);
    const scale = (event as PinchEvent).scale;
    gestureScale = Number.isFinite(scale) && scale > 0 ? scale : 1;
  };
  const gestureChange = (event: Event) => {
    consume(event);
    const scale = (event as PinchEvent).scale;
    if (!Number.isFinite(scale) || scale <= 0) return;
    // WebKit uses cumulative gesture scale. On touch devices, touch events own
    // the same gesture so a parallel GestureEvent must not apply it twice.
    if (!previousTouch) zoom(scale / (gestureScale ?? 1));
    gestureScale = scale;
  };
  const gestureEnd = (event: Event) => {
    consume(event);
    gestureScale = null;
  };
  const touchStart = (event: TouchEvent) => {
    previousTouch = touchPair(event.touches);
    // NGL still receives start/end for its one-finger rotation state.
  };
  const touchMove = (event: TouchEvent) => {
    const current = touchPair(event.touches);
    if (!current) {
      previousTouch = null;
      return;
    }
    consume(event);
    if (previousTouch && previousTouch.distance > 0) {
      zoom(current.distance / previousTouch.distance);
      // Allow the gesture center to move without losing pinch updates. NGL's
      // default heuristic discards slow pinches and treats drifting ones as pan.
      controls.pan(previousTouch.x - current.x, previousTouch.y - current.y);
    }
    previousTouch = current;
  };
  const touchEnd = (event: TouchEvent) => {
    previousTouch = touchPair(event.touches);
  };
  const reset = () => {
    previousTouch = null;
    gestureScale = null;
  };

  host.addEventListener('wheel', wheel, options);
  host.addEventListener('gesturestart', gestureStart, options);
  host.addEventListener('gesturechange', gestureChange, options);
  host.addEventListener('gestureend', gestureEnd, options);
  host.addEventListener('touchstart', touchStart, options);
  host.addEventListener('touchmove', touchMove, options);
  host.addEventListener('touchend', touchEnd, options);
  host.addEventListener('touchcancel', reset, options);
  window.addEventListener('blur', reset);

  return () => {
    host.removeEventListener('wheel', wheel, options);
    host.removeEventListener('gesturestart', gestureStart, options);
    host.removeEventListener('gesturechange', gestureChange, options);
    host.removeEventListener('gestureend', gestureEnd, options);
    host.removeEventListener('touchstart', touchStart, options);
    host.removeEventListener('touchmove', touchMove, options);
    host.removeEventListener('touchend', touchEnd, options);
    host.removeEventListener('touchcancel', reset, options);
    window.removeEventListener('blur', reset);
  };
}
