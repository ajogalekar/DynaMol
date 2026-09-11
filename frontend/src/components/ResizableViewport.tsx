import { useEffect, useRef, useState, type ReactNode } from 'react';
import { GripHorizontal } from 'lucide-react';

const storageKey = 'dynamol.viewer-height';
const minimum = 280;
const maximum = () => Math.max(minimum + 60, window.innerHeight - 230);
const clamp = (value: number) => Math.round(Math.max(minimum, Math.min(maximum(), value)));

export default function ResizableViewport({ children }: { children: ReactNode }) {
  const viewport = useRef<HTMLDivElement>(null);
  const drag = useRef<{ y: number; height: number } | null>(null);
  const [height, setHeight] = useState<number | null>(() => {
    try {
      const saved = localStorage.getItem(storageKey);
      return saved !== null && Number.isFinite(Number(saved)) ? clamp(Number(saved)) : null;
    } catch {
      return null;
    }
  });
  const [measured, setMeasured] = useState(400);
  const [dragging, setDragging] = useState(false);
  const [limit, setLimit] = useState(maximum);

  useEffect(() => {
    const update = () => {
      setLimit(maximum());
      setHeight((value) => (value === null ? null : clamp(value)));
    };
    window.addEventListener('resize', update);
    return () => window.removeEventListener('resize', update);
  }, []);
  useEffect(() => {
    const observer = new ResizeObserver(() => {
      if (viewport.current)
        setMeasured(Math.round(viewport.current.getBoundingClientRect().height));
    });
    if (viewport.current) observer.observe(viewport.current);
    return () => observer.disconnect();
  }, []);
  useEffect(() => {
    try {
      if (height === null) localStorage.removeItem(storageKey);
      else localStorage.setItem(storageKey, String(height));
    } catch {
      // Resizing remains available when the browser disables local storage.
    }
  }, [height]);

  const finish = () => {
    drag.current = null;
    setDragging(false);
  };
  return (
    <>
      <div
        ref={viewport}
        className="molecule-viewport"
        id="molecular-viewport"
        style={height === null ? undefined : { height, minHeight: height, flex: '0 0 auto' }}
      >
        {children}
      </div>
      <div
        className={`viewer-resize-handle${dragging ? ' dragging' : ''}`}
        role="separator"
        tabIndex={0}
        aria-label="Resize protein view"
        aria-orientation="horizontal"
        aria-controls="molecular-viewport"
        aria-valuemin={minimum}
        aria-valuemax={Math.max(limit, measured)}
        aria-valuenow={measured}
        aria-valuetext={`${measured} pixels high. Drag or use arrow keys; double-click or press Enter to reset.`}
        title="Drag to resize the protein view · double-click to reset"
        onDoubleClick={() => setHeight(null)}
        onPointerDown={(event) => {
          if (event.button !== 0) return;
          event.preventDefault();
          event.currentTarget.focus();
          event.currentTarget.setPointerCapture(event.pointerId);
          drag.current = {
            y: event.clientY,
            height: viewport.current?.getBoundingClientRect().height ?? measured,
          };
          setDragging(true);
        }}
        onPointerMove={(event) => {
          if (drag.current) setHeight(clamp(drag.current.height + event.clientY - drag.current.y));
        }}
        onPointerUp={(event) => {
          finish();
          if (event.currentTarget.hasPointerCapture(event.pointerId))
            event.currentTarget.releasePointerCapture(event.pointerId);
        }}
        onPointerCancel={finish}
        onLostPointerCapture={finish}
        onKeyDown={(event) => {
          const step = event.shiftKey ? 60 : 20;
          if (event.key === 'ArrowDown') setHeight(clamp(measured + step));
          else if (event.key === 'ArrowUp') setHeight(clamp(measured - step));
          else if (event.key === 'Home') setHeight(minimum);
          else if (event.key === 'End') setHeight(maximum());
          else if (event.key === 'Enter') setHeight(null);
          else return;
          event.preventDefault();
          event.stopPropagation();
        }}
      >
        <GripHorizontal size={18} aria-hidden="true" />
        <span>Drag to resize view</span>
      </div>
    </>
  );
}
