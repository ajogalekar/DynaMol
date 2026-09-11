import {
  forwardRef,
  useCallback,
  useEffect,
  useImperativeHandle,
  useLayoutEffect,
  useMemo,
  useRef,
  useState,
} from 'react';
import * as NGL from 'ngl';
import type { Dataset, Measurement, Representation, Visibility } from '../types';
import {
  atomGroup,
  atomIsVisible,
  atomSelection,
  createCoordinateInterpolator,
  visibleGroups,
  shortProteinFragments,
} from './viewerGeometry';
import './MolecularViewer.css';

export interface ViewerHandle {
  fit(): void;
  /** Factors greater than one zoom in; factors between zero and one zoom out. */
  zoom(factor: number): void;
  focus(indices: number[]): void;
  snapshot(): Promise<Blob | null>;
}

export interface MolecularViewerProps {
  dataset: Dataset | null;
  coordinates: Float32Array | null;
  frame: number;
  visibility: Visibility;
  representation: Representation;
  colorScheme: 'chain' | 'residue' | 'element';
  selectedAtoms: number[];
  measurements: Measurement[];
  picking: boolean;
  spin: boolean;
  onAtomPick(index: number): void;
  onReady(): void;
  onError(message: string): void;
}

// A small public-API boundary keeps the component independent of NGL's evolving
// generic representation parameter declarations. These methods are verified
// against nglviewer/ngl src/component and src/controls, not private WebGL state.
interface Vector {
  x: number;
  y: number;
  z: number;
}
interface AtomProxy extends Vector {
  index: number;
  atomname: string;
  resname: string;
  resno: number;
  chainname: string;
  element: string;
}
interface Signal<T extends (...args: never[]) => unknown = () => void> {
  add(fn: T): void;
  remove(fn: T): void;
}
interface RenderRepresentation {
  setSelection(selection: string): void;
  setVisibility(visible: boolean): void;
  setParameters(parameters: Record<string, unknown>): void;
  update(what: Record<string, boolean>): void;
}
interface StructureComponent {
  structure: {
    atomCount: number;
    updatePosition(positions: Float32Array, refresh?: boolean): void;
    getAtomProxy(index: number): AtomProxy;
  };
  addRepresentation(type: string, parameters: Record<string, unknown>): RenderRepresentation;
  removeRepresentation(representation: RenderRepresentation): void;
  updateRepresentations(what: Record<string, boolean>): void;
  autoView(selection?: string | number, duration?: number): void;
  getCenter(selection: string): Vector;
  getZoom(selection: string): number;
}
interface PickProxy {
  atom?: AtomProxy;
  closestBondAtom?: AtomProxy;
  position?: Vector;
  surface?: unknown;
  canvasPosition?: { x: number; y: number };
}
type PickSignal = {
  add(fn: (proxy?: PickProxy) => void): void;
  remove(fn: (proxy?: PickProxy) => void): void;
};
interface Stage {
  loadFile(file: Blob, options: Record<string, unknown>): Promise<StructureComponent>;
  removeComponent(component: StructureComponent): void;
  handleResize(): void;
  setSpin(spin: boolean): void;
  dispose(): void;
  makeImage(options: Record<string, unknown>): Promise<Blob>;
  signals: { clicked: PickSignal; hovered: PickSignal };
  mouseControls: { remove(trigger: string): void };
  viewer: {
    renderer?: unknown;
    width: number;
    height: number;
    requestRender(): void;
    signals: { rendered: Signal };
  };
  viewerControls: {
    getCameraDistance(): number;
    getPositionOnCanvas(position: Vector): { x: number; y: number };
    signals: { changed: Signal };
  };
  animationControls: {
    zoom(distance: number, duration?: number): void;
    zoomMove(center: Vector, distance: number, duration?: number): void;
  };
}

type RenderGroups = {
  polymer: RenderRepresentation;
  water: RenderRepresentation;
  ligands: RenderRepresentation;
  ions: RenderRepresentation;
  detail: RenderRepresentation;
  hydrogen: RenderRepresentation;
  selected: RenderRepresentation;
  selectedLabel: RenderRepresentation;
};

const COLOR_SCHEMES = { chain: 'chainname', residue: 'residueindex', element: 'element' };

function validateTopology(component: StructureComponent, dataset: Dataset) {
  if (
    component.structure.atomCount !== dataset.n_atoms ||
    dataset.atoms.length !== dataset.n_atoms
  ) {
    throw new Error(
      'Topology and trajectory atom counts differ. Reload matching topology and trajectory files.',
    );
  }
  // Check every atom: equal counts alone do not establish consistent ordering.
  // The server's canonical PDB can remap residue numbers/chain IDs for PDB limits,
  // so compare the atom and residue names that must be preserved by that writer.
  for (let index = 0; index < dataset.n_atoms; index++) {
    const proxy = component.structure.getAtomProxy(index);
    const atom = dataset.atoms[index];
    if (
      atom.index !== index ||
      // NGL's AtomMap uppercases atom names and elements (e.g. Cl -> CL).
      proxy.atomname.trim().toUpperCase() !== atom.name.trim().toUpperCase() ||
      proxy.resname.trim().toUpperCase() !== atom.residue.trim().slice(0, 3).toUpperCase() ||
      (atom.element !== 'X' && proxy.element.toUpperCase() !== atom.element.toUpperCase())
    ) {
      throw new Error(
        `Topology atom ordering differs at atom ${index + 1}. Display is stopped to prevent incorrect measurements.`,
      );
    }
  }
}

const MolecularViewer = forwardRef<ViewerHandle, MolecularViewerProps>(
  function MolecularViewer(props, ref) {
    const hostRef = useRef<HTMLDivElement>(null);
    const tooltipRef = useRef<HTMLDivElement>(null);
    const stageRef = useRef<Stage | null>(null);
    const componentRef = useRef<StructureComponent | null>(null);
    const groupsRef = useRef<RenderGroups | null>(null);
    const measurementRepresentations = useRef(new Map<string, RenderRepresentation>());
    const labelElements = useRef(new Map<string, HTMLDivElement>());
    const interpolatorRef = useRef<ReturnType<typeof createCoordinateInterpolator> | null>(null);
    const hoveredIndex = useRef<number | null>(null);
    const latest = useRef(props);
    latest.current = props;
    const [stageVersion, setStageVersion] = useState(0);
    const [componentVersion, setComponentVersion] = useState(0);
    const [status, setStatus] = useState<'empty' | 'loading' | 'ready' | 'error'>('empty');
    const [error, setError] = useState('');
    const activeDatasetId = useRef<string | null>(null);
    const representationType = useRef<Representation>('cartoon');
    const ribbonFallback = useMemo(
      () =>
        props.dataset && props.coordinates
          ? shortProteinFragments(props.dataset, props.coordinates)
          : new Set<number>(),
      [props.dataset?.id, props.coordinates],
    );

    const fail = useCallback((message: string) => {
      setError(message);
      setStatus('error');
      latest.current.onError(message);
    }, []);

    const updateLabels = useCallback(() => {
      const stage = stageRef.current;
      const component = componentRef.current;
      const { measurements, dataset } = latest.current;
      if (!stage || !component || !dataset || activeDatasetId.current !== dataset.id) return;
      const width = stage.viewer.width;
      const height = stage.viewer.height;
      for (const measurement of measurements) {
        const element = labelElements.current.get(measurement.id);
        if (!element || !measurement.atoms.length) continue;
        const center = { x: 0, y: 0, z: 0 };
        for (const index of measurement.atoms) {
          if (index < 0 || index >= dataset.n_atoms) continue;
          const atom = component.structure.getAtomProxy(index);
          center.x += atom.x / measurement.atoms.length;
          center.y += atom.y / measurement.atoms.length;
          center.z += atom.z / measurement.atoms.length;
        }
        const projected = stage.viewerControls.getPositionOnCanvas(center);
        const y = height - projected.y;
        const onScreen =
          Number.isFinite(projected.x) &&
          Number.isFinite(y) &&
          projected.x > 20 &&
          projected.x < width - 20 &&
          y > 20 &&
          y < height - 20;
        element.style.visibility = onScreen ? 'visible' : 'hidden';
        element.style.transform = `translate(${projected.x}px, ${y - 17}px) translate(-50%, -100%)`;
      }
    }, []);

    const focus = useCallback((indices: number[]) => {
      const component = componentRef.current;
      const stage = stageRef.current;
      const dataset = latest.current.dataset;
      if (!component || !stage || !dataset) return;
      const valid = indices.filter(
        (index) => Number.isInteger(index) && index >= 0 && index < dataset.n_atoms,
      );
      if (!valid.length) return;
      const selection = atomSelection(valid);
      // Keep a useful neighborhood around a single atom instead of zooming into
      // a near-zero bounding box. The camera orientation is preserved.
      stage.animationControls.zoomMove(
        component.getCenter(selection),
        -Math.max(24, Math.abs(component.getZoom(selection))),
        450,
      );
    }, []);

    useImperativeHandle(
      ref,
      () => ({
        fit() {
          const component = componentRef.current;
          const dataset = latest.current.dataset;
          if (!component || !dataset) return;
          const visible = dataset.atoms.filter((atom) =>
            atomIsVisible(atom, latest.current.visibility),
          );
          const solute = visible.filter((atom) => atomGroup(atom) !== 'water');
          component.autoView(
            atomSelection((solute.length ? solute : visible).map((atom) => atom.index)),
            500,
          );
        },
        zoom(factor) {
          const stage = stageRef.current;
          if (stage && Number.isFinite(factor) && factor > 0) {
            stage.animationControls.zoom(
              Math.max(3, stage.viewerControls.getCameraDistance() / factor),
              180,
            );
          }
        },
        focus,
        async snapshot() {
          const stage = stageRef.current;
          if (!stage || !componentRef.current) return null;
          try {
            const blob = await stage.makeImage({
              factor: 2,
              antialias: true,
              trim: false,
              transparent: false,
            });
            const host = hostRef.current;
            const hostBounds = host?.getBoundingClientRect();
            const labels = latest.current.measurements.flatMap((measurement) => {
              const element = labelElements.current.get(measurement.id);
              if (!element || element.style.visibility === 'hidden' || !hostBounds) return [];
              const rect = element.getBoundingClientRect();
              return [
                {
                  x: rect.left - hostBounds.left,
                  y: rect.top - hostBounds.top,
                  width: rect.width,
                  height: rect.height,
                  text: Array.from(element.children)
                    .map((child) => child.textContent)
                    .join(' '),
                  color: measurement.color,
                },
              ];
            });
            if (!labels.length || !hostBounds) return blob;
            // NGL exports the WebGL scene. Composite the saved-frame value labels
            // too, so a downloaded figure preserves the analysis visible on screen.
            const bitmap = await createImageBitmap(blob);
            const canvas = document.createElement('canvas');
            canvas.width = bitmap.width;
            canvas.height = bitmap.height;
            const context = canvas.getContext('2d');
            if (!context) {
              bitmap.close();
              return blob;
            }
            context.drawImage(bitmap, 0, 0);
            bitmap.close();
            context.scale(canvas.width / hostBounds.width, canvas.height / hostBounds.height);
            context.font = '600 11px sans-serif';
            context.textBaseline = 'middle';
            for (const label of labels) {
              context.fillStyle = '#0d1a28';
              context.strokeStyle = label.color;
              context.lineWidth = 0.5;
              context.beginPath();
              context.roundRect(label.x, label.y, label.width, label.height, 6);
              context.fill();
              context.stroke();
              context.fillStyle = label.color;
              context.fillText(label.text, label.x + 8, label.y + label.height / 2);
            }
            return await new Promise<Blob | null>((resolve) => canvas.toBlob(resolve, 'image/png'));
          } catch (cause) {
            latest.current.onError(
              `Could not create snapshot: ${cause instanceof Error ? cause.message : String(cause)}`,
            );
            return null;
          }
        },
      }),
      [focus],
    );

    useEffect(() => {
      const host = hostRef.current;
      if (!host) return;
      let stage: Stage;
      try {
        stage = new NGL.Stage(host, {
          backgroundColor: '#0b121b',
          cameraType: 'orthographic',
          quality: 'high',
          sampleLevel: 0,
          lightIntensity: 0.9,
          ambientIntensity: 0.45,
          clipNear: 0,
          clipFar: 100,
          fogNear: 75,
          fogFar: 100,
          tooltip: false,
          mousePreset: 'default',
          hoverTimeout: 70,
        }) as unknown as Stage;
        if (!stage.viewer?.renderer)
          throw new Error(
            'WebGL could not be initialized. Enable hardware acceleration in your browser, then reload.',
          );
        stageRef.current = stage;
        // The application's measurement state is authoritative. NGL's defaults
        // create a second measurement buffer on right/Ctrl-click and recenter
        // the scene on a normal click, both surprising while picking atoms.
        stage.mouseControls.remove('clickPick-*');
        stage.mouseControls.remove('drag-ctrl-shift-left');
        stage.mouseControls.remove('drag-ctrl-shift-right');
      } catch (cause) {
        fail(cause instanceof Error ? cause.message : 'Unable to initialize the molecular viewer.');
        return;
      }

      const atomFromPick = (proxy?: PickProxy): number | null => {
        const atom = proxy?.atom || proxy?.closestBondAtom;
        if (atom) return atom.index;
        // Surface pickers identify a triangle; resolve its nearest actual atom.
        const component = componentRef.current;
        const dataset = latest.current.dataset;
        if (!proxy?.surface || !proxy.position || !component || !dataset) return null;
        let index: number | null = null;
        let shortest = Infinity;
        for (const metadata of dataset.atoms) {
          if (!atomIsVisible(metadata, latest.current.visibility)) continue;
          const current = component.structure.getAtomProxy(metadata.index);
          const distance =
            (current.x - proxy.position.x) ** 2 +
            (current.y - proxy.position.y) ** 2 +
            (current.z - proxy.position.z) ** 2;
          if (distance < shortest) {
            shortest = distance;
            index = metadata.index;
          }
        }
        return index;
      };
      const onClick = (proxy?: PickProxy) => {
        const index = atomFromPick(proxy);
        if (latest.current.picking && index !== null) latest.current.onAtomPick(index);
      };
      const onHover = (proxy?: PickProxy) => {
        const index = atomFromPick(proxy);
        hoveredIndex.current = index;
        const tooltip = tooltipRef.current;
        const atom = index !== null ? latest.current.dataset?.atoms[index] : null;
        if (!tooltip) return;
        if (!atom || !proxy?.canvasPosition) {
          tooltip.style.display = 'none';
          return;
        }
        tooltip.textContent = `${atom.residue} ${atom.resid} · ${atom.name}${atom.chain ? ` · chain ${atom.chain}` : ''}${latest.current.picking ? ' — click to select' : ''}`;
        tooltip.style.display = 'block';
        tooltip.style.left = `${Math.max(12, Math.min(stage.viewer.width - 250, proxy.canvasPosition.x + 14))}px`;
        tooltip.style.top = `${Math.max(12, stage.viewer.height - proxy.canvasPosition.y - 40)}px`;
      };
      const onDoubleClick = () => {
        if (!latest.current.picking && hoveredIndex.current !== null) focus([hoveredIndex.current]);
      };
      const onContextLost = (event: Event) => {
        event.preventDefault();
        fail(
          'The graphics context was lost. Reload the page to restore the viewer. Your simulation keeps running in the background.',
        );
      };
      stage.signals.clicked.add(onClick);
      stage.signals.hovered.add(onHover);
      stage.viewer.signals.rendered.add(updateLabels);
      stage.viewerControls.signals.changed.add(updateLabels);
      host.addEventListener('dblclick', onDoubleClick);
      host.addEventListener('webglcontextlost', onContextLost, true);
      let previousExtent = Math.min(host.clientWidth, host.clientHeight);
      const resize = new ResizeObserver(() => {
        const extent = Math.min(host.clientWidth, host.clientHeight);
        const distance = stage.viewerControls.getCameraDistance();
        stage.handleResize();
        // Opening the studio reduces the canvas. Preserve the visible molecular
        // field instead of cropping a zoomed scene to the smaller viewport.
        if (previousExtent > 0 && extent > 0 && componentRef.current && extent !== previousExtent) {
          stage.animationControls.zoom(Math.max(3, (distance * previousExtent) / extent), 0);
        }
        previousExtent = extent;
        updateLabels();
      });
      resize.observe(host);
      setStageVersion((version) => version + 1);
      return () => {
        resize.disconnect();
        stage.signals.clicked.remove(onClick);
        stage.signals.hovered.remove(onHover);
        stage.viewer.signals.rendered.remove(updateLabels);
        stage.viewerControls.signals.changed.remove(updateLabels);
        host.removeEventListener('dblclick', onDoubleClick);
        host.removeEventListener('webglcontextlost', onContextLost, true);
        stageRef.current = null;
        componentRef.current = null;
        groupsRef.current = null;
        stage.dispose();
        // NGL releases GPU resources but leaves its canvas in the host. React's
        // StrictMode remount must not retain that blank canvas above the new one.
        host.replaceChildren();
      };
    }, [fail, focus, updateLabels]);

    useEffect(() => {
      const stage = stageRef.current;
      const dataset = props.dataset;
      // Stage creation updates stageVersion after this effect's first pass.
      // Wait for that update so a new scene starts only one topology request.
      if (!stage || stageVersion === 0) return;
      if (componentRef.current) stage.removeComponent(componentRef.current);
      componentRef.current = null;
      groupsRef.current = null;
      interpolatorRef.current = null;
      activeDatasetId.current = null;
      measurementRepresentations.current.clear();
      if (tooltipRef.current) tooltipRef.current.style.display = 'none';
      if (!dataset) {
        setStatus('empty');
        return;
      }

      const controller = new AbortController();
      let obsolete = false;
      setStatus('loading');
      setError('');
      (async () => {
        const response = await fetch(dataset.topology_url, { signal: controller.signal });
        if (!response.ok)
          throw new Error(`Could not load molecular topology (HTTP ${response.status}).`);
        const blob = await response.blob();
        if (obsolete) return;
        const component = await stage.loadFile(blob, {
          ext: 'pdb',
          defaultRepresentation: false,
          firstModelOnly: true,
          name: dataset.name,
        });
        if (obsolete) {
          if (stageRef.current === stage) stage.removeComponent(component);
          return;
        }
        try {
          validateTopology(component, dataset);
        } catch (cause) {
          stage.removeComponent(component);
          throw cause;
        }
        componentRef.current = component;
        activeDatasetId.current = dataset.id;
        const add = (type: string, parameters: Record<string, unknown>) =>
          component.addRepresentation(type, { sele: 'none', ...parameters });
        groupsRef.current = {
          polymer: add('cartoon', {
            colorScheme: 'residueindex',
            colorScale: ['#69ddd3', '#6cb4f5', '#a989e5'],
            quality: 'high',
            aspectRatio: 4,
            radiusScale: 0.82,
            roughness: 0.32,
            metalness: 0.08,
          }),
          water: add('ball+stick', {
            colorScheme: 'element',
            radiusScale: 0.35,
            aspectRatio: 1.5,
            opacity: 0.55,
          }),
          ligands: add('ball+stick', {
            colorScheme: 'element',
            colorValue: '#f5bf68',
            radiusScale: 0.9,
            aspectRatio: 1.8,
            multipleBond: 'symmetric',
          }),
          ions: add('spacefill', { colorScheme: 'element', radiusScale: 0.6 }),
          detail: add('ball+stick', {
            colorScheme: 'element',
            colorValue: '#86b7c9',
            radiusScale: 0.55,
            aspectRatio: 1.65,
            opacity: 0.9,
          }),
          hydrogen: add('ball+stick', {
            colorScheme: 'element',
            radiusScale: 0.3,
            aspectRatio: 1.5,
          }),
          selected: add('spacefill', {
            color: '#ffc778',
            radiusType: 'size',
            radiusSize: 0.55,
            opacity: 0.9,
            disablePicking: true,
          }),
          selectedLabel: add('label', {
            color: '#ffdfa8',
            labelType: 'atomname',
            labelGrouping: 'atom',
            radiusType: 'size',
            radiusSize: 1,
            fixedSize: true,
            fontFamily: 'sans-serif',
            fontWeight: 'bold',
            showBackground: true,
            backgroundColor: '#142333',
            backgroundOpacity: 0.85,
            backgroundMargin: 0.5,
            yOffset: 1,
            zOffset: 1.8,
            disablePicking: true,
          }),
        };
        representationType.current = 'cartoon';
        const currentCoordinates =
          latest.current.dataset?.id === dataset.id ? latest.current.coordinates : null;
        if (currentCoordinates) {
          interpolatorRef.current = createCoordinateInterpolator(dataset, currentCoordinates);
          component.structure.updatePosition(interpolatorRef.current(latest.current.frame), false);
          component.updateRepresentations({ position: true });
        }
        const polymer = dataset.atoms
          .filter((atom) => atomGroup(atom) === 'protein')
          .map((atom) => atom.index);
        const solute = dataset.atoms
          .filter((atom) => atomGroup(atom) !== 'water')
          .map((atom) => atom.index);
        component.autoView(
          atomSelection(
            dataset.solvation && latest.current.visibility.water
              ? dataset.atoms.map((atom) => atom.index)
              : polymer.length
                ? polymer
                : solute.length
                  ? solute
                  : dataset.atoms.map((atom) => atom.index),
          ),
          0,
        );
        // Give the opening composition a little breathing room.
        stage.animationControls.zoom(stage.viewerControls.getCameraDistance() * 0.92, 0);
        stage.setSpin(latest.current.spin);
        setComponentVersion((version) => version + 1);
        setStatus('ready');
        latest.current.onReady();
      })().catch((cause) => {
        if (!obsolete && !(cause instanceof DOMException && cause.name === 'AbortError')) {
          fail(cause instanceof Error ? cause.message : String(cause));
        }
      });
      return () => {
        obsolete = true;
        controller.abort();
      };
    }, [props.dataset?.id, props.dataset?.topology_url, stageVersion, fail]);

    useEffect(() => {
      const dataset = props.dataset;
      const component = componentRef.current;
      if (!dataset || !component || activeDatasetId.current !== dataset.id) return;
      if (!props.coordinates) {
        interpolatorRef.current = null;
        return;
      }
      try {
        interpolatorRef.current = createCoordinateInterpolator(dataset, props.coordinates);
        component.structure.updatePosition(interpolatorRef.current(latest.current.frame), false);
        component.updateRepresentations({ position: true });
        updateLabels();
      } catch (cause) {
        fail(cause instanceof Error ? cause.message : String(cause));
      }
    }, [props.coordinates, props.dataset?.id, componentVersion, fail, updateLabels]);

    useLayoutEffect(() => {
      const component = componentRef.current;
      const interpolate = interpolatorRef.current;
      if (!component || !interpolate || activeDatasetId.current !== props.dataset?.id) return;
      // Never reload the structure, recreate representations, or fit the camera
      // during playback. This updates their existing position buffers in place.
      component.structure.updatePosition(interpolate(props.frame), false);
      component.updateRepresentations({ position: true });
      updateLabels();
    }, [props.frame, props.dataset?.id, componentVersion, updateLabels]);

    useEffect(() => {
      const component = componentRef.current;
      const dataset = props.dataset;
      const groups = groupsRef.current;
      if (!component || !dataset || !groups || activeDatasetId.current !== dataset.id) return;
      const visible = visibleGroups(dataset, props.visibility);
      const polymerSelection = atomSelection(visible.protein);
      if (representationType.current !== props.representation) {
        component.removeRepresentation(groups.polymer);
        groups.polymer = component.addRepresentation(props.representation, {
          sele: polymerSelection,
          quality: 'high',
          radiusScale: props.representation === 'cartoon' ? 0.82 : 0.65,
          aspectRatio: props.representation === 'cartoon' ? 4 : 1.8,
          surfaceType: 'av',
          probeRadius: 1.4,
          useWorker: true,
          opacity: props.representation === 'surface' ? 0.85 : 1,
          roughness: 0.32,
          metalness: 0.08,
        });
        representationType.current = props.representation;
      }
      groups.polymer.setSelection(polymerSelection);
      groups.polymer.setParameters({
        colorScheme: COLOR_SCHEMES[props.colorScheme],
        colorScale:
          props.colorScheme === 'chain'
            ? ['#69ddd3', '#a989e5', '#f6bd76', '#6daaf8']
            : ['#69ddd3', '#6cb4f5', '#a989e5'],
        colorValue: '#7ed7cd',
      });
      groups.water.setSelection(atomSelection(visible.water));
      groups.ligands.setSelection(atomSelection(visible.ligands));
      groups.ions.setSelection(atomSelection(visible.ions));
      groups.detail.setSelection(
        atomSelection(
          props.picking && ['cartoon', 'surface'].includes(props.representation)
            ? visible.protein
            : props.representation === 'cartoon'
              ? visible.protein.filter((index) => ribbonFallback.has(index))
              : [],
        ),
      );
      // Keep the whole heavy-atom skeleton behind the selected hydrogens.
      // Selecting only H and its donor leaves disconnected side-chain fragments.
      // Short fragments already have a complete atomic trace in groups.detail.
      groups.hydrogen.setSelection(
        atomSelection(
          props.representation === 'cartoon' &&
            !props.picking &&
            props.visibility.hydrogens !== 'none'
            ? visible.protein.filter((index) => !ribbonFallback.has(index))
            : [],
        ),
      );
      stageRef.current?.viewer.requestRender();
    }, [
      props.visibility,
      props.representation,
      props.colorScheme,
      props.picking,
      props.dataset?.id,
      componentVersion,
      ribbonFallback,
    ]);

    useEffect(() => {
      const dataset = props.dataset;
      const groups = groupsRef.current;
      if (!groups || !dataset || activeDatasetId.current !== dataset.id) return;
      const valid = props.selectedAtoms.filter(
        (index) =>
          Number.isInteger(index) &&
          index >= 0 &&
          index < dataset.n_atoms &&
          atomIsVisible(dataset.atoms[index], props.visibility),
      );
      groups.selected.setSelection(atomSelection(valid));
      groups.selectedLabel.setSelection(atomSelection(valid));
    }, [props.selectedAtoms, props.visibility, props.dataset?.id, componentVersion]);

    useEffect(() => {
      const component = componentRef.current;
      const dataset = props.dataset;
      if (!component || !dataset || activeDatasetId.current !== dataset.id) return;
      const representations = measurementRepresentations.current;
      const currentIds = new Set(props.measurements.map((measurement) => measurement.id));
      for (const [id, representation] of representations) {
        if (!currentIds.has(id)) {
          component.removeRepresentation(representation);
          representations.delete(id);
        }
      }
      for (const measurement of props.measurements) {
        if (
          representations.has(measurement.id) ||
          measurement.atoms.some((index) => index < 0 || index >= dataset.n_atoms)
        )
          continue;
        const parameters: Record<string, unknown> = {
          color: measurement.color,
          labelVisible: false,
          linewidth: 2,
          lineOpacity: 0.8,
          lineVisible: true,
          vectorVisible: true,
          sectorVisible: true,
          sectorOpacity: 0.14,
          planeVisible: false,
          useCylinder: false,
          disablePicking: true,
        };
        let type = measurement.kind;
        if (measurement.kind === 'distance') parameters.atomPair = [measurement.atoms];
        if (measurement.kind === 'angle') parameters.atomTriple = [measurement.atoms];
        if (measurement.kind === 'dihedral') parameters.atomQuad = [measurement.atoms];
        if (measurement.kind === 'hbond') {
          type = 'distance';
          parameters.atomPair = [[measurement.atoms[0], measurement.atoms[2]]];
        }
        representations.set(measurement.id, component.addRepresentation(type, parameters));
      }
      updateLabels();
    }, [props.measurements, props.dataset?.id, componentVersion, updateLabels]);

    useEffect(() => {
      stageRef.current?.setSpin(props.spin);
    }, [props.spin, stageVersion]);

    const savedFrame = Math.max(
      0,
      Math.min((props.dataset?.n_frames || 1) - 1, Math.round(props.frame)),
    );
    return (
      <div className={`molecular-viewer${props.picking ? ' molecular-viewer--picking' : ''}`}>
        <div
          ref={hostRef}
          className="molecular-viewer__canvas"
          role="img"
          aria-label={
            props.dataset
              ? `Interactive 3D molecular structure of ${props.dataset.name}. Drag to rotate, scroll to zoom, right-drag to pan, double-click an atom to focus.`
              : 'Molecular structure viewer'
          }
        />
        <div ref={tooltipRef} className="molecular-viewer__tooltip" role="tooltip" />
        {status === 'loading' && (
          <div className="molecular-viewer__status" role="status">
            <span className="molecular-viewer__spinner" />
            <strong>Bringing molecules into view</strong>
            <span>Preparing the topology and rendering layers…</span>
          </div>
        )}
        {status === 'empty' && (
          <div className="molecular-viewer__status molecular-viewer__status--empty">
            <div className="molecular-viewer__orbit">
              <i />
              <i />
              <i />
            </div>
            <strong>A world in motion.</strong>
            <span>Open a structure or trajectory to begin exploring.</span>
          </div>
        )}
        {status === 'error' && (
          <div className="molecular-viewer__status molecular-viewer__status--error" role="alert">
            <strong>The viewer needs your attention</strong>
            <span>{error}</span>
          </div>
        )}
        {status === 'ready' &&
          props.measurements.map((measurement) => {
            const value = measurement.values[savedFrame];
            return (
              <div
                key={measurement.id}
                ref={(element) => {
                  if (element) labelElements.current.set(measurement.id, element);
                  else labelElements.current.delete(measurement.id);
                }}
                className="molecular-viewer__measurement"
                style={{ borderColor: `${measurement.color}70`, color: measurement.color }}
                title={`${measurement.label} · saved frame ${savedFrame + 1}. The numeric value comes from saved-frame analysis; connecting lines follow the displayed coordinates.`}
              >
                <span>{measurement.kind === 'hbond' ? 'H bond' : measurement.kind}</span>
                <strong>
                  {Number.isFinite(value)
                    ? value.toFixed(
                        measurement.unit === '°' || measurement.unit === 'degrees' ? 1 : 2,
                      )
                    : '—'}{' '}
                  {measurement.unit}
                </strong>
              </div>
            );
          })}
        {status === 'ready' && props.dataset?.has_unitcell && props.dataset.n_frames > 1 && (
          <div
            className="molecular-viewer__periodic-note"
            title="Periodic coordinates play at saved frames. Smooth interpolation requires verified unwrapped display coordinates to avoid artificial jumps across the unit cell."
          >
            Periodic trajectory · saved-frame playback
          </div>
        )}
      </div>
    );
  },
);

export default MolecularViewer;
