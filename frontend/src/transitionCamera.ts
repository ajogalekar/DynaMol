import type { Dataset } from './types';

/** Preserve scene scale/rotation and compensate for engine box recentering.
 * NGL's orientation translation is applied before scene rotation. Coordinates
 * and physical geometry stay unchanged; only the camera's center is shifted.
 */
export function transitionCamera(
  camera: number[] | null,
  source: Dataset,
  sourceCoordinates: Float32Array | null,
  sourceFrame: number,
  target: Dataset,
  targetCoordinates: Float32Array,
  targetFrame: number,
): number[] | null {
  if (!camera || !sourceCoordinates) return camera;
  function center(d: Dataset, xyz: Float32Array, frame: number) {
    const backbone = d.atoms.filter((a) => a.category === 'protein' && a.name === 'CA');
    const atoms = backbone.length
      ? backbone
      : d.atoms.filter((a) => !['water', 'ions'].includes(a.category) && a.element !== 'H');
    if (!atoms.length) return null;
    const result = [0, 0, 0];
    const offset = Math.max(0, Math.min(d.n_frames - 1, frame)) * d.n_atoms * 3;
    for (const atom of atoms)
      for (let axis = 0; axis < 3; axis++)
        result[axis] += xyz[offset + atom.index * 3 + axis] / atoms.length;
    return result;
  }
  const before = center(source, sourceCoordinates, sourceFrame);
  const after = center(target, targetCoordinates, targetFrame);
  if (!before || !after) return camera;
  const next = [...camera];
  for (let axis = 0; axis < 3; axis++) next[12 + axis] += before[axis] - after[axis];
  return next;
}
