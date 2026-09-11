import type { AtomInfo, Dataset, Visibility } from '../types';

export const atomSelection = (indices: number[]) =>
  indices.length ? `@${indices.join(',')}` : 'none';

export function atomGroup(atom: AtomInfo): 'protein' | 'water' | 'ligands' | 'ions' {
  const category = atom.category.toLowerCase();
  if (['protein', 'polymer', 'nucleic', 'nucleic acid', 'rna', 'dna'].includes(category))
    return 'protein';
  if (['water', 'waters', 'solvent'].includes(category)) return 'water';
  if (['ion', 'ions'].includes(category)) return 'ions';
  return 'ligands';
}

export function atomIsVisible(atom: AtomInfo, visibility: Visibility) {
  if (!visibility[atomGroup(atom)]) return false;
  if (['H', 'D', 'T'].includes(atom.element.toUpperCase())) {
    if (visibility.hydrogens === 'none') return false;
    if (visibility.hydrogens === 'polar' && atom.nonpolar_hydrogen) return false;
  }
  return true;
}

export function visibleGroups(dataset: Dataset, visibility: Visibility) {
  const groups = {
    protein: [] as number[],
    water: [] as number[],
    ligands: [] as number[],
    ions: [] as number[],
  };
  dataset.atoms.forEach((atom) => {
    if (atomIsVisible(atom, visibility)) groups[atomGroup(atom)].push(atom.index);
  });
  return groups;
}

/**
 * Interpolation belongs only to rendering. The input buffer remains immutable;
 * scientific measurements are obtained from the backend's saved coordinates.
 * Bonded molecules snap together if a large displacement suggests a PBC jump.
 */
export function createCoordinateInterpolator(dataset: Dataset, coordinates: Float32Array) {
  const stride = dataset.n_atoms * 3;
  const expectedLength = stride * dataset.n_frames;
  if (coordinates.length !== expectedLength) {
    throw new Error(
      `Coordinate count mismatch: expected ${expectedLength.toLocaleString()} values, received ${coordinates.length.toLocaleString()}.`,
    );
  }
  if (coordinates.some((value) => !Number.isFinite(value))) {
    throw new Error(
      'The trajectory contains non-finite coordinates and cannot be rendered safely.',
    );
  }
  const parents = Int32Array.from({ length: dataset.n_atoms }, (_, index) => index);
  const root = (index: number): number => {
    let current = index;
    while (parents[current] !== current) {
      parents[current] = parents[parents[current]];
      current = parents[current];
    }
    return current;
  };
  for (const [a, b] of dataset.bonds) {
    if (a >= 0 && b >= 0 && a < dataset.n_atoms && b < dataset.n_atoms) parents[root(b)] = root(a);
  }
  for (let atom = 0; atom < dataset.n_atoms; atom++) parents[atom] = root(atom);

  const output = new Float32Array(stride);
  const jumpGroups = new Uint8Array(dataset.n_atoms);
  let lastLowerFrame = -1;

  return (requestedFrame: number): Float32Array => {
    const bounded = Math.min(
      dataset.n_frames - 1,
      Math.max(0, Number.isFinite(requestedFrame) ? requestedFrame : 0),
    );
    // Unit-cell vectors and a verified unwrapping policy are not yet in the API.
    // Do not invent minimum-image paths for periodic coordinates.
    const frame = dataset.has_unitcell ? Math.round(bounded) : bounded;
    const lowerFrame = Math.floor(frame);
    const upperFrame = Math.min(lowerFrame + 1, dataset.n_frames - 1);
    const fraction = frame - lowerFrame;
    const lower = lowerFrame * stride;
    const upper = upperFrame * stride;

    if (!fraction || lowerFrame === upperFrame) {
      output.set(coordinates.subarray(lower, lower + stride));
      return output;
    }
    if (lastLowerFrame !== lowerFrame) {
      jumpGroups.fill(0);
      for (let atom = 0; atom < dataset.n_atoms; atom++) {
        const at = atom * 3;
        const dx = coordinates[upper + at] - coordinates[lower + at];
        const dy = coordinates[upper + at + 1] - coordinates[lower + at + 1];
        const dz = coordinates[upper + at + 2] - coordinates[lower + at + 2];
        if (dx * dx + dy * dy + dz * dz > 64) jumpGroups[parents[atom]] = 1;
      }
      lastLowerFrame = lowerFrame;
    }
    for (let atom = 0; atom < dataset.n_atoms; atom++) {
      const blend = jumpGroups[parents[atom]] ? (fraction < 0.5 ? 0 : 1) : fraction;
      for (let axis = 0; axis < 3; axis++) {
        const at = atom * 3 + axis;
        output[at] =
          coordinates[lower + at] + (coordinates[upper + at] - coordinates[lower + at]) * blend;
      }
    }
    return output;
  };
}
