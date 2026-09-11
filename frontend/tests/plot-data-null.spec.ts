import { test, expect } from '@playwright/test';
import { createPlotData } from '../src/components/PlotPanel';
import type { Measurement } from '../src/types';

const measurement = (overrides: Partial<Measurement>): Measurement => ({
  id: 'saved-frame-geometry',
  kind: 'distance',
  atoms: [0, 1],
  label: 'Observed geometry',
  color: '#74dfc5',
  unit: 'Å',
  values: [],
  times_ps: [],
  ...overrides,
});

test('missing geometry retains physical timestamps and breaks the trace without changing the mean', () => {
  const plot = createPlotData(measurement({ values: [2, null, 4, 6], times_ps: [0, 0.1, 0.7, 1] }));
  expect(plot).not.toBeNull();
  expect(plot!.mean).toBe(4);
  expect(plot!.lo).toBe(2);
  expect(plot!.hi).toBe(6);
  expect(plot!.segments).toHaveLength(2);
  expect(plot!.segments[0].trim().split(' ')).toHaveLength(1);
  expect(plot!.segments[1].trim().split(' ')).toHaveLength(2);
  expect((plot!.x(2) - plot!.x(0)) / (plot!.x(3) - plot!.x(0))).toBeCloseTo(0.7);
  expect(plot!.nearestIndex(0.1)).toBe(1);
  expect(plot!.segments.join(' ')).not.toMatch(/NaN|null|undefined/);
});

test('empty and undefined series remain distinct from finite torsions across the branch cut', () => {
  expect(createPlotData(measurement({}))).toBeNull();
  expect(createPlotData(measurement({ values: [null, null], times_ps: [0, 1] }))).toBeNull();
  const plot = createPlotData(
    measurement({
      kind: 'dihedral',
      unit: '°',
      values: [179, -179, null, 178],
      times_ps: [0, 1, 2, 3],
    }),
  );
  expect(plot!.segments).toHaveLength(3);
  expect(Math.abs(plot!.mean!)).toBeGreaterThan(178);
  expect(plot!.segments.join(' ')).not.toMatch(/NaN|null|undefined/);
});
