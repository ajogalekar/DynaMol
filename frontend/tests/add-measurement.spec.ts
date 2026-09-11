import { test } from './testWorkspace';
import { expect, type Page, type APIRequestContext } from '@playwright/test';

type Atom = { index: number; name: string; residue: string; resid: number; element: string };
async function pick(page: Page, atom: Atom) {
  await page
    .getByRole('textbox', { name: 'Find an atom', exact: true })
    .fill(`${atom.residue}${atom.resid} · ${atom.name}`);
  await page
    .locator('.atom-search-results')
    .getByRole('button', { name: new RegExp(`#${atom.index + 1} · ${atom.element}$`) })
    .click();
}
async function ready(page: Page, request: APIRequestContext) {
  const demo = await (await request.get('/api/datasets/demo')).json();
  const saved = await request.post('/api/workspace', {
    data: {
      state: { version: 1, dataset_id: demo.id, measurements: [], active_measurement: null },
    },
  });
  expect(saved.ok()).toBeTruthy();
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await expect(page.locator('.measurement-tab')).toHaveCount(0);
  return {
    demo,
    atoms: demo.atoms.filter((atom: Atom) => atom.name === 'CA').slice(0, 4) as Atom[],
  };
}
async function kind(page: Page, name: string) {
  await page
    .locator('.measure-types')
    .getByRole('button', { name: new RegExp(`${name}$`, 'i') })
    .click();
}

test('bottom Add automatically plots distance, angle and dihedral and restores saved plots', async ({
  page,
  request,
}, testInfo) => {
  const { atoms } = await ready(page, request);
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  const submitted: unknown[] = [];
  page.on('request', (req) => {
    if (req.url().endsWith('/measurements') && req.method() === 'POST')
      submitted.push(req.postDataJSON());
  });
  const results: Record<string, unknown>[] = [];
  for (const [index, [name, count]] of (
    [
      ['distance', 2],
      ['angle', 3],
      ['dihedral', 4],
    ] as const
  ).entries()) {
    await page
      .getByRole('button', { name: index ? 'Measure' : 'Add a measurement', exact: true })
      .click();
    await kind(page, name);
    await expect(page.locator('.measurement-draft')).toContainText(`0/${count} atoms selected`);
    await pick(page, atoms[0]);
    await expect(page.locator('.measurement-draft')).toContainText(`1/${count} atoms selected`);
    if (!index)
      await page.screenshot({ path: testInfo.outputPath('measurement-draft.png'), fullPage: true });
    const response = page.waitForResponse(
      (res) => res.url().endsWith('/measurements') && res.request().postDataJSON()?.kind === name,
    );
    for (const atom of atoms.slice(1, count)) await pick(page, atom);
    const measured = await (await response).json();
    results.push(measured);
    await expect(page.locator('.measurement-draft')).toHaveCount(0);
    await expect(page.locator('.measurement-tab')).toHaveCount(index + 1);
    await expect(page.locator('.measurement-tab.active')).toHaveCount(1);
    await expect(page.locator('.plot-readout strong')).toHaveText(
      `${measured.values[0].toFixed(2)}${measured.unit}`,
    );
    await expect(
      page.locator('.molecular-viewer__measurement').filter({ hasText: measured.unit }).first(),
    ).toBeVisible();
  }
  expect(submitted).toEqual(
    ['distance', 'angle', 'dihedral'].map((name, i) => ({
      kind: name,
      atoms: atoms.slice(0, i + 2).map((a) => a.index),
    })),
  );
  await expect
    .poll(
      async () => (await (await request.get('/api/workspace')).json()).state.measurements.length,
    )
    .toBe(3);
  const saved = (await (await request.get('/api/workspace')).json()).state;
  expect(
    saved.measurements.map(({ kind, atoms, values }: Record<string, unknown>) => ({
      kind,
      atoms,
      values,
    })),
  ).toEqual(results.map(({ kind, atoms, values }) => ({ kind, atoms, values })));
  expect(saved.active_measurement).toBe(saved.measurements[2].id);
  await page.screenshot({ path: testInfo.outputPath('three-measurements.png'), fullPage: true });
  await page.reload();
  await expect(page.locator('.measurement-tab')).toHaveCount(3);
  await expect(page.locator('.plot-readout strong')).toHaveText(
    `${(results[2].values as number[])[0].toFixed(2)}°`,
  );
  expect(errors).toEqual([]);
});

test('bottom draft shows a failed calculation and retries once successfully', async ({
  page,
  request,
}) => {
  const { atoms } = await ready(page, request);
  let calls = 0;
  await page.route('**/measurements', async (route) => {
    calls++;
    if (calls === 1)
      await route.fulfill({ status: 503, json: { detail: 'Calculation temporarily unavailable' } });
    else await route.continue();
  });
  await page.getByRole('button', { name: 'Add a measurement', exact: true }).click();
  for (const atom of atoms.slice(0, 2)) await pick(page, atom);
  await expect(page.locator('.measurement-draft [role="alert"]')).toContainText(
    'Calculation temporarily unavailable',
  );
  await expect(page.locator('.measurement-tab')).toHaveCount(0);
  expect(calls).toBe(1);
  await page.getByRole('button', { name: 'Retry adding measurement', exact: true }).click();
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  await expect(page.locator('.measurement-draft')).toHaveCount(0);
  expect(calls).toBe(2);
});

test('cancel and changed measurement types discard pending results', async ({ page, request }) => {
  const { atoms } = await ready(page, request);
  let release: () => void = () => {};
  let responseHeld = false;
  let delay = true;
  await page.route('**/measurements', async (route) => {
    if (!delay) return route.continue();
    delay = false;
    const response = await route.fetch();
    await new Promise<void>((resolve) => {
      release = resolve;
      responseHeld = true;
    });
    await route.fulfill({ response }).catch(() => {});
  });
  await page.getByRole('button', { name: 'Add a measurement', exact: true }).click();
  for (const atom of atoms.slice(0, 2)) await pick(page, atom);
  await expect(page.locator('.measurement-draft')).toContainText('Adding distance plot…');
  // Wait for the held real backend response, then cancel it through the UI.
  await expect.poll(() => responseHeld).toBe(true);
  await page.getByRole('button', { name: 'Cancel new measurement', exact: true }).click();
  release();
  await expect(page.locator('.measurement-draft')).toHaveCount(0);
  await expect(page.locator('.measurement-tab')).toHaveCount(0);

  delay = true;
  responseHeld = false;
  await page.getByRole('button', { name: 'Add a measurement', exact: true }).click();
  for (const atom of atoms.slice(0, 2)) await pick(page, atom);
  await expect(page.locator('.measurement-draft')).toContainText('Adding distance plot…');
  await expect.poll(() => responseHeld).toBe(true);
  await kind(page, 'angle');
  release();
  await expect(page.locator('.measurement-draft')).toContainText('Angle · 0/3 atoms selected');
  for (const atom of atoms.slice(0, 3)) await pick(page, atom);
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  await expect(page.locator('.plot-readout strong small')).toHaveText('°');
  await expect
    .poll(async () =>
      (await (await request.get('/api/workspace')).json()).state.measurements.map(
        (m: { kind: string }) => m.kind,
      ),
    )
    .toEqual(['angle']);
});

test('stopping or escaping an Add draft keeps subsequent inspection from creating a plot', async ({
  page,
  request,
}) => {
  const { atoms } = await ready(page, request);
  await page.getByRole('button', { name: 'Add a measurement', exact: true }).click();
  await pick(page, atoms[0]);
  await page.getByRole('button', { name: 'Stop picking', exact: true }).click();
  await expect(page.locator('.measurement-draft')).toHaveCount(0);
  await pick(page, atoms[1]);
  await expect(page.locator('.live-measurement__value')).toBeVisible();
  await expect(page.locator('.measurement-tab')).toHaveCount(0);
  await page.getByRole('button', { name: 'Add a measurement', exact: true }).click();
  await expect(page.locator('.measurement-draft')).toContainText('0/2 atoms selected');
  await page.keyboard.press('Escape');
  await expect(page.locator('.measurement-draft')).toHaveCount(0);
  await expect(page.locator('.atom-slot.filled')).toHaveCount(0);
});
