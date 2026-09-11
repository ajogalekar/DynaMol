import { test } from './testWorkspace';
import { expect, type Locator, type Page } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

const fixtures = fileURLToPath(new URL('../../docs/audit/preparation-fixtures/', import.meta.url));
const heightKey = 'dynamol.viewer-height';

async function ready(page: Page) {
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
}

async function studio(page: Page) {
  await page.goto('/');
  await ready(page);
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  return page.getByRole('dialog', { name: 'Set molecules in motion.' });
}

async function upload(page: Page, panel: Locator, filename: string) {
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/structures/upload') && r.request().method() === 'POST',
  );
  await panel.getByLabel('Upload simulation structure', { exact: true }).setInputFiles({
    name: filename,
    mimeType: 'chemical/x-pdb',
    buffer: await readFile(`${fixtures}/${filename}`),
  });
  const received = await response;
  expect(received.ok(), await received.text()).toBe(true);
  const dataset = await received.json();
  await expect(panel.locator('.source-summary strong')).toHaveText(dataset.name);
  await ready(page);
  return dataset;
}

async function onscreen(page: Page, target: Locator) {
  await expect(target).toBeVisible();
  const box = await target.boundingBox();
  const view = page.viewportSize()!;
  expect(box).not.toBeNull();
  expect(box!.x).toBeGreaterThanOrEqual(0);
  expect(box!.y).toBeGreaterThanOrEqual(0);
  expect(box!.x + box!.width).toBeLessThanOrEqual(view.width + 1);
  expect(box!.y + box!.height).toBeLessThanOrEqual(view.height + 1);
}

async function canvasFits(page: Page) {
  await expect
    .poll(async () => {
      const outer = await page.locator('#molecular-viewport').boundingBox();
      const canvas = await page.locator('.molecular-viewer__canvas canvas').boundingBox();
      if (!outer || !canvas) return false;
      return (
        canvas.width > 100 &&
        canvas.height > 100 &&
        canvas.x >= outer.x - 1 &&
        canvas.y >= outer.y - 1 &&
        canvas.x + canvas.width <= outer.x + outer.width + 1 &&
        canvas.y + canvas.height <= outer.y + outer.height + 1 &&
        Math.abs(canvas.height - outer.height) < 3
      );
    })
    .toBe(true);
  const outer = (await page.locator('#molecular-viewport').boundingBox())!;
  for (const button of await page.locator('.canvas-tools button').all()) {
    const box = (await button.boundingBox())!;
    expect(box.y).toBeGreaterThanOrEqual(outer.y);
    expect(box.y + box.height).toBeLessThanOrEqual(outer.y + outer.height);
  }
}

test('loop preference remains selectable on complete protein and still enables a supported known gap', async ({
  page,
}, info) => {
  const panel = await studio(page);
  await upload(page, panel, 'six_residues_intact.pdb');
  const repairs = panel.getByRole('region', { name: 'Missing atoms & residues' });
  const loop = repairs.getByRole('checkbox', {
    name: 'Build supported missing loops / residues',
    exact: true,
  });
  await expect(loop).toBeEnabled();
  await expect(repairs).toContainText('No missing protein sequence regions were detected.');
  await expect(repairs).toContainText('ligand atom repair is separate below');
  await loop.check();
  await expect(loop).toBeChecked();
  await expect(panel.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  await loop.uncheck();
  await upload(page, panel, 'six_residues_known_gap.pdb');
  await expect(repairs.locator('.missing-residues')).toContainText('ILE');
  await expect(repairs).toContainText('can build a starting model');
  await expect(loop).toBeEnabled();
  await loop.check();
  await expect(panel.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  await expect(
    panel.getByRole('button', { name: 'Preparation options', exact: true }),
  ).toHaveAttribute('aria-expanded', 'false');
  await repairs.scrollIntoViewIfNeeded();
  await info.attach('visible-loop-options', {
    body: await page.screenshot(),
    contentType: 'image/png',
  });
});

test('protein viewport drag and keyboard size persist, reset, and keep the canvas fitted', async ({
  page,
}, info) => {
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(e.message));
  await page.goto('/');
  await ready(page);
  const separator = page.getByRole('separator', { name: 'Resize protein view', exact: true });
  const viewport = page.locator('#molecular-viewport');
  const minimum = Number(await separator.getAttribute('aria-valuemin'));
  await separator.focus();
  await separator.press('Home');
  await expect(separator).toHaveAttribute('aria-valuenow', String(minimum));
  await canvasFits(page);
  const before = await viewport.boundingBox();
  const grip = await separator.boundingBox();
  await page.mouse.move(grip!.x + grip!.width / 2, grip!.y + grip!.height / 2);
  await page.mouse.down();
  await page.mouse.move(grip!.x + grip!.width / 2, grip!.y + grip!.height / 2 + 140, { steps: 10 });
  await page.mouse.up();
  await expect(separator).toHaveAttribute('aria-valuenow', String(minimum + 140));
  expect((await viewport.boundingBox())!.width).toBeCloseTo(before!.width, 0);
  await canvasFits(page);
  await separator.press('ArrowDown');
  await expect(separator).toHaveAttribute('aria-valuenow', String(minimum + 160));
  await separator.press('Shift+ArrowUp');
  await expect(separator).toHaveAttribute('aria-valuenow', String(minimum + 100));
  await expect
    .poll(() => page.evaluate((key) => localStorage.getItem(key), heightKey))
    .toBe(String(minimum + 100));
  await page.reload();
  await ready(page);
  await expect(separator).toHaveAttribute('aria-valuenow', String(minimum + 100));
  await canvasFits(page);
  await separator.focus();
  await separator.press('End');
  await expect(separator).toHaveAttribute('aria-valuenow', '770');
  await canvasFits(page);
  await page.setViewportSize({ width: 1280, height: 700 });
  await expect(separator).toHaveAttribute('aria-valuenow', '470');
  await canvasFits(page);
  await separator.press('Enter');
  await expect.poll(() => page.evaluate((key) => localStorage.getItem(key), heightKey)).toBeNull();
  await expect(viewport).not.toHaveAttribute('style', /height/);
  await canvasFits(page);
  await separator.press('Home');
  await expect(separator).toHaveAttribute('aria-valuenow', String(minimum));
  await separator.dblclick();
  await expect.poll(() => page.evaluate((key) => localStorage.getItem(key), heightKey)).toBeNull();
  await expect(viewport).not.toHaveAttribute('style', /height/);
  await canvasFits(page);
  await info.attach('resized-viewport', {
    body: await page.screenshot(),
    contentType: 'image/png',
  });
  expect(errors).toEqual([]);
});

for (const size of [
  { width: 1440, height: 1000 },
  { width: 1280, height: 600 },
]) {
  test(`full screen retains top actions and fitted protein view at ${size.width}x${size.height}`, async ({
    page,
  }, info) => {
    await page.setViewportSize(size);
    await page.goto('/');
    await ready(page);
    await page.getByRole('button', { name: 'Full screen', exact: true }).click();
    await expect
      .poll(() => page.evaluate(() => document.fullscreenElement?.classList.contains('app-shell')))
      .toBe(true);
    const open = page.getByRole('button', { name: 'Open files', exact: true });
    const simulation = page.getByRole('button', { name: 'New simulation', exact: true });
    await onscreen(page, open);
    await onscreen(page, simulation);
    await simulation.click();
    await expect(page.getByRole('dialog', { name: 'Set molecules in motion.' })).toBeVisible();
    await onscreen(page, open);
    await onscreen(page, page.getByRole('button', { name: 'Explore', exact: true }));
    await onscreen(page, page.getByRole('button', { name: 'Simulate', exact: true }));
    const separator = page.getByRole('separator', { name: 'Resize protein view', exact: true });
    await separator.focus();
    await separator.press('Home');
    await expect(separator).toHaveAttribute(
      'aria-valuenow',
      (await separator.getAttribute('aria-valuemin'))!,
    );
    await canvasFits(page);
    await info.attach('fullscreen-simulation-header', {
      body: await page.screenshot(),
      contentType: 'image/png',
    });
    await page.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
    await open.click();
    await expect(page.getByRole('dialog', { name: 'Bring your molecules.' })).toBeVisible();
    await page.getByRole('button', { name: 'Close import', exact: true }).click();
    await page.getByRole('button', { name: 'Exit full screen', exact: true }).click();
    await expect.poll(() => page.evaluate(() => document.fullscreenElement === null)).toBe(true);
    await onscreen(page, open);
    await onscreen(page, simulation);
  });
}

test('real incomplete ligands allow individual repair and removal with matching readiness and submission payload', async ({
  page,
  request,
}, info) => {
  const datasetId = '6ede70d3302d4645';
  const source = await request.get(`/api/datasets/${datasetId}`);
  expect(source.ok(), 'Copy the isolated 6A93 monomer fixture before this audit.').toBe(true);
  const original = await source.json();
  const saved = await request.post('/api/workspace', {
    data: {
      state: { version: 1, dataset_id: datasetId, measurements: [], active_measurement: null },
    },
  });
  expect(saved.ok(), await saved.text()).toBe(true);
  const submissions: Record<string, unknown>[] = [];
  // Only intercept submission: inspection and readiness use the real isolated API.
  // No preparation worker or native molecular dynamics is launched by this audit.
  await page.route('**/api/preparations', async (route) => {
    if (route.request().method() !== 'POST') return route.continue();
    submissions.push(route.request().postDataJSON());
    await route.fulfill({
      status: 400,
      contentType: 'application/json',
      body: JSON.stringify({
        detail: 'Audit captured preparation settings; computation was not started.',
      }),
    });
  });
  const panel = await studio(page);
  const loop = panel.getByRole('checkbox', {
    name: 'Build supported missing loops / residues',
    exact: true,
  });
  await expect(loop).toBeEnabled();
  await expect(panel.getByRole('region', { name: 'Missing atoms & residues' })).toContainText(
    'No missing protein sequence regions were detected.',
  );
  const firstKey = 'F:3004::1PE';
  const secondKey = 'H:3006::1PE';
  const first = panel.getByRole('combobox', {
    name: `Preparation choice ${firstKey}`,
    exact: true,
  });
  const second = panel.getByRole('combobox', {
    name: `Preparation choice ${secondKey}`,
    exact: true,
  });
  await expect(first).toBeVisible();
  await expect(second).toBeVisible();
  await expect(first).toHaveValue('');
  const prepare = panel.getByRole('button', { name: 'Prep complex', exact: true });
  await expect(prepare).toBeDisabled();
  await expect(panel.locator('[aria-label="Preparation readiness"]')).toContainText(
    'Before you continue',
  );
  await expect(panel.locator('.ion-preparation-summary')).toContainText('Zn+2');
  const inspection = page.waitForResponse(
    (r) =>
      r.url().endsWith(`/api/datasets/${datasetId}/inspection`) &&
      r.request().method() === 'POST' &&
      r.request().postDataJSON().ligand_actions[firstKey] === 'repair',
  );
  await panel.getByRole('button', { name: 'Repair missing ligand atoms', exact: true }).click();
  await expect(first).toHaveValue('repair');
  await expect(second).toHaveValue('repair');
  const inspected = await (await inspection).json();
  expect(
    inspected.ligands.filter((l: { selected_action?: string }) => l.selected_action === 'repair'),
  ).toHaveLength(2);
  await expect(prepare).toBeEnabled();
  const expectedActions = { [firstKey]: 'repair', [secondKey]: 'remove' };
  const readiness = page.waitForResponse(
    (r) =>
      r.url().endsWith(`/api/datasets/${datasetId}/readiness`) &&
      r.request().postDataJSON()?.mode === 'preparation' &&
      r.request().postDataJSON().settings.ligand_actions[secondKey] === 'remove',
  );
  await second.selectOption('remove');
  const readyResponse = await readiness;
  expect(readyResponse.request().postDataJSON().settings.ligand_actions).toEqual(expectedActions);
  expect((await readyResponse.json()).ready).toBe(true);
  await expect(prepare).toBeEnabled();
  await expect(
    panel
      .locator('.ligand-preparation-card')
      .filter({
        has: page.getByRole('combobox', { name: `Preparation choice ${secondKey}`, exact: true }),
      }),
  ).toContainText('Removal selected');
  const blocked = page.waitForResponse(
    (r) =>
      r.url().endsWith(`/api/datasets/${datasetId}/readiness`) &&
      r.request().postDataJSON()?.mode === 'preparation' &&
      !(firstKey in r.request().postDataJSON().settings.ligand_actions),
  );
  await first.selectOption('');
  expect((await (await blocked).json()).ready).toBe(false);
  await expect(prepare).toBeDisabled();
  await first.selectOption('repair');
  await expect(prepare).toBeEnabled();
  await prepare.click();
  await expect.poll(() => submissions.length).toBe(1);
  expect(submissions[0].dataset_id).toBe(datasetId);
  expect(submissions[0].ligand_actions).toEqual(expectedActions);
  expect(submissions[0].remove_heterogens).toBe(false);
  await expect(panel).toContainText(
    'Audit captured preparation settings; computation was not started.',
  );
  const after = await (await request.get(`/api/datasets/${datasetId}`)).json();
  expect(after).toEqual(original);
  await first.scrollIntoViewIfNeeded();
  await info.attach('individual-ligand-actions', {
    body: await page.screenshot(),
    contentType: 'image/png',
  });
});
