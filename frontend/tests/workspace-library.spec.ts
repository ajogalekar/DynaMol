import { expect, type Page, type APIRequestContext } from '@playwright/test';
import { test } from './testWorkspace';
import fs from 'node:fs/promises';

// Run against an isolated DYNAMOL_DATA_DIR server. These tests intentionally
// rename/archive/trash disposable uploads, and never operate on user datasets.
async function seed(request: APIRequestContext) {
  const demo = await (await request.get('/api/datasets/demo')).json();
  const response = await request.post('/api/workspace', {
    data: {
      state: { version: 1, dataset_id: demo.id, frame: 0, measurements: [], named_selections: [] },
    },
  });
  expect(response.ok(), await response.text()).toBeTruthy();
  return demo;
}
async function ready(page: Page) {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
}
async function pick(
  page: Page,
  atom: { index: number; name: string; residue: string; resid: number; element: string },
) {
  await page
    .getByRole('textbox', { name: 'Find an atom', exact: true })
    .fill(`${atom.residue}${atom.resid} · ${atom.name}`);
  await page
    .locator('.atom-search-results')
    .getByRole('button', { name: new RegExp(`#${atom.index + 1} · ${atom.element}$`) })
    .click();
}
async function current(request: APIRequestContext) {
  return (await (await request.get('/api/workspace')).json()).state;
}
async function projects(page: Page) {
  await page.getByRole('button', { name: 'Projects', exact: true }).click();
  await expect(page.locator('.workspace-library__busy')).toHaveCount(0);
}

test('workspace survives a reload with camera, frame, rendering, plots and named atom selections', async ({
  page,
  request,
}, info) => {
  const demo = await seed(request);
  const errors: string[] = [];
  page.on('pageerror', (e) => errors.push(e.message));
  await ready(page);
  const atoms = demo.atoms.filter((a: { name: string }) => a.name === 'CA').slice(0, 2);
  for (const atom of atoms) await pick(page, atom);
  await expect(page.locator('.live-measurement__value')).toBeVisible();
  await page.getByRole('button', { name: 'Plot over time', exact: true }).click();
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  for (const atom of atoms) await pick(page, atom);
  await page.getByRole('textbox', { name: 'Selection name', exact: true }).fill('Backbone anchors');
  await page.getByRole('button', { name: 'Save named selection', exact: true }).click();
  await page.locator('.trajectory-analysis > summary').click();
  await page
    .getByRole('combobox', { name: 'Structural analysis kind', exact: true })
    .selectOption('rmsf');
  await page
    .getByRole('combobox', { name: 'Analysis atom group', exact: true })
    .selectOption({ label: 'Backbone anchors (2 atoms)' });
  await page.getByRole('spinbutton', { name: 'Analysis reference frame', exact: true }).fill('2');
  await page.getByRole('button', { name: 'Sticks', exact: true }).click();
  await page
    .getByRole('combobox', { name: 'Color molecules by', exact: true })
    .selectOption('chain');
  await page.getByRole('button', { name: 'Polar only', exact: true }).click();
  await page.getByRole('slider', { name: 'Trajectory frame', exact: true }).fill('12');
  await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
  await expect
    .poll(async () => {
      const state = await current(request);
      return (
        state?.frame === 12 &&
        state?.representation === 'licorice' &&
        state?.named_selections.length === 1 &&
        state?.measurements.length === 1 &&
        state?.camera?.length === 16
      );
    })
    .toBeTruthy();
  await page.waitForTimeout(1000); // let the actual camera's 180 ms animation finish and autosave
  const before = await current(request);
  await page.reload();
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.getByRole('slider', { name: 'Trajectory frame', exact: true })).toHaveValue(
    '12',
  );
  await expect(page.getByRole('button', { name: 'Sticks', exact: true })).toHaveClass(/selected/);
  await expect(page.getByRole('button', { name: 'Polar only', exact: true })).toHaveClass(/active/);
  await expect(page.getByRole('combobox', { name: 'Color molecules by', exact: true })).toHaveValue(
    'chain',
  );
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  await expect(page.locator('.named-selections__row')).toContainText('Backbone anchors');
  await expect(page.locator('.trajectory-analysis')).toHaveAttribute('open', '');
  await expect(
    page.getByRole('combobox', { name: 'Structural analysis kind', exact: true }),
  ).toHaveValue('rmsf');
  await expect(
    page.getByRole('spinbutton', { name: 'Analysis reference frame', exact: true }),
  ).toHaveValue('2');
  await expect(page.getByRole('button', { name: 'Play trajectory', exact: true })).toBeVisible();
  await page.waitForTimeout(1000);
  const after = await current(request);
  expect(
    Math.max(...before.camera.map((value: number, i: number) => Math.abs(value - after.camera[i]))),
  ).toBeLessThan(0.0001);
  expect(after.measurements).toEqual(before.measurements);
  expect(after.selected_atoms).toEqual(before.selected_atoms);
  expect(after.analysis_settings).toEqual(before.analysis_settings);
  expect(errors).toEqual([]);
  await page.screenshot({ path: info.outputPath('restored-workspace.png'), fullPage: true });
  await fs.writeFile(
    info.outputPath('restored-workspace.json'),
    JSON.stringify({ before, after, errors }, null, 2),
  );
});

test('named project save, update, reopen, export and import preserve a reviewable scene', async ({
  page,
  request,
}, info) => {
  await seed(request);
  await ready(page);
  const name = `Workspace roundtrip ${Date.now()}`;
  await page.getByRole('slider', { name: 'Trajectory frame', exact: true }).fill('9');
  await projects(page);
  await page.getByRole('textbox', { name: 'New project name', exact: true }).fill(name);
  await page.getByRole('button', { name: 'Save new project', exact: true }).click();
  await expect(page.locator('.workspace-library__notice')).toContainText('Project saved');
  const row = page
    .locator('.workspace-library__row')
    .filter({ has: page.getByText(name, { exact: true }) });
  await expect(row).toHaveCount(1);
  const download = page.waitForEvent('download');
  await row.getByRole('button', { name: `Export project ${name}`, exact: true }).click();
  const zip = await download;
  const path = info.outputPath('DynaMol-roundtrip.zip');
  await zip.saveAs(path);
  await expect(page.locator('.workspace-library__busy')).toHaveCount(0);
  const exported = (await (await request.get('/api/projects')).json()).find(
    (p: { name: string }) => p.name === name,
  );
  await page.getByRole('button', { name: 'Close workspace library', exact: true }).click();
  await page.getByRole('slider', { name: 'Trajectory frame', exact: true }).fill('3');
  await projects(page);
  await row.getByRole('button', { name: 'Open', exact: true }).click();
  await expect(page.getByRole('slider', { name: 'Trajectory frame', exact: true })).toHaveValue(
    '9',
  );
  await page.getByRole('slider', { name: 'Trajectory frame', exact: true }).fill('4');
  await projects(page);
  await row.getByRole('button', { name: `Update project ${name}`, exact: true }).click();
  await page
    .locator('.workspace-library__confirm')
    .getByRole('button', { name: 'Confirm', exact: true })
    .click();
  await expect(page.locator('.workspace-library__notice')).toContainText('updated');
  expect((await (await request.get(`/api/projects/${exported.id}`)).json()).state.frame).toBe(4);
  await page.getByLabel('Import project backup file', { exact: true }).setInputFiles(path);
  await expect(page.locator('.workspace-library__notice')).toContainText('Imported');
  await expect(row).toHaveCount(2);
  const projectList = await (await request.get('/api/projects')).json();
  const imported = projectList.find(
    (p: { id: string; name: string }) => p.id !== exported.id && p.name === name,
  );
  expect((await (await request.get(`/api/projects/${imported.id}`)).json()).state.frame).toBe(9);
  await row
    .first()
    .getByRole('button', { name: `Remove project ${name}`, exact: true })
    .click();
  await page
    .locator('.workspace-library__confirm')
    .getByRole('button', { name: 'Confirm', exact: true })
    .click();
  await expect(row).toHaveCount(1);
});

test('molecule library search, rename, archive, recoverable trash and restore operate on a disposable upload', async ({
  page,
  request,
}) => {
  const demo = await seed(request);
  const name = `Disposable library ${Date.now()}`;
  const pdb = await (await request.get(`/api/datasets/${demo.id}/topology`)).body();
  const uploaded = await request.post('/api/datasets/upload', {
    multipart: {
      topology: { name: `${name}.pdb`, mimeType: 'chemical/x-pdb', buffer: pdb },
      stride: '1',
    },
  });
  expect(uploaded.ok()).toBeTruthy();
  const fixture = await uploaded.json();
  await ready(page);
  await page.getByTitle('Switch dataset', { exact: true }).click();
  await expect(page.locator('.workspace-library__busy')).toHaveCount(0);
  await page.getByRole('textbox', { name: 'Search molecule library', exact: true }).fill(name);
  let row = page
    .locator('.workspace-library__row')
    .filter({ has: page.getByText(name, { exact: true }) });
  await expect(row).toHaveCount(1);
  await row.getByRole('button', { name: `Rename ${name}`, exact: true }).click();
  const renamed = `${name} renamed`;
  await page.getByRole('textbox', { name: 'Rename molecule', exact: true }).fill(renamed);
  await page
    .locator('.workspace-library__confirm')
    .getByRole('button', { name: 'Confirm', exact: true })
    .click();
  row = page
    .locator('.workspace-library__row')
    .filter({ has: page.getByText(renamed, { exact: true }) });
  await expect(row).toHaveCount(1);
  await row.getByRole('button', { name: `Archive ${renamed}`, exact: true }).click();
  await expect(row).toHaveCount(0);
  await page
    .getByRole('combobox', { name: 'Library filter', exact: true })
    .selectOption('archived');
  await expect(row).toHaveCount(1);
  await row.getByRole('button', { name: `Unarchive ${renamed}`, exact: true }).click();
  await page.getByRole('combobox', { name: 'Library filter', exact: true }).selectOption('active');
  await expect(row).toHaveCount(1);
  await row.getByRole('button', { name: `Move ${renamed} to trash`, exact: true }).click();
  const confirm = page
    .locator('.workspace-library__confirm')
    .getByRole('button', { name: 'Confirm', exact: true });
  await expect(confirm).toBeDisabled();
  await page.getByRole('textbox', { name: 'Confirm dataset name', exact: true }).fill(renamed);
  await confirm.click();
  await expect(row).toHaveCount(0);
  expect((await request.get(`/api/datasets/${fixture.id}`)).status()).toBe(404);
  await page.getByRole('combobox', { name: 'Library filter', exact: true }).selectOption('trash');
  await expect(row).toHaveCount(1);
  await row.getByRole('button', { name: 'Restore', exact: true }).click();
  await expect(row).toHaveCount(0);
  const restored = await (await request.get(`/api/datasets/${fixture.id}`)).json();
  expect(restored.name).toBe(renamed);
  expect(restored.n_atoms).toBe(fixture.n_atoms);
  await expect(page.locator('.workspace-library__disk')).toContainText('free');
});

test('bad backup and protected source failures remain visible without losing the current workspace', async ({
  page,
  request,
}) => {
  const demo = await seed(request);
  await ready(page);
  await projects(page);
  await page.getByLabel('Import project backup file', { exact: true }).setInputFiles({
    name: 'invalid.zip',
    mimeType: 'application/zip',
    buffer: Buffer.from('not a zip'),
  });
  await expect(page.locator('.workspace-library__error')).toContainText('malformed or unreadable');
  await page.getByRole('button', { name: 'Molecule library', exact: true }).click();
  await page.getByRole('textbox', { name: 'Search molecule library', exact: true }).fill(demo.name);
  const row = page.locator('.workspace-library__row').filter({ hasText: demo.name });
  await row.getByRole('button', { name: `Move ${demo.name} to trash`, exact: true }).click();
  await page.getByRole('textbox', { name: 'Confirm dataset name', exact: true }).fill(demo.name);
  await page
    .locator('.workspace-library__confirm')
    .getByRole('button', { name: 'Confirm', exact: true })
    .click();
  await expect(page.locator('.workspace-library__error')).toContainText('bundled example');
  expect((await request.get('/api/datasets/demo')).ok()).toBeTruthy();
  await page.getByRole('button', { name: 'Close workspace library', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
});

test('a delayed project open cannot replace a newer molecule or close its simulation pane', async ({
  page,
  request,
}) => {
  const demo = await seed(request);
  const pdb = await (await request.get(`/api/datasets/${demo.id}/topology`)).body();
  const uploaded = await request.post('/api/datasets/upload', {
    multipart: {
      topology: {
        name: `Later molecule ${Date.now()}.pdb`,
        mimeType: 'chemical/x-pdb',
        buffer: pdb,
      },
      stride: '1',
    },
  });
  expect(uploaded.ok()).toBeTruthy();
  const later = await uploaded.json();
  const name = `Delayed project ${Date.now()}`;
  await request.post('/api/projects', {
    data: { name, state: { version: 1, dataset_id: demo.id, frame: 9 } },
  });
  await ready(page);
  let release!: () => void;
  let observed!: () => void;
  const gate = new Promise<void>((resolve) => {
    release = resolve;
  });
  const received = new Promise<void>((resolve) => {
    observed = resolve;
  });
  await page.route('**/api/datasets/demo', async (route) => {
    const response = await route.fetch();
    observed();
    await gate;
    await route.fulfill({ response });
  });
  try {
    await projects(page);
    await page
      .locator('.workspace-library__row')
      .filter({ has: page.getByText(name, { exact: true }) })
      .getByRole('button', { name: 'Open', exact: true })
      .click();
    await received;
    await page.getByRole('button', { name: 'Close workspace library', exact: true }).click();
    await page.getByTitle('Switch dataset', { exact: true }).click();
    await page
      .locator('.workspace-library__row')
      .filter({ hasText: later.name })
      .getByRole('button', { name: 'Open', exact: true })
      .click();
    await expect(page.locator('.structure-card h2')).toHaveText(later.name);
    await page.getByRole('button', { name: 'New simulation', exact: true }).click();
    release();
    await expect(page.getByRole('dialog', { name: 'Set molecules in motion.' })).toBeVisible();
    await page.waitForTimeout(500);
    await expect(page.locator('.structure-card h2')).toHaveText(later.name);
    await expect(page.getByRole('dialog', { name: 'Set molecules in motion.' })).toBeVisible();
  } finally {
    release();
  }
});
