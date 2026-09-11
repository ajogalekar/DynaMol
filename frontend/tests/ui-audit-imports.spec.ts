import { test } from './testWorkspace';
import { expect, type Page, type TestInfo } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

const path = (name: string) =>
  fileURLToPath(new URL(`../../docs/audit/ui-audit-fixtures/${name}`, import.meta.url));
async function ready(page: Page) {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
}
async function importReady(page: Page) {
  await ready(page);
  await page.getByRole('button', { name: 'Open files', exact: true }).click();
  return page.getByRole('dialog', { name: 'Bring your molecules.' });
}

for (const format of [
  'pdb',
  'cif',
  'gro',
  'h5',
  'xtc',
  'dcd',
  'trr',
  'nc',
  'xyz',
  'mdcrd',
  'lammpstrj',
]) {
  test(`audit import format: ${format} loads genuine coordinates, stride and explicit timing`, async ({
    page,
    request,
  }, info) => {
    const dialog = await importReady(page);
    const topologyOnly = ['pdb', 'cif', 'gro', 'h5'].includes(format);
    if (format === 'lammpstrj')
      await expect(dialog.locator('input[type=file]').nth(1)).toHaveAttribute(
        'accept',
        /\.lammpstrj/,
      );
    const name = `e2e-ui-audit-${format}-${Date.now()}`;
    await dialog
      .locator('input[type=file]')
      .first()
      .setInputFiles({
        name: `${name}.${topologyOnly ? format : 'pdb'}`,
        mimeType: 'application/octet-stream',
        buffer: await readFile(path(`short.${topologyOnly ? format : 'pdb'}`)),
      });
    if (!topologyOnly)
      await dialog
        .locator('input[type=file]')
        .nth(1)
        .setInputFiles(path(`short.${format}`));
    await dialog.getByRole('spinbutton', { name: /^Read every/ }).fill('2');
    await dialog.getByRole('spinbutton', { name: /^Time between original frames/ }).fill('0.25');
    const response = page.waitForResponse(
      (r) => r.url().endsWith('/api/datasets/upload') && r.request().method() === 'POST',
    );
    await dialog.getByRole('button', { name: 'Open in DynaMol', exact: true }).click();
    const network = await response;
    expect(network.ok(), await network.text()).toBe(true);
    const loaded = await network.json();
    expect(loaded.n_atoms).toBe(113);
    expect(loaded.n_frames).toBe(format === 'cif' ? 1 : 3);
    expect(loaded.times_ps).toEqual(format === 'cif' ? [0] : [0, 0.5, 1]);
    const reference = JSON.parse(await readFile(path('manifest.json'), 'utf8'));
    const measured = await request.post(`/api/datasets/${loaded.id}/measurements`, {
      data: { kind: 'distance', atoms: [0, 1] },
    });
    expect(measured.ok(), await measured.text()).toBe(true);
    const distances = (await measured.json()).values;
    distances.forEach((value: number, i: number) =>
      expect(Math.abs(value - reference.first_pair_distance_angstrom[i * 2])).toBeLessThan(0.03),
    );
    await expect(dialog).toHaveCount(0);
    await expect(page.locator('.structure-card h2')).toHaveText(name);
    await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
    await expect(page.getByRole('alert')).toHaveCount(0);
    await info.attach(`${format}-result.json`, {
      body: JSON.stringify(
        {
          id: loaded.id,
          name,
          n_atoms: loaded.n_atoms,
          n_frames: loaded.n_frames,
          times_ps: loaded.times_ps,
          warnings: loaded.warnings,
        },
        null,
        2,
      ),
      contentType: 'application/json',
    });
  });
}

async function scene(page: Page, info: TestInfo, name: string) {
  await page.mouse.move(2, 2);
  await page.waitForTimeout(450);
  const image = await page.locator('.molecular-viewer__canvas canvas').screenshot();
  await info.attach(name, { body: image, contentType: 'image/png' });
  return image;
}

test('audit: protein, ligand, water and ion visibility each changes a mixed molecular scene', async ({
  page,
}, info) => {
  const dialog = await importReady(page);
  await dialog
    .locator('input[type=file]')
    .first()
    .setInputFiles({
      name: `e2e-ui-audit-mixed-${Date.now()}.pdb`,
      mimeType: 'chemical/x-pdb',
      buffer: await readFile(path('mixed-groups.pdb')),
    });
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/datasets/upload') && r.request().method() === 'POST',
  );
  await dialog.getByRole('button', { name: 'Open in DynaMol', exact: true }).click();
  const loaded = await (await response).json();
  await expect(dialog).toHaveCount(0);
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  const counts = loaded.atoms.reduce((out: Record<string, number>, atom: { category: string }) => {
    out[atom.category] = (out[atom.category] || 0) + 1;
    return out;
  }, {});
  expect(counts.water).toBe(3);
  expect(counts.ligand ?? counts.ligands).toBe(3);
  expect(counts.ion ?? counts.ions).toBe(1);
  await page.getByRole('button', { name: 'Ball & stick', exact: true }).click();
  await page.getByRole('button', { name: 'Show water', exact: true }).click();
  await page.getByRole('button', { name: 'All', exact: true }).click();
  await page.getByRole('button', { name: 'Fit molecule', exact: true }).click();
  await page.waitForTimeout(600);
  for (const label of ['protein & nucleic acids', 'ligands', 'water', 'ions']) {
    const before = await scene(page, info, `all-before-${label}`);
    await page.getByRole('button', { name: `Hide ${label}`, exact: true }).click();
    await expect(page.getByRole('button', { name: `Show ${label}`, exact: true })).toHaveAttribute(
      'aria-pressed',
      'false',
    );
    const hidden = await scene(page, info, `hidden-${label}`);
    expect(hidden.equals(before), `${label} must affect visible molecular geometry.`).toBe(false);
    await page.getByRole('button', { name: `Show ${label}`, exact: true }).click();
    await expect(page.getByRole('button', { name: `Hide ${label}`, exact: true })).toHaveAttribute(
      'aria-pressed',
      'true',
    );
  }
  await expect(page.getByRole('alert')).toHaveCount(0);
});

test('audit: trajectory and topology mismatch is reported without replacing the current structure', async ({
  page,
  request,
}) => {
  const dialog = await importReady(page);
  const initial = await page.locator('.structure-card h2').innerText();
  const topology = await (await request.get('/api/datasets/demo/topology')).body();
  await dialog.locator('input[type=file]').first().setInputFiles({
    name: 'e2e-ui-audit-mismatch.pdb',
    mimeType: 'chemical/x-pdb',
    buffer: topology,
  });
  await dialog.locator('input[type=file]').nth(1).setInputFiles(path('short.xtc'));
  await dialog.getByRole('button', { name: 'Open in DynaMol', exact: true }).click();
  await expect(dialog.getByRole('alert')).toBeVisible();
  await expect(page.locator('.structure-card h2')).toHaveText(initial);
  await expect(dialog.getByRole('button', { name: 'Open in DynaMol', exact: true })).toBeEnabled();
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
});

test('audit: hidden atom groups and hydrogen modes also hide selection highlights', async ({
  page,
}, info) => {
  const dialog = await importReady(page);
  await dialog
    .locator('input[type=file]')
    .first()
    .setInputFiles({
      name: `e2e-ui-audit-hidden-selection-${Date.now()}.pdb`,
      mimeType: 'chemical/x-pdb',
      buffer: await readFile(path('mixed-groups.pdb')),
    });
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/datasets/upload') && r.request().method() === 'POST',
  );
  await dialog.getByRole('button', { name: 'Open in DynaMol', exact: true }).click();
  const dataset = await (await response).json();
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'Ball & stick', exact: true }).click();
  await page.getByRole('button', { name: 'Show water', exact: true }).click();
  type SelectedAtom = {
    index: number;
    category: string;
    element: string;
    name: string;
    residue: string;
    resid: number;
    nonpolar_hydrogen: boolean;
  };
  async function choose(atom: SelectedAtom) {
    await page
      .getByRole('textbox', { name: 'Find an atom', exact: true })
      .fill(`${atom.residue}${atom.resid} · ${atom.name}`);
    await page
      .locator('.atom-search-results')
      .getByRole('button', { name: new RegExp(`#${atom.index + 1} · ${atom.element}$`) })
      .click();
    await page.waitForTimeout(1200);
  }
  async function expectNoHiddenHighlight(label: string) {
    const hidden = await scene(page, info, `${label}-selected-hidden`);
    await page.getByRole('button', { name: 'Clear selection', exact: true }).click();
    const cleared = await scene(page, info, `${label}-selection-cleared`);
    expect(
      hidden.equals(cleared),
      `${label}: clearing already-hidden atom selection must not remove visible geometry.`,
    ).toBe(true);
  }
  for (const [category, label] of [
    ['protein', 'protein & nucleic acids'],
    ['ligands', 'ligands'],
    ['water', 'water'],
    ['ions', 'ions'],
  ]) {
    await choose(
      dataset.atoms.find(
        (atom: SelectedAtom) => atom.category === category && atom.element !== 'H',
      ),
    );
    await page.getByRole('button', { name: `Hide ${label}`, exact: true }).click();
    await expectNoHiddenHighlight(category);
    await page.getByRole('button', { name: `Show ${label}`, exact: true }).click();
  }
  for (const mode of ['Polar only', 'Hidden']) {
    await page.getByRole('button', { name: 'All', exact: true }).click();
    await choose(dataset.atoms.find((atom: SelectedAtom) => atom.nonpolar_hydrogen));
    await page.getByRole('button', { name: mode, exact: true }).click();
    await expectNoHiddenHighlight(mode);
  }
});

test('audit: polymer visibility includes a schematic nucleic-acid backbone', async ({
  page,
}, info) => {
  const dialog = await importReady(page);
  await dialog
    .locator('input[type=file]')
    .first()
    .setInputFiles({
      name: `e2e-ui-audit-nucleic-${Date.now()}.pdb`,
      mimeType: 'chemical/x-pdb',
      buffer: await readFile(path('schematic-nucleic.pdb')),
    });
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/datasets/upload') && r.request().method() === 'POST',
  );
  await dialog.getByRole('button', { name: 'Open in DynaMol', exact: true }).click();
  const dataset = await (await response).json();
  expect(dataset.n_atoms).toBe(24);
  expect(dataset.atoms.every((atom: { category: string }) => atom.category === 'nucleic')).toBe(
    true,
  );
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'Ball & stick', exact: true }).click();
  await page.waitForTimeout(600);
  const shown = await scene(page, info, 'schematic-nucleic-visible');
  await page.getByRole('button', { name: 'Hide protein & nucleic acids', exact: true }).click();
  const hidden = await scene(page, info, 'schematic-nucleic-hidden');
  expect(hidden.equals(shown)).toBe(false);
  await page.getByRole('button', { name: 'Show protein & nucleic acids', exact: true }).click();
  const restored = await scene(page, info, 'schematic-nucleic-restored');
  expect(restored.equals(shown)).toBe(true);
});
