import { test, expect, type Page, type TestInfo } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

const fixtures = fileURLToPath(new URL('../../docs/audit/preparation-fixtures/', import.meta.url));

async function openReady(page: Page) {
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await page.goto('/');
  // Opening Studio before initial dataset initialization finishes can close it.
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  return errors;
}

async function openStudio(page: Page) {
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const studio = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  await expect(studio).toBeVisible();
  return studio;
}

async function sourceLoaded(page: Page, name: string) {
  await expect(page.locator('.structure-card h2')).toHaveText(name);
  await expect(page.locator('.source-summary strong')).toHaveText(name);
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.getByRole('dialog', { name: 'Set molecules in motion.' })).toBeVisible();
  await expect(page.getByRole('button', { name: 'Play trajectory', exact: true })).toBeDisabled();
  await expect(page.getByRole('alert')).toHaveCount(0);
}

async function renderedCanvas(page: Page, testInfo: TestInfo, name: string) {
  const canvas = page.locator('.molecular-viewer__canvas canvas');
  await expect(canvas).toHaveCount(1);
  await expect(canvas).toBeVisible();
  const bounds = await canvas.boundingBox();
  expect(bounds?.width).toBeGreaterThan(200);
  expect(bounds?.height).toBeGreaterThan(200);
  const image = await canvas.screenshot();
  await testInfo.attach(name, { body: image, contentType: 'image/png' });
  const count = await page.evaluate(async (base64) => {
    const bitmap = await createImageBitmap(
      await (await fetch(`data:image/png;base64,${base64}`)).blob(),
    );
    const sample = document.createElement('canvas');
    sample.width = bitmap.width;
    sample.height = bitmap.height;
    const context = sample.getContext('2d')!;
    context.drawImage(bitmap, 0, 0);
    const pixels = context.getImageData(0, 0, sample.width, sample.height).data;
    let colored = 0;
    for (let i = 0; i < pixels.length; i += 4) {
      const maximum = Math.max(pixels[i], pixels[i + 1], pixels[i + 2]);
      const minimum = Math.min(pixels[i], pixels[i + 1], pixels[i + 2]);
      if (maximum > 75 && maximum - minimum > 35) colored++;
    }
    return colored;
  }, image.toString('base64'));
  expect(count, 'Imported atoms must produce actual colored molecular geometry.').toBeGreaterThan(
    500,
  );
  return image;
}

test('sample distance starts hidden; scene and individual controls hide labels while retaining the plot', async ({
  page,
}, testInfo) => {
  const errors = await openReady(page);
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  const eye = page.locator('.measurement-tab .measurement-eye');
  await expect(eye).toHaveAttribute('aria-pressed', 'false');
  await expect(
    page.getByRole('button', { name: 'Show measurements', exact: true }),
  ).toHaveAttribute('aria-pressed', 'false');
  const readout = await page.locator('.plot-readout strong').innerText();
  await page.mouse.move(20, 20);
  const hidden = await renderedCanvas(page, testInfo, 'distance-hidden-default');

  await page.getByRole('button', { name: 'Show measurements', exact: true }).click();
  await expect(eye).toHaveAttribute('aria-pressed', 'true');
  await expect(
    page.getByRole('button', { name: 'Hide measurements', exact: true }),
  ).toHaveAttribute('aria-pressed', 'true');
  await page.mouse.move(20, 20);
  const shown = await renderedCanvas(page, testInfo, 'distance-visible');
  expect(shown.equals(hidden), 'Showing the measurement must change the real WebGL canvas.').toBe(
    false,
  );

  await eye.click();
  await expect(eye).toHaveAttribute('aria-pressed', 'false');
  await expect(page.getByRole('button', { name: 'Show measurements', exact: true })).toBeVisible();
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  await expect(page.locator('.plot-readout strong')).toHaveText(readout);
  await expect(page.locator('.chart-wrap svg')).toBeVisible();
  await eye.click();
  await expect(eye).toHaveAttribute('aria-pressed', 'true');
  await page.getByRole('button', { name: 'Hide measurements', exact: true }).click();
  await expect(eye).toHaveAttribute('aria-pressed', 'false');
  await expect(page.locator('.plot-readout strong')).toHaveText(readout);
  expect(errors).toEqual([]);
});

test('Studio PDB upload stays beside the live view and exposes known missing loops as an explicit option', async ({
  page,
  request,
}, testInfo) => {
  const errors = await openReady(page);
  const studio = await openStudio(page);
  const name = `e2e-known-gap-${Date.now()}`;
  const buffer = await readFile(`${fixtures}/six_residues_known_gap.pdb`);
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/structures/upload') && r.request().method() === 'POST',
  );
  await studio
    .getByLabel('Upload simulation structure', { exact: true })
    .setInputFiles({ name: `${name}.pdb`, mimeType: 'chemical/x-pdb', buffer });
  const imported = await (await response).json();
  await sourceLoaded(page, name);
  await expect(studio.locator('.missing-residues')).toContainText('Missing sequence regions');
  await expect(studio.locator('.missing-residues')).toContainText('ILE');
  const build = studio.getByRole('checkbox', {
    name: 'Build supported missing loops / residues',
    exact: true,
  });
  await expect(build).toBeEnabled();
  await expect(build).not.toBeChecked();
  await build.check();
  await expect(build).toBeChecked();
  await build.uncheck();
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  await expect(studio.getByLabel('Preparation pH', { exact: true })).toHaveValue('7');
  const inspection = await (await request.get(`/api/datasets/${imported.id}/inspection`)).json();
  expect(inspection.has_sequence).toBe(true);
  expect(
    inspection.missing_residues.some(
      (r: { buildable: boolean; residues: string[] }) => r.buildable && r.residues.includes('ILE'),
    ),
  ).toBe(true);
  await renderedCanvas(page, testInfo, 'studio-known-loop-pdb');
  expect(errors).toEqual([]);
});

test('Studio MOL2 upload preserves explicit atoms, bonds, and source charges in the viewed molecule', async ({
  page,
}, testInfo) => {
  const errors = await openReady(page);
  const studio = await openStudio(page);
  const name = `e2e-ethanol-mol2-${Date.now()}`;
  const mol2 = `@<TRIPOS>MOLECULE\nethanol\n3 2 1 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 C1 0.000 0.000 0.000 C.3 1 LIG1 -0.1\n2 C2 1.500 0.000 0.000 C.3 1 LIG1 0.3\n3 O1 2.000 1.200 0.000 O.3 1 LIG1 -0.2\n@<TRIPOS>BOND\n1 1 2 1\n2 2 3 1\n@<TRIPOS>SUBSTRUCTURE\n1 LIG1 1 GROUP 0 A\n`;
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/structures/upload') && r.request().method() === 'POST',
  );
  await studio
    .getByLabel('Upload simulation structure', { exact: true })
    .setInputFiles({
      name: `${name}.mol2`,
      mimeType: 'chemical/x-mol2',
      buffer: Buffer.from(mol2),
    });
  const imported = await (await response).json();
  await sourceLoaded(page, name);
  expect(imported.n_atoms).toBe(3);
  expect(imported.bonds).toEqual([
    [0, 1],
    [1, 2],
  ]);
  expect(imported.atoms.map((a: { partial_charge: number }) => a.partial_charge)).toEqual([
    -0.1, 0.3, -0.2,
  ]);
  await expect(studio.locator('.inspection-summary')).toHaveText('Small-molecule structure');
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeDisabled();
  await renderedCanvas(page, testInfo, 'studio-ethanol-mol2');
  expect(errors).toEqual([]);
});

test('SMILES errors preserve the current scene, then a named 3D conformer loads without leaving Studio', async ({
  page,
}, testInfo) => {
  const errors = await openReady(page);
  const original = await page.locator('.structure-card h2').innerText();
  const studio = await openStudio(page);
  await studio.getByRole('tab', { name: 'SMILES', exact: true }).click();
  await studio
    .getByRole('textbox', { name: 'SMILES string', exact: true })
    .fill('not-a-valid-smiles');
  await studio.getByRole('button', { name: 'Build 3D structure & view', exact: true }).click();
  await expect(studio.getByRole('alert')).toContainText('SMILES');
  await expect(page.locator('.structure-card h2')).toHaveText(original);
  const name = `e2e-lactate-smiles-${Date.now()}`;
  await studio.getByRole('textbox', { name: 'SMILES molecule name', exact: true }).fill(name);
  await studio
    .getByRole('textbox', { name: 'SMILES string', exact: true })
    .fill('C[C@H](O)C(=O)[O-]');
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/structures/smiles') && r.request().method() === 'POST',
  );
  await studio.getByRole('button', { name: 'Build 3D structure & view', exact: true }).click();
  const imported = await (await response).json();
  await sourceLoaded(page, name);
  expect(imported.chemistry.total_formal_charge).toBe(-1);
  expect(imported.chemistry.canonical_isomeric_smiles).toContain('@');
  expect(imported.chemistry.conformer.method).toBe('RDKit ETKDGv3');
  expect(imported.warnings.join(' ')).toContain('Computed 3D conformer');
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeDisabled();
  await renderedCanvas(page, testInfo, 'studio-lactate-smiles');
  expect(errors).toEqual([]);
});

test('Studio fetch retrieves actual RCSB and PubChem structures into the visible scene', async ({
  page,
}, testInfo) => {
  test.setTimeout(120_000);
  const errors = await openReady(page);
  const studio = await openStudio(page);
  const fetched: { id: string; name: string; source: string }[] = [];
  await studio.getByRole('tab', { name: 'Fetch', exact: true }).click();
  await studio.getByRole('textbox', { name: 'PDB accession', exact: true }).fill('1UBQ');
  let response = page.waitForResponse(
    (r) => r.url().endsWith('/api/structures/fetch') && r.request().method() === 'POST',
  );
  await studio.getByRole('button', { name: 'Fetch', exact: true }).click();
  let network = await response;
  expect(network.ok(), await network.text()).toBe(true);
  const protein = await network.json();
  fetched.push({ id: protein.id, name: protein.name, source: protein.source });
  await sourceLoaded(page, protein.name);
  expect(protein.source).toBe('rcsb_pdb');
  expect(protein.n_atoms).toBeGreaterThan(600);
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  await renderedCanvas(page, testInfo, 'studio-rcsb-1ubq');
  await studio
    .getByRole('combobox', { name: 'Structure database', exact: true })
    .selectOption('pubchem');
  await studio.getByRole('textbox', { name: 'PubChem name or CID', exact: true }).fill('caffeine');
  response = page.waitForResponse(
    (r) => r.url().endsWith('/api/structures/fetch') && r.request().method() === 'POST',
  );
  await studio.getByRole('button', { name: 'Fetch', exact: true }).click();
  network = await response;
  expect(network.ok(), await network.text()).toBe(true);
  const molecule = await network.json();
  fetched.push({ id: molecule.id, name: molecule.name, source: molecule.source });
  await sourceLoaded(page, molecule.name);
  expect(molecule.source).toBe('pubchem');
  expect(molecule.n_atoms).toBe(24);
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeDisabled();
  await renderedCanvas(page, testInfo, 'studio-pubchem-caffeine');
  await testInfo.attach('fetched-datasets-for-cleanup', {
    body: JSON.stringify(fetched, null, 2),
    contentType: 'application/json',
  });
  console.log(
    'Preparation E2E fetched dataset IDs:',
    fetched.map((dataset) => dataset.id).join(', '),
  );
  expect(errors).toEqual([]);
});
