import { test } from './testWorkspace';
import { expect } from '@playwright/test';
import { fileURLToPath } from 'node:url';
import { readFile } from 'node:fs/promises';
import { execFileSync } from 'node:child_process';

const fixture = (name: string) =>
  fileURLToPath(new URL(`../../docs/audit/ui-audit-fixtures/${name}`, import.meta.url));

for (const format of ['cif', 'mol', 'sdf', 'smi', 'smiles']) {
  test(`audit Studio upload format: ${format} displays the imported structure`, async ({
    page,
  }, info) => {
    await page.goto('/');
    await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
    await page.getByRole('button', { name: 'New simulation', exact: true }).click();
    const dialog = page.getByRole('dialog', { name: 'Set molecules in motion.' });
    const name = `e2e-ui-audit-source-${format}-${Date.now()}`;
    const response = page.waitForResponse(
      (r) => r.url().endsWith('/api/structures/upload') && r.request().method() === 'POST',
    );
    await dialog.getByLabel('Upload simulation structure', { exact: true }).setInputFiles({
      name: `${name}.${format}`,
      mimeType: 'application/octet-stream',
      buffer: await readFile(fixture(`${format === 'cif' ? 'short' : 'ethanol'}.${format}`)),
    });
    const network = await response;
    expect(network.ok(), await network.text()).toBe(true);
    const result = await network.json();
    expect(result.n_atoms).toBe(format === 'cif' ? 113 : 9);
    await expect(page.locator('.structure-card h2')).toHaveText(name);
    await expect(dialog.locator('.source-summary strong')).toHaveText(name);
    await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
    await expect(page.getByRole('button', { name: 'Play trajectory', exact: true })).toBeDisabled();
    await info.attach(`${format}-source-result.json`, {
      body: JSON.stringify(
        { id: result.id, name: result.name, atoms: result.n_atoms, source: result.source },
        null,
        2,
      ),
      contentType: 'application/json',
    });
  });
}

test('audit: failed inspection retries against the real service without reloading the molecule', async ({
  page,
}) => {
  await page.route('**/api/datasets/demo/inspection?*', (route) => route.abort('failed'));
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  await expect(dialog.getByRole('button', { name: 'Retry inspection', exact: true })).toBeVisible();
  const original = await page.locator('.structure-card h2').innerText();
  await page.unroute('**/api/datasets/demo/inspection?*');
  await dialog.getByRole('button', { name: 'Retry inspection', exact: true }).click();
  await expect(dialog.locator('.inspection-summary')).toBeVisible();
  await expect(dialog.getByRole('button', { name: 'Retry inspection', exact: true })).toHaveCount(
    0,
  );
  await expect(page.locator('.structure-card h2')).toHaveText(original);
});

test('audit: ligand details, explicit-state SMILES and pH trigger real inspection updates', async ({
  page,
}, info) => {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  const source = page.waitForResponse(
    (r) => r.url().endsWith('/api/structures/upload') && r.request().method() === 'POST',
  );
  await dialog.getByLabel('Upload simulation structure', { exact: true }).setInputFiles({
    name: `e2e-ui-audit-override-${Date.now()}.pdb`,
    mimeType: 'chemical/x-pdb',
    buffer: await readFile(fixture('mixed-groups.pdb')),
  });
  const dataset = await (await source).json();
  const card = dialog.locator('.ligand-preparation-card').first();
  await expect(card).toBeVisible();
  await card.locator('summary').click();
  await expect(card.getByRole('textbox', { name: /^Ligand state override/ })).toBeVisible();
  const invalidResponse = page.waitForResponse(
    (r) =>
      r.url().includes(`/api/datasets/${dataset.id}/inspection`) &&
      r.request().method() === 'POST' &&
      Object.values(r.request().postDataJSON()?.ligand_overrides ?? {}).includes('not-a-smiles'),
  );
  await card.getByRole('textbox', { name: /^Ligand state override/ }).fill('not-a-smiles');
  const invalid = await invalidResponse;
  expect(invalid.ok(), await invalid.text()).toBe(true);
  expect((await invalid.json()).ligands[0].error).toBeTruthy();
  await expect(card).toContainText('Needs chemistry');
  const checked = page.waitForResponse(
    (r) =>
      r.url().includes(`/api/datasets/${dataset.id}/inspection`) &&
      r.request().method() === 'POST' &&
      Object.values(r.request().postDataJSON()?.ligand_overrides ?? {}).includes('CCO'),
  );
  await card.getByRole('textbox', { name: /^Ligand state override/ }).fill('CCO');
  const valid = await checked;
  expect(valid.ok(), await valid.text()).toBe(true);
  const result = await valid.json();
  expect(result.ligands[0].error ?? null).toBeNull();
  expect(result.ligands[0].formal_charge).toBe(0);
  await expect(card.locator('summary')).toContainText('0 charge');
  await expect(card.getByLabel(/^Selected SMILES/)).toBeVisible();
  const phCheck = page.waitForResponse(
    (r) =>
      r.url().includes(`/api/datasets/${dataset.id}/inspection`) &&
      r.request().method() === 'POST' &&
      r.request().postDataJSON()?.ph === 8,
  );
  await dialog.getByLabel('Preparation pH', { exact: true }).fill('8');
  expect((await phCheck).ok()).toBe(true);
  await dialog.getByRole('button', { name: 'Preparation options', exact: true }).click();
  await dialog
    .getByRole('checkbox', { name: 'Remove ligands and other non-protein residues', exact: true })
    .check();
  await expect(card.getByRole('textbox', { name: /^Ligand state override/ })).toBeDisabled();
  await expect(dialog.getByRole('button', { name: 'Prep protein', exact: true })).toBeVisible();
  await dialog
    .getByRole('checkbox', { name: 'Remove ligands and other non-protein residues', exact: true })
    .uncheck();
  await expect(dialog.getByRole('button', { name: 'Prep complex', exact: true })).toBeVisible();
  await card.locator('summary').click();
  await expect(card.getByRole('textbox', { name: /^Ligand state override/ })).not.toBeVisible();
  await info.attach('ligand-inspection.json', {
    body: JSON.stringify(result, null, 2),
    contentType: 'application/json',
  });
});

test('audit: validated complex result exposes prepared coordinates and parameter downloads', async ({
  page,
  request,
}, info) => {
  const datasetId = process.env.DYNAMOL_AUDIT_COMPLEX_DATASET;
  test.skip(
    !datasetId,
    'Set DYNAMOL_AUDIT_COMPLEX_DATASET to the completed native complex validation result.',
  );
  const response = await request.get(`/api/datasets/${datasetId}`);
  expect(response.ok(), await response.text()).toBe(true);
  const dataset = await response.json();
  expect(dataset.preparation?.ligand_parameters).toBeTruthy();
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByTitle('Switch dataset', { exact: true }).click();
  await page
    .locator('.workspace-library__row')
    .filter({ hasText: dataset.name })
    .first()
    .getByRole('button', { name: 'Open', exact: true })
    .click();
  await expect(page.locator('.structure-card h2')).toHaveText(dataset.name);
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  await expect(
    dialog.getByRole('combobox', { name: 'Solvent environment', exact: true }),
  ).toHaveValue('explicit');
  await expect(
    dialog
      .getByRole('combobox', { name: 'Solvent environment', exact: true })
      .locator('option[value="implicit"]'),
  ).toHaveJSProperty('disabled', true);
  const pdbEvent = page.waitForEvent('download');
  await dialog.getByRole('link', { name: 'Prepared PDB', exact: true }).click();
  const pdb = await readFile((await (await pdbEvent).path())!, 'utf8');
  expect(pdb.split('\n').filter((line) => /^(ATOM  |HETATM)/.test(line))).toHaveLength(
    dataset.n_atoms,
  );
  const parametersEvent = page.waitForEvent('download');
  await dialog.getByRole('link', { name: 'Parameters & preparation files', exact: true }).click();
  const download = await parametersEvent;
  const zipPath = info.outputPath('complex-parameters.zip');
  await download.saveAs(zipPath);
  const python = fileURLToPath(new URL('../../.venv/bin/python', import.meta.url));
  const entries = JSON.parse(
    execFileSync(
      python,
      [
        '-c',
        'import json,sys,zipfile;print(json.dumps(zipfile.ZipFile(sys.argv[1]).namelist()))',
        zipPath,
      ],
      { encoding: 'utf8' },
    ),
  );
  expect(entries.some((entry: string) => entry.endsWith('/ligand.xml'))).toBe(true);
  expect(entries.some((entry: string) => entry.endsWith('/charged.mol2'))).toBe(true);
  await info.attach('complex-download-manifest.json', {
    body: JSON.stringify(
      { dataset_id: dataset.id, name: dataset.name, n_atoms: dataset.n_atoms, entries },
      null,
      2,
    ),
    contentType: 'application/json',
  });
  await dialog.getByRole('button', { name: /GROMACS/ }).click();
  await expect(dialog.locator('.engine-compatibility')).toContainText(
    'DynaMol cannot yet transfer that state to GROMACS.',
  );
  await expect(
    dialog.getByRole('button', { name: 'Start simulation', exact: true }),
  ).toBeDisabled();
  await dialog.getByRole('button', { name: 'OpenMM', exact: true }).click();
});
