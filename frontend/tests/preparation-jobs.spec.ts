import { test } from './testWorkspace';
import { expect, type Page } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

test.use({ actionTimeout: 15_000 });

async function sceneReady(page: Page, name: string) {
  await expect(page.locator('.structure-card h2')).toHaveText(name);
  await expect(page.locator('.source-summary strong')).toHaveText(name);
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.getByRole('dialog', { name: 'Set molecules in motion.' })).toBeVisible();
}

test('Prep protein completes real repair in the background; explicit water appears and implicit restores the dry parent', async ({
  page,
  request,
}, testInfo) => {
  test.setTimeout(150_000);
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const studio = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  await expect(studio).toBeVisible();

  const name = `e2e-protein-prep-${Date.now()}`;
  const fixture = fileURLToPath(
    new URL('../../docs/audit/preparation-fixtures/six_residues_missing_atom.pdb', import.meta.url),
  );
  const importedResponse = page.waitForResponse(
    (r) => r.url().endsWith('/api/structures/upload') && r.request().method() === 'POST',
  );
  await studio
    .getByLabel('Upload simulation structure', { exact: true })
    .setInputFiles({
      name: `${name}.pdb`,
      mimeType: 'chemical/x-pdb',
      buffer: await readFile(fixture),
    });
  const source = await (await importedResponse).json();
  await sceneReady(page, name);
  const tracking = studio.getByRole('region', { name: 'Track during simulation', exact: true });
  await tracking.getByRole('button', { name: /^Track during simulation/ }).click();
  await tracking.getByRole('button', { name: 'Add measurement', exact: true }).click();
  for (const atom of source.atoms.filter((a: {name: string}) => a.name === 'CA').slice(0, 2)) {
    await tracking.getByRole('textbox', { name: 'Find an atom', exact: true }).fill(`${atom.residue}${atom.resid} · ${atom.name}`);
    await tracking.locator('.atom-search-results').getByRole('button', {name: new RegExp(`#${atom.index + 1} · ${atom.element}$`)}).click();
  }
  await expect(tracking.getByRole('checkbox')).toBeChecked();
  const assertTrackingPreserved = async (expected: {id: string; atoms: {index: number;name: string}[]}) => {
    await expect.poll(async () => (await (await request.get('/api/workspace')).json()).state.dataset_id).toBe(expected.id);
    const state = (await (await request.get('/api/workspace')).json()).state;
    expect(state.measurements).toHaveLength(1);
    expect(state.measurements[0].trackDuringRun).toBe(true);
    expect(state.measurements[0].atoms).toEqual(expected.atoms.filter((a) => a.name === 'CA').slice(0,2).map((a) => a.index));
  };
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  const inspection = await (await request.get(`/api/datasets/${source.id}/inspection`)).json();
  expect(
    inspection.missing_atoms.some((residue: { atoms: string[] }) => residue.atoms.includes('NZ')),
  ).toBe(true);
  expect(
    source.atoms.some(
      (atom: { name: string; resid: number }) => atom.name === 'NZ' && atom.resid === 6,
    ),
  ).toBe(false);
  await studio.getByRole('button', { name: 'Preparation options', exact: true }).click();
  await expect(studio.getByRole('checkbox', { name: /Refine side-chain rotamers/ })).toBeChecked();
  await expect(
    studio.getByRole('checkbox', { name: 'Add missing heavy atoms', exact: true }),
  ).toBeChecked();

  const startedPreparation = page.waitForResponse(
    (r) => r.url().endsWith('/api/preparations') && r.request().method() === 'POST',
  );
  await studio.getByRole('button', { name: 'Prep protein', exact: true }).click();
  const prepResponse = await startedPreparation;
  expect(prepResponse.status(), await prepResponse.text()).toBe(201);
  const preparationJob = await prepResponse.json();
  expect(preparationJob.config.dataset_id).toBe(source.id);
  expect(preparationJob.config.optimize_sidechains).toBe(true);
  const preparationCard = studio
    .locator('.job-card')
    .filter({ hasText: preparationJob.name })
    .first();
  await expect(preparationCard).toBeVisible();
  await expect(preparationCard.locator('.job-status')).toHaveText('completed', {
    timeout: 100_000,
  });
  await expect(preparationCard.locator('.job-progress b')).toHaveText('100%');
  const completedPreparation = await (await request.get(`/api/jobs/${preparationJob.id}`)).json();
  expect(completedPreparation.dataset_id).toBeTruthy();
  const prepared = await (
    await request.get(`/api/datasets/${completedPreparation.dataset_id}`)
  ).json();
  await sceneReady(page, prepared.name);
  await assertTrackingPreserved(prepared);
  await expect(studio.locator('.preparation-result')).toContainText('Preparation recorded · pH 7');
  await expect(page.getByRole('button', { name: 'Polar only', exact: true })).toHaveClass(/active/);
  expect(prepared.preparation.parent_dataset_id).toBe(source.id);
  expect(prepared.preparation.simulation_ready).toBe(true);
  expect(prepared.preparation.ph).toBe(7);
  expect(
    prepared.preparation.repaired_atoms.some((residue: { atoms: string[] }) =>
      residue.atoms.includes('NZ'),
    ),
  ).toBe(true);
  expect(
    prepared.atoms.some(
      (atom: { name: string; resid: number }) => atom.name === 'NZ' && atom.resid === 6,
    ),
  ).toBe(true);
  expect(
    prepared.atoms.filter((atom: { element: string }) => atom.element === 'H').length,
  ).toBeGreaterThan(0);
  expect(
    prepared.atoms.filter((atom: { category: string }) => atom.category === 'water'),
  ).toHaveLength(0);
  await expect(page.getByRole('alert')).toHaveCount(0);

  const startedSolvation = page.waitForResponse(
    (r) =>
      r.url().endsWith(`/api/datasets/${prepared.id}/solvate`) && r.request().method() === 'POST',
    { timeout: 30_000 },
  );
  await studio
    .getByRole('combobox', { name: 'Solvent environment', exact: true })
    .selectOption('explicit');
  const solventResponse = await startedSolvation;
  expect(solventResponse.status(), await solventResponse.text()).toBe(201);
  const solvationJob = await solventResponse.json();
  const solventCard = studio.locator('.job-card').filter({ hasText: solvationJob.name }).first();
  await expect(solventCard.locator('.job-status')).toHaveText('completed', { timeout: 40_000 });
  const completedSolvation = await (await request.get(`/api/jobs/${solvationJob.id}`)).json();
  const solvated = await (
    await request.get(`/api/datasets/${completedSolvation.dataset_id}`)
  ).json();
  await sceneReady(page, solvated.name);
  await assertTrackingPreserved(solvated);
  expect(solvated.solvation.parent_dataset_id).toBe(prepared.id);
  expect(solvated.solvation.water_model).toBe('tip3p');
  expect(solvated.solvation.equilibrated).toBe(false);
  expect(solvated.has_unitcell).toBe(true);
  const waterAtoms = solvated.atoms.filter(
    (atom: { category: string }) => atom.category === 'water',
  ).length;
  expect(waterAtoms).toBeGreaterThan(1000);
  expect(waterAtoms).toBe(solvated.solvation.water_atoms);
  expect(solvated.n_atoms).toBeGreaterThan(prepared.n_atoms);
  await expect(studio.locator('.solvent-preview-card')).toContainText(
    'Explicit water box is ready',
  );
  await expect(page.getByRole('button', { name: 'Hide water', exact: true })).toHaveAttribute(
    'aria-pressed',
    'true',
  );

  await page.mouse.move(20, 20);
  const canvas = page.locator('.molecular-viewer__canvas canvas');
  const withWater = await canvas.screenshot();
  await testInfo.attach('real-explicit-water-visible', {
    body: withWater,
    contentType: 'image/png',
  });
  await page.getByRole('button', { name: 'Hide water', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Show water', exact: true })).toHaveAttribute(
    'aria-pressed',
    'false',
  );
  await page.mouse.move(20, 20);
  const withoutWater = await canvas.screenshot();
  await testInfo.attach('same-system-water-hidden', {
    body: withoutWater,
    contentType: 'image/png',
  });
  const differingPixels = await page.evaluate(
    async ([left, right]) => {
      const arrays: Uint8ClampedArray[] = [];
      for (const image of [left, right]) {
        const bitmap = await createImageBitmap(
          await (await fetch(`data:image/png;base64,${image}`)).blob(),
        );
        const canvas = document.createElement('canvas');
        canvas.width = bitmap.width;
        canvas.height = bitmap.height;
        const context = canvas.getContext('2d')!;
        context.drawImage(bitmap, 0, 0);
        arrays.push(context.getImageData(0, 0, canvas.width, canvas.height).data);
      }
      let count = 0;
      for (let i = 0; i < arrays[0].length; i += 4) {
        const difference =
          Math.abs(arrays[0][i] - arrays[1][i]) +
          Math.abs(arrays[0][i + 1] - arrays[1][i + 1]) +
          Math.abs(arrays[0][i + 2] - arrays[1][i + 2]);
        if (difference > 40) count++;
      }
      return count;
    },
    [withWater.toString('base64'), withoutWater.toString('base64')],
  );
  expect(
    differingPixels,
    'Hiding actual water must remove visible molecular geometry.',
  ).toBeGreaterThan(500);
  await page.getByRole('button', { name: 'Show water', exact: true }).click();
  await studio
    .getByRole('combobox', { name: 'Solvent environment', exact: true })
    .selectOption('implicit');
  await sceneReady(page, prepared.name);
  await expect(
    studio.getByRole('combobox', { name: 'Solvent environment', exact: true }),
  ).toHaveValue('implicit');
  await expect(page.getByRole('button', { name: 'Show water', exact: true })).toHaveAttribute(
    'aria-pressed',
    'false',
  );
  await expect(page.locator('.source-summary small')).toContainText(
    `${prepared.n_atoms.toLocaleString()} atoms`,
  );
  await expect(page.getByRole('alert')).toHaveCount(0);
  const originalStillIntact = await (await request.get(`/api/datasets/${source.id}`)).json();
  expect(originalStillIntact.n_atoms).toBe(source.n_atoms);
  expect(originalStillIntact.preparation).toBeUndefined();
  const evidence = {
    source_dataset_id: source.id,
    prepared_dataset_id: prepared.id,
    solvated_dataset_id: solvated.id,
    preparation_job_id: preparationJob.id,
    solvation_job_id: solvationJob.id,
    prepared_atoms: prepared.n_atoms,
    water_atoms: waterAtoms,
    differing_water_pixels: differingPixels,
  };
  await testInfo.attach('real-preparation-solvation-lifecycle', {
    body: JSON.stringify(evidence, null, 2),
    contentType: 'application/json',
  });
  console.log('Preparation lifecycle evidence:', JSON.stringify(evidence));
  expect(errors).toEqual([]);
});
