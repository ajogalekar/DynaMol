import { test, expect } from '@playwright/test';

test('1UA2 distinguishes covalent TPO residues from ATP ligands and keeps genuine gap warnings', async ({
  page,
  request,
}) => {
  const id = process.env.DYNAMOL_MODIFIED_DATASET || 'f491cb784d7f4fbf';
  const response = await request.get(`/api/datasets/${id}`);
  expect(
    response.ok(),
    'The retained 1UA2 source dataset must be available for this regression.',
  ).toBe(true);
  const dataset = await response.json();
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByTitle('Switch dataset', { exact: true }).click();
  await page
    .locator('.library-popover button')
    .filter({
      has: page.locator('b', {
        hasText: new RegExp(`^${dataset.name.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}$`),
      }),
    })
    .click();
  await expect(page.locator('.structure-card h2')).toHaveText(dataset.name);
  await page.getByRole('button', { name: 'Simulate', exact: true }).click();
  const studio = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  await expect(studio.locator('.modified-residue-card')).toHaveCount(4);
  await expect(studio.locator('.ligand-preparation-card')).toHaveCount(4);
  for (let i = 0; i < 4; i++) {
    const modified = studio.locator('.modified-residue-card').nth(i);
    await expect(modified.locator('summary')).toContainText('TPO');
    await expect(modified.locator('summary')).toContainText('Protein template');
    await modified.locator('summary').click();
    await expect(modified).toContainText('phosaa14SB');
    await expect(modified).toContainText('Fixed residue charge: -2');
    await expect(
      studio.locator('.ligand-preparation-card').nth(i).locator('summary'),
    ).toContainText('ATP');
  }
  await expect(studio.locator('.inspection-summary')).toContainText('4 chain gaps');
  await expect(studio.locator('.modified-residue-list')).toContainText(
    'Ligand removal preserves modified protein residues',
  );
  await studio.getByRole('button', { name: 'Preparation options', exact: true }).click();
  await studio
    .getByRole('checkbox', { name: 'Remove ligands and other non-protein residues', exact: true })
    .check();
  await expect(studio.locator('.ligand-preparation-intro')).toContainText('removal selected');
  await expect(studio.locator('.modified-residue-card')).toHaveCount(4);
  await expect(
    studio.getByRole('checkbox', { name: 'Build supported missing loops / residues', exact: true }),
  ).toBeDisabled();
  await studio
    .getByRole('checkbox', { name: 'Remove ligands and other non-protein residues', exact: true })
    .uncheck();
  await expect(studio.getByRole('button', { name: 'Prep complex', exact: true })).toBeVisible();
});
