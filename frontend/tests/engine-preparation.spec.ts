import { test } from './testWorkspace';
import { expect, type Locator, type Page } from '@playwright/test';
import { readFile, writeFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

const proteinFixture = fileURLToPath(
  new URL('../../docs/audit/preparation-fixtures/six_residues_intact.pdb', import.meta.url),
);
const studioName = 'Set molecules in motion.';
const engineButton = (studio: Locator, engine: 'OpenMM' | 'GROMACS') =>
  studio.getByRole('button', { name: engine, exact: true });

function gate() {
  let release!: () => void;
  const promise = new Promise<void>((resolve) => {
    release = resolve;
  });
  return { promise, release };
}

async function openStudio(page: Page) {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  return page.getByRole('dialog', { name: studioName });
}

async function uploadProtein(page: Page, studio: Locator, suffix: string) {
  await studio.getByRole('tab', { name: 'Upload', exact: true }).click();
  const received = page.waitForResponse(
    (response) =>
      response.url().endsWith('/api/structures/upload') && response.request().method() === 'POST',
  );
  await studio.getByLabel('Upload simulation structure', { exact: true }).setInputFiles({
    name: `engine-first-${suffix}-${Date.now()}.pdb`,
    mimeType: 'chemical/x-pdb',
    buffer: await readFile(proteinFixture),
  });
  const response = await received;
  expect(response.ok(), await response.text()).toBe(true);
  const dataset = await response.json();
  await expect(studio.locator('.source-summary strong')).toHaveText(dataset.name);
  await expect(page.locator('.structure-card h2')).toHaveText(dataset.name);
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  return dataset;
}

async function expectGromacsWorkflow(studio: Locator) {
  await expect(engineButton(studio, 'GROMACS')).toHaveAttribute('aria-pressed', 'true');
  await expect(engineButton(studio, 'OpenMM')).toHaveAttribute('aria-pressed', 'false');
  const preparation = studio.getByLabel('GROMACS preparation', { exact: true });
  await expect(preparation).toContainText('Native setup at Start');
  await expect(preparation).toContainText('complete standard protein');
  await expect(preparation).toContainText('default protonation');
  await expect(preparation).toContainText('Existing hydrogen choices are not retained');
  await expect(studio.getByRole('button', { name: /^Prep (protein|complex)$/ })).toHaveCount(0);
  await expect(studio.getByLabel('Preparation pH', { exact: true })).toHaveCount(0);
  await expect(studio.getByLabel('Preparation readiness', { exact: true })).toHaveCount(0);
  await expect(studio.getByRole('button', { name: /^(Build water|Update box)$/ })).toHaveCount(0);
  const solvent = studio.getByRole('combobox', { name: 'Solvent environment', exact: true });
  await expect(solvent).toHaveValue('explicit');
  await expect(solvent).toBeDisabled();
  await expect(studio.locator('.solvent-preview-card')).toContainText('GROMACS');
}

for (const viewport of [
  { width: 1440, height: 1000 },
  { width: 547, height: 637 },
]) {
  test(`engine selection precedes structure preparation at ${viewport.width}px`, async ({
    page,
  }, info) => {
    await page.setViewportSize(viewport);
    const studio = await openStudio(page);
    const openmm = engineButton(studio, 'OpenMM');
    const gromacs = engineButton(studio, 'GROMACS');
    await expect(openmm).toHaveAttribute('aria-pressed', 'true');
    await expect(openmm).toBeInViewport({ ratio: 0.95 });
    await expect(gromacs).toBeInViewport({ ratio: 0.95 });
    await expect(openmm.locator('strong')).toBeInViewport({ ratio: 1 });
    await expect(gromacs.locator('strong')).toBeInViewport({ ratio: 1 });
    const positions = await studio.evaluate((element) => {
      const engines = element.querySelector('.engine-options')!;
      const workbench = element.querySelector('.structure-workbench')!;
      return {
        enginesBeforeWorkbench: !!(
          engines.compareDocumentPosition(workbench) & Node.DOCUMENT_POSITION_FOLLOWING
        ),
        enginesBottom: engines.getBoundingClientRect().bottom,
        structureTop: workbench.getBoundingClientRect().top,
      };
    });
    expect(positions.enginesBeforeWorkbench).toBe(true);
    expect(positions.enginesBottom).toBeLessThan(positions.structureTop);
    await expect(studio.getByLabel('Preparation pH', { exact: true })).toBeVisible();
    await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeVisible();
    await gromacs.click();
    await expectGromacsWorkflow(studio);
    await studio.locator('.drawer-scroll').evaluate((element) => {
      element.scrollTop = 0;
    });
    await expect(gromacs).toBeInViewport({ ratio: 0.95 });
    const screenshot = info.outputPath(`engine-first-${viewport.width}.png`);
    await page.screenshot({ path: screenshot });
    await info.attach(`engine-first-${viewport.width}.png`, {
      path: screenshot,
      contentType: 'image/png',
    });
  });
}

test('GROMACS retains source loading and its engine choice across reopening and raw import', async ({
  page,
}, info) => {
  const studio = await openStudio(page);
  let preparationPosts = 0;
  let solvationPosts = 0;
  page.on('request', (request) => {
    if (request.method() !== 'POST') return;
    if (request.url().endsWith('/api/preparations')) preparationPosts++;
    if (/\/api\/datasets\/[^/]+\/solvate$/.test(request.url())) solvationPosts++;
  });
  await engineButton(studio, 'GROMACS').click();
  await expectGromacsWorkflow(studio);
  await studio.getByRole('tab', { name: 'Fetch', exact: true }).click();
  await expect(studio.getByRole('textbox', { name: 'PDB accession', exact: true })).toBeEnabled();
  await studio.getByRole('tab', { name: 'SMILES', exact: true }).click();
  await expect(studio.getByRole('textbox', { name: 'SMILES string', exact: true })).toBeEnabled();
  await studio.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  await expectGromacsWorkflow(studio);
  const source = await uploadProtein(page, studio, 'gromacs');
  await expectGromacsWorkflow(studio);
  const readiness = studio.getByLabel('Simulation readiness', { exact: true });
  await expect(readiness).toContainText('Ready to run');
  await readiness.getByText('Model and preparation details', { exact: true }).click();
  await expect(readiness).toContainText('GROMACS · CPU · fixed-volume NVT · explicit solvent');
  await expect(studio.locator('.engine-preparation-note')).toContainText('Amber99SB-ILDN');
  await expect(readiness).toContainText('default protonation');
  await expect(readiness).not.toContainText('Requested/prepared pH');
  await expect(studio.getByRole('button', { name: 'Start simulation', exact: true })).toBeEnabled();
  await engineButton(studio, 'GROMACS').click();
  await expect(readiness).toContainText('Ready to run');
  await expect(studio.getByRole('button', { name: 'Start simulation', exact: true })).toBeEnabled();
  expect(preparationPosts).toBe(0);
  expect(solvationPosts).toBe(0);
  await engineButton(studio, 'OpenMM').click();
  await expect(engineButton(studio, 'OpenMM')).toHaveAttribute('aria-pressed', 'true');
  await expect(studio.getByLabel('Preparation pH', { exact: true })).toBeVisible();
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  await expect(studio.getByLabel('GROMACS preparation', { exact: true })).toHaveCount(0);
  const evidence = info.outputPath('engine-choice-real-source.json');
  await writeFile(
    evidence,
    JSON.stringify({ dataset_id: source.id, preparationPosts, solvationPosts }, null, 2),
  );
  await info.attach('engine-choice-real-source.json', {
    path: evidence,
    contentType: 'application/json',
  });
});

test('GROMACS submission sends native explicit settings and locks engine switching while pending', async ({
  page,
}) => {
  const studio = await openStudio(page);
  await engineButton(studio, 'GROMACS').click();
  await uploadProtein(page, studio, 'submission');
  const readiness = studio.getByLabel('Simulation readiness', { exact: true });
  await expect(readiness).toContainText('Ready to run');
  const reached = gate();
  const release = gate();
  let submitted: Record<string, unknown> | null = null;
  await page.route('**/api/jobs', async (route) => {
    if (route.request().method() !== 'POST') return route.continue();
    submitted = route.request().postDataJSON();
    reached.release();
    await release.promise;
    // A transport failure checks the submission lock without launching molecular dynamics.
    await route.fulfill({
      status: 503,
      contentType: 'application/json',
      body: JSON.stringify({ detail: 'Engine submission transport test' }),
    });
  });
  try {
    await studio.getByRole('button', { name: 'Start simulation', exact: true }).click();
    await reached.promise;
    await expect(engineButton(studio, 'OpenMM')).toBeDisabled();
    await expect(engineButton(studio, 'GROMACS')).toBeDisabled();
    expect(submitted).toMatchObject({ engine: 'gromacs', solvent: 'explicit' });
    release.release();
    await expect(studio).toContainText('Engine submission transport test');
    await expect(engineButton(studio, 'OpenMM')).toBeEnabled();
    await expect(engineButton(studio, 'GROMACS')).toHaveAttribute('aria-pressed', 'true');
  } finally {
    release.release();
    await page.unrouteAll({ behavior: 'wait' });
  }
});

test('real prepared and solvated OpenMM states block GROMACS with a direct recovery, and jobs lock switching', async ({
  page,
  request,
}, info) => {
  test.setTimeout(180_000);
  const studio = await openStudio(page);
  const source = await uploadProtein(page, studio, 'preparation');
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  const submitReached = gate();
  const submitRelease = gate();
  const pollRelease = gate();
  await page.route('**/api/preparations', async (route) => {
    submitReached.release();
    await submitRelease.promise;
    await route.continue();
  });
  await page.route('**/api/jobs', async (route) => {
    if (route.request().method() === 'GET') await pollRelease.promise;
    await route.continue();
  });
  try {
    const submitted = page.waitForResponse(
      (response) =>
        response.url().endsWith('/api/preparations') && response.request().method() === 'POST',
    );
    await studio.getByRole('button', { name: 'Prep protein', exact: true }).click();
    await submitReached.promise;
    await expect(engineButton(studio, 'GROMACS')).toBeDisabled();
    await expect(engineButton(studio, 'OpenMM')).toBeDisabled();
    submitRelease.release();
    const response = await submitted;
    expect(response.status(), await response.text()).toBe(201);
    const preparationJob = await response.json();
    await expect(studio.getByLabel('Protein preparation monitor', { exact: true })).toBeVisible();
    await expect(engineButton(studio, 'GROMACS')).toBeDisabled();
    await studio.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
    await page.getByRole('button', { name: 'New simulation', exact: true }).click();
    await expect(engineButton(studio, 'GROMACS')).toBeDisabled();
    pollRelease.release();
    await expect(studio.locator('.structure-job-monitor')).toContainText(
      'Prepared structure is in the viewer',
      { timeout: 100_000 },
    );
    await expect(engineButton(studio, 'GROMACS')).toBeEnabled();
    const completed = await (await request.get(`/api/jobs/${preparationJob.id}`)).json();
    expect(completed.status).toBe('completed');
    const prepared = await (await request.get(`/api/datasets/${completed.dataset_id}`)).json();
    expect(prepared.preparation.parent_dataset_id).toBe(source.id);
    expect(prepared.solvation).toBeFalsy();
    await engineButton(studio, 'GROMACS').click();
    await expectGromacsWorkflow(studio);
    const warning = studio.locator('.engine-compatibility');
    await expect(warning).toContainText('cannot yet transfer that state to GROMACS');
    await expect(
      studio.getByRole('button', { name: 'Start simulation', exact: true }),
    ).toBeDisabled();
    await warning
      .getByRole('button', { name: 'Use OpenMM for this structure', exact: true })
      .click();
    await expect(engineButton(studio, 'OpenMM')).toHaveAttribute('aria-pressed', 'true');
    await expect(warning).toHaveCount(0);
    await expect(studio.locator('.source-summary strong')).toHaveText(prepared.name);

    const solvating = page.waitForResponse(
      (result) =>
        result.url().endsWith(`/api/datasets/${prepared.id}/solvate`) &&
        result.request().method() === 'POST',
    );
    const waterButton = studio.getByRole('button', { name: 'Build water', exact: true });
    if (await waterButton.count()) await waterButton.click();
    else
      await studio
        .getByRole('combobox', { name: 'Solvent environment', exact: true })
        .selectOption('explicit');
    const solventResponse = await solvating;
    expect(solventResponse.status(), await solventResponse.text()).toBe(201);
    const solvationJob = await solventResponse.json();
    await expect(engineButton(studio, 'GROMACS')).toBeDisabled();
    await expect(studio.locator('.structure-job-monitor')).toContainText(
      'The explicit-water structure is in the viewer',
      { timeout: 100_000 },
    );
    await expect(engineButton(studio, 'GROMACS')).toBeEnabled();
    const solvatedJob = await (await request.get(`/api/jobs/${solvationJob.id}`)).json();
    expect(solvatedJob.status).toBe('completed');
    const solvated = await (await request.get(`/api/datasets/${solvatedJob.dataset_id}`)).json();
    expect(solvated.solvation.parent_dataset_id).toBe(prepared.id);
    expect(solvated.atoms.some((atom: { category: string }) => atom.category === 'water')).toBe(
      true,
    );
    await engineButton(studio, 'GROMACS').click();
    await expectGromacsWorkflow(studio);
    await expect(warning).toContainText('saved solvent box');
    await expect(
      studio.getByRole('button', { name: 'Start simulation', exact: true }),
    ).toBeDisabled();
    const evidence = info.outputPath('engine-preparation-real-results.json');
    await writeFile(
      evidence,
      JSON.stringify(
        {
          source_id: source.id,
          preparation_job: preparationJob.id,
          prepared_id: prepared.id,
          solvation_job: solvationJob.id,
          solvated_id: solvated.id,
        },
        null,
        2,
      ),
    );
    await info.attach('engine-preparation-real-results.json', {
      path: evidence,
      contentType: 'application/json',
    });
  } finally {
    submitRelease.release();
    pollRelease.release();
    await page.unrouteAll({ behavior: 'wait' });
  }
});
