import { test } from './testWorkspace';
import { expect, type APIRequestContext, type Page } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

test.use({ actionTimeout: 15_000 });

const fixture = fileURLToPath(
  new URL('../../docs/audit/preparation-fixtures/six_residues_missing_atom.pdb', import.meta.url),
);

function gate() {
  let release!: () => void;
  const promise = new Promise<void>((resolve) => (release = resolve));
  return { promise, release };
}

async function openNew(page: Page) {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  return page.getByRole('dialog', { name: 'Set molecules in motion.' });
}

async function sourceFile(suffix: string) {
  return {
    name: `new-setup-${suffix}-${Date.now()}.pdb`,
    mimeType: 'chemical/x-pdb',
    buffer: await readFile(fixture),
  };
}

async function prepareApi(request: APIRequestContext) {
  const file = await sourceFile('historical');
  const uploaded = await request.post('/api/structures/upload', { multipart: { file } });
  expect(uploaded.ok(), await uploaded.text()).toBe(true);
  const source = await uploaded.json();
  const submitted = await request.post('/api/preparations', {
    data: { dataset_id: source.id, name: `${source.name} prepared`, ph: 7.4 },
  });
  expect(submitted.status(), await submitted.text()).toBe(201);
  const job = await submitted.json();
  let completed = job;
  await expect
    .poll(
      async () => {
        completed = await (await request.get(`/api/jobs/${job.id}`)).json();
        return completed.status;
      },
      { timeout: 100_000, intervals: [250, 500, 1000] },
    )
    .toBe('completed');
  const prepared = await (await request.get(`/api/datasets/${completed.dataset_id}`)).json();
  return { source, job: completed, prepared };
}

test('New simulation clears a completed preparation panel and restores the entire setup without deleting provenance', async ({
  page,
  request,
}, testInfo) => {
  test.setTimeout(150_000);
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  const studio = await openNew(page);
  const uploaded = page.waitForResponse(
    (response) =>
      response.url().endsWith('/api/structures/upload') && response.request().method() === 'POST',
  );
  await studio
    .getByLabel('Upload simulation structure', { exact: true })
    .setInputFiles(await sourceFile('complete'));
  const source = await (await uploaded).json();
  await expect(studio.locator('.source-summary strong')).toHaveText(source.name);
  await studio.getByLabel('Preparation pH', { exact: true }).fill('7.4');
  const submitted = page.waitForResponse(
    (response) =>
      response.url().endsWith('/api/preparations') && response.request().method() === 'POST',
  );
  await studio.getByRole('button', { name: 'Prep protein', exact: true }).click();
  const job = await (await submitted).json();
  const monitor = studio.locator('.structure-job-monitor');
  await expect(monitor).toContainText('Protein preparation complete', { timeout: 100_000 });
  await expect(monitor).toContainText('Prepared structure is in the viewer');
  const completed = await (await request.get(`/api/jobs/${job.id}`)).json();
  const prepared = await (await request.get(`/api/datasets/${completed.dataset_id}`)).json();
  await expect(page.locator('.structure-card h2')).toHaveText(prepared.name);
  await expect(studio.locator('.preparation-result')).toContainText(
    'Preparation recorded · pH 7.4',
  );

  await studio.getByLabel('Simulation name', { exact: true }).fill('Discard this setup');
  await studio.getByRole('spinbutton', { name: 'Duration ps', exact: true }).fill('31');
  await studio.getByRole('spinbutton', { name: 'Temperature K', exact: true }).fill('321');
  await studio.getByRole('button', { name: 'Advanced controls Optional', exact: true }).click();
  await studio.getByLabel('Random seed', { exact: true }).fill('997');
  await studio.getByLabel('Minimize energy before equilibration', { exact: true }).uncheck();
  await studio.getByRole('button', { name: 'Preparation options', exact: true }).click();
  await studio.getByLabel('Add missing heavy atoms', { exact: true }).uncheck();
  await studio.getByLabel('Refine side-chain rotamers and clashes', { exact: true }).uncheck();
  await studio.getByLabel('Preparation pH', { exact: true }).fill('9.1');
  await studio.getByRole('tab', { name: 'SMILES', exact: true }).click();
  await studio.getByLabel('SMILES molecule name').fill('Discard this molecule entry');
  await studio.getByLabel('SMILES string').fill('not a smiles');
  await studio.getByRole('button', { name: 'Build 3D structure & view', exact: true }).click();
  await expect(studio.locator('.structure-workbench [role="alert"]')).toBeVisible();
  await studio.getByRole('button', { name: 'GROMACS', exact: true }).click();
  await monitor.getByRole('button', { name: 'View protein preparation log', exact: true }).click();
  await expect(studio.locator('.job-log')).toBeVisible();
  expect(
    await studio.locator('.drawer-scroll').evaluate((element) => element.scrollTop),
  ).toBeGreaterThan(0);

  await studio.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  await expect(monitor).toHaveCount(0);
  await expect(studio.getByRole('button', { name: 'GROMACS', exact: true })).toHaveAttribute(
    'aria-pressed',
    'true',
  );
  await expect(studio.locator('.source-summary strong')).toHaveText(prepared.name);
  await expect(studio.getByRole('tab', { name: 'Upload', exact: true })).toHaveAttribute(
    'aria-selected',
    'true',
  );
  await expect(studio.locator('.drawer-scroll')).toHaveJSProperty('scrollTop', 0);
  await expect(studio.getByLabel('Simulation name', { exact: true })).toHaveValue(
    'My molecular journey',
  );
  await expect(studio.getByRole('spinbutton', { name: 'Duration ps', exact: true })).toHaveValue(
    '10',
  );
  await expect(studio.getByRole('spinbutton', { name: 'Temperature K', exact: true })).toHaveValue(
    '300',
  );
  await expect(studio.locator('.advanced-fields')).toHaveCount(0);
  await expect(studio.locator('.job-log')).toHaveCount(0);
  await expect(studio.locator('.structure-workbench [role="alert"]')).toHaveCount(0);
  const history = studio.locator('.job-card').filter({ hasText: job.name });
  await expect(history.locator('.job-status')).toHaveText('completed');
  await expect(history.getByRole('button', { name: 'Open structure', exact: true })).toBeEnabled();
  await studio.getByRole('button', { name: 'OpenMM', exact: true }).click();
  await expect(studio.getByLabel('Preparation pH', { exact: true })).toHaveValue('7.4');
  await expect(studio.locator('.preparation-result')).toContainText(
    'Preparation recorded · pH 7.4',
  );
  await expect(
    studio.getByRole('button', { name: 'Preparation options', exact: true }),
  ).toHaveAttribute('aria-expanded', 'false');
  await studio.locator('.drawer-scroll').evaluate((element) => (element.scrollTop = 0));
  await testInfo.attach('new-simulation-fresh-panel', {
    body: await page.screenshot(),
    contentType: 'image/png',
  });
  await studio.getByRole('button', { name: 'Preparation options', exact: true }).click();
  await expect(studio.getByLabel('Add missing heavy atoms', { exact: true })).toBeChecked();
  await expect(
    studio.getByLabel('Refine side-chain rotamers and clashes', { exact: true }),
  ).toBeChecked();
  await studio.getByRole('button', { name: 'Advanced controls Optional', exact: true }).click();
  await expect(studio.getByLabel('Random seed', { exact: true })).toHaveValue('42');
  await expect(
    studio.getByLabel('Minimize energy before equilibration', { exact: true }),
  ).toBeChecked();
  await studio.getByRole('tab', { name: 'SMILES', exact: true }).click();
  await expect(studio.getByLabel('SMILES molecule name')).toHaveValue('');
  await expect(studio.getByLabel('SMILES string')).toHaveValue('');
  await expect(monitor).toHaveCount(0);
  expect((await (await request.get(`/api/datasets/${prepared.id}`)).json()).preparation).toEqual(
    prepared.preparation,
  );
  expect((await (await request.get(`/api/jobs/${job.id}`)).json()).status).toBe('completed');
  expect(errors).toEqual([]);
  await testInfo.attach('real-preparation-preserved-after-new-setup', {
    body: JSON.stringify({ source: source.id, job: job.id, prepared: prepared.id }, null, 2),
    contentType: 'application/json',
  });
});

test('a historical completed job arriving after New cannot restore the old preparation banner', async ({
  page,
  request,
}) => {
  test.setTimeout(150_000);
  const { job, prepared } = await prepareApi(request);
  const saved = await request.post('/api/workspace', {
    data: { state: { version: 1, dataset_id: prepared.id } },
  });
  expect(saved.ok()).toBe(true);
  const reached = gate();
  const released = gate();
  const delivered = gate();
  await page.route('**/api/jobs', async (route) => {
    const response = await route.fetch();
    reached.release();
    await released.promise;
    await route.fulfill({ response });
    delivered.release();
  });
  try {
    const studio = await openNew(page);
    await reached.promise;
    await expect(studio.locator('.structure-job-monitor')).toHaveCount(0);
    released.release();
    await delivered.promise;
    const history = studio.locator('.job-card').filter({ hasText: job.name });
    await expect(history.locator('.job-status')).toHaveText('completed');
    await expect(studio.locator('.source-summary strong')).toHaveText(prepared.name);
    await expect(studio.locator('.structure-job-monitor')).toHaveCount(0);
    // A subsequent real polling response must remain equally harmless.
    await page.waitForResponse((response) => response.url().endsWith('/api/jobs'));
    await expect(studio.locator('.structure-job-monitor')).toHaveCount(0);
  } finally {
    released.release();
    await page.unrouteAll({ behavior: 'wait' });
  }
});

for (const delayedStage of ['upload', 'coordinates'] as const) {
  test(`a late ${delayedStage} response from an earlier setup cannot replace the molecule after New`, async ({
    page,
    request,
  }) => {
    const studio = await openNew(page);
    const originalName = await page.locator('.structure-card h2').innerText();
    const original = await (await request.get('/api/workspace')).json();
    const reached = gate();
    const released = gate();
    const delivered = gate();
    let importedId = '';
    await page.route('**/api/structures/upload', async (route) => {
      const response = await route.fetch();
      importedId = (await response.json()).id;
      if (delayedStage === 'upload') {
        reached.release();
        await released.promise;
      }
      await route.fulfill({ response });
      if (delayedStage === 'upload') delivered.release();
    });
    if (delayedStage === 'coordinates') {
      await page.route('**/api/datasets/*/coordinates', async (route) => {
        const response = await route.fetch();
        reached.release();
        await released.promise;
        await route.fulfill({ response });
        delivered.release();
      });
    }
    try {
      await studio
        .getByLabel('Upload simulation structure', { exact: true })
        .setInputFiles(await sourceFile(`late-${delayedStage}`));
      await reached.promise;
      await expect(studio.getByText('Preparing your request…', { exact: true })).toBeVisible();
      await studio.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
      await page.getByRole('button', { name: 'New simulation', exact: true }).click();
      await expect(studio.getByText('Preparing your request…', { exact: true })).toHaveCount(0);
      await expect(studio.locator('.source-summary strong')).toHaveText(originalName);
      released.release();
      await delivered.promise;
      await page.waitForResponse((response) => response.url().endsWith('/api/jobs'));
      await expect(page.locator('.structure-card h2')).toHaveText(originalName);
      await expect(studio.locator('.source-summary strong')).toHaveText(originalName);
      await expect(page.locator('.canvas-loading')).toHaveCount(0);
      await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
      expect((await (await request.get('/api/workspace')).json()).state.dataset_id).toBe(
        original.state.dataset_id,
      );
      // The imported file remains on disk for deliberate future use; only stale scene mutation is rejected.
      expect(importedId).not.toBe(original.state.dataset_id);
      expect((await request.get(`/api/datasets/${importedId}`)).ok()).toBe(true);
    } finally {
      released.release();
      await page.unrouteAll({ behavior: 'wait' });
    }
  });
}

test('New ignores an older history load while keeping that result available for deliberate reopening', async ({
  page,
  request,
}) => {
  test.setTimeout(150_000);
  const { job, prepared } = await prepareApi(request);
  const studio = await openNew(page);
  const originalName = await page.locator('.structure-card h2').innerText();
  const reached = gate();
  const released = gate();
  const delivered = gate();
  const resultUrl = `**/api/datasets/${prepared.id}`;
  await page.route(resultUrl, async (route) => {
    const response = await route.fetch();
    reached.release();
    await released.promise;
    await route.fulfill({ response });
    delivered.release();
  });
  try {
    const history = studio.locator('.job-card').filter({ hasText: job.name });
    await history.getByRole('button', { name: 'Open structure', exact: true }).click();
    await reached.promise;
    await expect(history.getByRole('button', { name: 'Opening…', exact: true })).toBeDisabled();
    await expect(history.locator('.spin')).toBeVisible();
    await studio.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
    await page.getByRole('button', { name: 'New simulation', exact: true }).click();
    await expect(
      history.getByRole('button', { name: 'Open structure', exact: true }),
    ).toBeEnabled();
    released.release();
    await delivered.promise;
    await page.waitForResponse((response) => response.url().endsWith('/api/jobs'));
    await expect(page.locator('.structure-card h2')).toHaveText(originalName);
    await expect(studio.locator('.structure-job-monitor')).toHaveCount(0);
    await page.unroute(resultUrl);
    await history.getByRole('button', { name: 'Open structure', exact: true }).click();
    await expect(page.locator('.structure-card h2')).toHaveText(prepared.name);
    await expect(studio.locator('.source-summary strong')).toHaveText(prepared.name);
    await expect(page.locator('.canvas-loading')).toHaveCount(0);
    await expect(
      history.getByRole('button', { name: 'Open structure', exact: true }),
    ).toBeEnabled();
    await expect(studio.locator('.preparation-result')).toContainText(
      'Preparation recorded · pH 7.4',
    );
  } finally {
    released.release();
    await page.unrouteAll({ behavior: 'wait' });
  }
});
