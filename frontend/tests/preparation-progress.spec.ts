import { test } from './testWorkspace';
import { expect, type Page, type Locator } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

const fixture = fileURLToPath(
  new URL('../../docs/audit/preparation-fixtures/six_residues_missing_atom.pdb', import.meta.url),
);

function gate() {
  let release!: () => void;
  const promise = new Promise<void>((resolve) => {
    release = resolve;
  });
  return { promise, release };
}

async function importProtein(page: Page, suffix: string) {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const studio = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  const name = `e2e-prep-monitor-${suffix}-${Date.now()}`;
  const uploaded = page.waitForResponse(
    (response) =>
      response.url().endsWith('/api/structures/upload') && response.request().method() === 'POST',
  );
  await studio.getByLabel('Upload simulation structure', { exact: true }).setInputFiles({
    name: `${name}.pdb`,
    mimeType: 'chemical/x-pdb',
    buffer: await readFile(fixture),
  });
  const source = await (await uploaded).json();
  console.log('Preparation monitor source:', JSON.stringify({ id: source.id, name: source.name }));
  await expect(page.locator('.structure-card h2')).toHaveText(name);
  await expect(studio.locator('.source-summary strong')).toHaveText(name);
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  return { studio, source };
}

async function expectPinnedMonitor(page: Page, studio: Locator) {
  const monitor = studio.locator('.structure-job-monitor');
  const scroller = studio.locator('.drawer-scroll');
  await expect(monitor).toBeVisible();
  expect(await monitor.evaluate((element) => !!element.closest('.drawer-scroll'))).toBe(false);
  await scroller.evaluate((element) => {
    element.scrollTop = element.scrollHeight;
  });
  await expect(monitor).toBeInViewport({ ratio: 1 });
  const bounds = await monitor.boundingBox();
  expect(bounds).not.toBeNull();
  expect(bounds!.y).toBeGreaterThanOrEqual(0);
  expect(bounds!.y + bounds!.height).toBeLessThanOrEqual(page.viewportSize()!.height);
}

test('preparation has immediate feedback, survives reopening Studio, and announces the real viewed result', async ({
  page,
  request,
}, testInfo) => {
  test.setTimeout(150_000);
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  const { studio, source } = await importProtein(page, 'success');
  const postGate = gate();
  const postReached = gate();
  const pollingGate = gate();
  const resultGate = gate();
  const resultReached = gate();

  // Delay transport only: every job, stage, error, and result comes from the real API.
  await page.route('**/api/preparations', async (route) => {
    postReached.release();
    await postGate.promise;
    await route.continue();
  });
  // The small repair can finish in under a second. Hold polling while checking the actual
  // queued response so reopening/log assertions do not depend on machine speed.
  await page.route('**/api/jobs', async (route) => {
    if (route.request().method() === 'GET') await pollingGate.promise;
    await route.continue();
  });
  await page.route(/\/api\/datasets\/[a-f0-9]{16}$/, async (route) => {
    if (new URL(route.request().url()).pathname !== `/api/datasets/${source.id}`) {
      resultReached.release();
      await resultGate.promise;
    }
    await route.continue();
  });

  try {
    const submitted = page.waitForResponse(
      (response) =>
        response.url().endsWith('/api/preparations') && response.request().method() === 'POST',
    );
    await studio.getByRole('button', { name: 'Prep protein', exact: true }).click();
    await postReached.promise;
    const monitor = studio.locator('.structure-job-monitor');
    await expect(monitor).toContainText('Starting protein preparation…');
    const starting = studio.getByRole('button', { name: 'Starting preparation…', exact: true });
    await expect(starting).toBeDisabled();
    await expect(starting.locator('.spin')).toBeVisible();
    await expect(monitor.locator('.structure-job-monitor-spinner')).toBeVisible();

    await page.setViewportSize({ width: 547, height: 637 });
    await expectPinnedMonitor(page, studio);
    await testInfo.attach('preparation-starting-narrow-window', {
      body: await page.screenshot(),
      contentType: 'image/png',
    });
    await page.setViewportSize({ width: 1440, height: 1000 });
    postGate.release();
    const response = await submitted;
    expect(response.status(), await response.text()).toBe(201);
    const job = await response.json();
    console.log('Preparation monitor job:', job.id);
    expect(job.config.dataset_id).toBe(source.id);

    const progress = monitor.getByRole('progressbar', { name: 'Protein preparation progress' });
    await expect(progress).toBeVisible();
    await expect(
      studio.getByRole('button', { name: 'Preparing protein…', exact: true }),
    ).toBeDisabled();
    await expect(
      monitor.getByRole('button', { name: 'Cancel protein preparation', exact: true }),
    ).toBeEnabled();
    expect(job.status).toBe('queued');
    expect(job.progress).toBe(0);
    await expect(monitor).toContainText(job.stage);
    await expect(progress).toHaveAttribute('aria-valuemin', '0');
    await expect(progress).toHaveAttribute('aria-valuemax', '100');

    await studio.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
    await expect(studio).toHaveCount(0);
    await page.getByRole('button', { name: 'New simulation', exact: true }).click();
    await expect(monitor).toBeVisible();
    await expect(progress).toBeVisible();
    await expect(
      monitor.getByRole('button', { name: 'Cancel protein preparation', exact: true }),
    ).toBeEnabled();
    await monitor
      .getByRole('button', { name: 'View protein preparation log', exact: true })
      .click();
    const jobCard = studio.locator('.job-card').filter({ hasText: job.name }).first();
    await expect(jobCard).toBeInViewport();
    await expectPinnedMonitor(page, studio);

    pollingGate.release();
    await resultReached.promise;
    // The worker can finish before the dataset fetch/view load; do not claim it is shown early.
    await expect(monitor).not.toContainText('Prepared structure is in the viewer');
    resultGate.release();
    await expect(monitor).toContainText('Protein preparation complete', { timeout: 100_000 });
    await expect(monitor).toContainText('Prepared structure is in the viewer');
    await expect(progress).toHaveAttribute('aria-valuenow', '100');
    await expect(monitor).toContainText('100%');
    await expect(
      monitor.getByRole('button', { name: 'Cancel protein preparation', exact: true }),
    ).toHaveCount(0);
    const completed = await (await request.get(`/api/jobs/${job.id}`)).json();
    expect(completed.status).toBe('completed');
    const prepared = await (await request.get(`/api/datasets/${completed.dataset_id}`)).json();
    expect(prepared.preparation.parent_dataset_id).toBe(source.id);
    await expect(page.locator('.structure-card h2')).toHaveText(prepared.name);
    await expect(page.locator('.canvas-loading')).toHaveCount(0);
    await expect(studio.locator('.source-summary strong')).toHaveText(prepared.name);
    await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();

    await page.setViewportSize({ width: 547, height: 637 });
    await expectPinnedMonitor(page, studio);
    await testInfo.attach('preparation-complete-narrow-window', {
      body: await page.screenshot(),
      contentType: 'image/png',
    });
    await testInfo.attach('preparation-progress-real-job', {
      body: JSON.stringify(
        { source_id: source.id, job_id: job.id, prepared_id: prepared.id },
        null,
        2,
      ),
      contentType: 'application/json',
    });
    await monitor.getByRole('button', { name: 'Dismiss preparation status', exact: true }).click();
    await expect(monitor).toHaveCount(0);
    await expect(page.locator('.structure-card h2')).toHaveText(prepared.name);
    expect(errors).toEqual([]);
  } finally {
    postGate.release();
    pollingGate.release();
    resultGate.release();
    await page.unrouteAll({ behavior: 'wait' });
  }
});

test('preparation readiness blocks incomplete atoms before submission and clears when repair is restored', async ({
  page,
  request,
}, testInfo) => {
  const { studio, source } = await importProtein(page, 'validation');
  let submissions = 0;
  page.on('request', (response) => {
    if (response.url().endsWith('/api/preparations') && response.method() === 'POST') submissions++;
  });
  await studio.getByRole('button', { name: 'Preparation options', exact: true }).click();
  const repair = studio.getByRole('checkbox', { name: 'Add missing heavy atoms', exact: true });
  await repair.uncheck();
  const readiness = studio.getByLabel('Preparation readiness', { exact: true });
  await expect(readiness).toContainText('Before you continue');
  await expect(readiness).toContainText(
    'Missing heavy/terminal atoms prevent force-field preparation',
  );
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeDisabled();
  await expect(studio.locator('.structure-job-monitor')).toHaveCount(0);
  expect(submissions).toBe(0);
  const jobs = await (await request.get('/api/jobs')).json();
  expect(
    jobs.filter((job: { config: { dataset_id: string } }) => job.config.dataset_id === source.id),
  ).toEqual([]);
  await expect(page.locator('.structure-card h2')).toHaveText(source.name);
  await page.setViewportSize({ width: 547, height: 637 });
  await readiness.scrollIntoViewIfNeeded();
  await expect(readiness).toBeInViewport();
  await testInfo.attach('preparation-readiness-error-narrow-window', {
    body: await page.screenshot(),
    contentType: 'image/png',
  });
  await repair.check();
  await expect(readiness).toContainText('Ready to prepare');
  await expect(studio.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  expect(submissions).toBe(0);
});

test('a prepared topology display failure stops the spinner and retries the same real result', async ({
  page,
  request,
}, testInfo) => {
  test.setTimeout(150_000);
  const { studio, source } = await importProtein(page, 'viewer-retry');
  let failedTopologyUrl = '';
  let allowRetry = false;
  await page.route('**/api/datasets/*/topology', async (route) => {
    if (!allowRetry && new URL(route.request().url()).pathname !== source.topology_url) {
      // Keep the real preparation, metadata, coordinates, and PDB intact. Simulate a
      // transport outage for the first display attempt, then restore it for manual retry.
      failedTopologyUrl = route.request().url();
      await route.abort('failed');
      return;
    }
    await route.continue();
  });
  try {
    const submitted = page.waitForResponse(
      (response) =>
        response.url().endsWith('/api/preparations') && response.request().method() === 'POST',
    );
    await studio.getByRole('button', { name: 'Prep protein', exact: true }).click();
    const response = await submitted;
    expect(response.status(), await response.text()).toBe(201);
    const job = await response.json();
    console.log('Preparation monitor retry job:', job.id);
    const monitor = studio.locator('.structure-job-monitor');
    await expect(monitor).toContainText('Could not display the prepared structure:', {
      timeout: 100_000,
    });
    expect(failedTopologyUrl).toMatch(/\/api\/datasets\/[a-f0-9]{16}\/topology$/);
    await expect(monitor).toContainText('Failed to fetch');
    await expect(monitor).not.toContainText('Prepared structure is in the viewer');
    await expect(monitor.locator('.structure-job-monitor-spinner')).toHaveCount(0);
    const retry = monitor.getByRole('button', { name: 'Open result', exact: true });
    await expect(retry).toBeEnabled();
    await expect(
      monitor.getByRole('progressbar', { name: 'Protein preparation progress', exact: true }),
    ).toHaveAttribute('aria-valuenow', '100');

    const completed = await (await request.get(`/api/jobs/${job.id}`)).json();
    expect(completed.status).toBe('completed');
    const prepared = await (await request.get(`/api/datasets/${completed.dataset_id}`)).json();
    expect(prepared.preparation.parent_dataset_id).toBe(source.id);
    expect(new URL(failedTopologyUrl).pathname).toBe(prepared.topology_url);
    await testInfo.attach('preparation-viewer-fetch-failure', {
      body: await page.screenshot(),
      contentType: 'image/png',
    });

    const reloadedTopology = page.waitForResponse(
      (result) => result.url() === failedTopologyUrl && result.request().method() === 'GET',
    );
    allowRetry = true;
    await retry.click();
    expect((await reloadedTopology).status()).toBe(200);
    await expect(monitor).toContainText('Protein preparation complete');
    await expect(monitor).toContainText('Prepared structure is in the viewer');
    await expect(monitor).not.toContainText('Could not display the prepared structure:');
    await expect(monitor.locator('.structure-job-monitor-spinner')).toHaveCount(0);
    await expect(retry).toHaveCount(0);
    await expect(page.locator('.canvas-loading')).toHaveCount(0);
    await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
    await expect(page.locator('.structure-card h2')).toHaveText(prepared.name);
    await expect(studio.locator('.source-summary strong')).toHaveText(prepared.name);
    await testInfo.attach('preparation-viewer-retry-real-job', {
      body: JSON.stringify(
        {
          source_id: source.id,
          job_id: job.id,
          prepared_id: prepared.id,
          failed_topology: failedTopologyUrl,
        },
        null,
        2,
      ),
      contentType: 'application/json',
    });
  } finally {
    await page.unrouteAll({ behavior: 'wait' });
  }
});
