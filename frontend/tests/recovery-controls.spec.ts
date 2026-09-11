import { test, expect } from '@playwright/test';
import { stat } from 'node:fs/promises';

// Run only on an explicitly supplied isolated API/data directory. The input is
// an already minimized standard protein fixture, so this checks job controls
// without repeating lengthy molecular preparation or touching user projects.
const input = process.env.DYNAMOL_RECOVERY_TEST_INPUT;
test('Stop saves a checkpoint, Resume continues it, and terminal jobs retain Files', async ({
  page,
  request,
  baseURL,
}) => {
  test.skip(!input, 'Requires DYNAMOL_RECOVERY_TEST_INPUT in an isolated test server.');
  test.setTimeout(180_000);
  expect(process.env.DYNAMOL_TEST_ISOLATED).toBe('1');
  expect(['8765', '5173', '4173']).not.toContain(new URL(baseURL!).port);
  const initialized = await request.post('/api/workspace', {
    data: { state: { version: 1, dataset_id: input, measurements: [] } },
  });
  expect(initialized.ok()).toBe(true);

  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const readiness = page.getByLabel('Simulation readiness', { exact: true });
  await expect(readiness).toContainText('Ready to run');
  const response = await request.post('/api/jobs', {
    data: {
      dataset_id: input,
      engine: 'openmm',
      name: `Recovery controls ${Date.now()}`,
      duration_ps: 1,
      temperature_k: 300,
      timestep_fs: 2,
      report_interval: 50,
      friction_ps: 1,
      seed: 2026,
      solvent: 'implicit',
      minimize: false,
      equilibration_steps: 0,
      padding_nm: 1,
    },
  });
  expect(response.ok()).toBe(true);
  const job = await response.json();
  const card = page
    .locator('.job-card')
    .filter({ has: page.locator('strong', { hasText: job.name }) });
  try {
    await expect(card).toBeVisible();
    await expect(readiness).toContainText('Computation in progress');
    await expect
      .poll(
        async () =>
          (await (await request.get(`/api/jobs/${job.id}`)).json()).recovery?.checkpoint_saved,
        { timeout: 90_000 },
      )
      .toBe(true);
    await expect(card.locator('.run-diagnostics')).toBeVisible();
    await card.getByLabel('Diagnostic quantity', { exact: true }).selectOption('temperature_k');
    await expect(
      card.getByRole('img', { name: 'Temperature, Time (ps), K', exact: true }),
    ).toBeVisible();
    const diagnostics = await (await request.get(`/api/jobs/${job.id}/diagnostics`)).json();
    expect(diagnostics.available).toContain('temperature_k');
    expect(diagnostics.latest.temperature_k).toBeGreaterThan(0);
    await card.getByRole('button', { name: 'Stop', exact: true }).click();
    await expect
      .poll(async () => (await (await request.get(`/api/jobs/${job.id}`)).json()).status)
      .toBe('cancelled');
    await expect(card.getByRole('button', { name: 'Resume', exact: true })).toBeEnabled();
    await expect(card).toContainText('original configuration');
    const downloadPromise = page.waitForEvent('download');
    await card.getByRole('link', { name: 'Files', exact: true }).click();
    const download = await downloadPromise;
    expect(download.suggestedFilename()).toBe(`DynaMol-${job.id}.zip`);
    expect((await stat((await download.path())!)).size).toBeGreaterThan(0);
    const resumed = page.waitForResponse((response) =>
      response.url().endsWith(`/api/jobs/${job.id}/resume`),
    );
    await card.getByRole('button', { name: 'Resume', exact: true }).click();
    expect((await resumed).ok()).toBe(true);
    await expect
      .poll(async () => (await (await request.get(`/api/jobs/${job.id}`)).json()).status, {
        timeout: 90_000,
      })
      .toBe('completed');
    await expect(readiness).toContainText('Ready to run');
    await expect(page.getByRole('button', { name: 'Start simulation', exact: true })).toBeEnabled();
    await expect(card.getByRole('link', { name: 'Files', exact: true })).toBeVisible();
    await card.getByRole('button', { name: 'Open trajectory', exact: false }).click();
    await expect(page.locator('.structure-card h2')).toHaveText(job.name);
    const result = await (await request.get(`/api/jobs/${job.id}`)).json();
    expect(result.completed_steps).toBe(500);
    expect(result.restarts).toHaveLength(1);
    expect(errors).toEqual([]);
  } finally {
    await request.post(`/api/jobs/${job.id}/cancel`);
  }
});
