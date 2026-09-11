import { test } from './testWorkspace';
import { expect, type Page, type APIRequestContext } from '@playwright/test';
// Real native output, geometry, topology and persistence; only event delivery timing is replayed.
const fixtureId = process.env.DYNAMOL_TRACKING_JOB;
async function replay(page: Page, request: APIRequestContext) {
  const response = await request.get(`/api/jobs/${fixtureId}`);
  expect(response.ok(), 'Seed a verified native tracking result in this isolated API.').toBe(true);
  const job = await response.json();
  const source = await (await request.get(`/api/datasets/${job.config.dataset_id}`)).json();
  const output = await (await request.get(`/api/datasets/${job.dataset_id}`)).json();
  const snapshot = await (await request.get(`/api/jobs/${job.id}/measurements`)).json();
  const measurements = [];
  for (const definition of job.config.measurements) {
    const geometry = await (
      await request.post(`/api/datasets/${source.id}/measurement-preview`, { data: definition })
    ).json();
    measurements.push({ ...definition, ...geometry, trackDuringRun: true, visible: true });
  }
  expect(
    (
      await request.post('/api/workspace', {
        data: {
          state: {
            version: 1,
            dataset_id: source.id,
            measurements,
            active_measurement: measurements[0].id,
            camera: [80, 0, 0, 0, 0, 80, 0, 0, 0, 0, 80, 0, -20, -20, -20, 1],
          },
        },
      })
    ).ok(),
  ).toBe(true);
  let complete = false;
  await page.route('**/api/jobs', async (route) => {
    const response = await route.fetch();
    await route.fulfill({
      response,
      json: (await response.json()).map((j: { id: string }) =>
        j.id === job.id && !complete
          ? { ...j, status: 'running', stage: 'Running dynamics', progress: 30, dataset_id: null }
          : j,
      ),
    });
  });
  await page.route(`**/api/jobs/${job.id}/measurements`, async (route) => {
    const response = await route.fetch(),
      real = await response.json();
    await route.fulfill({
      response,
      json: complete
        ? real
        : {
            ...real,
            status: 'running',
            output_dataset_id: null,
            measurements: real.measurements.map(
              (m: { values: number[]; times_ps: number[]; frame_errors: unknown[] }) => ({
                ...m,
                values: m.values.slice(0, 5),
                times_ps: m.times_ps.slice(0, 5),
                frame_errors: m.frame_errors.slice(0, 5),
              }),
            ),
          },
    });
  });
  await page.addInitScript(
    ({ id, source }) => {
      if (!sessionStorage.getItem('tracking-fixture-initialized')) {
        localStorage.setItem(
          'dynamol.followed-run.v1',
          JSON.stringify({ id, source, automatic: true }),
        );
        sessionStorage.setItem('tracking-fixture-initialized', '1');
      }
    },
    { id: job.id, source: source.id },
  );
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.plot-heading')).toContainText('Plots update with saved frames');
  await expect(page.locator('.measurement-tab')).toHaveCount(3);
  return {
    job,
    source,
    output,
    snapshot,
    complete: () => {
      complete = true;
    },
  };
}

test.describe('recorded native run handoff', () => {
  test.skip(
    !fixtureId,
    'Set DYNAMOL_TRACKING_JOB to a verified native result in the isolated API.',
  );
  test('live plots stay across modes and complete into Explore with camera, indices and values preserved', async ({
    page,
    request,
  }, info) => {
    const errors: string[] = [];
    page.on('pageerror', (e) => errors.push(e.message));
    const f = await replay(page, request);
    await expect(page.locator('.plot-readout strong')).toHaveText(
      `${f.snapshot.measurements[0].values[4].toFixed(2)}Å`,
    );
    const simulate = page
      .getByRole('navigation', { name: 'Workspace', exact: true })
      .getByRole('button', { name: /^Simulate/ });
    await simulate.click();
    await expect(page.locator('.measurement-tab')).toHaveCount(3);
    await page.getByRole('button', { name: 'Explore', exact: true }).click();
    await expect(page.locator('.plot-heading')).toContainText('Plots update with saved frames');
    await simulate.click();
    f.complete();
    await expect(page.locator('.structure-card h2')).toHaveText(f.output.name);
    await expect(page.getByRole('button', { name: 'Explore', exact: true })).toHaveClass('active');
    await expect(page.locator('.plot-heading')).not.toContainText('Plots update with saved frames');
    await expect(page.locator('.measurement-tab')).toHaveCount(3);
    await expect
      .poll(async () => (await (await request.get('/api/workspace')).json()).state.dataset_id)
      .toBe(f.output.id);
    const saved = (await (await request.get('/api/workspace')).json()).state;
    expect(saved.camera.slice(0, 12)).toEqual([80, 0, 0, 0, 0, 80, 0, 0, 0, 0, 80, 0]);
    expect(saved.frame).toBe(f.output.n_frames - 1);
    for (const m of saved.measurements) {
      const native = f.snapshot.measurements.find((n: { id: string }) => n.id === m.id);
      expect(m.atoms).toEqual(native.output_atoms);
      expect(m.values).toEqual(native.values);
      expect(m.times_ps).toEqual(f.output.times_ps);
      expect(m.trackDuringRun).toBe(true);
    }
    await page.screenshot({ path: info.outputPath('completed-run-explore.png'), fullPage: true });
    await page.reload();
    await expect(page.locator('.measurement-tab')).toHaveCount(3);
    expect(errors).toEqual([]);
  });
  test('another scene keeps focus; explicitly opening the completed run restores its plots', async ({
    page,
    request,
  }) => {
    const f = await replay(page, request);
    await page.getByTitle('Switch dataset', { exact: true }).click();
    await page
      .locator('.workspace-library__row')
      .filter({ hasText: 'Ubiquitin · starting structure' })
      .getByRole('button', { name: 'Open', exact: true })
      .click();
    await expect(page.locator('.structure-card h2')).toHaveText('Ubiquitin · starting structure');
    await expect.poll(async () => (await (await request.get('/api/workspace')).json()).state.dataset_id).toBe('ubiquitin-start');
    const before = (await (await request.get('/api/workspace')).json()).state;
    f.complete();
    const notice = page.locator('.run-result-notice').filter({ hasText: 'is ready to explore' });
    await expect(notice).toBeVisible();
    await expect(page.locator('.structure-card h2')).toHaveText('Ubiquitin · starting structure');
    const after = (await (await request.get('/api/workspace')).json()).state;
    expect(after.camera).toEqual(before.camera);
    expect(after.measurements).toEqual(before.measurements);
    expect(after.visibility).toEqual(before.visibility);
    await notice.getByRole('button', { name: 'Open trajectory', exact: true }).click();
    await expect(page.locator('.structure-card h2')).toHaveText(f.output.name);
    await expect(page.locator('.measurement-tab')).toHaveCount(3);
    await expect(page.getByRole('button', { name: 'Explore', exact: true })).toHaveClass('active');
  });
  test('New simulation cancels automatic handoff and releases the old live plot', async ({
    page,
    request,
  }) => {
    const f = await replay(page, request);
    await page.getByRole('button', { name: 'New simulation', exact: true }).click();
    await expect(page.locator('.plot-heading')).not.toContainText('Plots update with saved frames');
    f.complete();
    await expect(
      page.locator('.run-result-notice').filter({ hasText: 'is ready to explore' }),
    ).toBeVisible();
    await expect(page.locator('.structure-card h2')).toHaveText(f.source.name);
    await expect(page.getByRole('dialog', { name: 'Set molecules in motion.' })).toBeVisible();
  });
  test('a completion already loading cannot overwrite a new setup', async ({ page, request }) => {
    const f = await replay(page, request);
    let release!: () => void;
    let pending = false;
    const gate = new Promise<void>((resolve) => { release = resolve; });
    await page.route(`**${f.output.coordinates_url}`, async (route) => {
      const response = await route.fetch();
      pending = true;
      await gate;
      await route.fulfill({ response });
    });
    f.complete();
    await expect.poll(() => pending).toBe(true);
    await page.getByRole('button', { name: 'New simulation', exact: true }).click();
    release();
    await expect(page.locator('.run-result-notice').filter({ hasText: 'is ready to explore' })).toBeVisible();
    await expect(page.locator('.structure-card h2')).toHaveText(f.source.name);
    await expect(page.getByRole('dialog', { name: 'Set molecules in motion.' })).toBeVisible();
  });

});
