import { test } from './testWorkspace';
import { expect, type Page, type APIRequestContext } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';

async function openStudio(page: Page, request: APIRequestContext) {
  const jobs = await (await request.get('/api/jobs')).json();
  expect(
    jobs.filter(
      (j: { status: string }) =>
        !['completed', 'failed', 'cancelled', 'interrupted'].includes(j.status),
    ),
    'Native audit must have exclusive ownership of the job queue.',
  ).toEqual([]);
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  return page.getByRole('dialog', { name: 'Set molecules in motion.' });
}

test('audit native: short GROMACS run, log show/hide, output ZIP and completed trajectory', async ({
  page,
  request,
}, info) => {
  test.setTimeout(180_000);
  const dialog = await openStudio(page, request);
  await dialog.getByRole('button', { name: /GROMACS/ }).click();
  await expect(dialog.getByRole('button', { name: /GROMACS/ })).toContainText('Installed');
  const name = `e2e-ui-audit-gromacs-${Date.now()}`;
  await dialog.getByLabel('Simulation name', { exact: true }).fill(name);
  await dialog.getByRole('spinbutton', { name: /^Duration/ }).fill('0.004');
  await dialog.getByRole('button', { name: /Advanced controls/ }).click();
  await dialog.getByRole('spinbutton', { name: /^Equilibration/ }).fill('0');
  await dialog.getByRole('spinbutton', { name: /^Save every/ }).fill('1');
  // Keep the normal minimization step for the newly solvated GROMACS system.
  await expect(
    dialog.getByLabel('Minimize energy before equilibration', { exact: true }),
  ).toBeChecked();
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/jobs') && r.request().method() === 'POST',
  );
  await dialog.getByRole('button', { name: 'Start simulation', exact: true }).click();
  const network = await response;
  expect(network.status(), await network.text()).toBe(201);
  const job = await network.json();
  const card = dialog.locator('.job-card').filter({ hasText: name });
  await expect(card).toBeVisible();
  await expect(card.locator('.job-log')).toBeVisible();
  await card.getByRole('button', { name: 'Hide log', exact: true }).click();
  await expect(card.locator('.job-log')).toHaveCount(0);
  await card.getByRole('button', { name: 'View log', exact: true }).click();
  await expect(card.locator('.job-log')).toBeVisible();
  await expect
    .poll(async () => (await (await request.get(`/api/jobs/${job.id}`)).json()).status, {
      timeout: 150_000,
    })
    .toMatch(/completed|failed|cancelled|interrupted/);
  await expect(card.locator('.job-status')).toHaveText('completed');
  await expect(card.locator('.job-progress b')).toHaveText('100%');
  const downloadEvent = page.waitForEvent('download');
  await card.getByRole('link', { name: 'Files', exact: true }).click();
  const download = await downloadEvent;
  const zipPath = info.outputPath('gromacs-output.zip');
  await download.saveAs(zipPath);
  expect((await readFile(zipPath)).subarray(0, 2).toString()).toBe('PK');
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
  expect(entries.some((entry: string) => entry.endsWith('.xtc'))).toBe(true);
  expect(entries).toContain('config.json');
  expect(entries).toContain('status.json');
  await info.attach('gromacs-output-files.json', {
    body: JSON.stringify(entries, null, 2),
    contentType: 'application/json',
  });
  await info.attach('gromacs-job.json', {
    body: await (await request.get(`/api/jobs/${job.id}`)).text(),
    contentType: 'application/json',
  });
  await card.getByRole('button', { name: 'Open trajectory', exact: false }).click();
  await expect(page.locator('.structure-card h2')).toHaveText(name);
  await expect(page.getByRole('button', { name: 'Play trajectory', exact: true })).toBeEnabled();
  await expect(page.getByRole('button', { name: 'Hide water', exact: true })).toHaveAttribute(
    'aria-pressed',
    'true',
  );
});

test('audit native: running OpenMM monitor reopens from footer and Stop cancels the worker', async ({
  page,
  request,
}, info) => {
  test.setTimeout(90_000);
  const dialog = await openStudio(page, request);
  const name = `e2e-ui-audit-cancel-${Date.now()}`;
  await dialog.getByLabel('Simulation name', { exact: true }).fill(name);
  await dialog.getByRole('spinbutton', { name: /^Duration/ }).fill('10');
  await dialog.getByRole('button', { name: /Advanced controls/ }).click();
  await dialog.getByRole('spinbutton', { name: /^Equilibration/ }).fill('0');
  await dialog.getByLabel('Minimize energy before equilibration', { exact: true }).uncheck();
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/jobs') && r.request().method() === 'POST',
  );
  await dialog.getByRole('button', { name: 'Start simulation', exact: true }).click();
  const network = await response;
  expect(network.status(), await network.text()).toBe(201);
  const job = await network.json();
  try {
    await dialog.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
    await page.getByRole('button', { name: /background job.*running · view progress/ }).click();
    const card = dialog.locator('.job-card').filter({ hasText: name });
    await expect(card.getByRole('button', { name: 'Stop', exact: true })).toBeVisible();
    await expect(
      dialog.getByRole('button', { name: 'Start simulation', exact: true }),
    ).toBeDisabled();
    await card.getByRole('button', { name: 'Stop', exact: true }).click();
    await expect(card.locator('.job-status')).toHaveText('cancelled');
    await expect(card.getByRole('button', { name: 'Stop', exact: true })).toHaveCount(0);
    await expect(
      dialog.getByRole('button', { name: 'Start simulation', exact: true }),
    ).toBeEnabled();
    await expect(
      page.getByRole('button', { name: /background job.*running · view progress/ }),
    ).toHaveCount(0);
    await info.attach('cancelled-job.json', {
      body: await (await request.get(`/api/jobs/${job.id}`)).text(),
      contentType: 'application/json',
    });
  } finally {
    const current = await (await request.get(`/api/jobs/${job.id}`)).json();
    if (!['completed', 'failed', 'cancelled', 'interrupted'].includes(current.status))
      await request.post(`/api/jobs/${job.id}/cancel`);
  }
});

test('audit native: preparation monitor Cancel stops a real preparation job', async ({
  page,
  request,
}, info) => {
  const dialog = await openStudio(page, request);
  const sourceName = `e2e-ui-audit-prep-cancel-${Date.now()}`;
  const topology = await (await request.get('/api/datasets/demo/topology')).body();
  await dialog
    .getByLabel('Upload simulation structure', { exact: true })
    .setInputFiles({ name: `${sourceName}.pdb`, mimeType: 'chemical/x-pdb', buffer: topology });
  await expect(dialog.locator('.source-summary strong')).toHaveText(sourceName);
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(dialog.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/preparations') && r.request().method() === 'POST',
  );
  await dialog.getByRole('button', { name: 'Prep protein', exact: true }).click();
  const network = await response;
  expect(network.status(), await network.text()).toBe(201);
  const job = await network.json();
  try {
    const monitor = dialog.getByRole('region', {
      name: 'Protein preparation monitor',
      exact: true,
    });
    await monitor.getByRole('button', { name: 'Cancel protein preparation', exact: true }).click();
    await expect
      .poll(async () => (await (await request.get(`/api/jobs/${job.id}`)).json()).status)
      .toBe('cancelled');
    await expect(monitor).toContainText(/cancelled/i);
    await expect(
      monitor.getByRole('button', { name: 'Cancel protein preparation', exact: true }),
    ).toHaveCount(0);
    await monitor
      .getByRole('button', { name: 'View protein preparation log', exact: true })
      .click();
    await expect(
      dialog.locator('.job-card').filter({ hasText: job.name }).first().locator('.job-log'),
    ).toBeVisible();
    await monitor.getByRole('button', { name: 'Dismiss preparation status', exact: true }).click();
    await expect(monitor).toHaveCount(0);
    await info.attach('cancelled-preparation.json', {
      body: await (await request.get(`/api/jobs/${job.id}`)).text(),
      contentType: 'application/json',
    });
  } finally {
    const current = await (await request.get(`/api/jobs/${job.id}`)).json();
    if (!['completed', 'failed', 'cancelled', 'interrupted'].includes(current.status))
      await request.post(`/api/jobs/${job.id}/cancel`);
  }
});

test('audit native: prepared PDB download, water update, water log and setup cancellation', async ({
  page,
  request,
}, info) => {
  test.setTimeout(120_000);
  await openStudio(page, request);
  // Establish a fresh, real prepared fixture. The prep button lifecycle itself is
  // covered in preparation-jobs.spec.ts; this case targets the result controls.
  const name = `e2e-ui-audit-water-controls-${Date.now()}`;
  const fixture = fileURLToPath(
    new URL('../../docs/audit/preparation-fixtures/six_residues_missing_atom.pdb', import.meta.url),
  );
  const imported = await request.post('/api/structures/upload', {
    multipart: {
      file: {
        name: `${name}-source.pdb`,
        mimeType: 'chemical/x-pdb',
        buffer: await readFile(fixture),
      },
    },
  });
  expect(imported.ok(), await imported.text()).toBe(true);
  const source = await imported.json();
  const preparedResponse = await request.post('/api/preparations', {
    data: { dataset_id: source.id, name, ph: 7, optimize_sidechains: false, seed: 2026 },
  });
  expect(preparedResponse.status(), await preparedResponse.text()).toBe(201);
  const job = await preparedResponse.json();
  await expect
    .poll(async () => (await (await request.get(`/api/jobs/${job.id}`)).json()).status, {
      timeout: 60_000,
    })
    .toMatch(/completed|failed|cancelled|interrupted/);
  const completed = await (await request.get(`/api/jobs/${job.id}`)).json();
  expect(completed.status, completed.error).toBe('completed');
  const prepared = await (await request.get(`/api/datasets/${completed.dataset_id}`)).json();
  await page.reload();
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByTitle('Switch dataset', { exact: true }).click();
  await page
    .locator('.workspace-library__row')
    .filter({ hasText: name })
    .first()
    .getByRole('button', { name: 'Open', exact: true })
    .click();
  await expect(page.locator('.structure-card h2')).toHaveText(name);
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  const downloaded = page.waitForEvent('download');
  await dialog.getByRole('link', { name: 'Prepared PDB', exact: true }).click();
  const pdb = await readFile((await (await downloaded).path())!, 'utf8');
  expect(pdb.split('\n').filter((line) => /^(ATOM  |HETATM)/.test(line))).toHaveLength(
    prepared.n_atoms,
  );
  await info.attach('prepared-download.pdb', { body: pdb, contentType: 'chemical/x-pdb' });
  const waterResponse = page.waitForResponse(
    (r) => r.url().endsWith('/solvate') && r.request().method() === 'POST',
  );
  await dialog
    .getByRole('combobox', { name: 'Solvent environment', exact: true })
    .selectOption('explicit');
  const waterJob = await (await waterResponse).json();
  await expect
    .poll(async () => (await (await request.get(`/api/jobs/${waterJob.id}`)).json()).status, {
      timeout: 60_000,
    })
    .toBe('completed');
  await expect(dialog.locator('.solvent-preview-card')).toContainText(
    'Explicit water is in the view',
  );
  let monitor = dialog.getByRole('region', { name: 'Explicit-water setup monitor', exact: true });
  await monitor.getByRole('button', { name: 'View explicit-water setup log', exact: true }).click();
  await expect(
    dialog.locator('.job-card').filter({ hasText: waterJob.name }).first().locator('.job-log'),
  ).toBeVisible();
  await monitor.getByRole('button', { name: 'Dismiss explicit-water status', exact: true }).click();
  await dialog.getByRole('button', { name: /Advanced controls/ }).click();
  await dialog.getByRole('spinbutton', { name: /^Box padding/ }).fill('1.1');
  const updateResponse = page.waitForResponse(
    (r) => r.url().endsWith('/solvate') && r.request().method() === 'POST',
  );
  await dialog.getByRole('button', { name: 'Update box', exact: true }).click();
  const updated = await (await updateResponse).json();
  expect(updated.config.padding_nm).toBe(1.1);
  await expect
    .poll(async () => (await (await request.get(`/api/jobs/${updated.id}`)).json()).status, {
      timeout: 60_000,
    })
    .toBe('completed');
  await expect(dialog.locator('.solvent-preview-card')).toContainText('1.1 nm padding');
  await dialog.getByRole('spinbutton', { name: /^Box padding/ }).fill('1.2');
  const cancelResponse = page.waitForResponse(
    (r) => r.url().endsWith('/solvate') && r.request().method() === 'POST',
  );
  await dialog.getByRole('button', { name: 'Update box', exact: true }).click();
  const cancelled = await (await cancelResponse).json();
  try {
    monitor = dialog.getByRole('region', { name: 'Explicit-water setup monitor', exact: true });
    await monitor.getByRole('button', { name: 'Cancel explicit-water setup', exact: true }).click();
    await expect
      .poll(async () => (await (await request.get(`/api/jobs/${cancelled.id}`)).json()).status)
      .toBe('cancelled');
    await expect(monitor).toContainText(/cancelled/i);
    await expect(dialog.locator('.solvent-preview-card')).toContainText('1.1 nm padding');
    await monitor
      .getByRole('button', { name: 'Dismiss explicit-water status', exact: true })
      .click();
    await info.attach('water-controls.json', {
      body: JSON.stringify(
        {
          source_dataset_id: source.id,
          prepared_job_id: job.id,
          initial_water_job_id: waterJob.id,
          updated_water_job_id: updated.id,
          cancelled_water_job_id: cancelled.id,
        },
        null,
        2,
      ),
      contentType: 'application/json',
    });
  } finally {
    const current = await (await request.get(`/api/jobs/${cancelled.id}`)).json();
    if (!['completed', 'failed', 'cancelled', 'interrupted'].includes(current.status))
      await request.post(`/api/jobs/${cancelled.id}/cancel`);
  }
});
