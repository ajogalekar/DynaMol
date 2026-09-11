import { expect } from '@playwright/test';
import { test } from './testWorkspace';
import { readFile, writeFile } from 'node:fs/promises';

test('real RMSD/RMSF plots, frame seeking, atom selection and CSV agree with saved-coordinate analysis', async ({
  page,
  request,
}, info) => {
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.locator('.trajectory-analysis > summary').click();
  const pane = page.locator('.trajectory-analysis');
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/structural-analysis') && r.request().method() === 'POST',
  );
  await pane.getByRole('button', { name: 'Calculate RMSD', exact: true }).click();
  const rmsdResponse = await response;
  expect(rmsdResponse.ok(), await rmsdResponse.text()).toBeTruthy();
  const rmsd = await rmsdResponse.json();
  expect(rmsd.values).toHaveLength(101);
  expect(rmsd.values[0]).toBeLessThan(1e-5);
  expect(rmsd.values.every(Number.isFinite)).toBeTruthy();
  await expect(pane.locator('.analysis-summary')).toContainText(
    `Mean ${rmsd.summary.mean.toFixed(3)} Å`,
  );
  const chart = pane.getByRole('img', { name: 'RMSD, Time (ps), Å', exact: true });
  await chart.scrollIntoViewIfNeeded();
  const bounds = await chart.boundingBox();
  await chart.click({
    position: { x: (bounds!.width * (65 + 715 * 0.75)) / 800, y: bounds!.height * 0.4 },
  });
  await expect(page.getByRole('slider', { name: 'Trajectory frame', exact: true })).toHaveValue(
    '75',
  );
  const downloadEvent = page.waitForEvent('download');
  await pane.getByRole('button', { name: 'Export CSV', exact: true }).click();
  const csvPath = info.outputPath('rmsd.csv');
  await (await downloadEvent).saveAs(csvPath);
  const lines = (await readFile(csvPath, 'utf8')).split('\n');
  expect(lines).toHaveLength(102);
  expect(lines[0]).toBe('"Time (ps)","label","rmsd_angstrom"');
  expect(Number(lines[76].split(',').at(-1)!.replaceAll('"', ''))).toBe(rmsd.values[75]);
  const recordDownload = page.waitForEvent('download');
  await pane.getByRole('button', { name: 'Analysis record', exact: true }).click();
  const recordPath = info.outputPath('rmsd-record.json');
  await (await recordDownload).saveAs(recordPath);
  const record = JSON.parse(await readFile(recordPath, 'utf8'));
  expect(record.analysis).toEqual(rmsd);
  expect(record.dataset.id).toBe('demo');
  expect(await pane.evaluate((el) => el.clientHeight >= el.scrollHeight - 1)).toBeTruthy();

  await pane
    .getByRole('combobox', { name: 'Structural analysis kind', exact: true })
    .selectOption('rmsf');
  await pane.getByRole('spinbutton', { name: 'Analysis first frame', exact: true }).fill('11');
  await pane.getByRole('spinbutton', { name: 'Analysis last frame', exact: true }).fill('91');
  await pane.getByRole('spinbutton', { name: 'Analysis stride', exact: true }).fill('2');
  const rmsfResponse = page.waitForResponse((r) => r.url().endsWith('/structural-analysis'));
  await pane.getByRole('button', { name: 'Calculate RMSF', exact: true }).click();
  const rmsf = await (await rmsfResponse).json();
  expect(rmsf.summary.frames).toBe(41);
  expect(rmsf.values).toHaveLength(76);
  await expect(pane.locator('.analysis-summary')).toContainText(
    `Mean ${rmsf.summary.mean.toFixed(3)} Å`,
  );
  const fluctuationChart = pane.getByRole('img', {
    name: 'RMSF, Selected residue, Å',
    exact: true,
  });
  await fluctuationChart.scrollIntoViewIfNeeded();
  const fluctuationBounds = await fluctuationChart.boundingBox();
  await fluctuationChart.click({ position: { x: (fluctuationBounds!.width * 65) / 800, y: 30 } });
  await expect
    .poll(async () => (await (await request.get('/api/workspace')).json()).state?.selected_atoms)
    .toEqual(rmsf.atom_groups[0]);
  await pane.screenshot({
    path: info.outputPath('structural-analysis-desktop.png'),
  });
  await writeFile(
    info.outputPath('structural-analysis.json'),
    JSON.stringify({ rmsd, rmsf, errors }, null, 2),
  );
  expect(errors).toEqual([]);
});

test('named ligand-sized selections support unfitted analysis and explain undefined fits', async ({
  page,
  request,
}) => {
  const saved = await request.post('/api/workspace', {
    data: {
      state: {
        version: 1,
        dataset_id: 'demo',
        named_selections: [{ id: 'pair', name: 'Local pair', atoms: [0, 1] }],
      },
    },
  });
  expect(saved.ok(), await saved.text()).toBeTruthy();
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.locator('.trajectory-analysis > summary').click();
  const pane = page.locator('.trajectory-analysis');
  await pane
    .getByRole('combobox', { name: 'Analysis atom group', exact: true })
    .selectOption({ label: 'Local pair (2 atoms)' });
  await pane
    .getByRole('combobox', { name: 'Analysis alignment', exact: true })
    .selectOption('selection');
  await pane.getByRole('button', { name: 'Calculate RMSD', exact: true }).click();
  await expect(pane.getByRole('alert')).toContainText('at least three non-collinear atoms');
  await pane
    .getByRole('combobox', { name: 'Analysis alignment', exact: true })
    .selectOption('none');
  await pane.getByRole('button', { name: 'Calculate RMSD', exact: true }).click();
  await expect(pane.locator('.analysis-summary')).toContainText('2 measured atoms');
  await expect(pane).toContainText('rigid translation and rotation');
  await pane.getByRole('spinbutton', { name: 'Analysis first frame', exact: true }).fill('20');
  await pane.getByRole('spinbutton', { name: 'Analysis last frame', exact: true }).fill('5');
  await pane.getByRole('button', { name: 'Calculate RMSD', exact: true }).click();
  await expect(pane.getByRole('alert')).toContainText('frame window inside');
});

test('readiness exposes resource limits before submission and recovers when settings are corrected', async ({
  page,
  request,
}, info) => {
  const before = await (await request.get('/api/jobs')).json();
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  const preparation = dialog.getByLabel('Preparation readiness', { exact: true });
  await expect(preparation).toContainText('Ready to prepare');
  await expect(dialog.getByRole('button', { name: 'Prep protein', exact: true })).toBeEnabled();
  await dialog.getByRole('button', { name: /GROMACS/ }).click();
  await dialog.getByRole('spinbutton', { name: /^Duration/ }).fill('18');
  await dialog.getByRole('button', { name: /Advanced controls/ }).click();
  await dialog.getByRole('spinbutton', { name: /^Save every/ }).fill('1');
  const readiness = dialog.getByLabel('Simulation readiness', { exact: true });
  await expect(readiness).toContainText('256 MiB coordinate limit');
  await expect(
    dialog.getByRole('button', { name: 'Start simulation', exact: true }),
  ).toBeDisabled();
  await dialog.getByRole('spinbutton', { name: /^Duration/ }).fill('0.01');
  await expect(readiness).toContainText('Ready to run');
  await expect(dialog.getByRole('button', { name: 'Start simulation', exact: true })).toBeEnabled();
  await readiness.locator('summary').click();
  await expect(readiness).toContainText('fixed-volume NVT');
  const rechecked = page.waitForResponse(
    (r) => r.url().endsWith('/readiness') && r.request().postDataJSON().mode === 'simulation',
  );
  await readiness.getByRole('button', { name: 'Recheck readiness', exact: true }).click();
  expect((await (await rechecked).json()).ready).toBeTruthy();
  await expect(readiness).toContainText('Ready to run');
  expect(await (await request.get('/api/jobs')).json()).toEqual(before);
  await readiness.screenshot({ path: info.outputPath('simulation-readiness.png') });
});

test('analysis remains usable in the narrow layout', async ({ page }, info) => {
  await page.setViewportSize({ width: 430, height: 932 });
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await page.locator('.trajectory-analysis > summary').click();
  await page.getByRole('button', { name: 'Calculate RMSD', exact: true }).click();
  await expect(page.locator('.analysis-summary')).toBeVisible();
  expect(
    await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth),
  ).toBeTruthy();
  await page
    .locator('.trajectory-analysis')
    .screenshot({ path: info.outputPath('structural-analysis-narrow.png') });
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  await expect(page.locator('.trajectory-analysis')).toBeHidden();
  await expect(page.locator('.molecular-viewer__canvas canvas')).toBeVisible();
});
