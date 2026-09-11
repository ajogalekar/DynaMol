import { test } from './testWorkspace';
import { expect, type Page } from '@playwright/test';

async function openStudio(page: Page) {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  return page.getByRole('dialog', { name: 'Set molecules in motion.' });
}

for (const viewport of [
  { width: 1440, height: 1000 },
  { width: 547, height: 637 },
]) {
  test(`water is a visible setup stage before run settings at ${viewport.width}px`, async ({
    page,
  }, info) => {
    await page.setViewportSize(viewport);
    const studio = await openStudio(page);
    const water = studio.getByRole('region', { name: 'Water & environment' });
    await expect(
      water.getByRole('combobox', { name: 'Solvent environment', exact: true }),
    ).toHaveValue('implicit');
    await expect(water.getByRole('button', { name: 'Build water', exact: true })).toBeVisible();
    await expect(water.getByRole('spinbutton', { name: /^Box padding/ })).toBeDisabled();
    const ordering = await studio.evaluate((element) => {
      const workbench = element.querySelector('.structure-workbench')!;
      const water = element.querySelector('.water-environment')!;
      const runName = element.querySelector('input[maxlength="100"]')!;
      return {
        afterPreparation: !!(
          workbench.compareDocumentPosition(water) & Node.DOCUMENT_POSITION_FOLLOWING
        ),
        beforeSettings: !!(
          water.compareDocumentPosition(runName) & Node.DOCUMENT_POSITION_FOLLOWING
        ),
      };
    });
    expect(ordering).toEqual({ afterPreparation: true, beforeSettings: true });
    if (viewport.width > 900) {
      await water.scrollIntoViewIfNeeded();
      await expect(water).toBeInViewport({ ratio: 0.9 });
    } else {
      for (const control of [
        water.getByRole('heading', { name: 'Water & environment' }),
        water.getByRole('combobox', { name: 'Solvent environment', exact: true }),
        water.getByRole('button', { name: 'Build water', exact: true }),
        water.getByRole('spinbutton', { name: /^Box padding/ }),
      ]) {
        await control.scrollIntoViewIfNeeded();
        await expect(control).toBeInViewport({ ratio: 0.95 });
      }
    }
    await page.screenshot({ path: info.outputPath(`water-stage-${viewport.width}.png`) });
    await studio.getByRole('button', { name: 'GROMACS', exact: true }).click();
    await expect(
      water.getByRole('combobox', { name: 'Solvent environment', exact: true }),
    ).toHaveValue('explicit');
    await expect(
      water.getByRole('combobox', { name: 'Solvent environment', exact: true }),
    ).toBeDisabled();
    await expect(water.getByRole('button', { name: 'Build water', exact: true })).toHaveCount(0);
    await expect(water).toContainText('built during native setup');
    await expect(water.getByRole('spinbutton', { name: /^Box padding/ })).toBeEnabled();
  });
}

test('tracking is optional, selects existing geometry, and stays with the scene after reopening', async ({
  page,
  request,
}) => {
  const studio = await openStudio(page);
  const tracking = studio.getByRole('region', { name: 'Track during simulation', exact: true });
  const expand = tracking.getByRole('button', { name: /^Track during simulation/ });
  await expect(expand).toHaveAttribute('aria-expanded', 'false');
  await expect(tracking.getByRole('button', { name: 'Add measurement', exact: true })).toHaveCount(
    0,
  );
  const requests: string[] = [];
  page.on('request', (request) => {
    if (request.method() === 'POST' && request.url().endsWith('/api/jobs'))
      requests.push(request.url());
  });
  await expand.click();
  const selected = tracking.getByRole('checkbox');
  await expect(selected).toHaveCount(1);
  await expect(selected).not.toBeChecked();
  await selected.check();
  await expect(expand).toContainText('1 selected');
  await expect
    .poll(async () => {
      const workspace = await (await request.get('/api/workspace')).json();
      return workspace.state.measurements[0]?.trackDuringRun;
    })
    .toBe(true);
  await expand.click();
  await expect(expand).toHaveAttribute('aria-expanded', 'false');
  await studio.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  await expect(expand).toContainText('1 selected');
  await expand.click();
  await expect(selected).toBeChecked();
  await selected.uncheck();
  await expect(expand).toContainText('Optional');
  await studio.getByRole('button', { name: 'GROMACS', exact: true }).click();
  await expect(tracking).toContainText('Track heavy-atom geometry here');
  expect(requests).toEqual([]);
});

test('embedded measurements track distance and angle without submitting the simulation and survive mode switches', async ({
  page,
  request,
}, info) => {
  const demo = await (await request.get('/api/datasets/demo')).json();
  const atoms = demo.atoms.filter((atom: { name: string }) => atom.name === 'CA').slice(0, 3);
  const reset = await request.post('/api/workspace', {
    data: {
      state: { version: 1, dataset_id: demo.id, measurements: [], active_measurement: null },
    },
  });
  expect(reset.ok()).toBe(true);
  const studio = await openStudio(page);
  const tracking = studio.getByRole('region', { name: 'Track during simulation', exact: true });
  await tracking.getByRole('button', { name: /^Track during simulation/ }).click();
  let starts = 0;
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  page.on('request', (request) => {
    if (request.method() === 'POST' && request.url().endsWith('/api/jobs')) starts++;
  });
  for (const [kind, count] of [
    ['distance', 2],
    ['angle', 3],
  ] as const) {
    await tracking.getByRole('button', { name: 'Add measurement', exact: true }).click();
    const editor = tracking.locator('.tracking-editor');
    await expect(editor).toBeVisible();
    await editor
      .locator('.measure-types')
      .getByRole('button', { name: new RegExp(`${kind}$`, 'i') })
      .click();
    const search = editor.getByRole('textbox', { name: 'Find an atom', exact: true });
    await search.fill('nothing-matches');
    await search.press('Enter');
    await expect(
      studio.getByRole('button', { name: 'Start simulation', exact: true }),
    ).toBeDisabled();
    const response = page.waitForResponse(
      (res) =>
        res.url().endsWith('/measurements') &&
        res.request().method() === 'POST' &&
        res.request().postDataJSON()?.kind === kind,
    );
    for (const atom of atoms.slice(0, count)) {
      await search.fill(`${atom.residue}${atom.resid} · ${atom.name}`);
      await editor
        .locator('.atom-search-results')
        .getByRole('button', { name: new RegExp(`#${atom.index + 1} · ${atom.element}$`) })
        .click();
    }
    const measured = await (await response).json();
    await expect(editor).toHaveCount(0);
    const row = tracking
      .locator('.tracking-measurement')
      .filter({ has: page.locator('small', { hasText: kind }) });
    await expect(row.getByRole('checkbox')).toBeChecked();
    await expect(page.locator('.plot-readout strong')).toHaveText(
      `${measured.values[0].toFixed(2)}${measured.unit}`,
    );
  }
  await expect(page.locator('.measurement-tab')).toHaveCount(2);
  await expect(tracking.locator('.tracking-measurement input:checked')).toHaveCount(2);
  const plotBefore = await page.locator('.plot-readout strong').innerText();
  await page.getByRole('button', { name: 'Explore', exact: true }).click();
  await expect(studio).toHaveCount(0);
  await expect(page.locator('.measurement-tab')).toHaveCount(2);
  await expect(page.locator('.plot-readout strong')).toHaveText(plotBefore);
  await page.getByRole('button', { name: 'Simulate', exact: true }).click();
  await tracking.getByRole('button', { name: /^Track during simulation/ }).click();
  await expect(tracking.locator('.tracking-measurement input:checked')).toHaveCount(2);
  await expect
    .poll(async () => {
      const workspace = await (await request.get('/api/workspace')).json();
      return workspace.state.measurements.filter(
        (measurement: { trackDuringRun?: boolean }) => measurement.trackDuringRun,
      ).length;
    })
    .toBe(2);
  await page.screenshot({ path: info.outputPath('embedded-simulation-measurements.png') });
  await page.reload();
  await expect(page.locator('.measurement-tab')).toHaveCount(2);
  await expect(page.locator('.plot-readout strong')).toHaveText(plotBefore);
  await page.getByRole('button', { name: 'Simulate', exact: true }).click();
  await tracking.getByRole('button', { name: /^Track during simulation/ }).click();
  await expect(tracking.locator('.tracking-measurement input:checked')).toHaveCount(2);
  expect(starts).toBe(0);
  expect(errors).toEqual([]);
});
