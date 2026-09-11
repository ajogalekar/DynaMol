import { test, expect, type Page, type APIRequestContext, type TestInfo } from '@playwright/test';
import { readFile } from 'node:fs/promises';

type Atom = { index: number; name: string; residue: string; resid: number; element: string };
type Dataset = {
  id: string;
  name: string;
  n_frames: number;
  n_atoms: number;
  atoms: Atom[];
  bonds: [number, number][];
};

async function realDemo(request: APIRequestContext): Promise<Dataset> {
  const response = await request.get('/api/datasets/demo');
  expect(response.ok(), 'Generate the bundled real demo before running E2E tests.').toBeTruthy();
  return response.json();
}

async function openReady(page: Page) {
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await expect(page.getByRole('alert')).toHaveCount(0);
  return errors;
}

async function pickAtom(page: Page, atom: Atom) {
  await page
    .getByRole('textbox', { name: 'Find an atom', exact: true })
    .fill(`${atom.residue}${atom.resid} · ${atom.name}`);
  await page
    .locator('.atom-search-results')
    .getByRole('button', { name: new RegExp(`#${atom.index + 1} · ${atom.element}$`) })
    .click();
}

async function chooseDataset(page: Page, name: string) {
  await page.getByTitle('Switch dataset', { exact: true }).click();
  await page
    .locator('.library-popover')
    .getByRole('button')
    .filter({ hasText: name })
    .first()
    .click();
  await expect(page.locator('.structure-card h2')).toHaveText(name);
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
}

async function assertRenderedScene(page: Page, testInfo: TestInfo, name: string) {
  await expect(page.locator('.molecular-viewer__canvas canvas')).toHaveCount(1);
  const canvas = page.locator('.molecular-viewer__canvas canvas');
  await expect(canvas).toBeVisible();
  const bounds = await canvas.boundingBox();
  expect(bounds?.width).toBeGreaterThan(200);
  expect(bounds?.height).toBeGreaterThan(200);
  const image = await canvas.screenshot();
  await testInfo.attach(name, { body: image, contentType: 'image/png' });
  // Inspect actual screenshot pixels, not component state: a blank WebGL canvas must fail.
  const coloredPixels = await page.evaluate(async (base64) => {
    const bitmap = await createImageBitmap(
      await (await fetch(`data:image/png;base64,${base64}`)).blob(),
    );
    const sample = document.createElement('canvas');
    sample.width = bitmap.width;
    sample.height = bitmap.height;
    const context = sample.getContext('2d')!;
    context.drawImage(bitmap, 0, 0);
    const pixels = context.getImageData(0, 0, sample.width, sample.height).data;
    let count = 0;
    for (let i = 0; i < pixels.length; i += 4) {
      const maximum = Math.max(pixels[i], pixels[i + 1], pixels[i + 2]);
      const minimum = Math.min(pixels[i], pixels[i + 1], pixels[i + 2]);
      if (maximum > 75 && maximum - minimum > 35) count++;
    }
    return count;
  }, image.toString('base64'));
  expect(
    coloredPixels,
    'The real molecular canvas must contain colored molecular geometry.',
  ).toBeGreaterThan(500);
  return image;
}

test('desktop renders molecular geometry; playback, seeking, camera and visibility respond', async ({
  page,
}, testInfo) => {
  const errors = await openReady(page);
  const original = await assertRenderedScene(page, testInfo, 'desktop-ribbons');
  const slider = page.getByRole('slider', { name: 'Trajectory frame', exact: true });
  await page.getByRole('button', { name: 'Play trajectory', exact: true }).click();
  await expect.poll(async () => Number(await slider.inputValue())).toBeGreaterThan(2);
  await page.getByRole('button', { name: 'Pause trajectory', exact: true }).click();
  const paused = await slider.inputValue();
  await page.waitForTimeout(300);
  await expect(slider).toHaveValue(paused);
  await slider.fill('50');
  await expect(slider).toHaveValue('50');
  await expect(page.locator('.time-display b')).toHaveText('1.00');
  await page.getByRole('button', { name: 'Next frame', exact: true }).click();
  await expect(slider).toHaveValue('51');
  await page.getByRole('button', { name: 'First frame', exact: true }).click();
  await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
  await page.waitForTimeout(300);
  const zoomed = await page.locator('.molecular-viewer__canvas canvas').first().screenshot();
  expect(zoomed.equals(original)).toBeFalsy();
  await page.getByRole('button', { name: 'Zoom out', exact: true }).click();
  await page.getByRole('button', { name: 'Fit molecule', exact: true }).click();
  await page.getByRole('button', { name: 'Hide protein & nucleic acids', exact: true }).click();
  await expect(
    page.getByRole('button', { name: 'Show protein & nucleic acids', exact: true }),
  ).toHaveAttribute('aria-pressed', 'false');
  await page.getByRole('button', { name: 'Show protein & nucleic acids', exact: true }).click();
  await page.getByRole('button', { name: 'Show water', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Hide water', exact: true })).toHaveAttribute(
    'aria-pressed',
    'true',
  );
  await page.getByRole('button', { name: 'Hide water', exact: true }).click();
  await page.getByRole('button', { name: 'Polar only', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Polar only', exact: true })).toHaveClass(/active/);
  await page.getByRole('button', { name: 'All', exact: true }).click();
  await page.getByRole('button', { name: 'Ball & stick', exact: true }).click();
  await page.getByRole('combobox', { name: 'Color molecules by' }).selectOption('element');
  await assertRenderedScene(page, testInfo, 'desktop-ball-stick');
  await page.getByRole('textbox', { name: 'Find a residue', exact: true }).fill('35');
  await page.getByRole('button', { name: 'Focus residue', exact: true }).click();
  await expect(page.getByRole('status')).toContainText('Focused on 35');
  const snapshot = page.waitForEvent('download');
  await page.getByRole('button', { name: 'Save snapshot', exact: true }).click();
  const download = await snapshot;
  expect(download.suggestedFilename()).toMatch(/\.png$/);
  expect((await readFile((await download.path())!)).length).toBeGreaterThan(5000);
  expect(errors).toEqual([]);
});

test('atom search produces real distance and hydrogen-bond plots with CSV export', async ({
  page,
  request,
}) => {
  const demo = await realDemo(request);
  await openReady(page);
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  const alphaCarbons = demo.atoms.filter((atom) => atom.name === 'CA');
  await pickAtom(page, alphaCarbons[0]);
  await pickAtom(page, alphaCarbons[1]);
  const distanceResponse = page.waitForResponse(
    (response) =>
      response.url().endsWith('/measurements') &&
      response.request().postDataJSON()?.atoms?.[0] === alphaCarbons[0].index,
  );
  await page.getByRole('button', { name: 'Plot over time', exact: true }).click();
  const distance = await (await distanceResponse).json();
  await expect(page.locator('.measurement-tab')).toHaveCount(2);
  await expect(page.locator('.plot-readout strong')).toContainText(distance.values[0].toFixed(2));
  const exportEvent = page.waitForEvent('download');
  await page.getByRole('button', { name: 'Export measurement CSV', exact: true }).click();
  const exported = await exportEvent;
  const lines = (await readFile((await exported.path())!, 'utf8')).trim().split('\n');
  expect(lines[0]).toBe('time_ps,distance_angstrom');
  expect(lines.length).toBe(demo.n_frames + 1);
  expect(Number(lines[1].split(',')[1])).toBeCloseTo(distance.values[0], 6);
  const donor = demo.atoms.find(
    (atom) =>
      atom.element === 'N' &&
      demo.bonds.some(
        ([a, b]) =>
          (a === atom.index && demo.atoms[b].element === 'H') ||
          (b === atom.index && demo.atoms[a].element === 'H'),
      ),
  )!;
  const bond = demo.bonds.find(
    ([a, b]) =>
      (a === donor.index && demo.atoms[b].element === 'H') ||
      (b === donor.index && demo.atoms[a].element === 'H'),
  )!;
  const hydrogen = demo.atoms[bond[0] === donor.index ? bond[1] : bond[0]];
  const acceptor = demo.atoms.find((atom) => atom.element === 'O')!;
  await page.getByRole('button', { name: /Hydrogen bond$/ }).click();
  for (const atom of [donor, hydrogen, acceptor]) await pickAtom(page, atom);
  const hbondResponse = page.waitForResponse(
    (response) =>
      response.url().endsWith('/measurements') &&
      response.request().postDataJSON()?.kind === 'hbond',
  );
  await page.getByRole('button', { name: 'Plot over time', exact: true }).click();
  const hbond = await (await hbondResponse).json();
  await expect(page.locator('.measurement-tab')).toHaveCount(3);
  await expect(page.locator('.occupancy')).toHaveText(
    `${(hbond.occupancy * 100).toFixed(1)}% geometric occupancy`,
  );
  expect(hbond.angle_values).toHaveLength(demo.n_frames);
  await expect(page.getByRole('alert')).toHaveCount(0);
});

test('import reports invalid data, then loads a genuine topology with correct single-frame behavior', async ({
  page,
  request,
}) => {
  const demo = await realDemo(request);
  const topology = await (await request.get('/api/datasets/demo/topology')).body();
  await openReady(page);
  await page.getByRole('button', { name: 'Open files', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Bring your molecules.' });
  await dialog
    .locator('input[type=file]')
    .first()
    .setInputFiles({
      name: 'invalid-e2e.pdb',
      mimeType: 'chemical/x-pdb',
      buffer: Buffer.from('This is not a molecular coordinate file.\n'),
    });
  await dialog.getByRole('button', { name: 'Open in DynaMol', exact: true }).click();
  await expect(dialog.getByRole('alert')).toBeVisible();
  const name = `e2e-structure-${Date.now()}`;
  await dialog
    .locator('input[type=file]')
    .first()
    .setInputFiles({ name: `${name}.pdb`, mimeType: 'chemical/x-pdb', buffer: topology });
  await dialog.getByRole('button', { name: 'Open in DynaMol', exact: true }).click();
  await expect(dialog).toHaveCount(0);
  await expect(page.locator('.structure-card h2')).toHaveText(name);
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.getByRole('button', { name: 'Play trajectory', exact: true })).toBeDisabled();
  await expect(page.locator('.timeline-meta')).toContainText('1 saved frames');
  await expect(page.locator('.structure-card')).toContainText(
    `${demo.n_atoms.toLocaleString()} atoms`,
  );
  await expect(page.locator('.measurement-tab')).toHaveCount(0);
  await expect(page.getByRole('alert')).toHaveCount(0);
});

test('a delayed real measurement cannot contaminate a newly selected dataset', async ({
  page,
  request,
}) => {
  const demo = await realDemo(request);
  const list = (await (await request.get('/api/datasets')).json()) as Dataset[];
  const target = list.find((dataset) => dataset.id === 'ubiquitin-start')!;
  expect(target).toBeTruthy();
  await openReady(page);
  const atoms = demo.atoms.filter((atom) => atom.name === 'CA').slice(0, 2);
  for (const atom of atoms) await pickAtom(page, atom);
  let release!: () => void, received!: () => void;
  const gate = new Promise<void>((resolve) => {
    release = resolve;
  });
  const observed = new Promise<void>((resolve) => {
    received = resolve;
  });
  // No synthetic response: hold a genuine API response to reproduce network scheduling.
  await page.route('**/api/datasets/demo/measurements', async (route) => {
    const response = await route.fetch();
    if (route.request().postDataJSON()?.atoms?.[0] === atoms[0].index) {
      received();
      await gate;
    }
    await route.fulfill({ response });
  });
  await page.getByRole('button', { name: 'Plot over time', exact: true }).click();
  await observed;
  await chooseDataset(page, target.name);
  release();
  await expect(page.getByRole('button', { name: 'Plot over time', exact: true })).not.toContainText(
    'Loading',
  );
  await page.waitForTimeout(300);
  await expect(page.locator('.measurement-tab')).toHaveCount(0);
  await expect(page.locator('.structure-card h2')).toHaveText(target.name);
});

test('a slower earlier dataset lookup cannot overwrite the latest library selection', async ({
  page,
  request,
}) => {
  const list = (await (await request.get('/api/datasets')).json()) as Dataset[];
  const older = list.find((dataset) => dataset.id === 'ubiquitin-start')!;
  const topology = await (await request.get('/api/datasets/demo/topology')).body();
  const fixture = await request.post('/api/datasets/upload', {
    multipart: {
      topology: {
        name: `e2e-race-${Date.now()}.pdb`,
        mimeType: 'chemical/x-pdb',
        buffer: topology,
      },
    },
  });
  expect(fixture.ok()).toBeTruthy();
  const latest = (await fixture.json()) as Dataset;
  expect(older).toBeTruthy();
  await openReady(page);
  let release!: () => void, received!: () => void;
  const gate = new Promise<void>((resolve) => {
    release = resolve;
  });
  const observed = new Promise<void>((resolve) => {
    received = resolve;
  });
  await page.route(`**/api/datasets/${older.id}`, async (route) => {
    const response = await route.fetch();
    received();
    await gate;
    await route.fulfill({ response });
  });
  await page.getByTitle('Switch dataset', { exact: true }).click();
  await page
    .locator('.library-popover')
    .getByRole('button')
    .filter({ hasText: older.name })
    .first()
    .click();
  await observed;
  await page
    .locator('.library-popover')
    .getByRole('button')
    .filter({ hasText: latest.name })
    .first()
    .click();
  await expect(page.locator('.structure-card h2')).toHaveText(latest.name);
  release();
  await page.waitForTimeout(400);
  await expect(page.locator('.structure-card h2')).toHaveText(latest.name);
});

test('simulation studio launches real OpenMM dynamics, reports 100%, and opens its output', async ({
  page,
  request,
}) => {
  test.setTimeout(120_000);
  const originalDemo = await realDemo(request);
  await openReady(page);
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  await expect(dialog.getByRole('button', { name: /OpenMM/ })).toContainText('Installed');
  const name = `E2E OpenMM ${Date.now()}`;
  await dialog.getByLabel('Simulation name', { exact: true }).fill(name);
  await dialog.getByLabel('Duration', { exact: false }).fill('0.02');
  await dialog.getByRole('button', { name: /Advanced controls/ }).click();
  await dialog.getByRole('spinbutton', { name: /^Equilibration/ }).fill('0');
  await dialog.getByLabel('Save every', { exact: false }).fill('5');
  await dialog.getByLabel('Minimize energy before equilibration', { exact: true }).uncheck();
  const started = page.waitForResponse(
    (response) => response.url().endsWith('/api/jobs') && response.request().method() === 'POST',
  );
  await dialog.getByRole('button', { name: 'Start simulation', exact: true }).click();
  const jobResponse = await started;
  expect(jobResponse.status()).toBe(201);
  const job = await jobResponse.json();
  expect(job.config.dataset_id).toBe('demo');
  expect(job.config.duration_ps).toBe(0.02);
  expect(job.config.minimize).toBe(false);
  const card = dialog.locator('.job-card').filter({ hasText: name });
  await expect(card).toBeVisible();
  await expect(card.locator('.job-status')).toHaveText('completed', { timeout: 90_000 });
  await expect(card.locator('.job-progress')).toContainText('10 / 10 steps');
  await expect(card.locator('.job-progress b')).toHaveText('100%');
  await expect(card.locator('.progress-track > div')).toHaveAttribute('style', /width: 100%/);
  await card.getByRole('button', { name: 'Open trajectory', exact: false }).click();
  await expect(dialog).toBeVisible();
  await expect(page.locator('.structure-card h2')).toHaveText(name);
  await expect(page.getByRole('button', { name: 'Play trajectory', exact: true })).toBeEnabled();
  await expect(page.locator('.timeline-meta')).toContainText('3 saved frames');
  const unchangedDemo = await realDemo(request);
  expect(unchangedDemo.n_frames).toBe(originalDemo.n_frames);
  expect(unchangedDemo.n_atoms).toBe(originalDemo.n_atoms);
});

test('responsive workspace retains a usable scene, playback, imports and simulation controls', async ({
  page,
}, testInfo) => {
  await page.setViewportSize({ width: 820, height: 1180 });
  const errors = await openReady(page);
  await assertRenderedScene(page, testInfo, 'tablet-scene');
  const slider = page.getByRole('slider', { name: 'Trajectory frame', exact: true });
  await page.getByRole('button', { name: 'Play trajectory', exact: true }).click();
  await expect.poll(async () => Number(await slider.inputValue())).toBeGreaterThan(1);
  await page.getByRole('button', { name: 'Pause trajectory', exact: true }).click();
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  await expect(page.getByRole('dialog', { name: 'Set molecules in motion.' })).toBeVisible();
  const sceneBounds = await page.locator('.molecular-viewer__canvas canvas').boundingBox();
  const studioBounds = await page
    .getByRole('dialog', { name: 'Set molecules in motion.' })
    .boundingBox();
  expect(sceneBounds!.height).toBeGreaterThan(150);
  expect(studioBounds!.y).toBeGreaterThan(sceneBounds!.y + sceneBounds!.height);
  expect(studioBounds!.y + studioBounds!.height).toBeLessThanOrEqual(1181);
  await page.screenshot({ path: testInfo.outputPath('tablet-studio.png') });
  await page.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
  await page.setViewportSize({ width: 390, height: 844 });
  await expect(page.getByRole('button', { name: 'Open files', exact: true })).toBeVisible();
  const overflow = await page.evaluate(
    () => document.documentElement.scrollWidth - window.innerWidth,
  );
  expect(
    overflow,
    'Mobile viewport should not force horizontal page scrolling.',
  ).toBeLessThanOrEqual(2);
  await page.getByRole('button', { name: 'Open files', exact: true }).click();
  await expect(page.getByRole('dialog', { name: 'Bring your molecules.' })).toBeVisible();
  await expect(page.getByRole('button', { name: 'Close import', exact: true })).toBeVisible();
  await page.screenshot({ path: testInfo.outputPath('phone-import.png'), fullPage: true });
  expect(errors).toEqual([]);
});
