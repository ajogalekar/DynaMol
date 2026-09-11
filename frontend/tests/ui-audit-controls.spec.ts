import { test, expect, type Page, type TestInfo } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';

type Atom = { index: number; name: string; residue: string; resid: number; element: string };
const fixture = (name: string) =>
  fileURLToPath(new URL(`../../docs/audit/preparation-fixtures/${name}`, import.meta.url));

async function ready(page: Page) {
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  return errors;
}
async function idleKeyboard(page: Page) {
  await page.locator('.structure-card h2').click();
}
async function studio(page: Page) {
  await page.getByRole('button', { name: 'Simulate', exact: true }).click();
  const dialog = page.getByRole('dialog', { name: 'Set molecules in motion.' });
  await expect(dialog).toBeVisible();
  return dialog;
}
async function pixels(page: Page, info?: TestInfo, name = 'scene') {
  await page.waitForTimeout(350);
  await page.mouse.move(3, 3);
  const image = await page.locator('.molecular-viewer__canvas canvas').screenshot();
  if (info) await info.attach(name, { body: image, contentType: 'image/png' });
  return image;
}
async function pick(page: Page, atom: Atom) {
  await page
    .getByRole('textbox', { name: 'Find an atom', exact: true })
    .fill(`${atom.residue}${atom.resid} · ${atom.name}`);
  await page
    .locator('.atom-search-results')
    .getByRole('button', { name: new RegExp(`#${atom.index + 1} · ${atom.element}$`) })
    .click();
}

test('audit: navigation, structure details, library, guide links and modal focus controls', async ({
  page,
}, info) => {
  const errors = await ready(page);
  await page.getByRole('button', { name: 'Structure details', exact: true }).click();
  await expect(page.locator('.dataset-details')).toContainText('Source:');
  await page.getByRole('button', { name: 'Structure details', exact: true }).click();
  await expect(page.locator('.dataset-details')).toHaveCount(0);
  await page.getByTitle('Switch dataset', { exact: true }).click();
  await expect(page.locator('.library-popover')).toBeVisible();
  await idleKeyboard(page);
  await page.keyboard.press('Escape');
  await expect(page.locator('.library-popover')).toHaveCount(0);
  await studio(page);
  await page.getByRole('button', { name: 'Explore', exact: true }).click();
  await expect(page.getByRole('dialog')).toHaveCount(0);
  await page.getByRole('button', { name: 'New simulation', exact: true }).click();
  await page.getByRole('button', { name: 'Close simulation studio', exact: true }).click();
  await page.getByRole('button', { name: 'Open quick guide', exact: true }).click();
  const guide = page.getByRole('dialog', { name: 'Meet DynaMol.' });
  await expect(guide).toBeVisible();
  await expect(guide.getByRole('button', { name: 'Close quick guide', exact: true })).toBeFocused();
  await page.keyboard.press('Shift+Tab');
  await expect(guide.getByRole('link', { name: 'MDTraj', exact: true })).toBeFocused();
  await page.keyboard.press('Tab');
  await expect(guide.getByRole('button', { name: 'Close quick guide', exact: true })).toBeFocused();
  const links = await guide.getByRole('link').evaluateAll((elements) =>
    elements.map((element) => ({
      text: element.textContent?.trim(),
      href: element.getAttribute('href'),
      target: element.getAttribute('target'),
    })),
  );
  expect(links).toHaveLength(4);
  expect(links.every((link) => link.href?.startsWith('https://') && link.target === '_blank')).toBe(
    true,
  );
  for (const link of links) {
    const opened = page.waitForEvent('popup');
    await guide.getByRole('link', { name: link.text!, exact: true }).click();
    const popup = await opened;
    await expect.poll(() => popup.url()).toBe(link.href);
    await popup.close();
  }
  await info.attach('guide-links.json', {
    body: JSON.stringify(links, null, 2),
    contentType: 'application/json',
  });
  await guide.getByRole('button', { name: 'Close quick guide', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Open quick guide', exact: true })).toBeFocused();
  await page.getByRole('button', { name: /Keyboard shortcuts/ }).click();
  await page.keyboard.press('Escape');
  await expect(guide).toHaveCount(0);
  await page.getByRole('button', { name: 'Meet DynaMol', exact: true }).click();
  await expect(guide).toBeVisible();
  await guide.getByRole('button', { name: 'Close quick guide', exact: true }).click();
  await page.getByRole('link', { name: 'DynaMol home', exact: true }).click();
  await expect(page.locator('.structure-card h2')).toContainText('Ubiquitin');
  expect(errors).toEqual([]);
});

test('audit: every representation and color choice changes the actual rendered geometry', async ({
  page,
}, info) => {
  const errors = await ready(page);
  const images: Buffer[] = [];
  for (const label of ['Ribbons', 'Ball & stick', 'Sticks', 'Surface']) {
    await page.getByRole('button', { name: label, exact: true }).click();
    await expect(page.getByRole('button', { name: label, exact: true })).toHaveClass(/selected/);
    images.push(await pixels(page, info, label));
  }
  for (let i = 1; i < images.length; i++)
    expect(images[i].equals(images[i - 1]), 'Representation must change the actual canvas.').toBe(
      false,
    );
  await page.getByRole('button', { name: 'Ball & stick', exact: true }).click();
  const colors: Buffer[] = [];
  for (const value of ['residue', 'chain', 'element']) {
    await page
      .getByRole('combobox', { name: 'Color molecules by', exact: true })
      .selectOption(value);
    await expect(
      page.getByRole('combobox', { name: 'Color molecules by', exact: true }),
    ).toHaveValue(value);
    colors.push(await pixels(page, info, `color-${value}`));
  }
  expect(colors[0].equals(colors[1])).toBe(false);
  expect(colors[1].equals(colors[2])).toBe(false);
  expect(errors).toEqual([]);
});

test('audit: mouse camera controls, auto rotation, focus validation and full screen', async ({
  page,
}, info) => {
  const errors = await ready(page);
  const canvas = page.locator('.molecular-viewer__canvas canvas');
  const box = (await canvas.boundingBox())!;
  const x = box.x + box.width * 0.5,
    y = box.y + box.height * 0.5;
  let before = await pixels(page);
  for (const button of ['left', 'right'] as const) {
    await page.mouse.move(x, y);
    await page.mouse.down({ button });
    await page.mouse.move(x + 50, y + 30, { steps: 8 });
    await page.mouse.up({ button });
    const after = await pixels(page, info, `camera-${button === 'left' ? 'rotate' : 'pan'}`);
    expect(after.equals(before)).toBe(false);
    before = after;
  }
  await page.mouse.move(x, y);
  await page.mouse.wheel(0, -180);
  const wheeled = await pixels(page);
  expect(wheeled.equals(before)).toBe(false);
  await page.getByRole('button', { name: 'Auto rotate', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Auto rotate', exact: true })).toHaveAttribute(
    'aria-pressed',
    'true',
  );
  const spinning = await pixels(page);
  expect((await pixels(page)).equals(spinning)).toBe(false);
  await page.getByRole('button', { name: 'Auto rotate', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Auto rotate', exact: true })).toHaveAttribute(
    'aria-pressed',
    'false',
  );
  await page.waitForTimeout(1200);
  const stopped = await pixels(page);
  const settled = await pixels(page);
  const changedPixels = await page.evaluate(
    async (images) => {
      const arrays: Uint8ClampedArray[] = [];
      for (const image of images) {
        const bitmap = await createImageBitmap(
          await (await fetch(`data:image/png;base64,${image}`)).blob(),
        );
        const canvas = document.createElement('canvas');
        canvas.width = bitmap.width;
        canvas.height = bitmap.height;
        const context = canvas.getContext('2d')!;
        context.drawImage(bitmap, 0, 0);
        arrays.push(context.getImageData(0, 0, canvas.width, canvas.height).data);
      }
      let different = 0;
      for (let i = 0; i < arrays[0].length; i += 4)
        if ([0, 1, 2].some((c) => Math.abs(arrays[0][i + c] - arrays[1][i + c]) > 8)) different++;
      return different;
    },
    [stopped.toString('base64'), settled.toString('base64')],
  );
  // NGL progressively refines edge antialiasing after camera motion stops.
  expect(changedPixels, 'After stopping rotation the rendered molecule must settle.').toBeLessThan(
    200,
  );
  for (const query of ['GLY35', 'A:35', '999999']) {
    await page.getByRole('textbox', { name: 'Find a residue', exact: true }).fill(query);
    await page.getByRole('button', { name: 'Focus residue', exact: true }).click();
    await expect(page.getByRole('status')).toContainText(
      query === '999999' ? 'No matching residue' : `Focused on ${query}`,
    );
    await page.getByRole('button', { name: 'Dismiss notification', exact: true }).click();
    await expect(page.getByRole('status')).toHaveCount(0);
  }
  await page.getByRole('button', { name: 'Full screen', exact: true }).click();
  await expect.poll(() => page.evaluate(() => !!document.fullscreenElement)).toBe(true);
  await page.getByRole('button', { name: 'Exit full screen', exact: true }).click();
  await expect.poll(() => page.evaluate(() => !!document.fullscreenElement)).toBe(false);
  expect(errors).toEqual([]);
});

test('audit: timeline boundaries, every speed, loop behavior and keyboard shortcuts', async ({
  page,
  request,
}) => {
  const errors = await ready(page);
  const demo = await (await request.get('/api/datasets/demo')).json();
  const slider = page.getByRole('slider', { name: 'Trajectory frame', exact: true });
  await page.getByRole('button', { name: 'Last frame', exact: true }).click();
  await expect(slider).toHaveValue(String(demo.n_frames - 1));
  await page.getByRole('button', { name: 'Next frame', exact: true }).click();
  await expect(slider).toHaveValue(String(demo.n_frames - 1));
  await page.getByRole('button', { name: 'Previous frame', exact: true }).click();
  await expect(slider).toHaveValue(String(demo.n_frames - 2));
  await page.getByRole('button', { name: 'First frame', exact: true }).click();
  await page.getByRole('button', { name: 'Previous frame', exact: true }).click();
  await expect(slider).toHaveValue('0');
  for (const speed of ['0.25', '0.5', '1', '2', '4']) {
    await page.getByRole('combobox', { name: 'Playback speed', exact: true }).selectOption(speed);
    await expect(page.getByRole('combobox', { name: 'Playback speed', exact: true })).toHaveValue(
      speed,
    );
    await page.getByRole('button', { name: 'Play trajectory', exact: true }).click();
    const initial = Number(await slider.inputValue());
    await expect.poll(async () => Number(await slider.inputValue())).toBeGreaterThan(initial);
    await page.getByRole('button', { name: 'Pause trajectory', exact: true }).click();
  }
  await page.getByRole('button', { name: 'Loop playback', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Loop playback', exact: true })).toHaveAttribute(
    'aria-pressed',
    'false',
  );
  await slider.fill(String(demo.n_frames - 2));
  await page.getByRole('button', { name: 'Play trajectory', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Play trajectory', exact: true })).toBeVisible();
  await expect(slider).toHaveValue(String(demo.n_frames - 1));
  await page.getByRole('button', { name: 'Loop playback', exact: true }).click();
  await slider.fill(String(demo.n_frames - 2));
  await page.getByRole('button', { name: 'Play trajectory', exact: true }).click();
  await expect.poll(async () => Number(await slider.inputValue())).toBeLessThan(30);
  await page.getByRole('button', { name: 'Pause trajectory', exact: true }).click();
  await page.getByRole('button', { name: 'First frame', exact: true }).click();
  await idleKeyboard(page);
  await page.keyboard.press('ArrowRight');
  await expect(slider).toHaveValue('1');
  await page.keyboard.press('ArrowLeft');
  await expect(slider).toHaveValue('0');
  await page.keyboard.press('Space');
  await expect(page.getByRole('button', { name: 'Pause trajectory', exact: true })).toBeVisible();
  await page.keyboard.press('Space');
  await expect(page.getByRole('button', { name: 'Play trajectory', exact: true })).toBeVisible();
  await page.keyboard.press('m');
  await expect(page.getByRole('button', { name: 'Stop picking', exact: true })).toBeVisible();
  await page.keyboard.press('Escape');
  await expect(page.getByRole('button', { name: 'Stop picking', exact: true })).toHaveCount(0);
  await page.keyboard.press('?');
  await expect(page.getByRole('dialog', { name: 'Meet DynaMol.' })).toBeVisible();
  await page.keyboard.press('Escape');
  await idleKeyboard(page);
  await page.keyboard.press('f');
  await page.getByRole('textbox', { name: 'Find a residue', exact: true }).fill('1');
  const paused = await slider.inputValue();
  await page.keyboard.press('ArrowRight');
  await expect(slider).toHaveValue(paused);
  expect(errors).toEqual([]);
});

test('audit: angle and dihedral, atom removal, plot hover and seeking, tab selection and deletion', async ({
  page,
  request,
}, info) => {
  const errors = await ready(page);
  const demo = await (await request.get('/api/datasets/demo')).json();
  const atoms: Atom[] = demo.atoms.filter((atom: Atom) => atom.name === 'CA').slice(0, 4);
  await page.getByRole('textbox', { name: 'Find an atom', exact: true }).fill('zzzz-no-atom');
  await expect(page.locator('.atom-search-results')).toHaveText('No atoms found.');
  await page.getByRole('button', { name: 'Clear atom search', exact: true }).click();
  await expect(page.getByRole('textbox', { name: 'Find an atom', exact: true })).toHaveValue('');
  await pick(page, atoms[0]);
  await page
    .getByRole('button', { name: `Remove atom ${atoms[0].index + 1}`, exact: true })
    .click();
  await expect(page.locator('.atom-slot.filled')).toHaveCount(0);
  await pick(page, atoms[0]);
  await page.getByRole('button', { name: 'Clear selection', exact: true }).click();
  await expect(page.locator('.atom-slot.filled')).toHaveCount(0);
  await page.getByRole('button', { name: 'Pick atoms in the scene', exact: true }).click();
  await page.getByRole('button', { name: 'Stop picking', exact: true }).click();
  for (const [kind, count] of [
    ['angle', 3],
    ['dihedral', 4],
  ] as const) {
    await page
      .locator('.measure-types')
      .getByRole('button', { name: new RegExp(`${kind}$`, 'i') })
      .click();
    for (const atom of atoms.slice(0, count)) await pick(page, atom);
    const response = page.waitForResponse(
      (r) => r.url().endsWith('/measurements') && r.request().postDataJSON()?.kind === kind,
    );
    await page.getByRole('button', { name: 'Plot over time', exact: true }).click();
    const result = await (await response).json();
    expect(result.kind).toBe(kind);
    expect(result.values).toHaveLength(demo.n_frames);
    expect(result.values.every(Number.isFinite)).toBe(true);
    await expect(page.locator('.plot-readout strong')).toContainText(result.values[0].toFixed(2));
    const exported = page.waitForEvent('download');
    await page.getByRole('button', { name: 'Export measurement CSV', exact: true }).click();
    const csv = await readFile((await (await exported).path())!, 'utf8');
    expect(csv.split('\n')[0]).toBe(`time_ps,${kind}_degrees`);
    await info.attach(`${kind}.csv`, { body: csv, contentType: 'text/csv' });
  }
  await expect(page.locator('.measurement-tab')).toHaveCount(3);
  const chart = page.locator('.chart-wrap svg');
  const bounds = (await chart.boundingBox())!;
  await page.mouse.move(bounds.x + bounds.width * 0.51, bounds.y + 30);
  await expect(page.locator('.plot-readout')).toContainText('At cursor');
  await chart.click({ position: { x: bounds.width * 0.51, y: 30 } });
  const frame = Number(
    await page.getByRole('slider', { name: 'Trajectory frame', exact: true }).inputValue(),
  );
  expect(frame).toBeGreaterThan(demo.n_frames * 0.4);
  expect(frame).toBeLessThan(demo.n_frames * 0.6);
  await page.mouse.move(3, 3);
  await expect(page.locator('.plot-readout')).toContainText('Current frame');
  await page.locator('.measurement-tab').first().getByRole('button').first().click();
  await expect(page.locator('.measurement-tab').first()).toHaveClass(/active/);
  while (await page.locator('.measurement-tab').count())
    await page
      .locator('.measurement-tab')
      .first()
      .getByRole('button', { name: /^Remove / })
      .click();
  await expect(
    page.getByRole('button', { name: 'Export measurement CSV', exact: true }),
  ).toBeDisabled();
  await expect(page.getByRole('button', { name: 'Show measurements', exact: true })).toBeDisabled();
  await page.getByRole('button', { name: 'Close measurement inspector', exact: true }).click();
  await expect(page.locator('.measurement-panel')).toHaveCount(0);
  await page.getByRole('button', { name: 'Toggle measurement inspector', exact: true }).click();
  await expect(page.locator('.measurement-panel')).toBeVisible();
  await page.getByRole('button', { name: 'Close measurement inspector', exact: true }).click();
  await page.getByRole('button', { name: 'Add a measurement', exact: false }).click();
  await expect(page.locator('.measurement-panel')).toBeVisible();
  await page.getByRole('button', { name: 'Close measurement inspector', exact: true }).click();
  await page.getByRole('button', { name: 'Measure', exact: true }).click();
  await expect(page.locator('.measurement-panel')).toBeVisible();
  expect(errors).toEqual([]);
});

test('audit: Escape exits picking after the pick button receives focus', async ({ page }) => {
  await ready(page);
  await page.getByRole('button', { name: 'Pick atoms in the scene', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Stop picking', exact: true })).toBeVisible();
  await page.keyboard.press('Escape');
  await expect(page.getByRole('button', { name: 'Stop picking', exact: true })).toHaveCount(0);
});

test('audit: keyboard picking pauses a playing trajectory', async ({ page }) => {
  await ready(page);
  await page.getByRole('button', { name: 'Play trajectory', exact: true }).click();
  await idleKeyboard(page);
  await page.keyboard.press('m');
  await expect(page.getByRole('button', { name: 'Stop picking', exact: true })).toBeVisible();
  await expect(page.getByRole('button', { name: 'Play trajectory', exact: true })).toBeVisible();
});

test('audit: real canvas atom hover, click selection, deselection and double-click focus', async ({
  page,
}, info) => {
  const errors = await ready(page);
  await page.getByRole('button', { name: 'Ball & stick', exact: true }).click();
  await page.getByRole('button', { name: 'Pick atoms in the scene', exact: true }).click();
  const image = await pixels(page);
  const candidates = await page.evaluate(async (base64) => {
    const bitmap = await createImageBitmap(
      await (await fetch(`data:image/png;base64,${base64}`)).blob(),
    );
    const canvas = document.createElement('canvas');
    canvas.width = bitmap.width;
    canvas.height = bitmap.height;
    const context = canvas.getContext('2d')!;
    context.drawImage(bitmap, 0, 0);
    const data = context.getImageData(0, 0, canvas.width, canvas.height).data;
    const points: { x: number; y: number; radius: number }[] = [];
    for (let y = 100; y < canvas.height - 100; y += 12)
      for (let x = 100; x < canvas.width - 100; x += 12) {
        const i = (y * canvas.width + x) * 4;
        if (
          Math.max(data[i], data[i + 1], data[i + 2]) > 130 &&
          Math.max(data[i], data[i + 1], data[i + 2]) -
            Math.min(data[i], data[i + 1], data[i + 2]) >
            40
        )
          points.push({ x, y, radius: (x - canvas.width / 2) ** 2 + (y - canvas.height / 2) ** 2 });
      }
    return points.sort((a, b) => a.radius - b.radius).slice(0, 40);
  }, image.toString('base64'));
  const box = (await page.locator('.molecular-viewer__canvas canvas').boundingBox())!;
  let picked: { x: number; y: number } | undefined;
  for (const point of candidates) {
    await page.mouse.move(box.x + point.x, box.y + point.y);
    await page.waitForTimeout(100);
    if (await page.getByRole('tooltip').isVisible()) {
      picked = point;
      break;
    }
  }
  expect(picked, 'At least one visible atom must be hoverable.').toBeTruthy();
  await expect(page.getByRole('tooltip')).toContainText('click to select');
  await page.mouse.click(box.x + picked!.x, box.y + picked!.y);
  await expect(page.locator('.atom-slot.filled')).toHaveCount(1);
  await page.mouse.click(box.x + picked!.x, box.y + picked!.y);
  await expect(page.locator('.atom-slot.filled')).toHaveCount(0);
  await page.getByRole('button', { name: 'Stop picking', exact: true }).click();
  await page.mouse.move(box.x + picked!.x, box.y + picked!.y);
  await page.waitForTimeout(150);
  await expect(page.getByRole('tooltip')).toBeVisible();
  const before = await page.locator('.molecular-viewer__canvas canvas').screenshot();
  await page.mouse.dblclick(box.x + picked!.x, box.y + picked!.y);
  await page.waitForTimeout(700);
  expect((await pixels(page, info, 'double-click-focused-atom')).equals(before)).toBe(false);
  expect(errors).toEqual([]);
});

test('audit: invalid hydrogen-bond selection reports a dismissible error and keeps the workspace', async ({
  page,
  request,
}) => {
  await ready(page);
  const dataset = await (await request.get('/api/datasets/demo')).json();
  await page
    .locator('.measure-types')
    .getByRole('button', { name: /Hydrogen bond$/ })
    .click();
  for (const atom of dataset.atoms.filter((a: Atom) => a.name === 'CA').slice(0, 3))
    await pick(page, atom);
  await page.getByRole('button', { name: 'Plot over time', exact: true }).click();
  await expect(page.getByRole('alert')).toBeVisible();
  await expect(page.getByRole('alert')).toContainText(/hydrogen|donor/i);
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  await page.getByRole('button', { name: 'Dismiss error', exact: true }).click();
  await expect(page.getByRole('alert')).toHaveCount(0);
  await page.getByRole('button', { name: 'Clear selection', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Plot over time', exact: true })).toBeDisabled();
});

test('audit: import chooser actions, cancel, drag and drop, form bounds and timing override', async ({
  page,
  request,
}) => {
  const errors = await ready(page);
  const topology = await (await request.get('/api/datasets/demo/topology')).body();
  const open = () => page.getByRole('button', { name: 'Open files', exact: true }).click();
  await open();
  const dialog = page.getByRole('dialog', { name: 'Bring your molecules.' });
  await expect(dialog.getByRole('button', { name: 'Open in DynaMol', exact: true })).toBeDisabled();
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toHaveCount(0);
  await open();
  await dialog.getByRole('button', { name: 'Close import', exact: true }).click();
  await open();
  const chooserPromise = page.waitForEvent('filechooser');
  await dialog.getByRole('button', { name: /Choose structure/ }).click();
  await (
    await chooserPromise
  ).setFiles({ name: 'e2e-ui-audit-selected.pdb', mimeType: 'chemical/x-pdb', buffer: topology });
  await expect(dialog.getByRole('button', { name: /e2e-ui-audit-selected.pdb/ })).toBeVisible();
  const trajectoryChooser = page.waitForEvent('filechooser');
  await dialog.getByRole('button', { name: /Add trajectory/ }).click();
  await (
    await trajectoryChooser
  ).setFiles({ name: 'e2e-ui-audit-trajectory.pdb', mimeType: 'chemical/x-pdb', buffer: topology });
  await expect(dialog.getByRole('button', { name: /e2e-ui-audit-trajectory.pdb/ })).toBeVisible();
  await dialog.getByRole('spinbutton', { name: /^Read every/ }).fill('0');
  await dialog.getByRole('button', { name: 'Open in DynaMol', exact: true }).click();
  expect(
    await dialog
      .getByRole('spinbutton', { name: /^Read every/ })
      .evaluate((element: HTMLInputElement) => element.validity.rangeUnderflow),
  ).toBe(true);
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await open();
  const name = `e2e-ui-audit-dropped-${Date.now()}`;
  const files = await page.evaluateHandle(
    ({ text, name }) => {
      const data = new DataTransfer();
      data.items.add(new File([text], `${name}.pdb`, { type: 'chemical/x-pdb' }));
      return data;
    },
    { text: topology.toString(), name },
  );
  await dialog.locator('.dropzone').dispatchEvent('drop', { dataTransfer: files });
  await files.dispose();
  await expect(dialog.getByRole('button', { name: new RegExp(`${name}.pdb`) })).toBeVisible();
  await dialog.getByRole('spinbutton', { name: /^Read every/ }).fill('2');
  await dialog.getByRole('spinbutton', { name: /^Time between original frames/ }).fill('0.5');
  const response = page.waitForResponse(
    (r) => r.url().endsWith('/api/datasets/upload') && r.request().method() === 'POST',
  );
  await dialog.getByRole('button', { name: 'Open in DynaMol', exact: true }).click();
  const loaded = await (await response).json();
  expect(loaded.n_frames).toBe(1);
  expect(loaded.n_atoms).toBeGreaterThan(1000);
  await expect(dialog).toHaveCount(0);
  await expect(page.locator('.structure-card h2')).toHaveText(name);
  expect(errors).toEqual([]);
});

test('audit: Studio source chooser, invalid fetch, all default and advanced configuration controls', async ({
  page,
}) => {
  const errors = await ready(page);
  const dialog = await studio(page);
  await expect(dialog.getByLabel('Simulation name', { exact: true })).toHaveValue(
    'My molecular journey',
  );
  await expect(dialog.getByLabel('Duration', { exact: false })).toHaveValue('10');
  await expect(dialog.getByLabel('Temperature', { exact: false })).toHaveValue('300');
  await expect(dialog.getByLabel('Solvent environment')).toHaveValue('implicit');
  await dialog.getByRole('button', { name: /Advanced controls/ }).click();
  const fields = [
    ['Integration step', '1'],
    ['Save every', '25'],
    ['Random seed', '2026'],
    ['Friction', '2'],
    ['Equilibration', '0'],
  ] as const;
  for (const [label, value] of fields) {
    await dialog.getByRole('spinbutton', { name: new RegExp(`^${label}`) }).fill(value);
    await expect(dialog.getByRole('spinbutton', { name: new RegExp(`^${label}`) })).toHaveValue(
      value,
    );
  }
  await expect(dialog.getByLabel('Box padding', { exact: false })).toBeDisabled();
  await dialog.getByLabel('Minimize energy before equilibration', { exact: true }).uncheck();
  await expect(
    dialog.getByLabel('Minimize energy before equilibration', { exact: true }),
  ).not.toBeChecked();
  await dialog.getByLabel('Minimize energy before equilibration', { exact: true }).check();
  await dialog.getByRole('button', { name: /GROMACS/ }).click();
  await expect(dialog.getByLabel('Solvent environment')).toHaveValue('explicit');
  await expect(dialog.getByLabel('Solvent environment')).toBeDisabled();
  await dialog.getByLabel('Box padding', { exact: false }).fill('1.2');
  await expect(dialog.getByLabel('Box padding', { exact: false })).toHaveValue('1.2');
  await dialog.getByRole('button', { name: /OpenMM/ }).click();
  await dialog.getByLabel('Solvent environment').selectOption('implicit');
  await expect(dialog.getByLabel('Box padding', { exact: false })).toBeDisabled();
  await dialog.getByRole('button', { name: /Advanced controls/ }).click();
  await expect(dialog.locator('.advanced-fields')).toHaveCount(0);
  await dialog.getByRole('button', { name: 'Preparation options', exact: true }).click();
  for (const label of [
    'Add missing heavy atoms',
    'Refine side-chain rotamers and clashes',
    'Remove existing waters before preparation',
    'Remove ligands and other non-protein residues',
  ]) {
    const checkbox = dialog.getByRole('checkbox', { name: new RegExp(`^${label}`) });
    const initial = await checkbox.isChecked();
    await checkbox.setChecked(!initial);
    expect(await checkbox.isChecked()).toBe(!initial);
    await checkbox.setChecked(initial);
  }
  await dialog.getByLabel('Preparation pH', { exact: true }).fill('6.5');
  await expect(dialog.getByLabel('Preparation pH', { exact: true })).toHaveValue('6.5');
  await dialog.getByRole('button', { name: 'Preparation options', exact: true }).click();
  await dialog.getByRole('tab', { name: 'Fetch', exact: true }).click();
  await expect(dialog.getByRole('button', { name: 'Fetch', exact: true })).toBeDisabled();
  await dialog.getByLabel('PDB accession', { exact: true }).fill('invalid-code');
  await dialog.getByRole('button', { name: 'Fetch', exact: true }).click();
  await expect(dialog.getByRole('alert')).toContainText('PDB');
  await dialog.getByRole('tab', { name: 'SMILES', exact: true }).click();
  await expect(
    dialog.getByRole('button', { name: 'Build 3D structure & view', exact: true }),
  ).toBeDisabled();
  await dialog.getByRole('tab', { name: 'Upload', exact: true }).click();
  const chooser = page.waitForEvent('filechooser');
  await dialog.getByRole('button', { name: /Drop a structure, or choose a file/ }).click();
  await (
    await chooser
  ).setFiles({
    name: `e2e-ui-audit-chooser-${Date.now()}.pdb`,
    mimeType: 'chemical/x-pdb',
    buffer: await readFile(fixture('six_residues_missing_atom.pdb')),
  });
  await expect(page.locator('.source-summary strong')).toContainText('e2e-ui-audit-chooser');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  expect(errors).toEqual([]);
});
