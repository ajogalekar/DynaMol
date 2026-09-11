import { test, expect, type Page } from '@playwright/test';

type Atom = { index: number; name: string; residue: string; resid: number; element: string };
async function pick(page: Page, atom: Atom) {
  await page
    .getByRole('textbox', { name: 'Find an atom', exact: true })
    .fill(`${atom.residue}${atom.resid} · ${atom.name}`);
  await page
    .locator('.atom-search-results')
    .getByRole('button', { name: new RegExp(`#${atom.index + 1} · ${atom.element}$`) })
    .click();
}
async function ready(page: Page) {
  await page.goto('/');
  await expect(page.getByRole('button', { name: 'Save snapshot', exact: true })).toBeEnabled();
  await expect(page.locator('.canvas-loading')).toHaveCount(0);
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
}

test('live measurements appear before Plot, update with frames, and reset with selection', async ({
  page,
  request,
}) => {
  await ready(page);
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  const demo = await (await request.get('/api/datasets/demo')).json();
  const atoms: Atom[] = demo.atoms.filter((atom: Atom) => atom.name === 'CA').slice(0, 4);
  const live = page.locator('.live-measurement__value');
  const slider = page.getByRole('slider', { name: 'Trajectory frame', exact: true });
  for (const [kind, count] of [
    ['distance', 2],
    ['angle', 3],
    ['dihedral', 4],
  ] as const) {
    await page
      .locator('.measure-types')
      .getByRole('button', { name: new RegExp(`${kind}$`, 'i') })
      .click();
    const response = page.waitForResponse(
      (res) =>
        res.url().endsWith('/measurement-preview') && res.request().postDataJSON()?.kind === kind,
    );
    for (const atom of atoms.slice(0, count)) await pick(page, atom);
    const data = await (await response).json();
    const format = (frame: number) =>
      `${data.values[frame].toFixed(kind === 'distance' ? 2 : 1)} ${data.unit}`;
    await expect(live).toHaveText(format(0));
    await expect(page.locator('.molecular-viewer__measurement--preview')).toContainText(format(0));
    await expect(page.locator('.measurement-tab')).toHaveCount(1);
    await slider.fill(String(demo.n_frames - 1));
    await expect(live).toHaveText(format(demo.n_frames - 1));
    await expect(page.locator('.live-measurement__heading')).toContainText(
      `Frame ${demo.n_frames}`,
    );
    await page
      .getByRole('button', { name: `Remove atom ${atoms[count - 1].index + 1}`, exact: true })
      .click();
    await expect(live).toHaveCount(0);
    await expect(page.locator('.molecular-viewer__measurement--preview')).toHaveCount(0);
    await slider.fill('0');
  }
  await page
    .locator('.measure-types')
    .getByRole('button', { name: /distance$/i })
    .click();
  for (const atom of atoms.slice(0, 2)) await pick(page, atom);
  await expect(live).toBeVisible();
  const responses: string[] = [];
  page.on('request', (req) => {
    if (req.url().endsWith('/measurement-preview')) responses.push(req.url());
  });
  await page.getByRole('button', { name: 'Play trajectory', exact: true }).click();
  await expect.poll(async () => Number(await slider.inputValue())).toBeGreaterThan(4);
  await page.getByRole('button', { name: 'Pause trajectory', exact: true }).click();
  await expect(page.locator('.live-measurement__heading')).toContainText(
    `Frame ${Math.round(Number(await slider.inputValue())) + 1}`,
  );
  expect(responses).toHaveLength(0);
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  await page.getByRole('button', { name: 'Plot over time', exact: true }).click();
  await expect(page.locator('.measurement-tab')).toHaveCount(2);
  await expect(live).toHaveCount(0);
  await expect(page.locator('.molecular-viewer__measurement--preview')).toHaveCount(0);
  expect(errors).toEqual([]);
});

test('live readout rejects invalid hydrogen-bond selection immediately and handles a new selection', async ({
  page,
  request,
}) => {
  await ready(page);
  const demo = await (await request.get('/api/datasets/demo')).json();
  const atoms: Atom[] = demo.atoms.filter((atom: Atom) => atom.name === 'CA').slice(0, 3);
  await page
    .locator('.measure-types')
    .getByRole('button', { name: /Hydrogen bond$/ })
    .click();
  for (const atom of atoms) await pick(page, atom);
  await expect(page.locator('.live-measurement__error')).toContainText(
    'The second atom must be hydrogen',
  );
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
  await page.getByRole('button', { name: 'Clear selection', exact: true }).click();
  await expect(page.locator('.live-measurement')).toHaveCount(0);
  await page.route('**/measurement-preview', async (route) => {
    await new Promise((resolve) => setTimeout(resolve, 350));
    await route.continue();
  });
  await page
    .locator('.measure-types')
    .getByRole('button', { name: /Distance$/ })
    .click();
  for (const atom of atoms.slice(0, 2)) await pick(page, atom);
  await page.getByRole('button', { name: 'Clear selection', exact: true }).click();
  await expect(page.locator('.live-measurement')).toHaveCount(0);
  for (const atom of atoms.slice(1, 3)) await pick(page, atom);
  const expected = await (
    await request.post(`/api/datasets/${demo.id}/measurement-preview`, {
      data: { kind: 'distance', atoms: atoms.slice(1, 3).map((atom) => atom.index) },
    })
  ).json();
  await expect(page.locator('.live-measurement__value')).toHaveText(
    `${expected.values[0].toFixed(2)} Å`,
  );
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
});

test('live scene labels follow molecular visibility while inspector retains selected geometry', async ({
  page,
  request,
}) => {
  await ready(page);
  const demo = await (await request.get('/api/datasets/demo')).json();
  const atoms: Atom[] = demo.atoms.filter((atom: Atom) => atom.name === 'CA').slice(0, 2);
  for (const atom of atoms) await pick(page, atom);
  const label = page.locator('.molecular-viewer__measurement--preview');
  await expect(label).toBeVisible();
  await page.getByRole('button', { name: 'Hide protein & nucleic acids', exact: true }).click();
  await expect(label).toHaveCount(0);
  await expect(page.locator('.live-measurement__value')).toBeVisible();
  await page.getByRole('button', { name: 'Show protein & nucleic acids', exact: true }).click();
  await expect(label).toBeVisible();
  await page.getByRole('button', { name: 'Close measurement inspector', exact: true }).click();
  await expect(label).toBeVisible();
  await expect(page.locator('.measurement-tab')).toHaveCount(1);
});
